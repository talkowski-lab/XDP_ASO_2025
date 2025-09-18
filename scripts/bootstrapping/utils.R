################################################
# Dependencies
################################################
suppressMessages(library(BSDA))
suppressMessages(library(Hmisc))
suppressMessages(library(DESeq2))
suppressMessages(library(MASS))
suppressMessages(library(aod))
suppressMessages(library(genefilter))
suppressMessages(library(future.apply))
suppressMessages(library(progressr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(lme4))
suppressMessages(library(matrixStats))
suppressMessages(library(pbkrtest))
suppressMessages(library(sva))
suppressMessages(library(tidyverse))

##############
# Important Directories
MAIN_DIR <- '/data/talkowski/Samples/XDP/ASO_dSVA/sid_results'
DEGRESULTS_DIR <- sprintf('%s/DEGResults', MAIN_DIR)
FROZEN_DIR <- '/data/talkowski/Samples/XDP/ASO_dSVA/Results_draft1'
# Gene IDs/lists
TAF1_ENSEMBLID <- 'ENSG00000147133'
# Likely artifacts/errors since they have insignificant and massive LFCs often
VOLCANO_EXCLUDE_GENES <- 
    c(
        'ENSG00000259785',
        'ENSG00000224201',
        'ENSG00000185304',
        'ENSG00000235649'
    )

##############
# Identify SVs using SVA and produce corrected counts (DESeq2)
run_sva_from_dds <- function(
    dds, 
    design_var){
    # run SVA to compute surrogate variables (SVs) for downstream batch correction
    # adds SVs as varaibles to the DESeq2 object and updates the formula with SVs
    svseq <- svaseq(counts(dds,
                           normalized=TRUE
                    ),
                    model.matrix(formula(sprintf("~ %s", design_var)),
                                 colData(dds)
                    ),
                    model.matrix(~ 1, 
                                 colData(dds)
                    ) 
             )  
    # renames rows/cols of esitmated SV matrix
    SV <- svseq$sv
    colnames(SV) <- paste0("SV", 1:svseq$n.sv)
    rownames(SV) <- colnames(dds)
    colData(dds) <- cbind(colData(dds), SV)
    # update design matrix for DESeq to account for all SVs
    design(dds) <- formula(sprintf("~ %s + %s", 
                                   design_var,
                                   paste("SV", 
                                         1:ncol(SV), 
                                         sep="", 
                                         collapse=" + "
                                   )
                           )
    )
    return(dds)
}

get_corrected_counts_from_dds <- function(
    dds,
    design_var='Condition'){
    # Given a DESeq2 object with SVs, return SVA-corrected normalized counts
    dds <- DESeq(dds)
    # SV variable names
    sv_names <- grep('^SV[0-9]+', colnames(colData(dds)), value=T)
    # Design matrix for Condition including all SVs in the formula
    mod <- as.matrix(model.matrix(design(dds), data=colData(dds))[, sv_names])
    # beta coefficients for each SV as estimate by DESeq2 GLM
    beta <- as.matrix(coef(dds)[, sv_names])
    # Model estimates betas for b_geno + b_s1 * SV1 + b_s2 * SV2 + ... + b_sn SVn
    # here we subtract out the effects of each SV estimated by the model
    # This should provide us with "SV corrected" counts
    # the log2 transforms the raw counts to same scale the SVs are estimated on?
    clean_counts <- log2(counts(dds, normalized=T) + 1) - t(mod %*% t(beta))
    # Deal with -Inf and such caused by talking the log2(x+1) of counts above
    # clean_counts <- log2(2**clean_counts + 1)
    return(clean_counts)
}

run_sva_from_raw <- function(
    raw_counts,
    metadata,
    design_var='Condition'){
    # run SVA to compute surrogate variables (SVs) for downstream batch correction 
    # takes a count matrix + meta and return a corrected counts + meta + SVs 
    svseq_dds <- 
        metadata %>%
        mutate(!!design_var:=as.factor(get(design_var))) %>%
        DESeqDataSetFromMatrix(
            countData=raw_counts,
            colData=.,
            design=formula(sprintf("~ %s", design_var))
        ) %>%  
        estimateSizeFactors() %>%
        run_sva_from_dds(design_var=design_var)
    # get corrected counts using SVs
    corrected_counts <- 
        get_corrected_counts_from_dds(
            svseq_dds,
            design_var
        )
    # return everything
    return(
        list(
            corrected_counts=corrected_counts,
            metadata=svseq_dds %>% colData() %>% as_tibble()
        )
    )
}
# Identify SVs using SVA and produce corrected counts (GLM)
run_glm_sva_from_raw <- function(
    counts,
    sample_metadata,
    design_var='Condition',
    return_all=FALSE){
    # Run SVA to compute SVs + correct counts with them using GLM estimates of
    # their beta + residuals
    # Takes a count matrix + meta and return a corrected counts + meta + SVs 
    glm_models <- 
        sample_metadata %>%
        mutate(!!design_var:=as.factor(get(design_var))) %>%
        buildDEmodel(
            counts,
            .,
            method="glm",
            estimate_offset=TRUE,
            apply_sva=TRUE,
            test_terms=c(design_var),
            distro_family="NB"
        )
    if (return_all) {
        return(glm_models)
    } else {
        # return SV corrected counts and metadata + SVs
        corrected_counts <- glm_models@deSV
        colnames(corrected_counts) <- glm_models@colData$SampleID
        return(
            list(
                corrected_counts=corrected_counts,
                metadata=glm_models@colData
            )
        )
    }
}

##############
# DESeq2 functions
make_contrast <- function(meta, analysis_params){
    # create a contrast using the appropruate condition as a baseline (denominator)
    if (!(is.null(analysis_params$control_label))){
        contrast <- c(analysis_params$design_var, 
                      analysis_params$condition_label,
                      analysis_params$control_label
        )
    } else {
        if ('NO.CON' %in% meta$Condition){
            control_label <- 'NO.CON'
        } else if (any(grepl('CON', unique(meta$Condition)))){
            control_label <- grep('CON', unique(meta$Condition), value=TRUE)[[1]]
        } else if ('NO.XDP' %in% meta$Condition){
            control_label <- 'NO.XDP'
        } else {
            control_label <- unique(meta$Condition)[[1]]
        }
        condition_label <- setdiff(unique(meta$Condition), control_label)
        contrast <- c(analysis_params$design_var, 
                      condition_label,
                      control_label
        )
    }
    return(contrast)
}
deseq_wrapper <- function(counts, meta, analysis_params,
                          toSave=FALSE, output_name=NA){
    # Wrapper to run DESeq2 for DEG analysis on the input
    # Expects meta to have column named 'Condition' to do the comparison along
    contrast <- make_contrast(meta, analysis_params)
    message(sprintf("Doing DE with the contrast %s", paste(contrast, collapse=", ")))
    message(table(meta[, analysis_params$design_var]))
    # create dds object from input data
    dds <- DESeqDataSetFromMatrix(countData=counts,
                                  colData=meta,
                                  design=formula(sprintf("~ %s", analysis_params$design_var))
    )  
    colData(dds)[, analysis_params$design_var] <- factor(colData(dds)[, analysis_params$design_var])
    dds <- estimateSizeFactors(dds)  
    dds <- run_sva(dds, design_var=analysis_params$design_var) # add SVs to the deseq object
    message(as.character(design(dds)))  # should be like `~ Condition + SV1 + SV2 + SV3`
    dds <- DESeq(dds) # run deseq and return the results
    res <- results(dds, 
                   alpha=analysis_params$alpha, 
                   pAdjustMethod="fdr",
                   independentFiltering=FALSE, 
                   contrast=contrast
    )
    message(summary(res))
    if (toSave){
        saveRDS(counts, sprintf('%s_DESeq2.counts.rds', output_name))
        saveRDS(meta, sprintf('%s_DESeq2.meta.rds', output_name))
        saveRDS(dds, sprintf('%s_DESeq2.dds.rds', output_name))
    }
    return(res[complete.cases(res),])
}

convert_deseq_to_tibble <- function(deseq_obj){
    # convert a DESeq results object to a tibble
    deseq_obj <- 
        deseq_obj %>% 
        as.data.frame() %>% 
        rownames_to_column(var='EnsemblID') %>% 
        as_tibble()
    return(deseq_obj)
}

##############
# Object class for GLM results for a set of genes 
setClass(
    "GLMres",
    representation(
        y="matrix",
        colData="data.frame",  # sample metadata
        design="formula",      # formula for the GLM
        beta="matrix",         # Estimated beta coefficients
        stderr="matrix",       # estimated std err of coefficients
        pval="matrix",         # pvalue for each gene
        fdr="matrix",          # FDR for each gene
        vcov="data.frame",     # variance-covariance matrix 
        resid="matrix",        # GLM residuals
        deSV="matrix",         # matrix of SVs estimated
        status="matrix",       # binary matrix of each genes GLM convergence status
        args="list",           # input arguments
        models="list"          # lme4 model objects for each gene
    )
)

setMethod(
    f="show",
    signature="GLMres",
    definition=
        function(object){
            print("class: GLMres")
            print("slots: y, colData, design, beta, stderr, pval, fdr, vcov, resid, deSV, status, args, models")
            print(sprintf("deSV dim: %s %s", nrow(object@deSV), ncol(object@deSV)))
        }
)

##############
# GLM functions
fit_single_model <- function(
    x, 
    meta_data=NULL,
    approach=NULL,
    fit_formula=NULL,
    offset_str=NULL,
    regression_caller=NULL,
    distro=NULL,
    defaultlink=NULL){
    # Get sample labels and model formula to fit
    meta_data[["Y"]] <- x
    fit_formula_final <- paste0("Y", fit_formula)
    # Run the appropriate model on the input i.e. raw RNASeq counts + SVs
    res <- tryCatch({
             if (is.null(offset_str)){
               if (approach == 0){
                 # non-NB distributions without offset called by glm,
                 # or non-gaussian and non-NB distribution without offset 
                 # called by glmer
                 fit_formula_final %>%
                     formula() %>%
                     regression_caller(data=meta_data, 
                                       family=distro
                     )
               } else if ((approach == 1) || (approach == 2)){
                 # NB distribution without offset called by glm.nb,
                 # or gaussian distribution without offset called by lmer,
                 # or NB distribution without offset called by glmer.nb
                 fit_formula_final %>%
                     formula() %>%
                     regression_caller(data=meta_data)
               }
             } else {
               if (approach == 0){
                 # non-NB distributions with offset called by glm,
                 # or non-gaussian distribution with offset called by glmer
                 fit_formula_final %>%
                     formula() %>%
                     regression_caller(data=meta_data, 
                                       family=distro,
                                       offset=defaultlink(meta_data[[offset_str]])
                     )
               } else if (approach == 1){
                 # Gaussian distribution with offset called by lmer,
                 # or NB distribution with offset called by glmer.nb
                 fit_formula_final %>%
                     formula() %>%
                     regression_caller(data=meta_data, 
                                       offset=defaultlink(meta_data[[offset_str]])
                     )
               } else if (approach == 2){
                 # NB distribution without offset called by glm.nb
                 fit_formula_final %>%
                     paste0(" + offset(log(", offset_str, "))") %>%
                     formula() %>%
                     regression_caller(data=meta_data)
               }
             }
           }, 
           error=function(e){ return(NULL) }
    )
    return(res)
}

merge_vcov_from_model_list <- function(model_list){
    # Merge Variance-CoVariance matrices for each gene model
    tt=lapply(model_list, vcov)
    tt1=lapply(names(tt), 
               FUN=function(x){
                   tmp=as.data.frame(as.matrix(tt[[x]]))
                   rownames(tmp)=NULL 
                   tmp[["name"]]=x 
                   return(tmp)
               }
    ) 
    tt2=do.call(rbind, tt1)
    return(tt2)
}

get_fdr <- function(
    res,
    beta=NULL,
    pval=NULL,
    fdr=NULL,
    term,
    class,
    keyword=NULL){
    # Reorganize DEG data from the model into a single matrix 
    if (!is.null(beta)){
        fc <- read.table(beta, header=T, row.names=1, sep="\t", check.names=F)
        pv <- read.table(pval, header=T, row.names=1, sep="\t", check.names=F)
        bh <- read.table(fdr,  header=T, row.names=1, sep="\t", check.names=F)
    } else {
        fc <- as.data.frame(res@beta)
        pv <- as.data.frame(res@pval)
        bh <- as.data.frame(res@fdr)
    }
    new0 <- fc[,grep(term, colnames(fc)), drop=F]
    colnames(new0) <- paste0("LFC.",colnames(fc)[grep(term, colnames(fc))])
    new1 <- pv[,grep(term, colnames(pv)), drop=F]
    colnames(new1) <- paste0("pval.",colnames(pv)[grep(term, colnames(pv))])
    new2 <- bh[,grep(term, colnames(bh)), drop=F]
    colnames(new2) <- paste0("FDR.",colnames(bh)[grep(term, colnames(bh))])
    if (!is.null(keyword)){
        new0 <- new0[keyword, , drop=F]
        new1 <- new1[keyword, , drop=F]
        new2 <- new2[keyword, , drop=F]
    }
    new <- as.data.frame(new2)
    new <- cbind(new, new0, new1)
    #colnames(new) <- paste0("V",1:ncol(new))
    new$name <- as.vector(rownames(new))
    # new$class <- class
    rownames(new) <- NULL
    return(new)
}

build_DE_DESeq2_model <- function(
    irm,
    irc,
    fml,
    estimate_offset,
    offset_str,
    apply_sva,
    n_sv,
    test_terms,
    remove_terms,
    random_terms){
    message("Constructing DESeq2 object...")
    mode(irm) <- "integer"
    ddssva <- DESeqDataSetFromMatrix(irm, irc, ~1)
    # set offset variables
    if (estimate_offset){
        message("Estimating size factors...")
        ddssva <- estimateSizeFactors(ddssva)
        offset_str <- "buildDEmodel.sizeFactor"
        colData(ddssva)[[offset_str]] <- sizeFactors(ddssva)
    } else {
        if (is.null(offset_str)){
            message("No size factor information provided. Setting to 1...")
            offset_str <- "buildDEmodel.sizeFactor"
            colData(ddssva)[[offset_str]] <- rep(1, nrow(irc))
            sizeFactors(ddssva) <- colData(ddssva)[[offset_str]]
        } else {
            if (all(offset_str %in% colnames(sampleTable)) == FALSE){
                stop("'offset_str' must be found in the column names of 'sampleTable'.")
            }
            message("Using '", offset_str, "' column in the design matrix as size factors...")
            sizeFactors(ddssva) <- colData(ddssva)[[offset_str]]
        }
    }
    # Apply SVA to counts if specified
    if (apply_sva){
        message("Identifying surrogate variables...")
# set null model to include covariate to remove
        if (is.null(remove_terms)){
            null_model_str="1"
        } else {
            null_model_str=paste(c(remove_terms), collapse=" + ")
        }
        mod <- model.matrix(formula(fml),
                            colData(ddssva)
        )
        message("svaseq FULL model: ", formula(fml))
        mod0 <- model.matrix(formula(paste0("~",
                                            null_model_str
                                     )
                             ), 
                             colData(ddssva)
        )
        message("svaseq null model: ", formula(paste0("~",null_model_str)))
        svseq <- svaseq(counts(ddssva,
                               normalized=TRUE
                        ), 
                        mod,
                        mod0,
                        n_sv
        )
        for (i in 1:dim(svseq$sv)[2]){
            ddssva[[paste0("SV",i)]] <- svseq$sv[,i]
            fml <- paste0(fml, " + SV", i)
            sv_terms <- c(sv_terms, paste0("SV", i))
        }
    }
    # DESeq2 based approach for gene expression modelling and exon splicing
    message("Initiating modeling...")
    design(ddssva) <- formula(fml)
    message("Fitting ", method, " models with NB distribution: ", fml, " ...")
    ddssva <- DESeq(ddssva, modelMatrixType="standard")
    if (length(c(remove_terms, sv_terms)) > 0){
        message("Correcting against unwanted effects: ", 
                paste(c(remove_terms, sv_terms), collapse=", "),
                " ..."
        )
        m <- as.data.frame(model.matrix(formula(fml), 
                                        data=colData(ddssva)
                           )
        )
        for (keyword in remove_terms){
            remove_cols=c(remove_cols, which(grepl(paste0("^", keyword), colnames(m))))
        }
        for (svword in sv_terms){
            remove_cols=c(remove_cols, which(colnames(m) == svword))
        }
        X <- as.matrix(m[,remove_cols])
        beta <- as.matrix(coef(ddssva)[, remove_cols])
        y <- log2(counts(ddssva, normalized=T))
        cleany <- y - t(X %*% t(beta))
    } else {
        cleany <- log2(counts(ddssva, normalized=T))
    }
    message("Collecting results...")
    ddssva@metadata$deSV <- cleany
    message("Done.")
    message("Initiating modeling...")
    return(ddssva)
}

buildDEmodel <- function(
    feature_matrix,
    sampleTable,
    method="DESeq2",
    estimate_offset=TRUE,
    offset_str=NULL,
    apply_sva=TRUE,
    n_sv=NULL,
    test_terms=NULL,
    remove_terms=NULL,
    random_terms=NULL,
    distro_family="poisson"){
    # define object for output format for GLM/GLMM methods
    setClass("GLMres",
             representation(y="matrix",
                            colData="data.frame",
                            design="formula",
                            beta="matrix",
                            stderr="matrix",
                            pval="matrix",
                            fdr="matrix",
                            vcov="data.frame",
                            resid="matrix",
                            deSV="matrix",
                            status="data.frame",
                            args="list",
                            models="list"
             )
    )
    setMethod(f="show",
              signature="GLMres",
              definition=function(object){
                            print("class: GLMres")
                            print("slots: y, colData, design, beta, stderr, pval, fdr, vcov, resid, deSV, status, args, models")
                            print(sprintf("deSV dim: %s %s", 
                                          nrow(object@deSV), 
                                          ncol(object@deSV)
                                  )
                            )
                         }
    )
    # validate input terms
    if (all(test_terms %in% colnames(sampleTable)) == FALSE){
        stop("'test_terms' must be found in the column names of 'sampleTable'.")
    }
    if (all(remove_terms %in% colnames(sampleTable)) == FALSE){
        stop("'remove_terms' must be found in the column names of 'sampleTable'.")
    }
    if (all(random_terms %in% colnames(sampleTable)) == FALSE){
        stop("'random_terms' must be found in the column names of 'sampleTable'.")
    }
    if (all(remove_terms %in% test_terms) == FALSE){
        stop("'remove_terms' must be a subset of 'test_terms'.")
    }
    # set up inputs for the model
    irm <- as.matrix(feature_matrix)   # RNASeq counts
    irc <- as.data.frame(sampleTable)  # Sample metadata
    fml <- paste0("~", 
                  paste(c(test_terms), 
                        collapse=" + "
                  )
    )  # format model formula
    sv_terms <- c()
    remove_cols <- c()
    # Select approach and build model
    message(method, " approach selected.")
    if ((method == "deseq2") || (method == "DESeq2")){
        build_DE_DESeq2_model(irm,
                              irc,
                              fml,
                              estimate_offset,
                              offset_str,
                              apply_sva,
                              n_sv,
                              test_terms,
                              remove_terms,
                              random_terms
        )
    } else { # if not DESeq2 approach
        # Estimating GLM model offsets 
        if (estimate_offset){
            message("Estimating size factors...")
            offset_str <- "buildDEmodel.sizeFactor"
            irc[[offset_str]] <- estimateSizeFactorsForMatrix(irm)
        } else {
            if (is.null(offset_str)){
                message("No size factor information provided. Will NOT use 'offset' option during modeling.")
            } else {
                if (all(offset_str %in% colnames(sampleTable)) == FALSE){
                    stop("'offset_str' must be found in the column names of 'sampleTable'.")
                }
                message("Using '", offset_str, "' column in the design matrix as size factors...")
            }
        }
        # Choose distribution and set GLM caller accordingly
        if ((distro_family == "poisson") || 
            (distro_family == "Poisson")){
            if (!all(sapply(irm, `%%`, 1) == 0)){
                stop(distro_family, " distribution is selected. The 'feature_matrix' must be all integers.")
            }
            distro_family <- "poisson"
            defaultlink <- log
            reverselink <- exp
            distro <- poisson
            if (is.null(offset_str)){
                stop(distro_family, " distribution is selected, but no size factor information provided. Either use 'estimate_offset=TRUE' or provide 'offset_str'.")
            }
            irm_norm <- t(t(irm) / as.vector(irc[[offset_str]]))
            irm_norm_1p <- irm_norm + 1
        } else if ((distro_family == "nb") || 
                   (distro_family == "NB") || 
                   (distro_family == "negative.binomial")){
            if (!all(sapply(irm, `%%`, 1) == 0)){
                stop(distro_family, " distribution is selected. The 'feature_matrix' must be all integers.")
            }
            distro_family <- "NB"
            defaultlink <- log
            reverselink <- exp
            distro <- NULL
            if (is.null(offset_str)){
                stop(distro_family, " distribution is selected, but no size factor information provided. Either use 'estimate_offset=TRUE' or provide 'offset_str'.")
            }
            irm_norm <- t(t(irm)/as.vector(irc[[offset_str]]))
            irm_norm_1p <- irm_norm + 1
        } else if ((distro_family == "gaussian") || 
                   (distro_family == "Gaussian")){
            if (all(sapply(irm, `%%`, 1) == 0)){
                warning(distro_family, " distribution is selected, but the values in 'feature_matrix' are all integers. Consider using Poisson distribution instead?")
            }
            distro_family <- "gaussian"
            defaultlink <- identity
            reverselink <- identity
            distro <- gaussian
            irm_norm <- irm
            irm_norm_1p <- irm + 1
            if (!is.null(offset_str)){
                warning(distro_family, " distribution is selected. The provided 'offset_str' option will be ignored.")
                offset_str <- NULL
            }
        }
        # Apply SVA to raw counts to estimate SVs for each sample
        sv_model_str <- ""
        if (apply_sva){
            if (is.null(remove_terms)){
                null_model_str="1"
            } else {
                null_model_str=paste(c(remove_terms), collapse=" + ")
            }
            mod <- model.matrix(formula(fml), irc)
            message("svaseq FULL model: ", formula(fml))
            mod0 <- model.matrix(formula(paste0("~",null_model_str)), irc)
            message("svaseq null model: ", formula(paste0("~",null_model_str)))
            svseq <- svaseq(as.matrix(irm_norm), mod, mod0, n_sv)
            for (i in 1:dim(svseq$sv)[2]){
                irc[[paste0("SV",i)]] <- svseq$sv[,i] 
                sv_model_str <- paste0(sv_model_str," + SV",i) 
                sv_terms <- c(sv_terms, paste0("SV",i))
            }
        }
        # Set up model formula and params to be called for each gene
        sp <- 0
        if ((method == "glm") || (method == "GLM")){
            fml_fit <- paste0(fml, sv_model_str)
            fml_fixed <- fml_fit
            if ((distro_family == "poisson") || (distro_family == "gaussian")){
                regression_caller <- glm
            } else if (distro_family == "NB"){
                regression_caller <- glm.nb
                sp <- 2
            }
        # If using a Mixture model (i.e. including random terms)
        } else if ((method == "glmm") || 
                   (method == "GLMM") || 
                   (method == "lme4") || 
                   (method == "glmer") || 
                   (method == "lmer") || 
                   (method == "mixed")){
            if (is.null(random_terms) || 
               (nchar(random_terms) == 0)){
                stop("'glmm' is selected so that 'random_terms' cannot be empty or NULL.")
            }
            main_terms <- test_terms[which(!(test_terms %in% random_terms))]
            main_model_st <- paste0("~", paste(c(main_terms), collapse=" + "))
            random_model_str <- ""
            for (r_t in random_terms){
                random_model_str <- paste0(random_model_str, " + (1|", r_t, ")")
            }
            fml_fit <- paste0(main_model_str, sv_model_str, random_model_str)
            fml_fixed <- paste0(main_model_str, sv_model_str)
            if (distro_family == "poisson"){
                regression_caller <- glmer
            } else if (distro_family == "NB"){
                regression_caller <- glmer.nb
                sp=1
            } else if (distro_family == "gaussian"){
                regression_caller <- lmer
                sp <- 1
            }
        }
        message(sprintf('Fitting %s models with %s distribution: %s ...',
                        method, distro_family, fml_fit
                )
        )
        # Fit models for each gene with specified params in parallel
        message(c('cores:', parallel::detectCores()))
        plan(multisession, workers=parallel::detectCores())
        models0 <- future_apply(irm, 
                                1,
                                fit_single_model,
                                meta_data=irc,   # sample metadata
                                approach=sp,  # 2 for GLM NB
                                fit_formula=fml_fit,
                                offset_str=offset_str, # "buildDEmodel.sizeFactor"
                                regression_caller=regression_caller, # glm.nb
                                distro=distro,  # NB
                                defaultlink=defaultlink  # ln
                   )
        plan(sequential)
        names(models0) <- rownames(irm)
        # Filter results where models (genes) failed to be built/estimated
        model_failed <- lapply(models0, 
                               function(model) {
                                   if (is.null(model)){
                                       return(TRUE)
                                   } else {
                                       if (method == 'glmm'){
                                           if (model@optinfo$message == "parameter values converged to within tolerance"){
                                               return(FALSE)
                                           } else {
                                               return(TRUE)
                                           }
                                       } else {
                                           return(!model$converged)
                                       }
                                   }
                               }
                        ) %>% unlist()
        models <- models0[which(!model_failed)]
        # Determine which models produced warnings (e.g. singularity, lack of convergence etc.)
        model_warned <- sapply(models0, 
                               function(x){
                                 if (method == 'glmm'){
                                       if (is.null(x)){
                                           warned <- FALSE  # failed, not a warning
                                       } else {
                                           warned <- (length(x@optinfo$conv$lme4$messages) > 0)
                                       }
                                   } else {
                                       warned <- FALSE
                                   }
                                 return(warned)
                              }
        )
        # Create binary matrix to record model success/failure/warnings 
        model_status <- tibble(gene=names(models0),
                               warned=as.integer(model_warned),
                               failed=as.integer(model_failed)
                        ) %>% 
                        mutate(ok=as.integer(!warned & !failed)) %>%
                        as.data.frame() 
        message(sprintf('%s out of %s entries successfully modeled (%s with warnings).',
                        length(models), nrow(irm), length(which(model_warned))
                )
        )
        message(sprintf('Models failed to converge for the following genes: %s',
                        paste(names(models0[which(model_failed)]), 
                              sep='', 
                              collapse=','
                        )
                )
        )
        # Extract all parameter values for each model (1 model per gene)
        message("Extracting model parameters...")
        # Function to call if there is an error in extracting params from models
        err_fnc <- function(err){
            message("PERMUTATION FAILED TO GET PARAMS")
            message(paste(sampleTable$SampleID, sep='', collapse=','))
            message("Failed with the following error")
            message(err)
            return(NA)
        }
        # Estimted Beta value for each parameter for each modes
        beta_all <- tryCatch(t(sapply(models,
                                      function(x){ coef(summary(x))[,1] }
                               )
                             ), 
                             error=err_fnc
                    )
        rownames(beta_all)=names(models)
        # Std. Err of each model
        stderr_all <- tryCatch(t(sapply(models,
                                        function(x){ coef(summary(x))[,2] }
                                 )
                               ), 
                               error=err_fnc
                      )
        rownames(stderr_all)=names(models)
        # Model Residuals
        resid_all <- tryCatch(t(sapply(models, 
                                       resid, 
                                       type="response"
                                )
                              ), 
                              error=err_fnc
                     )
        rownames(resid_all)=names(models)
        colnames(resid_all)=colnames(irm)
        # Raw pvalues
        p_all <- tryCatch(t(sapply(models, 
                                   function(x){ coef(summary(x))[,4] }
                            )
                          ),
                          error=err_fnc
                 )
        rownames(p_all)=names(models)
        # Adjusted pvalues
        q_all <- tryCatch(apply(as.matrix(p_all), 
                                2, 
                                p.adjust, 
                                method="BH"
                          ),
                          error=err_fnc
                 )
        if (nrow(p_all) == 1){ 
            q_all <- rbind(c(), q_all) 
        }
        rownames(q_all)=names(models)
        # Variance Co-variance matrix
        vcov_all <- tryCatch(merge_vcov_from_model_list(models), 
                             error=err_fnc
                    )
        y <- as.matrix(defaultlink(irm_norm[names(models), , drop=FALSE]))
        message('Collected all params')
        # Provide "corrected" counts if any co-variates/SVs are in the formula
        if ((length(c(remove_terms, sv_terms)) > 0) || (method != "glm")){
            message(sprintf("Correcting against unwanted effects: %s...",
                            paste(c(remove_terms, 
                                    sv_terms
                                  ),
                                  collapse=", "
                            )
                    )
            )
            m <- as.data.frame(model.matrix(formula(fml_fixed), data=irc))
            for (keyword in remove_terms){
                remove_cols <- c(remove_cols, 
                                 which(grepl(paste0("^", keyword), colnames(m)))
                               )
            }
            for (svword in sv_terms){
                remove_cols <- c(remove_cols, which(colnames(m) == svword))
            }
            remove_cols <- unique(remove_cols)
            #X=m[, remove_cols]
            #X=as.matrix(X)
            #beta_sub=beta_all[, remove_cols]
            #beta_sub=as.matrix(beta_sub)
            #cleany=y - t(as.matrix(X) %*% t(beta_sub))
            m[, remove_cols] <- 0
            cleanY <- reverselink(as.matrix(beta_all) %*% t(as.matrix(m))) + resid_all
            cleanY[which(cleanY < 0)] <- 0
            cleany <- defaultlink(cleanY)
        } else {
            cleany <- y
        }
        ## Package output into a single object to return
        message("Collecting results...")
        rownames(cleany) <- names(models)
        glm_res <- new("GLMres",
				       y=irm,
				       colData=irc,
				       design=formula(fml_fit),
				       beta=as.matrix(beta_all),
				       stderr=as.matrix(stderr_all),
				       pval=as.matrix(p_all),
				       fdr=as.matrix(q_all),
				       vcov=vcov_all,
				       resid=as.matrix(resid_all),
				       deSV=as.matrix(cleany),
				       status=model_status,
				       models=models0,
				       args=list(method=method,
                                 estimate_offset=estimate_offset,
                                 offset_str=offset_str,
                                 apply_sva=apply_sva,
                                 n_sv=n_sv,
                                 test_terms=test_terms,
                                 remove_terms=remove_terms,
                                 random_terms=random_terms,
                                 distro_family=distro_family
                       )
                   )
        message("Done.")
        return(glm_res)
    }
}

summarize_glm <- function(
    model_data,
    results,
    title_str,
    outfile){
    # Summarize GLM model outputs for convinience of inspecting model data after being computed
    summary_results <- list()
    summary_results$tested <- sum(colSums(model_data@status))
    summary_results <- c(summary_results, 
                         # status of each gene's test
                         colSums(model_data@status) %>% as.list(),
                         # DEG results
                         results %>%
                         as_tibble() %>%
                         setNames(c('fdr', 'lfc', 'pval', 'EnsemblID')) %>%
                         dplyr::summarize(nominal_up_genes=sum(pval < 0.05 & lfc > 0),
                                          nominal_down_genes=sum(pval < 0.05 & lfc < 0),
                                          fdr_up_genes=sum(fdr < 0.1 & lfc > 0),
                                          fdr_down_genes=sum(fdr < 0.1 & lfc < 0)
                         ) %>% 
                         as.list()
    )
    names(summary_results) <- paste0('genes_', names(summary_results))
    summary_results$n_sv <- length(grep('^SV[0-9]+$', colnames(model_data@colData)))
    summary_results$design <- Reduce(paste, deparse(model_data@design))
    write.table(as.data.frame(summary_results), 
                outfile,
                quote=FALSE,
                sep='\t',
                row.names=FALSE
    )
    title_str <- sprintf('%s (samples=%s, genes=%s)', 
                     title_str, 
                     ncol(model_data@y), 
                     nrow(model_data@y)
    )
    rnaPCA(model_data@y, 
           model_data@colData, 
           intgroup='Condition',
           title_str=sprintf('%s (preSVA)', title_str),
           outfile=sprintf('%s_PCA_preSVA.pdf', outfile)
    )
    colnames(model_data@deSV) <- colnames(model_data@y)
    rnaPCA(model_data@deSV, 
           model_data@colData, 
           intgroup='Condition',
           title_str=sprintf('%s (postSVA)', title_str),
           outfile=sprintf('%s_PCA_postSVA.pdf', outfile)
    )
    colnames(results) <- c('padj', 'lfc', 'pvalue', 'EnsemblID')
    volcanoPlotForGLM(results,
                      title_str=title_str,
                      outfile=sprintf('%s_volcanoPlot.pdf', outfile)
    )
}

##############
# GLM DEG Analysis Wrappers/Utilities 
level_metadata <- function(
    comparison_vector,
    meta,
    control_label=NA,
    condition_label=NA){
    # Extract the correct samples to use for the specified comparison from the full metadata
    # and create a column 'Condition' that is a factor with the correct levels for
    # any differential expression analysis
    # comparison should be a vector of c(treatment1, genotype1, treatment2, genotype2) e.g.
    # c('NO', 'CON', 'NO', 'XDP') or c('NO', 'CON', 'NO', 'CombinedXDP') 
    if (grepl('CombinedXDP', comparison_vector[2])){
        meta1 <- meta %>% 
                 filter(Treatment == comparison_vector[1] & Genotype == 'XDP') %>%
                 mutate(geno=sprintf('Combined%s', Genotype)) 
    } else if (grepl('CombinedCON', comparison_vector[2])){
        meta1 <- meta %>% 
                 filter(Treatment == comparison_vector[1] & 
                        (geno == 'XDP_dSVA' | geno == 'CON')
                 ) %>%
                 mutate(geno='CombinedCON')
    } else {
        meta1 <- meta %>% filter(Treatment == comparison_vector[1] & geno == comparison_vector[2]) 
    }
    if (grepl('CombinedXDP', comparison_vector[4])){
        meta2 <- meta %>% 
                 filter(Treatment == comparison_vector[3] & Genotype == 'XDP') %>%
                 mutate(geno=sprintf('Combined%s', Genotype)) 
    } else if (grepl('CombinedCON', comparison_vector[4])){
        meta2 <- meta %>% 
                 filter(Treatment == comparison_vector[3] &
                        (geno == 'XDP_dSVA' | geno == 'CON')
                 ) %>%
                 mutate(geno='CombinedCON')
    } else {
        meta2 <- meta %>% filter(Treatment == comparison_vector[3] & geno == comparison_vector[4]) 
    }
    # combine sample data for the 2 groups
    comparison_meta <- bind_rows(meta1, meta2) %>%
                       unite(Condition,
                             c(Treatment, geno),
                             sep='.',
                             remove=FALSE
                       ) 
    # set factor levels of condition col appropriately
    if (is.na(control_label) | is.na(condition_label)){
        if ('NO.CON' %in% comparison_meta$Condition){
            control_label <- 'NO.CON'
        } else if (any(grepl('CON', unique(comparison_meta$Condition)))){
            control_label <- grep('CON', unique(comparison_meta$Condition), value=TRUE)[[1]]
        } else if ('NO.XDP' %in% comparison_meta$Condition){
            control_label <- 'NO.XDP'
        } else {
            control_label <- unique(comparison_meta$Condition)[[1]]
        }
        condition_label <- setdiff(unique(meta$Condition), control_label)
    }
    comparison_meta <- comparison_meta %>%
                       mutate(Condition=factor(Condition, 
                                               levels=c(control_label,
                                                        condition_label
                                               )
                                        )
                       )
    return(comparison_meta)
}

combine_results <- function(
    model_data,
    test_terms){
    # Function to massage GLM model data into a tsv friendly shortcut
    model_results <- lapply(test_terms,
                            function(term){
                                get_fdr(model_data, term=term) %>% as_tibble() 
                            }
                     ) %>% 
                     bind_rows() %>%
                     pivot_longer(!name,
                                  names_to='metric',
                                  values_to='values'
                     ) %>%
                     filter(!(is.na(values))) %>%
                     separate(metric, 
                              c('metric', 'Variable'),
                              sep='\\.'
                     ) %>%
                     pivot_wider(names_from=metric,
                                 values_from=values
                     ) %>%
                     rename(padj = FDR,
                            lfc = LFC,
                            pvalue = pval,
                            EnsemblID = name
                     )
    return(model_results)
}

evaluation_wrapper <- function(
    method,
    counts,
    meta,
    analysis_params, 
    title_str=NA,
    output_file=NA){
    # wrapper for various differential expression analyses methods
    if (method == 'deseq'){
        # Just use DESeq for standard DEG analysis and save the resutls
        results <- deseq_wrapper(counts,
                                 meta,
                                 analysis_params
        )
        if (!(is.na(output_file))){
            results <- convert_deseq_to_tibble(results)
            write_tsv(results, output_file)
        } else {
            return(results)
        }
    } else if (method == 'glm'){
        model_data <- buildDEmodel(counts,
                                   meta,
                                   estimate_offset=TRUE,
                                   apply_sva=TRUE,
                                   method='glm',
                                   distro_family="nb",
                                   test_terms=c(analysis_params$design_var),
                                   remove_terms=NULL,
                                   random_terms=NULL,
        )
        model_results <- get_fdr(model_data, term=analysis_params$design_var)
        if (!(is.na(output_file))){
            write_tsv(as_tibble(model_results), output_file)
        } else {
            return(list(results=model_results, model_data=model_data))
        }
    } else if (method == 'glm_interact2'){
        # Use a GLM with an interaction term between Treatment and Genotype 
        # So the formula is ~ geno + Interaction
        meta <- meta %>%
                mutate(Interaction=case_when(Treatment == 'NO' ~ 0,
                                             Treatment != 'NO' & geno == 'CON' ~ 1,
                                             Treatment != 'NO' & geno != 'CON' ~ 2
                       )
                ) %>% 
                mutate(Interaction=as.factor(Interaction)) %>%
                droplevels()
        test_terms <- c('geno', 'Interaction')
        model_data <- buildDEmodel(counts,
                                   meta,
                                   estimate_offset=TRUE,
                                   apply_sva=TRUE,
                                   method='glm',
                                   distro_family="nb",
                                   test_terms=test_terms,
                                   remove_terms=NULL,
                                   random_terms=NULL
                      )
        model_results <- combine_results(model_data, test_terms)
        return(model_results)
    } else if (method == 'glm_interact'){
        # Use a GLM with an interaction term between Treatment and Genotype 
        # so the formula is ~ geno + Treatment + Interaction
        meta <- meta %>%
                mutate(Interaction=ifelse(Treatment != 'NO' & geno != 'CON',
                                          1,
                                          0
                       )
                ) %>% 
                mutate(Interaction=as.factor(Interaction)) %>%
                droplevels()
        test_terms <- c('geno', 
                        'Treatment',
                        'Interaction'
                      )
        model_data <- buildDEmodel(counts,
                                   meta,
                                   estimate_offset=TRUE,
                                   apply_sva=TRUE,
                                   method='glm',
                                   distro_family="nb",
                                   test_terms=test_terms,
                                   remove_terms=NULL,
                                   random_terms=NULL
                      )
        model_results <- combine_results(model_data, test_terms)
        return(model_results)
    } else if (method == 'glmm'){
        # Build mixture models, including 'MIN' (i.e. sample origin lines) as a 
        # random term to account for samples grown from the same cells
        model_data <- buildDEmodel(counts,
                                   meta,
                                   estimate_offset=TRUE,
                                   apply_sva=TRUE,
                                   method='glmm',
                                   distro_family="nb",
                                   test_terms=c(analysis_params$design_var),
                                   remove_terms=NULL,
                                   random_terms=c('MIN'),
        )
        model_results <- get_fdr(model_data, term=analysis_params$design_var)
        if (!(is.na(output_file))){
            write_tsv(as_tibble(model_results), output_file)
        } else {
            return(list(results=model_results, model_data=model_data))
        }
    } else {
        quit(save="no")
    }
}

