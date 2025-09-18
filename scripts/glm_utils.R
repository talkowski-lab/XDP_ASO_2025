################################################
# Dependencies
################################################
library(DESeq2)
library(sva)
library(MASS)
library(parallel)
library(future)
library(future.apply)
library(furrr)
library(tidyverse)

################################################
# Filter low expr genes by median CPM
################################################
filter_by_counts <- function(
    mat,
    conditionlist,
    cut=0.1,
    pct=0.5){
    hit = false  # if unchanged, gene does not pass cpm filter
    for (i in conditionlist){
        # number of sample in condition {i}
        samplenum = length(i)  
        # initiation: set number of sample that expresses the current gene to 0. e.g. no sample expressed this gene
        expressed = 0          
        for (x in i){
        # if sample {x} expressed this gene according to cutoff, increment counter
            if (mat[x] >= cut){expressed = expressed + 1} 
        }
        # if enough samples in condition {i} expressed this gene, keep it for downstream analysis
        if ((expressed/samplenum) >= pct){ hit = true; break }
    }
    return(hit)
}

cpm_filter_counts <- function(
    counts,
    meta,
    condition_var,
    cpm_cutoff=0.1, 
    samples_threshold=0.5){
    # keep any gene that has cpm >= cpm_cutoff in at least 
    # samples_threshold fraction of samples from any condition
    cpms <- 1000000 * t(t(counts) / meta$uniquely_mapped_reads)
    meta[[condition_var]] = as.character(meta[[condition_var]])
    message(paste(unique(meta[[condition_var]])))
    # for each condition (e.g. for genotype is con, xdp) get all samples that
    # meet that condition
    conditionlist <- list()
    for (condition in unique(meta[[condition_var]])){
        conditionlist[[condition]] <- which(meta[, condition_var] == condition)
    }
    # remove genes that do not pass the cpm filter using samples in any condition
    expressedgenes <- 
        apply(
            cpms, 
            1, 
            filter_by_counts, 
            conditionlist=conditionlist, 
            cut=cpm_cutoff, 
            pct=samples_threshold
    )
    message(sprintf("%s genes retained after filtering", sum(expressedgenes)))
    # returns list of genes to keep after cpm filtering
    return(expressedgenes)
}

################################################
# Object class definitions
################################################
setClass(
    "glmres",
    representation(
        y="matrix",           # raw input count matrix
        coldata="data.frame", # sample metadata
        design="formula",     # formula for the glm
        beta="matrix",        # estimated beta coefficients
        stderr="matrix",      # estimated std err of coefficients
        pval="matrix",        # pvalue for each gene
        fdr="matrix",         # fdr for each gene
        vcov="data.frame",    # variance-covariance matrix 
        resid="matrix",       # glm residuals
        desv="matrix",        # matrix of svs estimated
        status="matrix",      # binary matrix of gene convergence
        # models="list"         # lme4 model objects for each gene
        args="list"           # input arguments
    )
)

setMethod(
    f="show",
    signature="glmres",
    definition=
        function(object){
            print("class: GLMres")
            print("slots: y, colData, design, beta, stderr, pval, fdr, vcov, resid, deSV, status, args") #, models")
            print(
                sprintf(
                    "deSV dim: %s %s", 
                    nrow(object@deSV), 
                    ncol(object@deSV)
                )
            )
        }
)

################################################
# Functions to build GLM models
################################################
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
    res <- 
        tryCatch(
            {
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
    tt <- lapply(model_list, vcov)
    tt1 <- 
        lapply(
            names(tt), 
            FUN=
                function(x){
                    tmp=as.data.frame(as.matrix(tt[[x]]))
                    rownames(tmp)=NULL 
                    tmp[["name"]]=x 
                    return(tmp)
                }
        ) 
    tt2=do.call(rbind, tt1)
    return(tt2)
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
    setClass(
        "GLMres",
        representation(
            y="matrix",
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
            # models="list",
            args="list"
        )
    )
    setMethod(
        f="show",
        signature="GLMres",
        definition=
            function(object){
                print("class: GLMres")
                print("slots: y, colData, design, beta, stderr, pval, fdr, vcov, resid, deSV, status, args") #, models")
                print(sprintf("deSV dim: %s %s", nrow(object@deSV), ncol(object@deSV))) 
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
				       # models=models0,  # mostly large, redudant info 
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

################################################
# Compute DEG results for all contrasts from GLM 
################################################
get_fdr <- function(
    res,
    term){
    # Reorganize DEG data from list of model objects into a results table
    # Get FC resutls
    fc <- as.data.frame(res@beta)
    new0 <- fc[,grep(term, colnames(fc)), drop=F]
    colnames(new0) <- paste0("LFC.",colnames(fc)[grep(term, colnames(fc))])
    # Get p-values
    pv <- as.data.frame(res@pval)
    new1 <- pv[,grep(term, colnames(pv)), drop=F]
    colnames(new1) <- paste0("pval.",colnames(pv)[grep(term, colnames(pv))])
    # get BH corrected p-values
    bh <- as.data.frame(res@fdr)
    new2 <- bh[,grep(term, colnames(bh)), drop=F]
    colnames(new2) <- paste0("FDR.",colnames(bh)[grep(term, colnames(bh))])
    # Bind everything together into a single results table
    new <- as.data.frame(new2)
    new <- cbind(new, new0, new1)
    new$name <- as.vector(rownames(new))
    rownames(new) <- NULL
    return(new)
}

get_intercept_contrasts <- function(
    glm_models,
    design_var){
    intercept_level <- 
        glm_models@colData[[design_var]] %>% 
        levels() %>%
        paste0(design_var, .) %>%
        {.[1]}
    deg_results <- 
        get_fdr(
            glm_models, 
            term=design_var
        ) %>%
        as_tibble() %>%
        # pivot so data for each level of design_var in a separate row
        pivot_longer(
            !name,
            names_to='tmp',
            values_to='value'
        ) %>%
        separate_wider_delim(
            tmp,
            delim=fixed('.'),
            names=c('metric', 'level'),
            too_many='merge',
        ) %>%
        pivot_wider(
            names_from=metric,
            values_from=value
        ) %>%
        setNames(
            c(
                'EnsemblID',
                'Contrast',
                'fdr',
                'lfc',
                'pvalue'
            )
        ) %>%
        rowwise() %>%
        mutate(Contrast=sprintf('%s vs %s', Contrast, intercept_level)) 
    return(deg_results)
}

get_specific_contrast <- function(
    num_level,
    denom_level,
    contrast_vec, 
    glm_models){
    # Given the GLM results and a vector specifying a contrast compute the
    # LFC and p-value for all genes and return as a table
    # Calculate new betas (LFCs) for each gene
    ctr_betas <- glm_models@beta %*% as.vector(t(contrast_vec))
    # Compute variance for each gene for this contrast and flatten to vector
    nvars <- ncol(glm_models@beta)  # number of variables in the glm models
    gene_names <- rownames(glm_models@beta)  # all gene names in the analysis
    # Vector of variance of beta for each gene
    ctr_vars <- 
        sapply(
            gene_names,
            function(gene_name){
                # Get Var-CoVar matrix for this gene
                i <- which(gene_names == gene_name)
                gene_vcv <- 
                    glm_models@vcov %>%
                    {.[((i - 1) * nvars + 1) : (i * nvars), 1:nvars]} %>% 
                    as.matrix()
                # Get variance of Beta for this gene in this contrast
                gene_var <- 
                    t(as.vector(contrast_vec)) %*% 
                    tcrossprod(gene_vcv, t(as.vector(contrast_vec)))
                return(gene_var)
            }
        ) %>% 
        as.matrix()
    # Compute p-value from z-value (i.e. beta / std.err)
    ctr_pvalues <- 2 * pnorm(q=abs(ctr_betas / sqrt(ctr_vars)), lower.tail=FALSE)
    results_df <- 
        merge(
            ctr_betas,
            ctr_pvalues,
            by='row.names'
        ) %>%
        setNames(
            c(
                'EnsemblID',
                'lfc',
                'pvalue'
            )
        ) %>%
        add_column(Contrast=glue('{num_level} vs {denom_level}'))
    return(results_df)
}

get_all_contrasts <- function(
    glm_models,
    design_var){
    # Given the data for all the contrasts including the intercept, we want to 
    # compute all contrasts between non-intercept levels of the desgin_var
    # e.g. variable TYPE has possible values A B C D with A being the intercept,
    # Then here we would calculate the contrasts B vs C, B vs D and C vs D
    # and return all those results in a single table
    # using the data we already have of A vs B, A vs C etc...

    intercept_results <- get_intercept_contrasts(glm_models, design_var)
    # Get all levels of design_var that are not the intercept
    # i.e. the first factor level of design_var, which must be a factor
    all_variables <- colnames(glm_models@beta)
    all_levels <- 
        glm_models@colData[[design_var]] %>% 
        levels() %>%
        paste0(design_var, .)
    intercept_level <- all_levels[1]
    all_other_levels <- all_levels[2:length(all_levels)]
    # No other combinations to test (only 2 levels including intercept
    if (length(all_other_levels) < 2){ return(intercept_results) }
    # Get all untested combinations of the non-intercept levels to compare
    # and produce a contrast vector for each contrast
    all_other_contrasts <- 
        combn(all_other_levels, 2)  %>%
        t() %>%
        as_tibble() %>%
        setNames(c('denom_level', 'num_level')) %>%
        rowwise() %>%
        mutate(
            contrast_vec=
                list(
                    grepl(  num_level, all_variables, fixed=TRUE) -
                    grepl(denom_level, all_variables, fixed=TRUE)
                )
        )
    # Get LFCs and pvlaues for all the specified contrasts above
    all_other_results <- 
        all_other_contrasts %>%
        # for each contrast compute DEG results for all genes and bind together
        # rowwise() %>%
        # future_pmap(get_specific_contrast,
        purrr::pmap(
            get_specific_contrast,
            glm_models=glm_models,
            .progress=TRUE
        ) %>%
        bind_rows() %>%
        # adjust pvalue of all genes together per contrast
        group_by(Contrast) %>%  
        mutate(fdr=p.adjust(pvalue, method='BH')) %>%
        ungroup() 
    # Bind intercept contrasts and non-intercept contrasts into 1 table
    all_results <- bind_rows(intercept_results, all_other_results)
    return(all_results)
}

################################################
# Wrappers for convienience
################################################
run_glm <- function(
    metadata,
    raw_counts,
    design_var){
    # NOTE: sample_metadata should only contain samples of interest 
    # Subset full raw counts to only samples of interest
    counts <- raw_counts[, metadata$SampleID]
    # Filter low expression genes per group, defined by design_var
    genes_to_keep <- 
        CPM_filter_counts(
            counts,
            metadata, 
            condition_var=design_var,
            cpm_cutoff=0.1, 
            samples_threshold=0.5
        )
    counts <- counts[genes_to_keep, ]
    # Build GLM model for each gene and save the results
    glm_models <- 
        buildDEmodel(
            counts,
            metadata,
            test_terms=design_var,
            estimate_offset=TRUE,
            apply_sva=TRUE,
            method='glm',
            distro_family="nb",
            remove_terms=NULL,
            random_terms=NULL
        )
    return(glm_models)
}

glm_wrapper <- function(
    metadata,
    counts,
    design_var,
    name_only=FALSE,
    output_name=NA,
    output_dir=GLM_RESULTS_DIR){
    # NOTE: sample_metadata should ONLY contain samples being compared
    # Run the GLM model on the given set of samples and save 
    # 1. The full GLM data for each gene (gzip compressed)
    # 2. DEG results table 
    # 3. SV Adjusted SV counts
    # Name results files by which samples are being analyzed
    if (is.na(output_name)) {
        output_name <- 
            metadata %>%
            {.[[design_var]]} %>% 
            levels() %>%
            paste(collapse='+')
    } 
    print(output_dir)
    print(output_name)
    print(table(metadata[[design_var]]))
    if (name_only){ return(output_name) }
    # Produce a GLM model for each gene, predicting design_var with expression
    glm_models <- 
        run_glm(
            metadata,
            counts,
            design_var
        )
    saveRDS(
        glm_models, 
        compress='gzip',
        sprintf('%s/%s-GLMData.rds', output_dir, output_name)
    )
    # Summarize DEG results for all possible contrasts in design_var
    # deg_results <- get_intercept_contrasts(glm_models, design_var) %>%
    glm_models %>%
        get_all_contrasts(design_var) %>%
        # add gene names to EnsemblIDs
        left_join(
            gene_metadata %>%
            filter(type == 'gene') %>%
            dplyr::select(gene_id, gene_name),
            by=join_by(EnsemblID == gene_id),
        ) %>%
        relocate(EnsemblID, gene_name, Contrast) %>%
        write_tsv(sprintf('%s/%s-DEGResults.tsv', output_dir, output_name))
    # Save SV adjusted counts as tsvs for easy loading
    glm_models@deSV %>%
        as.data.frame() %>%
        setNames(glm_models@colData$SampleID) %>% 
        rownames_to_column(var='EnsemblID') %>%
        as_tibble() %>%
        write_tsv(sprintf('%s/%s-SVAdjustedCounts.tsv', output_dir, output_name))
}

################################################
# Load Results
################################################
load_deg_results <- function(){
    list.files(
        GLM_RESULTS_DIR,
        pattern='*-DEGResults.tsv',
        full.names=TRUE
    ) %>% 
    tibble(file=.) %>%
    mutate(comparison=basename(file)) %>%
    mutate(comparison=basename(str_remove(file, '-DEGResults.tsv'))) %>% 
    mutate(
        results=
            pmap(
                .l=list(file),
                .f=read_tsv,
                show_col_types=FALSE
            )
    ) %>% 
    dplyr::select(-c(file)) %>% 
    unnest(results) %>% 
    dplyr::rename('contrast'=Contrast) %>% 
    mutate(contrast=str_remove_all(contrast, 'Condition')) %>% 
    separate_wider_delim(
        contrast,
        delim=' vs ',
        names=c('numerator', 'denominator')
    )
}

load_corrected_counts_list <- function(){
    c(
        # Comparisons for all Controls vs Untreated XDP
        list(
             'cXDP'=
                file.path(GLM_RESULTS_DIR, 'NO.CON+NO.XDP-SVAdjustedCounts.tsv') %>% 
                read_tsv(show_col_types=FALSE),
            
             'XDP_naive'=
                file.path(GLM_RESULTS_DIR, 'NO.CON+NO.XDP_naive-SVAdjustedCounts.tsv') %>% 
                read_tsv(show_col_types=FALSE),
                
             'XDP_non_edit'=
                file.path(GLM_RESULTS_DIR, 'NO.CON+NO.XDP_non_edit-SVAdjustedCounts.tsv') %>% 
                read_tsv(show_col_types=FALSE)
        ),
        # CRISPR edited XDP (dSVA) vs Untreated XDP and vs Untreated Control
        list(
            'dSVA'=
                file.path(GLM_RESULTS_DIR, 'NO.CON+NO.XDP+NO.dSVA-SVAdjustedCounts.tsv') %>% 
                read_tsv(show_col_types=FALSE)
        ),
        # Comparisons for vs Untreated XDP vs ASO Treated XDP
        sapply(
            ASO_TREATMENTS,
            function(treatment) {
                file.path(
                    GLM_RESULTS_DIR,
                    glue('NO.CON+NO.XDP+{treatment}.XDP-SVAdjustedCounts.tsv')
                ) %>% 
                read_tsv(show_col_types=FALSE) 
            },
           simplify=FALSE,
           USE.NAMES=TRUE
        ) %>% 
        setNames(paste('ASO.', names(.), sep=''))
    )
}

load_corrected_counts_table <- function(
    outfile=CORRECTED_COUNTS_FILE,
    gene_metadata=NULL){
    if (file.exists(outfile)) {
        return(read_tsv(outfile, show_col_types=FALSE))
    } 
    mkdir(dirname(outfile))
    load_corrected_counts_list() %>% 
    # {.[grep('ASO|dSVA|XDP', names(.), value=TRUE)]} %>%
    {
        sapply(
            names(.),
            function(name){
                .[[name]] %>% 
                pivot_longer(
                    !EnsemblID,
                    names_to='SampleID',
                    values_to='counts'
                )
            },
            simplify=FALSE
        )
    } %>%
    bind_rows(.id='origin') %>%
    left_join(
        sample_metadata %>%
            mutate(isTreated=ifelse(Treatment == 'NO', FALSE, TRUE)),
        by='SampleID',
        multiple='all',
    ) %>% 
    mutate(
        Genotype=ifelse(Genotype == 'dSVA', 'XDP', Genotype),
        Treatment=ifelse(Treatment == 'XDP', 'Untreated', Treatment),
        Category=paste(Treatment, Genotype, sep='.')
    ) %>%
    left_join(
        gene_metadata %>%
            filter(type == 'gene') %>%
            dplyr::select(gene_id, gene_name),
        by=join_by(EnsemblID == gene_id)
    ) %>% 
    dplyr::select(!isTreated) %T>%
    write_tsv(outfile)
}

