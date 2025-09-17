################################################
# Dependencies
################################################
library(furrr)
library(tidyverse)

################################################
# Data Read/writing utils
################################################
check_cached_results <- function(
    results_file,
    force_redo=FALSE,
    return_data=TRUE,
    show_col_types=FALSE,
    results_fnc,
    ...){
    # Now check if the results file exists and load it
    tic()
    if (is.null(results_file)) {
        message("No results file, will just return data")
        results <- results_fnc(...)
        return_data <- TRUE
    } else {
        # Set read/write functions based on filetype
        output_filetype <- results_file %>% str_extract('\\.[^\\.]*$') 
        if (output_filetype == '.rds') {
            load_fnc <- readRDS
            save_fnc <- saveRDS
        } else if (output_filetype %in% c('.txt', '.tsv')) {
            load_fnc <- partial(read_tsv, show_col_types=show_col_types)
            save_fnc <- write_tsv
        } else {
            stop(glue('Invalid file extesion: {extension}'))
        }
        if (file.exists(results_file) & !force_redo) {
            message(glue('{results_file} exists, not recomputing results'))
            if (return_data) {
                message('Loading cached results...')
                results <- load_fnc(results_file)
            }
        # or force recompute+cache the data and return it
        } else {
            if (file.exists(results_file) & force_redo) {
                message(glue('{results_file} exists, recomputing results anyways'))
            } else {
                message(glue('No cached results, generating: {results_file}'))
            }
            dir.create(dirname(results_file), recursive=TRUE, showWarnings=FALSE)
            results <- results_fnc(...) %T>% save_fnc(results_file)
            # If results dont exist or force_redo is TRUE compute + cache results
            # Assumes save_fnc is of fomr save_fnc(result_object, filename)
        }
    }
    # Or dont save the results, just return the results
    toc()
    if (return_data) {
        return(results)
    } else { 
        return(invisible(NULL))
    }
}

save_output <- function(df, output_filename){
    output_file <- file.path(OUTPUT_DIR, output_filename)
    write_tsv(df, output_file)
    glue('Result is saved to {output_file}') %>% message()
}

mkdir <- 
    partial(
        dir.create,
        mode="775",
        recursive=TRUE,
        showWarnings=FALSE
    )
################################################
# Load Data
################################################
load_counts <- function(){
    # load raw counts for rnaseq data
    READ_COUNTS_FILE %>% 
    read.table(check.names=FALSE) 
}

load_meta <- function(filter_outliers=TRUE){
    # load sample metadata
    SAMPLE_METADATA_FILE %>% 
    read_tsv(show_col_types=FALSE) %>% 
    # exclude pre-selected samples for technical reasons
    filter(Exclude == 'NO') %>%
    # create the condition column
    unite(
        Condition, 
        c(Treatment, geno),
        sep='.',
        remove=FALSE
    ) %>%
    # create original patient lvl id for each sample
    mutate(PatientID=substr(ParentLine, 1, 5)) %>%
    # reorder levels for formulas later 
    mutate(
        geno=
            factor(
                geno,
                levels=
                    c(
                        'CON',
                        'XDP',
                        'XDP_dSVA',
                        'XDP_non_edit'
                    )
            ),
        Treatment=
            factor(
                Treatment,
                levels=
                    c(
                        'NO',
                        '307',
                        '876',
                        '877',
                        '879',
                        '880',
                        '131',
                        '881'
                    )
           )
    ) %>%
    # consistent outlier samples identified via pca clustering
    {
        if (filter_outliers) {
            filter(
                .,
                !SampleID %in% 
                    c(
                        'NO_ASO_XDP_33363C_n_15',
                        'NO_ASO_XDP_non_edit_35833A_ISO_SVA_1B5'
                    )
            )
        } else {
            .
        }
    }
}

load_gtf <- function(){
    # get annotation of genes and genomic positions
    # '/data/talkowski/samples/microdeldup/15q/celllines/data/ref/homo_sapiens.grch37.75.ercc.gtf' %>%
    #     rtracklayer::import() %>%
    #     as_tibble() 
    GENE_METADATA_FILE %>% 
    read_tsv(show_col_types=FALSE)
}

make_bootstrap_results <- function(
    outfile=BOOTSTRAP_RESULTS_FILE,
    gene_metadata=NA){
    if (file.exists(outfile)){
        bootstrap_results <- outfile %>% read_tsv(show_col_types=FALSE) 
    } else {
        # load bootstrap summary data, produced by evaluate_and_summarize.r 
        # each input file represents a summary of 200 bootstraps, so
        # these are all read in and the summary data is combined
        # gene_metadata <- ifelse(is.na(gene_metadata), load_gtf(), gene_metadata)
        bootstrap_results <- 
            BOOTSTRAP_RESULTS_DIR %>% 
            list.files(
                pattern='summary-glm-original-NO.CON.NO.CombinedXDP.tsv.[0-9.]*',
                full.names=true
            ) %>%
            {
                lapply(
                    seq_along(.),
                    function(idx) {
                        .[[idx]] %>%
                        read_tsv(progress=false) %>% 
                        add_column(index=idx)
                    }
                )
            } %>% 
            bind_rows() %>%
            dplyr::select(
                starts_with('frac'), 
                'ensemblid',
                'index',
                'bootstraps'
            ) %>%
            pivot_longer(
                starts_with('frac_'),
                names_to='expression_pattern',
                values_to='fraction_of_bootstraps',
                names_prefix='frac_'
            ) %>%
            group_by(ensemblid, expression_pattern) %>%
            mutate(raw_value=fraction_of_bootstraps * bootstraps) %>%
            summarize(
                bootstraps=sum(bootstraps),
                raw_value=sum(raw_value)
            ) %>%
            mutate(fraction_of_bootstraps=raw_value / bootstraps) %>%
            ungroup() %>%
            write_tsv(bootstrap_results, outfile)
    }
    bootstrap_results %>%
    mutate(
        expression_pattern=
            factor(
                expression_pattern,
                levels=
                    c(
                        'up',
                        'nominal', 
                        'nominal_and_up',
                        'nominal_and_down',
                        'down',
                        'fdr',
                        'fdr_and_up',
                        'fdr_and_down'
                    )
            )
    )
}

load_threshold_bootstrap_results <- function(
    bootstrap_results,
    cpm_thresh=0.75,
    deg_threshes=c(0.5, 0.6, 0.7, 0.8)){
    # count genes per threshold and significance level
    # bootstrap_results <- 
    #     ifelse(is.null(bootstrap_results),
    #            load_bootstrap_results(),
    #            bootstrap_results
    #     )
    bootstrap_results %>%
    # filter for genes passing cpm > 75% of the time
    filter(bootstraps / max(bootstraps) >= cpm_thresh) %>%
    # merge directionally distinct groups together
    mutate(expression_pattern=as.character(expression_pattern)) %>%
    filter(grepl('_and_', expression_pattern)) %>%
    mutate(sig.lvl=gsub('_and_.*$', '', expression_pattern)) %>%
    # check if each gene passes each deg threshold
    add_column(threshes=rep(list(deg_threshes), nrow(.))) %>%
    rowwise() %>%
    mutate(deg_threshes=list(threshes[(fraction_of_bootstraps > threshes)])) %>%
    unnest(deg_threshes) %>%
    # make binary columns (gene passes or not) per threshold 
    add_column(tmp=1) %>%
    pivot_wider(
        names_from=deg_threshes,
        names_prefix='deg_thresh-',
        values_from=tmp,
        values_fill=0
    ) %>%
    # count(sig.lvl, expression_pattern)
    # dplyr::arrange(desc(ensemblid), sig.lvl) %>%
    # select(ensemblid, sig.lvl, starts_with('deg_thresh-')) %>% head(10)
    # pivot all binary column into 2 cols for group_by 
    pivot_longer(
        starts_with('deg_thresh-'),
        names_to='deg_thresh',
        names_prefix='deg_thresh-',
        values_to='bin'
    ) %>%
    # filter genes that dont pass each specified threshold
    filter(bin != 0) %>%
    # count genes per group and get list of ensemblid
    group_by(sig.lvl, deg_thresh) %>%
    summarize(
        n_genes=sum(bin),
        genes=list(ensemblid)
    ) %>%
    ungroup() %>%
    add_column(cpm_thresh=cpm_thresh)
}

load_taf1_consistent_genes <- function(
    outfile=taf1_consistent_genes_file,
    bootstrap_min=0.75){
    # get genes whose expression patterns are as consistent as taf1 
    if (file.exists(outfile)){
        consistent_genes <- 
            outfile %>%
            read_tsv() 
    } else {
        # expression patterns to evalue consistency of
        significance_lvls <- 
            c('nominal_and_up',
             'nominal_and_down',
             'fdr_and_up',
             'fdr_and_down'
            )
        # bootstrap results
        bootstrap_results <- load_bootstrap_results()
        bootstrap_total <- max(bootstrap_results$bootstraps)
        # actual threshods of consistency for taf1 for each lvl
        taf1_threshes <- 
            bootstrap_results %>%
            filter(ensemblid == taf1_ensemblid) %>%
            select(expression_pattern, fraction_of_bootstraps) %>%
            pivot_wider(
                names_from=expression_pattern,
                values_from=fraction_of_bootstraps
            ) %>%
            select(all_of(significance_lvls)) %>%
            as.list()
        # get genes more consistent than taf1 per pattern
        consistent_genes <- 
            lapply(
                1:length(significance_lvls),
                function(i){
                    sig_lvl <- sig_lvls[i]
                    thresh <- taf1_threshes[[gsub('up', 'down', sig_lvl)]]
                    bootstrap_results %>%
                    filter(bootstraps / max(.$bootstraps) >= bootstrap_min) %>%
                    filter(expression_pattern == sig_lvl) %>%
                    filter(fraction_of_bootstraps >= thresh)
                }
            ) %>% 
            bind_rows()
        write_tsv(outfile, consistent_genes)
    }
    return(consistent_genes)
}

