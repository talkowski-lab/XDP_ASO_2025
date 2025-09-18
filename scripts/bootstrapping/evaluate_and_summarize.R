################################################
# Dependencies
################################################
library(here)
here::i_am('scripts/bootstrapping/evaluate_and_summarize.R')
BASE_DIR <- here()
SCRIPTS_DIR <- here('scripts')
source(file.path(SCRIPTS_DIR, 'bootstrapping/utils.R'))
source(file.path(SCRIPTS_DIR, 'locations.R'))
source(file.path(SCRIPTS_DIR, 'constants.R'))
source(file.path(SCRIPTS_DIR, 'data_utils.R'))
library(glue)
library(purrr)
library(tidyverse)
library(furrr)

################################################
# Functions
################################################
analyse_bootstrap <- function(
    samples,
    swapped_labels,
    index,
    counts,
    meta,
    comparison,
    method,
    analysis_params, 
    p=NA){
    # subset full data for only this bootstrap's specified sample count and meta data
    comparison_vector <- str_split(comparison, fixed('.'), simplify=TRUE)
    bootstrap_samples <- as.vector(str_split(samples, ',', simplify=TRUE))
    bootstrap_meta <- 
        level_metadata(
            comparison_vector,
            meta, 
            control_label=paste(comparison_vector[1:2], collapse='.'),
            condition_label=paste(comparison_vector[3:4], collapse='.')
        ) %>% 
        filter(SampleID %in% bootstrap_samples)
    # Get the counts and permute the labels if specified (not used)
    bootstrap_counts <- counts[, bootstrap_meta$SampleID]
    if (analysis_params$swap){  # if want to swap labels
        swapped_labels <- as.vector(str_split(swapped_labels, ',', simplify=TRUE))
        colnames(bootstrap_counts) <- swapped_labels
        bootstrap_meta$SampleID <- swapped_labels
    }
    # filter out low expression genes for this bootstrap
    genes_to_keep <- 
        CPM_filter_counts(
            bootstrap_counts,
            bootstrap_meta, 
            condition_var='Condition'
        )
    bootstrap_counts <- bootstrap_counts[genes_to_keep, ]
    # preform DEG analysis on the bootstrap samples with the specified parameters
    results <- 
        evaluation_wrapper(
            method, 
            bootstrap_counts, 
            bootstrap_meta, 
            analysis_params,
            title_str=comparison,
            output_file=NA
        )
    if (!(is.na(results))){  
        results <- as_tibble(results$results) 
    } 
    if (!(is.na(p))) p()  # increment progress bar
    return(results)
}

fix_column_names <- function(bootstrap_results, method){
    if (method == 'glm' | method == 'glm'){
        colnames(bootstrap_results) <- c('padj', 'lfc', 'pvalue', 'EnsemblID')
    } else if (method == 'deseq'){
        colnames(bootstrap_results) <- c('EnsemblID', 'baseMean', 'lfc', 'lfcSE', 'stat', 'pvalue', 'padj')
    }
    return(bootstrap_results)
}

summarize_bootstap_results <- function(bootstrap_results, alpha_nom, alpha_fdr){
    # summarize the DEG results for a chunk of analyzed bootstraps
    # so we can avoid caching all DEG results and just keep summaries
    # basically count frequencies/values for each gene across all bootstraps
    # as the columns 
    de_frequencies <- bootstrap_results %>%
                      group_by(EnsemblID) %>%
               dplyr::summarize(bootstraps=n(),
                                mean_lfc=mean(lfc),
                                var_lfc=var(lfc),
                                raw_up=sum(lfc > 0),
                                frac_up=raw_up / n(),
                                frac_down=sum(lfc < 0) / n(),
                                raw_nominal=sum(pvalue <= alpha_nom),
                                frac_nominal=raw_nominal / n(),
                                raw_fdr=sum(padj <= alpha_fdr),
                                frac_fdr=raw_fdr / n(),
                                frac_nominal_and_up=sum(lfc > 0 & pvalue <= alpha_nom) / n(),
                                frac_nominal_and_down=sum(lfc < 0 & pvalue <= alpha_nom) / n(),
                                frac_fdr_and_up=sum(lfc > 0 & padj <= alpha_fdr) / n(),
                                frac_fdr_and_down=sum(lfc < 0 & padj <= alpha_fdr) / n()
                      ) 
    return(de_frequencies)
}

main <- function(counts, meta, bootstraps_file, method, analysis_params){
    # set up output file
    output_dir <- dirname(bootstraps_file)
    swapped <- c('original', 'swapped')[as.numeric(analysis_params$swap)+1]
    comparison <- dirname(bootstraps_file) %>% basename()
    bootstrap_index <- str_split(bootstraps_file, fixed('.')) %>% 
                       unlist() %>% 
                       last() #%>% as.integer()
    output_file <- sprintf("%s/summary-%s-%s-%s.tsv.%s",
                           output_dir,
                           method,
                           swapped,
                           comparison,
                           bootstrap_index
    )
    # read in previously generated bootstraps 
    bootstraps <- read_tsv(bootstraps_file, 
                           col_names=c('samples', 'swapped_labels')
                  ) %>%
                  tibble::rowid_to_column('index')
    # split bootstraps into chunks
    num_chunks <- (nrow(bootstraps) %/% analysis_params$chunk_size)
    if (nrow(bootstraps) %% analysis_params$chunk_size != 0){
        num_chunks <- num_chunks + 1
    }
    for (chunk in 1:num_chunks){
        # for a single chunk, get all the bootstraps (n=200)
        # chunking so R can handle loading/computing these bootstraps all at once,
        # saving them without crashing
        chunk_start <- 1 + (chunk - 1) * analysis_params$chunk_size 
        chunk_end <- min(chunk * analysis_params$chunk_size, nrow(bootstraps))
        chunk_output_file <- sprintf('%s.%s-%s',
                                     output_file,
                                     chunk_start,
                                     chunk_end
                             )
        if (file.exists(chunk_output_file)){ next }
        message(sprintf("Analyzing DEGs for bootstraps chunk %s of %s...",
                        chunk,
                        num_chunks
                )
        )
        # Now compute DEG results for this chunk in parallel 
        plan(multisession, workers=analysis_params$workers)
        with_progress({
            bootstrap_results <- future_pmap(bootstraps[chunk_start:chunk_end, ],
                                             analyse_bootstrap,
                                             counts=counts,
                                             meta=meta,
                                             comparison=comparison,
                                             method=method,
                                             analysis_params=analysis_params,
                                             .progress=TRUE,
                                             .options=furrr_options(stdout=TRUE,
                                                                    seed=9
                                             )
            )
        })
        plan(sequential)
        # filter our NA's i.e. failed permutations
        bootstrap_results <- bootstrap_results[!is.na(bootstrap_results)]
        message(sprintf("%s of %s permutations passed in chunk %s", 
                        length(bootstrap_results),
                        analysis_params$chunk_size,
                        chunk
                )
        )
        message(sprintf("Summarizing results for bootstraps chunk %s of %s...",
                        chunk,
                        num_chunks
                )
        )
        # Now summarize the full DEG results into 1 table and save it to a tsv file
        de_frequencies <- bootstrap_results %>%
                          bind_rows() %>%
                          fix_column_names(., method) %>%
                          summarize_bootstap_results(.,
                                                     alpha_nom=analysis_params$alpha_nom,
                                                     alpha_fdr=analysis_params$alpha_fdr
                          )
        write_tsv(de_frequencies, chunk_output_file)
        message(sprintf("finished summarizing bootstrap chunk %s of %s",
                        chunk,
                        num_chunks
                )
        )
        message(sprintf("results saved to %s", chunk_output_file))
    }
    message(sprintf('Finished, summarized %s bootstraps from set %s',
                    nrow(bootstraps),
                    chunk
            )
    ) 
}

################################################
# Compute + summarize reuslts for this chunk of bootstraps
################################################
counts <- load_counts()
meta <- load_meta()
analysis_params <- list()
args <- commandArgs(TRUE)
bootstraps_file <- args[1]                         # specified set of bootstraps to analyze
method <- args[2]                                  # DEG anlysis method
analysis_params$workers <- as.integer(args[3])     # workers for parallelizing
analysis_params$chunk_size <- as.integer(args[4])  # permutations to compute before summarizing 
analysis_params$swap <- 0                          # dont swap samples labels
analysis_params$design_var <- 'Condition'          # column variable for DE analysis
analysis_params$alpha_nom <- 0.05                  # alpha threshold for nominal significance 
analysis_params$alpha_fdr <- 0.1                   # alpha threshold for FDR significance
message(
    paste(
        c('input', 'method', names(analysis_params)),
        c(bootstraps_file, method, unlist(analysis_params)),
        sep=': ',
        collapse='\n'
    )
)
main(counts, meta, bootstraps_file, method, analysis_params)
# Author: Siddharth Reed, Talkowski Lab
# To compute DEG results an summaries for a specified set of bootstraps (produced by permutations.R)
# Example Usage 
# $ ./evaluate_and_summarize.R 
#   '/data/talkowski/Samples/XDP/ASO_dSVA/sid_results/BootstrapResults/NO.CON.NO.CombinedXDP/bootstraps-NO.CON.NO.CombinedXDP.tsv.001' 
#   'glm'
#   16
#   200
