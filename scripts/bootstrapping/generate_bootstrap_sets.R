################################################
# Dependencies
################################################
library(here)
here::i_am('scripts/bootstrapping/generate_bootstrap_sets.R')
BASE_DIR <- here()
SCRIPTS_DIR <- here('scripts')
source(file.path(SCRIPTS_DIR, 'locations.R'))
source(file.path(SCRIPTS_DIR, 'constants.R'))
source(file.path(SCRIPTS_DIR, 'glm_utils.R'))
library(tidyverse)
OUTPUT_DIR <- file.path(BOOTSTRAP_RESULTS_DIR, 'NO.CON.NO.CombinedXDP')
mkdir(OUTPUT_DIR)

################################################
# Generate random subset of samples to use as bootstraps
################################################
generate_bootstraps <- function(
    sample_metadata, 
    group_col,
    sampleID_col,
    n_bootstraps,
    chunk_size,
    sample_size,
    seed){ 
    # group_col='Condition'; sampleID_col='SampleID'; n_bootstraps=550; chunk_size=100; sample_size=6; seed=9
    set.seed(seed) 
    # Groups of samples 
    groups <- 
        sample_metadata[[group_col]] %>% 
        as.character() %>%
        unique() %>%
        sort()
    if (length(groups) != 2) { stop('Comparison must have 2 groups only') }
    # name template for output files
    comparison <- groups %>% paste(collapse='.')
    # Samples in group 1
    group1_samples <- 
        sample_metadata %>%
        filter(get({{group_col}}) == groups[1]) %>%
        pull(get({{sampleID_col}}))
    # Samples in group 2
    group2_samples <- 
        sample_metadata %>%
        filter(get({{group_col}}) == groups[2]) %>%
        pull(get({{sampleID_col}}))
    # Produce bootstraps of samples
    # Chunk output so there are chunk_size bootstraps per file
    total=0; chunk=1
    while(total < n_bootstraps){
        # Create a list of chunk_size bootstraps, write to a single file
        lapply(
            1:chunk_size,
            function(j){
                # break if we hit the total number of bootstraps
                if (j > n_bootstraps - (i * chunk_size)) { break }
                # For a single bootstrap pick random samples per group
                all_samples <- 
                c(
                    sample(group1_samples, sample_size),
                    sample(group2_samples, sample_size)
                )
                # Return single row dataframe representing 1 bootstrap
                bootstrap <- 
                    tibble(
                        samples=paste(all_samples, collapse=','),
                        swapped_labels=paste(sample(all_samples), collapse=',')
                    )
            }
        ) %>%
        bind_rows() %>%
        write_tsv(
            sprintf(
                "%s/bootstraps-%s.tsv.%s",
                OUTPUT_DIR, comparison, chunk
            )
        )
        total <- chunk * chunk_size 
        chunk <- chunk + 1
    }
}

################################################
# Main Call
################################################
# args
args <- commandArgs(TRUE)
n_bootstraps <- as.integer(args[1])  # total number of bootstraps to generate
chunk_size <- as.integer(args[2])    # How many bootstraps to write per output file
sample_size <- as.integer(args[3])   # n/2 samples to use per bootstraps
seed <- as.integer(args[4])          # random seed
# Main
load_meta() %>% 
    mutate(Condition=paste(Treatment, Genotype, sep='.')) %>%
    filter(Condition %in% c('NO.CON', 'NO.XDP')) %>%
    generate_bootstraps(
        group_col='Condition',
        sampleID_col='SampleID',
        n_bootstraps=n_bootstraps,
        chunk_size=chunk_size,
        sample_size=sample_size,
        seed=seed
    )
# Author: Siddharth Reed, Talkowski Lab
# Script to generate bootstraps for DEG analyses (sets of sample names)
# each bootstrap is just a set of sample names chosen, they are saved in .tsv
# files and used later to compute DEG results
# General Usage 
# $ Rscript ./permutations.R 
    # ${number of permutations} 
    # ${permutations per file} 
    # ${n_samples to use for each condition in each permutation} 
    # ${random seed}
# Example
# Rscript ./generate_bootstrap_sets.R 100000 1000 6 9
