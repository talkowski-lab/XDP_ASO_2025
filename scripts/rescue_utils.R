################################################
# Dependencies
################################################
library(glue)
library(magrittr)
library(tidyverse)
# GLM_RESULTS_DIR <- sprintf('%s/rm2Outliers_GLMResults', FROZEN_DIR)
# RESCUE_RESULTS_DIR <- sprintf('%s/rm2Outliers_GLMPairComparison', FROZEN_DIR)
# BASELINE_XDP_FILE <- sprintf('%s/NO.CON+NO.XDP-LFC+Var.tsv', RESCUE_RESULTS_DIR)
# ALL_CORRECTED_COUNTS <- sprintf('%s/all.adjusted.counts.tsv', RESCUE_RESULTS_DIR)
# CLEAN_RESCUE_DATA_TSV <- sprintf('%s/all.rescue.data.tsv', RESCUE_RESULTS_DIR)

################################################
# Functions
################################################
load_xdp_results <- function(output_file=BASELINE_XDP_FILE){
    # Load the baseline DEG results + add lfc variance estimated from models
    if (file.exists(output_file)){
        return(read_tsv(output_file))
    }
    # Load full model results to get confidence/uncertainty from models objects
    xdp_glms <- 
        sprintf('%s/NO.CON+NO.XDP-GLMData.rds', GLM_RESULTS_DIR) %>% 
        readRDS()
    # Load DEG results 
    # Add Variance of LFCs to DEG results 
    sprintf('%s/NO.CON+NO.XDP-DEGResults.tsv', GLM_RESULTS_DIR) %>%
        read_tsv() %>%
        # dplyr::select(!Contrast) %>%
        # Get variance of the LFC estimate (Genotype Beta value)
        left_join(xdp_glms@vcov %>%
                  as_tibble() %>%
                  dplyr::rename('EnsemblID'=name) %>%
                  # Turn symmetric matrix into longform table of pairs
                  add_column(beta1=rep(colnames(dplyr::select(., !EnsemblID)),
                                       length(unique(.$EnsemblID))
                             )
                  ) %>%
                  pivot_longer(!c(EnsemblID, beta1),
                               names_to='beta2',
                               values_to='var'
                  ) %>% 
                  # Keep only variance estimated for the Genotype param (no SVs
                  filter(beta1 == 'ConditionNO.XDP',
                         beta2 == 'ConditionNO.XDP'
                  ) %>% 
                  dplyr::select(!starts_with('beta')),
                  by='EnsemblID'
        ) %T>%
        write_tsv(output_file)
}

load_corrected_counts <- function(
    outfile=ALL_CORRECTED_COUNTS,
    sample_metadata=NULL,
    gene_metadata=NULL){
    # Load cached results if they exist
    if (file.exists(outfile)) { return(read_tsv(outfile)) } 
    # Load adjusted counts for Untreated samples
    nonaso_sv_counts <- 
        c(
          # 'NO.CON+NO.XDP_naive-SVAdjustedCounts.tsv',
          # 'NO.CON+NO.XDP_non_edit-SVAdjustedCounts.tsv',
          'NO.CON+NO.XDP+NO.dSVA-SVAdjustedCounts.tsv',
          'NO.CON+NO.XDP-SVAdjustedCounts.tsv'
        ) %>%
        lapply(function(filename){
                   sprintf('%s/%s', GLM_RESULTS_DIR, filename) %>%
                       read_tsv()
               }
        ) %>%
        setNames(c('dSVA', 'XDP'))
        # setNames(c('dSVA', 'XDP', 'XDP_naive', 'XDP_non_edit'))
    # Load adjusted counts for Treated samples
    aso_sv_counts <- 
        ASO_TREATMENTS %>%
        sapply(function(treatment){
                   sprintf('%s/NO.CON+NO.XDP+%s.XDP-SVAdjustedCounts.tsv',
                           GLM_RESULTS_DIR, treatment
                   ) %>% 
                   read_tsv() 
               },
               simplify=FALSE,
               USE.NAMES=TRUE
        ) %>%
        setNames(paste('ASO', names(.), sep=''))
    # Combine all counts into a single tidy matrix
    c(aso_sv_counts, nonaso_sv_counts) %>%
        # {.[grep('ASO|dSVA|XDP', names(.), value=T)]} %>%
        {sapply(names(.),
                function(name){
                    .[[name]] %>% 
                    pivot_longer(!EnsemblID,
                                 names_to='SampleID',
                                 values_to='counts'
                    ) %>%
                    add_column(Treatment=gsub('\\.', '', name))
                },
                simplify=FALSE
         )
        } %>%
        bind_rows() %>%
        left_join(sample_metadata %>%
                  mutate(isTreated=(Treatment != 'NO')) %>%
                  select(SampleID, Genotype, isTreated),
                  multiple='all',
                  by='SampleID'
        ) %>% 
        left_join(gene_metadata %>%
                  filter(type == 'gene') %>%
                  dplyr::select(gene_id, gene_name),
                  by=join_by(EnsemblID == gene_id)
        # ) %>% 
        # filter(Treatment == 'XDP' | 
        #        isTreated | 
        #        (Treatment == 'dSVA' & Genotype == 'dSVA')
        # ) %>%
        # select(!isTreated) %>%
        # mutate(Genotype=ifelse(Genotype == 'dSVA', 'XDP', Genotype),
        #        Treatment=ifelse(Treatment == 'XDP', 'Untreated', Treatment),
        #        Category=paste(Treatment, Genotype, sep='.')
        ) %T>%
        write_tsv(outfile)
}

################################################
# Get LFCs from all sample pairs
################################################
get_dsva_sample_pairs <- function(
    corrected_counts,
    sample_metadata){
    long_count_data <- 
        corrected_counts %>%
        # long form dataframe
        pivot_longer(!EnsemblID,
                     names_to='SampleID',
                     values_to='SVACounts'
        ) %>%
        # Join Patient/Clone ID to gene counts
        left_join(sample_metadata,
                  multiple='all',
                  by='SampleID'
        ) %>%
        dplyr::select(EnsemblID,
                      Genotype,
                      PatientID,
                      SVACounts
        ) %>%
        filter(Genotype != 'CON') %>%
        # Add index for each sepearte sample that can be compared
        # i.e. same gene + PatientID + Genotype
        group_by(EnsemblID, PatientID, Genotype) %>%
        mutate(GenotypeSampleID=sequence(n())) %>%
        pivot_wider(names_from='Genotype',
                    values_from='SVACounts',
                    names_prefix='SVACounts_'
        ) 
    # Identify all dSVA/XDP Sample pairs that come from the same Patient
    all_sample_pairs <- 
        long_count_data %>%
        dplyr::select(!starts_with('SVACounts_')) %>%
        # Generate all pairs of comparable samples
        group_by(EnsemblID, PatientID) %>%
        crossing(GenotypeSampleID, GenotypeSampleID) %>%
        dplyr::select(!GenotypeSampleID) %>%
        # Remove redundant sample pairs
        distinct() %>%
        rename_with(~ gsub('...', '', ., fixed=TRUE)) %>%
        # Get XDP corrected counts for first sample of each pair
        left_join(long_count_data %>%
                  dplyr::select(!SVACounts_dSVA),
                  by=join_by(EnsemblID,
                             PatientID,
                             GenotypeSampleID3 == GenotypeSampleID
                  )
        ) %>%
        # Get dSVA corrected counts for the second sample of each pair
        left_join(long_count_data %>%
                  dplyr::select(!SVACounts_XDP),
                  by=join_by(EnsemblID,
                             PatientID,
                             GenotypeSampleID2 == GenotypeSampleID
                  )
        ) %>%
        # Remove pairs where XDP or dSVA counts are missing
        filter(if_all(starts_with('SVACounts_'), ~ !is.na(.))) %>%
        # rename columns for downstream convenience
        dplyr::rename('SVACounts_Treated'=SVACounts_dSVA,
                      'SVACounts_Untreated'=SVACounts_XDP,
        ) %>%
        add_column(Treatment='dSVA')
}

get_aso_sample_pairs <- function(
    corrected_counts,
    sample_metadata,
    treatment){
    sample_pairs <- 
        corrected_counts %>%
        pivot_longer(!EnsemblID,
                     names_to='SampleID',
                     values_to='SVACounts'
        ) %>%
        left_join(sample_metadata,
                  by='SampleID',
                  multiple='all'
        ) %>%
        dplyr::select(EnsemblID,
                      Treatment,
                      PatientID,
                      Clone,
                      geno,
                      SVACounts
        ) %>%
        # Remove Control samples 
        filter(geno != 'CON') %>%
        # rename columns to be treatment agnostic
        mutate(Treatment=case_when(Treatment == treatment ~ 'Treated',
                                   Treatment == 'NO' ~ 'Untreated'
               )
        ) %>%
        # Pivot Treated/Untreated Counts to their own columns
        pivot_wider(names_from='Treatment',
                    values_from='SVACounts',
                    names_prefix='SVACounts_'
        ) %>%
        # mostly length 1 lists() but still need to unnest
        unnest(starts_with('SVACounts_')) %>%
        # Remove pairs where Treated or Untreated counts are missing
        filter(if_all(starts_with('SVACounts_'),
                      ~ !is.na(.)
               )
        ) %>% 
        add_column(Treatment=treatment)
    return(sample_pairs)
}

compute_CIs <- function(
    aso_pair_df,
    dsva_pair_df,
    baseline_df){
    # Compute 95% Confidence Intervals of LFCs across sample pairs per comparison 
    # e.g. NO.XDP vs 880.XDP ~ 34 Sample Pairs LFCs -> 1 mean LFC + 95% CI

    # First compute the number of sample pairs to use for CI calculation, 
    # since for dSVA the number of pairs would be highly inflated due to 
    # no sample-specific pairing (only patient lvl) for dSVAs
    aso_pair_df <- aso_pair_df %>% mutate(Treatment=as.character(Treatment))
    dsva_pair_df <- dsva_pair_df %>% mutate(Treatment=as.character(Treatment))
    n_samples_df <- 
        aso_pair_df %>%
        group_by(
            Treatment,
            EnsemblID
        ) %>% 
        # Count samples pairs per treatment
        summarize(
            n_samples=n(),
            .groups='drop'
        ) %>%
        # Use larger number of samples from XDP/dSVA groups
        bind_rows(
            dsva_pair_df %>% 
            group_by(Treatment, EnsemblID, PatientID) %>% 
            summarize(n_samples=max(GenotypeSampleID2, GenotypeSampleID3)) %>%
            summarize(n_samples=sum(n_samples), .groups='drop')
        ) %>%
        distinct(Treatment, n_samples)
    # Now combine the aso and dsva pairs dfs together
    bind_rows(
        dsva_pair_df %>% 
        dplyr::select(-c(GenotypeSampleID2, GenotypeSampleID3,)),
        aso_pair_df %>% 
        dplyr::select(-c(Clone, geno))
    ) %>%
    # filter genes with -ve corrected counts to get LFC
    filter(if_all(starts_with('SVACounts_'), ~ . >= 0)) %>%
    # filter(if_all(starts_with('SVACounts_'), ~ . != -Inf)) %>%
    # Compute LFC for each sample pair (counts on log scale)
    mutate(lfc=SVACounts_Treated - SVACounts_Untreated) %>% #{.} -> all_pair_df
    # filter genes with negative counts to get LFC
    group_by(
        Treatment, 
        EnsemblID
    ) %>%
    # Calculate mean/var across sample pairs for each gene, each treatment
    summarize(
        var=var(lfc),
        lfc=mean(lfc),
        .groups='drop'
    ) %>%
    # Bind sample sizes for each treatment 
    left_join(
        n_samples_df,
        multiple='all',
        by= join_by(Treatment)
    ) %>%
    # Join Baseline LFCs from GLM 
    bind_rows(
        baseline_df %>%
        dplyr::select(EnsemblID, lfc, var) %>%
        add_column(
            Treatment='Control',
            n_samples=1  # since this is directly from GLM
        )
    ) %>%
    # Calculate confidence intervals for each comparison
    mutate(
        CI95.lower=lfc - qnorm(0.975) * (sqrt(var) / sqrt(n_samples)),
        CI95.upper=lfc + qnorm(0.975) * (sqrt(var) / sqrt(n_samples))
    ) 
}

################################################
# Compute statistical test use to define gene rescues
################################################
make_naive_rescue_results <- function(inverse_basline_df){
    # Get all the precomputed DESeq2 results for the relecant contrasts
    tibble(
        v1=c('dSVA', ASO_TREATMENTS),
        v2=c('NO.dSVA', paste0(ASO_TREATMENTS, '.XDP'))
    ) %>%
    mutate(
        filepath=
            file.path(
                GLM_RESULTS_DIR,
                glue('NO.CON+NO.XDP+{v2}-DEGResults.tsv')
            )
    ) %>%
    dplyr::select(v1, filepath) %>% 
    deframe() %>% 
    sapply(
        read_tsv,
        show_col_types=FALSE,
        simplify=FALSE
    ) %>% 
    bind_rows(.id='Treatment') %>%
    mutate(Contrast=str_remove_all(Contrast, 'Condition')) %>% 
    separate_wider_delim(
        Contrast,
        delim=' vs ',
        names=c('numerator', 'denominator')
    ) %>% 
    filter(
        denominator == 'NO.CON',
        numerator != 'NO.XDP'
    ) %>%
    bind_rows(
        inverse_baseline_df %>%
        dplyr::select(-c(var)) %>%
        add_column(
            Treatment='Control',
            numerator='NO.CON',
            denominator='NO.XDP'
        ) 
    ) %>%
    dplyr::rename('padj'=fdr) %>%
    mutate(isRescued=pvalue > 0.05)
}

ztest_vs_XDP <- function(
    lfc,
    var,
    n_samples,
    mu=0,  # H_0 -> no difference 
    conf.level=0.95,
    ...){
    # Taken from the z.test source code
    # 1-sample 2-sided z test to see if the mean LFC for a comparison is 0
    stderr <- sqrt(var) / sqrt(n_samples)
    zobs <- (lfc - mu) / stderr
    pvalue <- 2 * pnorm(-abs(zobs))
    alpha <- 1 - conf.level
    CI.upper <- mu + zobs * stderr + qnorm((1 - alpha / 2)) * stderr
    CI.lower <- mu + zobs * stderr - qnorm((1 - alpha / 2)) * stderr 
    list(
        statistic=zobs,
        stderr=stderr,
        CI.lower=CI.lower,
        CI.upper=CI.upper,
        pvalue=pvalue
    ) %>%
    c(., ...) %>%
    # c(., ..., lfc=lfc, var=var, n_samples=n_samples) %>%
    as_tibble_row()
}

ztest_vs_CON <- function(
    lfc.vsTreated, 
    var.vsTreated,
    n_samples.vsTreated,
    lfc.vsCON,
    var.vsCON,
    n_samples.vsCON,
    mu=0,  # H_0 -> no difference 
    conf.level=0.95,
    ...){
    # Taken from the z.test source code, but skipping the mean/std calculations
    # 2-sample 2-sided z test to see if the difference in mean LFC 
    # between the samples is statistically signigicant for this gene
    stderr <- 
        sqrt(
            (var.vsTreated / n_samples.vsTreated) + 
            (var.vsCON     / n_samples.vsCON)
        )
    zobs <- (lfc.vsTreated - lfc.vsCON - mu) / stderr
    pvalue <- 2 * pnorm(-abs(zobs))
    alpha <- 1 - conf.level 
    CI.lower <- mu + zobs * stderr - qnorm((1 - alpha / 2)) * stderr
    CI.upper <- mu + zobs * stderr + qnorm((1 - alpha / 2)) * stderr
    list(
        statistic=zobs,
        stderr=stderr,
        CI.lower=CI.lower,
        CI.upper=CI.upper,
        pvalue=pvalue
    ) %>%
    c(., ...) %>%
    as_tibble_row()
}

get_ztest_results <- function(
    ci_df,
    mu.Baseline=0,
    mu.Treatment=0,
    conf.level=0.95){
    # 2 Z-Tests for determining Rescue status of genes 
    # Test 1: Sig. test comparing LFC of NO.XDP vs ASO.XDP (gene is affected)
    results_vs_XDP <- 
        ci_df %>%
        dplyr::select(-c(starts_with('CI95'))) %>%
        filter(Treatment != 'Control') %>%
        future_pmap(
            ztest_vs_XDP,
            conf.level=conf.level,
            mu=mu.Treatment,
            .progress=TRUE
        ) %>%
        bind_rows() %>%
        mutate(
            numerator=glue('{Treatment}.XDP'),
            denominator='NO.XDP'
        )
    # Test 2: N.S test comparing LFCs of NO.CON vs NO.XDP && NO.XDP vs ASO.XDP 
    #         (ASO XPD and Untreated CON are ~ equally different from NO.XDP)
    results_vs_CON <- 
        ci_df %>% 
        filter(Treatment != 'Control') %>%
        dplyr::select(-c(starts_with('CI95'))) %>%
        # Join Baseline LFCs from GLM 
        left_join(
            ci_df %>%
            filter(Treatment == 'Control') %>%
            dplyr::select(EnsemblID, lfc, var, n_samples),
            suffix=c('.vsTreated', '.vsCON'),
            multiple='all',
            by='EnsemblID'
        ) %>%
        # Compute Z.Test on each row comparing Baseline vs Treatment LFC
        future_pmap(
            ztest_vs_CON,
            conf.level=conf.level,
            mu=mu.Baseline,
            .progress=TRUE
        ) %>%
        bind_rows() %>%
        mutate(
            numerator=glue('{Treatment}.XDP'),
            denominator='NO.CON'
        )
    # Join test results together with LFC stats
    ztest_rescue_data <-
        ci_df %>% 
        filter(Treatment != 'Control') %>%
        left_join(
            bind_rows(
                results_vs_XDP,
                results_vs_CON
            ),
            by=
                join_by(
                    Treatment,
                    EnsemblID
                )
        ) %>%
        # Adjust p-values per treatment (each gene is an indpt. test)
        group_by(
            Treatment,
            numerator,
            denominator
        ) %>% 
        mutate(padj=p.adjust(pvalue, method='BH')) %>%
        ungroup()
}

