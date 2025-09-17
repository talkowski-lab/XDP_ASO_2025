################################################
# Dependencies
################################################
library(here)
here::i_am('scripts/paired_analysis.R')
BASE_DIR <- here()
source(file.path(BASE_DIR, 'scripts/locations.R'))
source(file.path(BASE_DIR, 'scripts/constants.R'))
source(file.path(BASE_DIR, 'scripts/glm_utils.R'))
source(file.path(BASE_DIR, 'scripts/rescue_utils.R'))
library(grid)
library(gtable)
library(cowplot)
library(ComplexUpset)
INPUT_DIR <- GLM_RESULTS_DIR
OUTPUT_DIR <- FIGURES_DIR
OUTPUT_DATA_DIR <- file.path(FIGURES_DIR, 'rescue.data')
mkdir(OUTPUT_DATA_DIR)

################################################
# Load input data
################################################
sample_metadata <- load_meta()
# EnsemblID to gene name mappings
gene_metadata <- load_gtf()
# Baseline Untreated XDP vs CON LFCs + Variance
baseline_df <- load_xdp_results(output_file=BASELINE_DF_FILE)

################################################
# Identify pairs of samples to compare for dSVA vs. XDP
################################################
dsva_pair_df <- 
    sprintf(
        '%s/NO.CON+NO.XDP+NO.dSVA-SVAdjustedCounts.tsv', 
        INPUT_DIR
    ) %>%
    read_tsv() %>%
    # Join corrected counts to sample data including PatientIDs
    get_dsva_sample_pairs(sample_metadata) %T>%
    write_tsv(file.path(OUTPUT_DATA_DIR, 'dSVA_SamplePairs.tsv'))

################################################
# Identify pairs of samples to compare for ASO Treated XDP vs. Untreated XDP
################################################
aso3_pair_df <- 
    ASO_TREATMENTS %>%
    lapply(function(treatment){
             # Load corrected counts
             sprintf('%s/NO.CON+NO.XDP+%s.XDP-SVAdjustedCounts.tsv',
                     INPUT_DIR, treatment
             ) %>%
             read_tsv(show_col_types=FALSE) %>%
             # Get all Treated/Untreated matched sample pairs 
             get_aso_sample_pairs(sample_metadata,
                                  treatment
             )
       }
    ) %>% 
    bind_rows() %T>%
    write_tsv(file.path(OUTPUT_DIR, 'ASO_lvl3_SamplePairs.tsv'))

################################################
# Combined dSVA and ASO treatment pairs 
################################################
all3_pair_df <-
    combine_sample_pairs(dsva_pair_df,
                         aso3_pair_df
    ) %>% 
    # Add gene names 
    left_join(gene_metadata %>%
            filter(type == 'gene') %>%
            dplyr::select(gene_id, gene_name),
            by=join_by(EnsemblID == gene_id)
    ) %T>%
    write_tsv(file.path(OUTPUT_DATA_DIR, 'All_lvl3_SamplePairs.tsv'))

################################################
# Compute 95% CI for the mean LFCs estimated from the paired treatment data
################################################
ci3_df <- 
    compute_CIs(all3_pair_df,
                dsva_pair_df,
                aso3_pair_df,
                # Flip LFC values from  XDP / Control -> Control / XDP
                # since for all treatments the LFC represents Treated / XDP
                # This way all comparisons have the same baseline (denom): XDP
                baseline_df %>% mutate(lfc=-lfc)
    ) %>%
    # Add gene names 
    left_join(gene_metadata %>%
            filter(type == 'gene') %>%
            dplyr::select(gene_id, gene_name),
            by=join_by(EnsemblID == gene_id)
    ) %T>%
    write_tsv(file.path(OUTPUT_DATA_DIR, 'All_lvl3_LFC+CIs.tsv'))

################################################
# Check how many pairs had negative counts for at least 1 samples
################################################
ndp=nrow(dsva_pairs); ntp=nrow(aso3_pair_df); nap=nrow(all3_pair_df) 
sprintf('Total treated pairs: %s 
Pairs filtered out:  %s 
Fraction retained:   %.5f', 
ndp + ntp, 
ndp + ntp - nap, 
nap / (ndp + ntp)
) %>% 
message()

################################################
# 1. Naive test for determining Rescue status of genes 
################################################
#    - N.S test comparing NO.CON vs ASO.XDP (Treated XDP returns to Control lvl)

# Get DEG results comparing NO.CON vs ASO.XDP and any N.S. => rescued to normal
# Baseline in these results is Untreated CON
# Define "Rescued" genes as having a pvalue > 0.05 i.e.
# Untreated Control is not significantly different from Treated XDP
naive_rescue_data <- 
    c('NO.dSVA', paste0(ASO_TREATMENTS, '.XDP')) %>%
    lapply(
        function(treatment){
            file.path(INPUT_DIR, glue('NO.CON+NO.XDP+{treatment}-DEGResults.tsv')) %>% 
            read_tsv(show_col_types=FALSE) %>%
            filter(Contrast == sprintf('Condition%s vs ConditionNO.CON', treatment)) %>% 
            add_column(Treatment=treatment)
       }
    ) %>%
    bind_rows() %>%
    dplyr::select(!Contrast) %>%
    # Define Rescues 
    mutate(
        Treatment=gsub('NO\\.|\\.XDP', '', Treatment),
        isRescued=
            case_when(
                pvalue > 0.05 ~ 'Rescued',
                TRUE ~ 'Not Rescued'
           )
    ) %T>%
    write_tsv(file.path(OUTPUT_DATA_DIR, 'NaiveTest_RescueData.tsv'))

