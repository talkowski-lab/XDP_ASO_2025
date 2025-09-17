################################################
# Dependencies
################################################
library(here)
here::i_am('scripts/generate_GLMResults.R')
BASE_DIR <- here()
SCRIPT_DIR <- here('scripts')
source(file.path(SCRIPT_DIR, 'constants.R')
source(file.path(SCRIPT_DIR, 'locations.R')
source(file.path(SCRIPT_DIR, 'glm_utils.R')
OUTPUT_DIR <- GLM_RESULTS_DIR 

################################################
# Load all raw data to build GLMs SD
################################################
# Load sample metadata
sample_metadata <- 
    load_meta() %>% 
    mutate(Condition=paste(Treatment, geno, sep='.'))
# Load all raw counts from STAR
raw_counts <- load_counts()
# Filter out samples meant to be excluded
raw_counts <- raw_counts[, sample_metadata$SampleID]
# EnsemblID to gene name mappings
gene_metadata <- load_gtf()

################################################
# NO.CON + NO.XDP_naive
################################################
sample_metadata %>% 
    filter(Condition %in% c('NO.CON', 'NO.XDP')) %>%
    mutate(
        Condition=
            case_when(
                geno == 'XDP' ~ 'NO.XDP_naive',
                TRUE ~ Condition
            ) %>%
            factor(levels=c('NO.CON', 'NO.XDP_naive'))
    ) %>% 
    glm_wrapper(
        raw_counts,
        # name_only=TRUE,
        'Condition',
        output_dir=OUTPUT_DIR
    )

################################################
# NO.CON + NO.XDP_non_edit
################################################
sample_metadata %>% 
    filter(Condition %in% c('NO.CON', 'NO.XDP_non_edit')) %>%
    mutate(Condition=factor(Condition, levels=c('NO.CON', 'NO.XDP_non_edit'))) %>% 
    glm_wrapper(
        raw_counts,
        # name_only=TRUE,
        'Condition',
        output_dir=OUTPUT_DIR
    )

################################################
# NO.CON + NO.XDP 
################################################
sample_metadata %>% 
    mutate(Condition=paste(Treatment, Genotype, sep='.')) %>%
    filter(Condition %in% c('NO.CON', 'NO.XDP')) %>%
    mutate(Condition=factor(Condition, levels=c('NO.CON', 'NO.XDP'))) %>% 
    glm_wrapper(
        raw_counts,
        # name_only=TRUE,
        'Condition',
        output_dir=OUTPUT_DIR
    )

################################################
# NO.CON + NO.dSVA
################################################
sample_metadata %>% 
    mutate(Condition=paste(Treatment, Genotype, sep='.')) %>%
    filter(Condition %in% c('NO.CON', 'NO.dSVA')) %>%
    mutate(Condition=factor(Condition, levels=c('NO.CON', 'NO.dSVA'))) %>% 
    glm_wrapper(
        raw_counts,
        # name_only=TRUE,
        'Condition',
        output_dir=OUTPUT_DIR
    )

################################################
# NO.CON + NO.XDP + NO.dSVA (for paired analysis)
################################################
sample_metadata %>% 
    mutate(Condition=paste(Treatment, Genotype, sep='.')) %>%
    filter(Condition %in% c('NO.CON', 'NO.XDP', 'NO.dSVA')) %>%
    mutate(Condition=factor(Condition, levels=c('NO.CON', 'NO.XDP', 'NO.dSVA'))) %>%
    glm_wrapper(
        raw_counts,
        # name_only=TRUE,
        'Condition',
        output_dir=OUTPUT_DIR
    )

################################################
# Generate results for 3 levels and 5 levels in GLM for all ASOs
################################################
ASO_TREATMENTS %>% 
lapply(
    function(treatment){
    # 5 level results
    conditions <- 
        c(
            'NO.XDP',
            'NO.XDP_non_edit',
            sprintf('%s.XDP', treatment),
            sprintf('%s.XDP_non_edit', treatment)
        )
    sample_metadata %>%
    # Create dummy condition variable 
    mutate(Condition=paste(Treatment, geno , sep='.')) %>%
    filter(Condition %in% conditions) %>%
    mutate(Condition=factor(Condition, levels=conditions)) %>%
    glm_wrapper(
        raw_counts,
        'Condition',
        output_name=sprintf('%s-CON.Intercept5', treatment),
        output_dir=file.path(MISC_RESULTS_DIR, '5vs3levels_GLMResults')
    )
    # NO.CON + NO.XDP + ASO.XDP for all ASOs (for paired analysis)
    # 3 level results
    conditions <- 
        c(
            'NO.CON',
            'NO.XDP',
            sprintf('%s.XDP', treatment)
        )
    sample_metadata %>%
    # Create dummy condition variable 
    mutate(Condition=paste(Treatment, Genotype, sep='.')) %>%
    filter(Condition %in% conditions) %>% 
    mutate(Condition=factor(Condition, levels=conditions)) %>%
    glm_wrapper(
        raw_counts,
        'Condition',
        # output_name=sprintf('%s-CON.Intercept3', treatment),
        output_dir=OUTPUT_DIR
    )
    }
)

################################################
# Save sessionInfo to txt file
################################################
writeLines(
    capture.output(sessionInfo()),
    file.path(
        OUTPUT_DIR,
        sprintf('SessionInfo_%s.txt', format(Sys.time(), "%Y-%m-%d-%H:%M"))
    )
)

