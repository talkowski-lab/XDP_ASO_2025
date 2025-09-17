################################################
# Dependencies
################################################
library(here)
here::i_am('scripts/outlier_analysis.R')
BASE_DIR <- here()
SCRIPTS_DIR <- here('scripts')
source(file.path(SCRIPTS_DIR, 'locations.R'))
source(file.path(SCRIPTS_DIR, 'constants.R'))
source(file.path(SCRIPTS_DIR, 'data_utils.R'))
source(file.path(SCRIPTS_DIR, 'glm_utils.R'))
source(file.path(SCRIPTS_DIR, 'plot_utils.R'))
library(glue)
library(purrr)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(corrplot)
INPUT_DIR <- GLM_RESULTS_DIR
# OUTPUT_DIR <- OUTLIER_ANALYSIS_DIR
OUTPUT_DIR <- FIGURES_DIR
mkdir(file.path(OUTPUT_DIR, 'PCA'))
mkdir(file.path(OUTPUT_DIR, 'results'))

################################################
# Functions
################################################
is_outside_ellipse <- function(
    PCvalue.PC1, PCvalue.PC2,
    centroid.PC1, centroid.PC2,
    sd.PC1, sd.PC2){
    # Test whether the point (PCValue.PC1, PCValue.PC2) is outside
    # an ellipse @ (centroid.PC1, centroid.PC2) w/ radiii sd.PC1, sd.PC2
    # Used for defining outliers in PCA space
    xdiff <- ((PCvalue.PC1 - centroid.PC1)**2) / sd.PC1**2
    ydiff <- ((PCvalue.PC2 - centroid.PC2)**2) / sd.PC2**2
    return((xdiff + ydiff) > 1)
}

get_outlier_from_pca <- function(
    counts,
    metadata,
    sampleColumn='SampleID',
    groupColumn='Condition',
    ntop=500){
    # Project points into PCA space and find outliers points that are 
    # outside of 2 SDs of the centroid of the group the point belongs to
    # Assumes the column `Condition` contains the sample grouping and
    # that    the column `SampleID` contains IDs unique to each sample 

    # First filter all genes to the ntop most variable genes
    samples <- metadata[[sampleColumn]]
    mat <- as.matrix(counts[, samples]) 
    # filter genes to only keep the ntop most variable genes
    genes_to_keep <- 
        mat %>%
        rowVars() %>%
        order(decreasing=TRUE) %>%
        {.[seq_len(min(ntop, length(.)))]}
    # create dataframe with intgroups and PC loadings for plotting
    pca_results <- prcomp(t(mat[genes_to_keep, ]))
    # percent_var <- 100 * (pca_results$sdev^2 / sum(pca_results$sdev^2))
    # Reshape PCa results into a longform tibble
    pc_df <- 
        pca_results$x %>%
        as.data.frame() %>%
        rownames_to_column(var=sampleColumn) %>%
        as_tibble() %>%
        left_join(
            metadata %>%
            dplyr::select(all_of(c(sampleColumn, groupColumn))),
            by=sampleColumn
        ) %>%
        pivot_longer(
            starts_with('PC'),
            names_to='PC',
            values_to='PCvalue'
        ) 
    # Calculate centroids + radii per group of samples
    centroids <- 
        pc_df %>%
        group_by(across(all_of(c(groupColumn, 'PC')))) %>%
        summarize(
            centroid=mean(PCvalue),
            sd=sd(PCvalue) * 2
        ) %>% 
        ungroup()
    # Create df where each row is a sample + its centroid info
    # so we can test each row for outlier status
    ellipse_test_df <- 
        pc_df %>%
        filter(PC %in% c('PC1', 'PC2')) %>%
        left_join(
            centroids,
            multiple='all',
            by=c(groupColumn, 'PC')
        ) %>%
        pivot_wider(
            names_from=PC,
            values_from=c(PCvalue, centroid, sd),
            names_sep='.'
        ) %>%
        mutate(
            isOutlier=
                pmap_lgl(
                    dplyr::select(
                        .,
                        !all_of(
                            c(
                                sampleColumn,
                                groupColumn
                            )
                        )
                    ),
                    is_outside_ellipse
                )
        ) 
    # Reformat ellipse test df to be easier to parse visually
    ellipse_test_df %>%
    mutate(across(where(is.factor), as.character)) %>%
    pivot_longer(
        !all_of(
            c(
                sampleColumn,
                groupColumn,
                'isOutlier'
            )
        ),
        names_to='tmp',
        values_to='value'
    ) %>%
    separate_wider_delim(
        tmp,
        delim='.',
        names=c('metric', 'PC')
    ) %>%
    pivot_wider(
        names_from=metric,
        values_from=value
    ) 
}

################################################
# Load all count + sample metadata
################################################
# Load sample metadata
sample_metadata <- load_meta(filter_outliers=FALSE)
# Load all raw counts from STAR
raw_counts <- load_counts()
# Filter out samples meant to be excluded
raw_counts <- raw_counts[, sample_metadata$SampleID]
# Load correctd counts for everything
all_sv_counts <- load_corrected_counts_list()
sapply(all_sv_counts, dim) %>% t() %T>% {colnames(.) <- c('n.genes', 'n.samples')}

################################################
# PCA of Samples on CPMs 
################################################
lapply(
    names(all_sv_counts),
    function(sample_set){
       # Plot PCA for all sample sets from raw counts (not SVAdjusted)
       samples <- colnames(all_sv_counts[[sample_set]]) 
       samples <- samples[!grepl('EnsemblID|name', samples, perl=T)]
       # Prep metadata
        meta <- 
            sample_metadata %>%
            # keep relevant samples
            filter(SampleID %in% samples) %>%
            # shorter SampleID for labelling
            mutate(labels=gsub(' ', '_', paste(PatientID, Clone))) %>%
            # Variable for grouping samples
            mutate(Condition=paste(Treatment, geno, sep='.')) %>%
            mutate(Condition=as.factor(Condition))
        # Get CPMs
        genes_to_keep <- all_sv_counts[[sample_set]]$EnsemblID
        counts <- raw_counts[genes_to_keep, samples]
        cpms <- 1e6 * t(t(counts) / meta$Uniquely_Mapped_Reads)
        pdf(file.path(OUTPUT_DIR, 'PCA', glue('{sample_set}-CPM-Sample.PCA.pdf')))
        pca_plot <- 
            rnaPCA(
                cpms, 
                meta,
                groupColumn='Condition',
                title_str=paste('Raw', sample_set),
                labelColumn='labels',
                xlim=NULL, ylim=NULL
            )
        print(pca_plot)
        dev.off()
        return(invisible(NULL))
   }
)

################################################
# PCA of Samples on adjusted counts 
################################################
lapply(
    names(all_sv_counts),
    function(sample_set){
        counts <- all_sv_counts[[sample_set]]
        meta <- 
            sample_metadata %>%
            # keep relevant samples
            filter(SampleID %in% colnames(counts)) %>%
            # shorter SampleID for labelling
            mutate(labels=gsub(' ', '_', paste(PatientID, Clone))) %>%
            # Variable for grouping samples
            mutate(Condition=paste(Treatment, geno, sep='.')) %>%
            mutate(Condition=as.factor(Condition))
        # Plot non ASO counts
        pdf(file.path(OUTPUT_DIR, 'PCA', glue('{sample_set}-Adj-Sample.PCA.pdf')))
        pca_plot <- 
            rnaPCA(
                counts, 
                meta,
                title_str=sample_set,
                groupColumn='Condition',
                labelColumn='labels',
            )
        print(pca_plot)
        dev.off()
        return(invisible(NULL))
    }
)

################################################
# Test corr of NO.CON, NO.XDP SVAdjustedCounts across Treatments
################################################
# Load adjust counts for each treatment into easy matrix
common_sv_counts <- 
    sapply(
        names(all_sv_counts) %>% grep('XDP_', ., invert=TRUE, value=TRUE),
        function(name){
            all_sv_counts[[name]] %>% 
            pivot_longer(
                !EnsemblID,
                names_to='SampleID',
                values_to=as.character(glue('{name}-SVCount'))
            ) 
        },
        simplify=FALSE
    ) %>%
    purrr::reduce(
        dplyr::inner_join, 
        by=c('EnsemblID', 'SampleID')
    )
# Get correlations matrices (treatment x treatment) for each individual sample 
corr_mats <- 
    common_sv_counts %>%
    filter(if_all(ends_with('-SVCount'), ~ .x != -Inf)) %>%
    split(f=.$SampleID) %>%
    sapply(
        function(df) { 
            df %>%
            dplyr::select(ends_with('-SVCount')) %>%
            dplyr::rename_with(~ str_remove(.x, '-SVCount')) %>% 
            cor()
        },
        simplify=FALSE,
        USE.NAMES=TRUE
    )
# Plot Avg Treatment-Treatment Correlation matrix across XDP + CON samples
pdf(file.path(OUTPUT_DIR, 'XDP+CON_Sample.SVCount.Corr.pdf'))
corr_mats %>% 
    Reduce('+', .) %>%
    {. / length(corr_mats)} %>%
    corrplot(
        method='number',
        # method='color',
        type='upper',
        is.corr=TRUE,
        col.lim=c(0, 1),
        number.digits=3
    )
dev.off()

################################################
# Outlier detection from PCA space
################################################
sapply(
    names(all_sv_counts) %>% grep('XDP_', ., invert=TRUE, value=TRUE),
    function(sample_set){
        # keep relevant samples
        counts <- 
            all_sv_counts[[sample_set]]
        metadata <- 
            sample_metadata %>%
            # mutate(Condition=glue('{Treatment}.{Genotype}')) %>% 
            filter(SampleID %in% colnames(counts))
        # Detect outlier samples by distance in PCA space
        get_outlier_from_pca(
            counts=counts,
            metadata=metadata,
            sampleColumn='SampleID',
            groupColumn='Condition',
            ntop=500  # use top 500 most variable genes for PCA
        )
    },
    simplify=FALSE
) %>% 
bind_rows(.id='SVCountSet') %>%
write_tsv(file.path(OUTPUT_DIR, 'results', 'Sample.Outlier.Status.tsv'))

################################################
# Check if outliers overlap across sample sets
################################################
all_outliers_df <- 
    file.path(OUTPUT_DIR, 'results', 'Sample.Outlier.Status.tsv') %>%
    read_tsv(show_col_types=FALSE) %>% 
    pivot_wider(
        names_from=PC, 
        values_from=c(PCvalue, centroid, sd),
        names_sep='.'
    ) 
# Count outliers per comparison
all_outliers_df %>% 
    filter(isOutlier) %>%
    dplyr::count(SVCountSet, Condition) %>% 
    arrange(n) %>%
    print(n=Inf)
# Count outliers by PatientID
all_outliers_df %>% 
    filter(isOutlier) %>% 
    dplyr::select(-c(isOutlier, ends_with(c('.PC1', '.PC2')))) %>%
    left_join(
        sample_metadata %>%
        dplyr::select(SampleID, PatientID),
        by='SampleID'
    ) %>%
    add_count(SampleID) %>%
    arrange(desc(n), PatientID, SampleID, Condition) %>% 
    relocate(PatientID, SampleID, Condition) %>%
    {.}

