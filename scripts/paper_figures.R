################################################
# Dependencies
################################################
library(here)
# here::i_am('scripts/paper_figures.R')
here::i_am('2025-09-16_code/scripts/paper_figures.R')
BASE_DIR <- here()
# SCRIPTS_DIR <- here('scripts')
SCRIPTS_DIR <- here('2025-09-16_code/scripts')
source(file.path(SCRIPTS_DIR, 'locations.R'))
source(file.path(SCRIPTS_DIR, 'constants.R'))
source(file.path(SCRIPTS_DIR, 'data_utils.R'))
source(file.path(SCRIPTS_DIR, 'enrichment_utils.R'))
source(file.path(SCRIPTS_DIR, 'glm_utils.R'))
source(file.path(SCRIPTS_DIR, 'plot_utils.R'))
library(ggplot2)
library(ggpubr)
library(ggrepel)
# library(ggforce)
library(corrplot)
library(see)  # for  Okabe-Ito color scheme
# library(RColorBrewer)
library(glue)
library(magrittr)
OUTPUT_DIR <- FIGURES_DIR

################################################
# Load all Data that is reference multiple times
################################################
# Load sample metadata
sample_metadata <- load_meta()
# EnsemblID to gene name mappings
gene_metadata <- load_gtf()
# XDP vs CON Bootstrapping results
bootstrap_results <- make_bootstrap_results()
# # Combine all adjusted counts 
all_sv_counts <- load_corrected_counts_list()
sapply(all_sv_counts, dim) %>% t() %T>% {colnames(.) <- c('n.genes', 'n.samples')}
counts_df <- load_corrected_counts_table(gene_metadata=gene_metadata)

################################################
# PCA of Samples from Top genes
################################################
# Color naive and non_edit samples together
pca_wrapper(
    all_sv_counts[['cXDP']],
    sample_metadata,
    x_var='Genotype',
    outfile=file.path(OUTPUT_DIR, 'PCA-Genotype-NO.CON+NO.XDP.pdf'),
    title_str='Sample PCA'
)
# Color naive and non_edit samples separately
pca_wrapper(
    all_sv_counts[['cXDP']],
    sample_metadata,
    x_var='geno',
    outfile=file.path(OUTPUT_DIR, 'PCA-geno-NO.CON+NO.XDP.pdf'),
    title_str='Sample PCA'
)
# NO.CON vs XDP unedited only 
pca_wrapper(
    all_sv_counts[['XDP_non_edit']],
    sample_metadata,
    x_var='geno',
    outfile=file.path(OUTPUT_DIR, 'PCA-geno-NO.CON+NO.XDP_non_edit.pdf'),
    title_str='Sample PCA'
)
# NO.CON vs XDP unexposed
sample_metadata %>%
    mutate(geno=case_when(geno == 'XDP' ~ 'XDP_unexposed', TRUE ~ geno)) %>%
    pca_wrapper(
        all_sv_counts[['XDP_naive']],
        .,
        x_var='geno',
        outfile=file.path(OUTPUT_DIR, 'PCA-geno-NO.CON+NO.XDP_unexposed.pdf'),
        title_str='Sample PCA'
    )

################################################
# Box Plots of Adjusted Counts
################################################
BOX_DIR <- file.path(OUTPUT_DIR, 'BoxPlots')
# XDP_naive vs control
# counts_df %>% dplyr::count(origin)
counts_df %>%
    filter(origin == 'XDP_naive') %>% 
    plot_ExprBoxPlot(
        sample_metadata %>% mutate(geno=gsub('XDP$', "XDP_naive", geno)),
        x_var='geno',
        gene_names=c('TAF1'),
        gene_metadata=gene_metadata
    ) %>% 
    ggsave(
        file.path(BOX_DIR, 'TAF1-CON+XDP_naive-SVAdjustedCountBoxPlot.pdf'),
        ., width=6, height=6, units='in'
    )
# XDP_non_edit vs control
counts_df %>%
    filter(origin == 'XDP_non_edit') %>% 
    plot_ExprBoxPlot(
        sample_metadata,
        x_var='geno',
        gene_names=c('TAF1'),
        gene_metadata=gene_metadata
    ) %>% 
    ggsave(
        file.path(BOX_DIR, 'TAF1-CON+XDP_non_edit-SVAdjustedCountBoxPlot.pdf'),
        ., width=6, height=6, units='in'
    )
# CombinedXDP 3 colors
counts_df %>%
    filter(origin == 'cXDP') %>% 
    plot_ExprBoxPlot(
        sample_metadata %>% mutate(geno=gsub('XDP$', "XDP_naive", geno)),
        x_var='geno',
        gene_names=c('TAF1'),
        gene_metadata=gene_metadata
    ) %>% 
    ggsave(
        file.path(BOX_DIR, 'TAF1-CON+XDP-geno-SVAdjustedCountBoxPlot.pdf'),
        ., width=6, height=6, units='in'
    )
# CombinedXDP 2 colors
counts_df %>% 
    filter(origin == 'cXDP') %>% 
    # filter(Category %in% c('Untreated.CON', 'Untreated.XDP'))
    plot_ExprBoxPlot(
        sample_metadata,
        x_var='Genotype',
        gene_names=c('TAF1'),
        gene_metadata=gene_metadata
    ) %>% 
    ggsave(
        file.path(BOX_DIR, 'TAF1-CON+XDP-Genotype-SVAdjustedCountBoxPlot.pdf'),
        ., width=6, height=6, units='in'
    )

################################################
# Scatter plot of XDP_non_edit vs XDP_naive LFCs 
################################################
# Load DEG resutls
XDP_naive_results <- 
    file.path(GLM_RESULTS_DIR, 'NO.CON+NO.XDP_naive-DEGResults.tsv') %>%
    read_tsv(show_col_types=FALSE) %>%
    dplyr::select(!Contrast)
XDP_non_edit_results <- 
    file.path(GLM_RESULTS_DIR, 'NO.CON+NO.XDP_non_edit-DEGResults.tsv') %>%
    read_tsv(show_col_types=FALSE) %>%
    dplyr::select(!Contrast)
# Massage data for plotting
plot_df <- 
    XDP_naive_results %>%
    left_join(
        XDP_non_edit_results,
        suffix=c('.XDP_Unexposed', '.XDP_Unedited'),
        by='EnsemblID'
    ) %>%
    mutate(
        isDEG=
            case_when(
                pvalue.XDP_Unexposed < 0.05 & pvalue.XDP_Unedited < 0.05 ~ 'Both',
                pvalue.XDP_Unexposed < 0.05                              ~ 'Naive XDP',
                pvalue.XDP_Unedited  < 0.05                              ~ 'Unedited XDP',
                TRUE                                                     ~ 'Neither'

           )
    ) %>%
    filter(isDEG != 'Neither') %>%
    mutate(isDEG=factor(isDEG)) %>%
    filter(abs(lfc.XDP_Unedited) < 30 & abs(lfc.XDP_Unexposed) < 30)
# Plot 
g <- 
    ggscatter(
        plot_df,
        x='lfc.XDP_Unedited',
        y='lfc.XDP_Unexposed',
        color='isDEG',
        size=2,
        alpha=0.6
    ) +
    stat_cor(
        aes(color=isDEG),
        method='pearson',
        show.legend=FALSE
    ) + 
    geom_vline(xintercept=0, color='black') +
    geom_hline(yintercept=0, color='black') +
    geom_abline(intercept=0, slope=1, color='black') +
    xlim(-2, 2) + ylim(-2, 2) +
    labs(
        x='LFC of Unedited XDP DEGs',
        y='LFC of Naive XDP DEGs',
     ) + 
     scale_color_okabeito() +
     theme(legend.position='top') + 
     add_ggtheme()
ggsave(
    file.path(OUTPUT_DIR, 'Unedited.vs.Unexposed_ScatterPlot.pdf'),
    g, height=7, width=7, units='in'
)

################################################
# Bootstrapping figures 
################################################
# Plot only FDR_and_down distributions over all bootstraps
bootstrap_results %>%
    filter(expression_pattern == 'fdr_and_down') %>%
    plot_BootstrapDists() %>%
    ggsave(
        file.path(OUTPUT_DIR, 'FDR.Down_BootstrapDistributions.pdf'),
        ., width=6, height=6, units='in'
    )
# Plot only nominal_and_down distributions over all bootstraps
bootstrap_results %>%
    filter(expression_pattern == 'nominal_and_down') %>%
    plot_BootstrapDists(label.y=5000) %>%
    ggsave(
        file.path(OUTPUT_DIR, 'Nominal.Down_BootstrapDistributions.pdf'),
        ., width=6, height=6, units='in'
    )
# Plot expression pattern distributions over all bootstraps
bootstrap_results %>%
    plot_BootstrapDists(strip.size=8) %>%
    ggsave(
        file.path(OUTPUT_DIR, 'All_BootstrapDistributions.pdf'),
        ., width=9, height=5, units='in'
    )
# Plot distributions for multiple sub-samples for each pattern
plot_BootstrapSubSamples(
    n_subsamples=5,
    chunks_per_subsample=2,
) %>%
ggsave(
    sprintf("%s/All_BootstrapSubSampleDistributions.pdf", OUTPUT_DIR),
    ., width=9, height=5, units='in'
)

################################################
# Correlation Heatmap of Untreated XDP Adjusted Counts
################################################
XDP.Sample.List <- 
    sample_metadata %>% 
    filter(Genotype == 'XDP', Treatment == 'NO') %>% 
    pull(SampleID)
plot_df <- 
    # Pivot adjusted counts for easy joining
    sapply(
        names(all_sv_counts) %>% grep('(ASO|dSVA|XDP$)', ., value=TRUE),
        function(name){
            all_sv_counts[[name]] %>% 
            pivot_longer(
                !EnsemblID,
                names_to='SampleID',
                values_to=name
           ) 
       },
       simplify=FALSE
    ) %>%
    # Only keep samples found in every comparison i.e. Untreated XDP + CON samples
    purrr::reduce(
        dplyr::inner_join, 
        by=c('EnsemblID', 'SampleID')
    ) %>%
    # -Infs screw up the correlation
    filter(if_all(!c(EnsemblID, SampleID), ~ .x != -Inf)) %>%
    # Keep only Untreated XDP samples
    filter(SampleID %in% XDP.Sample.List) %>% 
    # Split into list of dataframes, 1 per sample
    split(f=.$SampleID) %>%
    # Calculate Treatment-Treatment correlations for each sample across genes
    sapply(
        function(df){
            df %>% 
            dplyr::select(!c(EnsemblID, SampleID)) %>%
            cor() %>%
            as.data.frame() %>%
            rownames_to_column(var='Treatment1') %>%
            pivot_longer(
                !Treatment1, 
                names_to='Treatment2',
                values_to='pearson'
            )
        },
        simplify=FALSE,
        USE.NAMES=TRUE
    ) %>%
    bind_rows(.id='SampleID') %>%
    # Get mean Treatment-Treatment count correlation across samples
    group_by(Treatment1, Treatment2) %>%
    summarize(mean_pearson=mean(pearson)) %>%
    # Remove redundant pairs for plotting heatmap
    rowwise() %>%
    mutate(uniq_pair=sort(c(Treatment1, Treatment2)) %>% paste(collapse='.')) %>%
    group_by(uniq_pair) %>%
    distinct(uniq_pair, .keep_all=T) %>%
    # Format correlation digits
    mutate(label=sprintf('%0.3f', mean_pearson))
# Plot heatmap
g <- 
    ggplot(
        plot_df,
        aes(
            Treatment1,
            Treatment2,
            fill=mean_pearson,
        )
     ) +
     geom_tile() +
     geom_text(aes(label=label), size=4, color='white') +
     theme(axis.line=element_blank()) +
     add_ggtheme()
ggsave(
    file.path(OUTPUT_DIR, 'AdjustedCountCorr_NO.XDP.pdf'),
    g, width=7, height=5
)

