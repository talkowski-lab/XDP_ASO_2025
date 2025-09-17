################################################
# Dependencies
################################################
library(ggplot2)
library(magrittr)
library(tidyverse)

################################################
# Utils
################################################
add_ggtheme <- function(){
    # Create uniform theme for figures
    # From Dadi Gao
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank(), 
          panel.background=element_blank(),
          axis.line=element_line(color="black",
                                 linewidth=1/2.13
          ),
          axis.ticks=element_line(color="black",
                                  linewidth=1/2.13
          ),
          axis.ticks.length=unit(3, "pt"),
          axis.title=element_text(family="sans",
                                  face="bold",
                                  color="black",
                                  size=10
          ),
          axis.text=element_text(family="sans",
                                 face="bold",
                                 color="black",
                                 size=8
          )
    )
}

################################################
# Various Plots
################################################
plot_ExprBoxPlot <- function(
    count_matrix,
    sample_metadata,
    gene_metadata=gene_metadata,
    x_var='geno',
    gene_names=c('TAF1'),
    ncol=2){
    # Reshape data for plotting
    plot_df <- 
        count_matrix %>%
        filter(gene_name %in% gene_names) %>%
        droplevels()
    # Make box plot per gene
    g <- 
        ggplot(
            plot_df,
            aes(
                x=.data[[x_var]],
                y=counts
            )
        ) + 
        geom_boxplot() +
        stat_compare_means(
            aes(group=.data[[x_var]]),
            method='t.test',
            comparisons=
                plot_df[[x_var]] %>% 
                unique() %>%
                as.character() %>%
                combn(2, simplify=FALSE)
        ) +
        facet_wrap(
            ~ gene_name,
            ncol=ncol
        ) +
        theme_bw() +
        add_ggtheme()
    return(g)
}

plot_SimpleExprBoxPlot <- function(
    count_matrix,
    x_var='geno',
    ncol=2){
    # Reshape data for plotting
    plot_df <- 
        count_matrix %>%
        droplevels()
    # Make box plot per gene
    g <- 
        ggplot(
            plot_df,
            aes(
                x=.data[[x_var]],
                y=counts
            )
        ) + 
        geom_boxplot() +
        facet_wrap(
            ~ gene_name,
            ncol=ncol
        ) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
        add_ggtheme()
    return(g)
}

pca_wrapper <- function(
    count_matrix,
    sample_metadata,
    x_var,
    outfile,
    height=5,
    width=5,
    ...){
    pca_plot <- 
        sample_metadata %>%
        filter(SampleID %in% colnames(count_matrix)) %>%
        mutate(labels=SampleID) %>%
        # mutate(labels=paste(PatientID, Clone) %>% gsub(' ', '_', .)) %>%
        # mutate(Condition=paste(Treatment, get(x_var), sep='.')) %>%
        # mutate(Condition=as.factor(Condition)) %>%
        rnaPCA(
            count_matrix,
            .,
            groupColumn=x_var,
            ...
        )
    ggsave(
        outfile, 
        pca_plot + scale_color_okabeito() + add_ggtheme(),
        height=height, width=width, units='in'
    )
}

rnaPCA <- function (
    object,
    sampleTable,
    groupColumn,
    sampleColumn='SampleID',
    labelColumn='labels',
    show_centroids=FALSE, show_sds=FALSE,
    ntop=500, PC1=1, PC2=2,
    title_str=NA,
    alpha=0.5, pointSize=3, labelSize=2,
    xlim=c(-15, 15), ylim=c(-15, 15)){
    # object is the raw counts, should be a dataframe
    # sampleTable is the metadata, should be a tibble
    # make sure metadata and counts have matching samples
    samples <- sampleTable[[sampleColumn]]
    mat <- as.matrix(object[, samples]) 
    # filter genes to only keep the ntop most variable genes
    genes_to_keep <- 
        mat %>%
        rowVars() %>%
        order(decreasing=TRUE) %>%
        {.[seq_len(min(ntop, length(.)))]}
    pca_results <- prcomp(t(mat[genes_to_keep, samples]))
    percent_var <- 100 * (pca_results$sdev^2 / sum(pca_results$sdev^2))
    # create dataframe with intgroups and PC loadings for plotting
    plot_df <- 
        pca_results$x %>%
        as.data.frame() %>%
        rownames_to_column(var=sampleColumn) %>%
        as_tibble() %>%
        left_join(
            sampleTable %>%
            dplyr::select(
                all_of(
                    c(
                        sampleColumn,
                        groupColumn,
                        labelColumn
                    )
                )
            ),
            by=sampleColumn
        )
    # Calculate centroids per group
    centroids <- 
        plot_df %>%
        dplyr::select(!all_of(labelColumn)) %>%
        pivot_longer(
            !all_of(c(groupColumn, sampleColumn)),
            names_to='PC',
            values_to='value'
        ) %>%
        group_by(pick(all_of(c(groupColumn, 'PC')))) %>%
        summarize(
            centroid=mean(value),
            sd=sd(value) * 2
        ) %>%
        ungroup()
    # Set plot title
    title_str <- 
        sprintf(
            '%s (Top %s genes) (n_samples=%s) (n_genes=%s)',
            ifelse(is.na(title_str), 'PCA of Samples', title_str),
            ntop,
            ncol(object),
            nrow(object)
        )
    # Plot PCA
    p <- 
        ggplot(
            plot_df, 
            aes(
                x=.data[[paste0("PC", PC1)]], 
                y=.data[[paste0("PC", PC2)]], 
                color=.data[[groupColumn]]
            )
        ) + 
        # Plot individual samples in PCA space
        geom_point(
           size=pointSize, 
           alpha=alpha
        ) + 
        # Plot individual sample names
        geom_text_repel(
           aes(label=.data[[labelColumn]]),
           size=labelSize,
           max.overlaps=nrow(sampleTable)
        ) +
        # Fix xlim/ylin scales to more easily compare figures
        coord_cartesian(xlim=xlim, ylim=ylim) +
        # Add title and variance explained by plotted PCs
        labs(
           title=title_str,
           x=sprintf('PC%s (%.2f%% variance)', PC1, percent_var[PC1]),
           y=sprintf('PC%s (%.2f%% variance)', PC2, percent_var[PC2])
        ) +
        geom_hline(yintercept=0) +
        geom_vline(xintercept=0) +
        # theme_bw() +
        theme(
            plot.title=element_text(size=10),
            legend.position='top'
        )
    # Draw centroids per group
    if (show_centroids) {
        centroid.df <- 
            centroids %>% 
            filter(PC %in% 
                c(
                    paste0("PC", PC1), 
                    paste0("PC", PC2)
                )
            ) %>%
            dplyr::select(
                all_of(
                    c(
                        groupColumn,
                        'PC',
                        'centroid'
                    )
                )
            ) %>%
            pivot_wider(
                names_from=PC,
                values_from=centroid
            )
        p <- 
            p +
            geom_point(
                data=centroid.df,
                aes(
                    x=.data[[paste0("PC", PC1)]], 
                    y=.data[[paste0("PC", PC2)]]
                ),
                shape=5,
                stroke=2,
                size=pointSize * 2
            )
    }
    # Draw circles around centrods w radius == 2 * sd per level
    if (show_sds) {
        centroid.df <- 
            centroids %>% 
            filter(PC %in% 
                c(
                    paste0("PC", PC1), 
                    paste0("PC", PC2)
                )
            ) %>%
            pivot_longer(
                c(centroid, sd),
                names_to='metric',
                values_to='value'
            ) %>%
            mutate(metric=paste(metric, PC, sep='.')) %>%
            dplyr::select(!PC) %>%
            pivot_wider(
                names_from=metric,
                values_from=value
            )
        p <- 
            p +
            geom_ellipse(
                data=centroid.df,
                aes(x0=centroid.PC1, 
                  y0=centroid.PC2,
                  a=sd.PC1,
                  b=sd.PC2,
                  angle=0,
                  color=.data[[groupColumn]]
                ),
                inherit.aes=FALSE,
                alpha=0.1
            )
    }
    return(p)
}

plot_BootstrapDists <- function(
    bootstrap_results, 
    strip.size=20,
    label.size=3,
    label.y=8000,
    nrow=2){
    plot_df <- 
        bootstrap_results 
    taf1_data <- 
        plot_df %>%
        filter(EnsemblID == TAF1_ENSEMBLID) %>%
        mutate(taf1_label=paste('TAF1', 
                                round(fraction_of_bootstraps, 3), 
                                sep=' '
               )
        )
    # Plot the distributions per expression pattern
    g <- ggplot(plot_df, 
                aes(x=fraction_of_bootstraps)
         ) + 
         geom_histogram(aes(x=fraction_of_bootstraps, 
                            y=after_stat(count)
                            # fill=expression_pattern
                        ), 
                        show.legend=FALSE,
                        alpha=0.8
         ) +
         geom_vline(aes(xintercept=fraction_of_bootstraps),
                    data=taf1_data,
                    size=0.75,
                    linetype='dashed'
         ) +
         geom_text(aes(x=0.1,
                       y=label.y,
                       label=taf1_label
                   ),
                   angle=0,
                   data=taf1_data,
                   hjust=-0.25,
                   size=label.size
         ) +
         facet_wrap(~ expression_pattern,
                    nrow=nrow,
                    scales='fixed'
         ) +
         labs(title=sprintf('Expression Patterns in Bootstraped Data (n=%s)',
                            max(bootstrap_results$bootstraps)
              ),
              y='Number of Genes',
              x='Fraction of All Permutations'
         ) +
         theme_bw() + 
         theme(legend.position='top',
               strip.text=element_text(size=strip.size,
                                       face='bold'
               ),
         ) + add_ggtheme()
    return(g)
}

merge_BootstrapSubSamples <- function(file_list){
    combined_results <- 
        lapply(1:length(file_list),
               function(idx) {
                   read_tsv(file_list[[idx]], 
                            progress=FALSE,
                            show_col_types=FALSE
                   ) %>% 
                   add_column(index=idx)
               }
        ) %>% 
        bind_rows() %>%
        dplyr::select(
            starts_with('frac'), 
            'EnsemblID',
            'index',
            'bootstraps'
        ) %>%
        pivot_longer(starts_with('frac_'),
                     names_to='expression_pattern',
                     values_to='fraction_of_bootstraps',
                     names_prefix='frac_'
        ) %>%
        group_by(EnsemblID,
                 expression_pattern
        ) %>%
        mutate(raw_value=fraction_of_bootstraps * bootstraps) %>%
        summarize(bootstraps=sum(bootstraps),
                  raw_value=sum(raw_value)
        ) %>%
        mutate(fraction_of_bootstraps=raw_value / bootstraps) %>%
        ungroup() 
    return(combined_results)
}

plot_BootstrapSubSamples <- function(
    n_subsamples=5,
    chunks_per_subsample=5,
    seed=9,
    CHUNK_SIZE=200){
    # Each file contains the summarized info for $CHUNK_SIZE permutations, so
    # the minimum possible n_samples is $CHUNK_SIZE
    set.seed(seed)
    plot_df <- 
        # Get list of all bootstrap results (200 bootstraps summarized per file)
        file.path(BOOTSTRAP_RESULTS_DIR, "NO.CON.NO.CombinedXDP") %>% 
        list.files(
            pattern='summary-.*NO.CON.NO.CombinedXDP.tsv\\.*',
            full.names=TRUE
        ) %>%
        # Randomly sample chunks of results based on arguments
        sample(n_subsamples * chunks_per_subsample) %>%
        # Combine sets of files into individual subsamples
        {
            mapply(
                function(i_start, chunk_num){
                    i_end <- i_start + chunks_per_subsample - 1
                    df <- 
                        merge_BootstrapSubSamples(.[i_start:i_end]) %>%
                        add_column(subset=chunk_num)
                    return(df)
                },
                seq.int(1, length(.), by=chunks_per_subsample),
                seq.int(1, length(.) / chunks_per_subsample, 1),
                SIMPLIFY=FALSE
            )
        } %>%
        bind_rows() %>%
        mutate(
            subset=paste0('subset-', subset),
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
    # Get data for drawing lines for TAF1
    taf1_data <- 
        plot_df %>%
        filter(EnsemblID == TAF1_ENSEMBLID) %>%
        mutate(taf1_label=sprintf('TAF1: %.3f', fraction_of_bootstraps))
    # Make the plot
    g <- 
        ggplot(
            plot_df,
            aes(x=fraction_of_bootstraps)
        ) +
        geom_density(
           aes(
               fraction_of_bootstraps, 
               color=subset,
               fill=subset
           ), 
           alpha=0.6
        ) +
        geom_vline(
           data=taf1_data,
           aes(
               xintercept=fraction_of_bootstraps,
               color=subset
           ),
           size=0.6,
           linetype='dashed'
        ) +
        geom_text(
           data=taf1_data,
           aes(
               x=0.1, 
               y=
                   rep(
                       seq(0.75, 1.25, length.out=n_subsamples) * 10, 
                       length(unique(expression_pattern))
                   ),
               label=taf1_label,
               color=subset,
           ),
           hjust=-0.25,
           size=2.5
        ) +
        facet_wrap(
           ~ expression_pattern, 
           nrow=2,
           # scales='free_y'
           scales='fixed'
        ) +
        labs(
           title=
               sprintf(
                   'Bootstrap Subsampling (bootstraps per subset=%s)',
                   chunks_per_subsample * CHUNK_SIZE
               ),
           x='Fraction of All Permutations',
           y='Density'
        ) +
        theme_bw() + 
        theme(
            legend.position='top',
            # legend.key.size=unit(0.75, 'in'),
            legend.text=element_text(size=8),
            legend.title=element_blank(),
            strip.text=element_text(size=8, face='bold')
        ) + 
        add_ggtheme()
    return(g)
}

