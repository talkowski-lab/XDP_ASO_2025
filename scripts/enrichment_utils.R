################################################
# Dependencies
################################################
library(tidyverse)
library(broom)
library(furrr)
library(future)
source('./locations.R')

################################################
# Load + merge all gene-functional group annotations
################################################
load_syngo_annotations <- function(){
    file.path(
        ANNOTATION_DATA_DIR,
        "SynGO_ontologies.txt",
    ) %>% 
    read_tsv(
        colnames=
            c(
                "GO_term_id",
                "GO_domain",
                "GO_term_name",
                "GO_parent_term_ID",
                "hgnc_id",
                "hgnc_symbol"
            )
    ) %>%
    mutate(id=paste(GO_term_id,
                    GO_domain,
                    GO_term_name,
                    sep="_"
           )
    ) %>%
    dplyr::select(id, hgnc_symbol) %>%
    setNames(c('term', 'genes')) %>%
    filter(genes != "") %>%
    add_column(Ontology='SYNGO') 
}

load_genelist_annotations <- function(){
    sprintf("%s/genelist_mef2c.txt", ANNOTATION_DATA_DIR) %>%
    read_tsv() %>%
    setNames(c('term', 'n_genes', 'genes')) %>%
    add_column(Ontology='GENELISTS')
}

load_hpo_annotations <- function(){
    sprintf("%s/HPO_09_12_2020.txt", ANNOTATION_DATA_DIR) %>%
    read_tsv() %>%
    setNames(c('to_rm', 'term', 'n_genes', 'genes')) %>%
    dplyr::select(!to_rm) %>%
    add_column(Ontology='HPO')
}

load_ndd_annotations <- function(){
    sprintf('%s/NDG.txt', ANNOTATION_DATA_DIR) %>%
    read_tsv(col_names=FALSE) %>%
    rename('NDG'=X1) %>%
    pivot_longer(
        NDG,
        names_to='Ontology',
        values_to='genes'
    ) %>%
    group_by(Ontology) %>%
    summarize(genes=paste(genes, collapse=';')) %>%
    add_column(term='NDG List')
}

load_TAF1_annotations <- function(){
    sprintf('%s/TAF1_targets.txt', ANNOTATION_DATA_DIR) %>%
    read_tsv(col_names=FALSE) %>%
    rename('TAF1'=X1) %>%
    pivot_longer(
        TAF1,
        names_to='Ontology',
        values_to='genes'
    ) %>%
    group_by(Ontology) %>%
    summarize(genes=paste(genes, collapse=';')) %>%
    add_column(term='TAF1 Targets')
}

load_all_annotations <- function(outfile=COMBINED_ANNOTATION_FILE){
    if (file.exists(outfile)){
        ontologies <- read_tsv(outfile) 
    } else {
        sprintf('Creating new Gene Set Annotation File %s', outfile) %>% message()
        # Input is named list of ontology lists filenames (gmt format)
        # Final output is a dataframe with 4 columns 
        ##  Ontology: which file/ontology the term comes from
        ##  term: the name of the term/functional group
        ##  genes: , delimited string of gene names in the group 
        ##  n_genes: number of genes in the group
        ontologies <- 
            list(
                "C2.CP"="c2.cp.v2022.1.Hs.symbols.gmt",
                "C3.TFT"="c3.tft.v2022.1.Hs.symbols.gmt",
                 "C5.GOBP"="c5.go.bp.v2022.1.Hs.symbols.gmt",
                 "C5.GOCC"="c5.go.cc.v2022.1.Hs.symbols.gmt",
                 "C5.GOMF"="c5.go.mf.v2022.1.Hs.symbols.gmt",
                "C5.HPO"="c5.hpo.v2022.1.Hs.symbols.gmt"
            ) %>%
            lapply(
                function(x){
                    sprintf('%s/%s', ANNOTATION_DATA_DIR, x) %>% 
                    read.gmt() 
                }
            ) %>%
            t() %>%
            as_tibble() %>%
            pivot_longer(
                everything(),
                names_to='Ontology',
                values_to='term_list'
            ) %>%
            unnest_wider(term_list) %>%
            mutate(
                across(
                    !Ontology, 
                    function(genes){
                        case_when(
                            lengths(genes) != 0 ~ paste(sort(unlist(genes)), collapse=','),
                            TRUE ~ NA
                      )
                    }
                )
            ) %>%
            pivot_longer(
                !Ontology,
                names_to='term',
                values_to='genes'
            ) %>% 
            filter(!is.na(genes)) %>%
            bind_rows(
                load_syngo_annotations(),
                load_genelist_annotations(),
                load_hpo_annotations(),
                load_ndd_annotations(),
                load_TAF1_annotations()
            ) %>% 
            mutate(genes=gsub(';', ',', genes)) %>%
            mutate(n_genes=str_count(genes, pattern=',') + 1) 
        write_tsv(ontologies, outfile)
    }
    # Convert ',' delimited string of gene names to list-column
    return(ontologies %>% mutate(genes=strsplit(genes, ',')))
}

format_all_annotations <- function(functional.annotations.df){
    # Conviniently rename columns, generate background sets for each database (Ontology)
    functional.annotations.df %>% 
    group_by(Ontology) %>%
    mutate(background.lst.genes=genes %>% unlist() %>% c() %>% unique() %>% list()) %>%
    ungroup() %>% 
    dplyr::rename(
        'functional.name'=term,
        'functional.lst.genes'=genes,
        'functional.n.genes'=n_genes
    ) %>%
    dplyr::select(-c(functional.n.genes)) %>% 
    rowwise() %>%
    mutate(functional.lst.genes=list(unique(functional.lst.genes))) %>% 
    ungroup()
}

################################################
# Functions
################################################
load_symbol_mapping <- function(){
    # Get annotation of genes and genomic positions
    GENE_ID_MAPPING_FILE %>% 
    rtracklayer::import() %>%
    as_tibble() %>%
    filter(type == 'gene')
}

check_mapping <- function(
    df, 
    gene_list,
    silent){
    # Background set and how many EnsemblID are matched to symbols 
    size_background <- 
        df %>% 
        pull(symbol) %>% unique() %>% length()
    size_background_mapped <- 
        df %>% 
        filter(!is.na(symbol)) %>% 
        pull(symbol) %>% unique() %>% length()
    # Gene set of interest set and how many EnsemblID are matched to symbols 
    size_hits <- 
        df %>% 
        filter(is_hit == 1) %>% 
        pull(symbol) %>% unique() %>% length()
    size_hits_mapped <- 
        df %>% 
        filter(is_hit == 1) %>% 
        filter(!is.na(symbol)) %>% 
        pull(symbol) %>% unique() %>% length()
    # Gene in the Functional Group testing againts and how many EnsemblID are matched to symbols 
    size_term <- length(gene_list)
    size_term_mapped <- 
        df %>% 
        filter(symbol %in% gene_list) %>% 
        pull(symbol) %>% unique() %>% length()
    if (!silent){
        # Backgroun set size
        sprintf('Background size: %s', size_background)
        sprintf(
            '%s of %s background EnsemblIDs mapped to symbols',
            size_background_mapped,
            size_background
        ) %>% 
        message()
        # Hit list size
        sprintf('Hit list size: %s', size_hits)
        sprintf(
            '%s of %s hits EnsemblIDs mapped to symbols',
            size_hits_mapped,
            size_hits
        ) %>% 
        message()
        sprintf('Term list size: %s', size_term)
        sprintf('%s of %s term EnsemblIDs mapped to symbols',
                size_term_mapped,
                size_term
        ) %>% 
        message()
    }
    return(size_term_mapped)
}
# Do a single fisher test
do_fisher_test <- function(
    query_genes,
    term_genes,
    background_genes){
    # For statistical reference on enrichment score see
    # https://metascape.org/blog/?p=122
    # message(lengths(c(query_genes, term_genes, background_genes)))
    # Check rate of EnsemblIDs mapping to Entrez IDs from the input 
    # term_ids_mapped <- check_mapping(input_data, term_genes, silent=silent)
    # List of term genes that are also in the background
    observed_term_genes  <- intersect(term_genes,  background_genes) 
    observed_query_genes <- intersect(query_genes, background_genes) 
    non_background_genes <- union(observed_query_genes, observed_term_genes)
    # Gene IS a hit AND IS in the functional group
    is_term_is_hit <- length(intersect(observed_query_genes, observed_term_genes))
    # Gene IS NOT a hit AND IS in the functional group
    is_term_not_hit <- length(observed_term_genes) - is_term_is_hit
    # Gene IS a hit AND IS NOT in the functional group
    not_term_is_hit <- length(observed_query_genes) - is_term_is_hit
    # Gene IS NOT a hit AND IS NOT in the functional group
    # not_term_not_hit <- length(background_genes) - length(observed_query_genes) - is_term_not_hit 
    not_term_not_hit <- length(background_genes) - length(non_background_genes)
    # Contigency table for fisher test
    contingency <- 
        matrix(
            c(
                is_term_is_hit,  is_term_not_hit,
                not_term_is_hit, not_term_not_hit
            ), 
            byrow=TRUE,
            ncol=2
        )
    # Fisher test results, will unpack later
    fisher_result <- fisher.test(contingency, alternative='greater')
    # Calculate enrichment Z score (from Serkan's original code)
    # observed_term_genes <- intersect(background, term_genes) 
    n_background <- length(background_genes)
    n_query <- length(query_genes)
    expected <- length(observed_term_genes) * n_query / n_background
    variance_term_1 <- (n_background - length(observed_term_genes)) / n_background
    variance_term_2 <- (n_background - n_query) / (n_background - 1)
    std_dev <- sqrt(expected * variance_term_1 * variance_term_2)
    enrichment_zscore <- (is_term_is_hit - expected) / std_dev
    # return data as a single row tibble 
    list(
        is_term_is_hit=is_term_is_hit,
        is_term_not_hit=is_term_not_hit,
        not_term_is_hit=not_term_is_hit,
        not_term_not_hit=not_term_not_hit,
        # observed_term_genes=length(observed_term_genes),
        term_and_hit_genes=
            paste(
                intersect(
                    query_genes,
                    observed_term_genes
                ), 
                collapse=','
            ),
         enrichment_zscore=enrichment_zscore
    ) %>%
    as_tibble() %>%
    # unpack fisher test object into the dataframe neatly
    bind_cols(tidy(fisher_result))
}
# convienience wrapper
compute_fisher_test <- function(
    functional.lst.genes,
    background.lst.genes,
    query.lst.genes,
    ...){
    do_fisher_test(
        query_genes=query.lst.genes,
        term_genes=functional.lst.genes,
        background_genes=background.lst.genes
    )
}
# Run fisher tests for all query gene sets against all functional gene groups
compute_all_fisher_tests <- function(
    functional.annotations.df,
    query.df,
    nrows=Inf,
    ...){
    functional.annotations.df %>% 
    head(nrows) %>% 
    cross_join(query.df) %>% 
    mutate( 
        fisher.results=
            pmap(
                .l=.,
                .f=compute_fisher_test,
                .progress=TRUE
            )
    ) %>% 
    unnest(fisher.results) %>% 
    rowwise() %>% 
    mutate(
        across(
            ends_with('.lst.genes'), 
            .names="{sub('.lst.genes', '.n.genes', col)}",
            length
        )
    ) %>% 
    ungroup() %>% 
    dplyr::select(
        -c(
            ends_with('.lst.genes'),
            method,
            alternative
        )
    ) %>%
    group_by(
        query.name,
        Ontology
    ) %>% 
    mutate(fdr=p.adjust(p.value, method='BH')) %>%
    arrange(fdr)
}

################################################
# Same logic, different columns names
################################################
format_functional_groups <- function(functional.annotations.df){
    functional.annotations.df
    dplyr::rename(
        'term_name'=term,
        'term_genes'=genes
    ) %>% 
    rowwise() %>% 
    mutate(
        term_genes=
            term_genes %>%
            strsplit(',') %>%
            unlist() %>% 
            unique() %>%
            list()
    ) %>% 
    ungroup() %>% 
    group_by(Ontology) %>%
    mutate(db_genes=term_genes %>% unlist() %>% c() %>% unique() %>% list()) %>%
    ungroup()
}

run_enrichments <- function(
    term_df,
    input_df,
    silent=FALSE){
    # Get all genes present in each Ontology, used for subsetting query/bg
    # to avoid deflating significance of enrichment tests i.e.
    # avoid counting query/bg genes that we cannot detect because they
    # are not present anywhere in the specified Ontology DB
    term_df %>%
    cross_join(input_df) %>%
    # head(30) %>% 
    # Only consider genes present in at least 1 functional term in the entire Ontology
    rowwise() %>% 
    mutate(
        across(
            c(
                background_genes,
                hit_genes
            ),
            ~ list(intersect(.x, db_genes))
        )
    ) %>% 
    ungroup() %>% 
    # Compute all enrichments in parallel
    mutate(
        enrichment.results=
            # future_pmap(
            pmap(
                .l=.,
                .f=do_fisher_test,
                .progress=TRUE
            )
    ) %>%
    unnest(enrichment.results) %>%
    # Do multiple testing correction per ontology
    group_by(Ontology) %>%
    mutate(padj=p.adjust(p.value, method='BH')) %>%
    ungroup() %>%
    arrange(padj, p.value) %>% 
    # Clean upt columns
    rowwise() %>% 
    mutate(
        across(
            c(hit_genes, term_genes, db_genes, background_genes),
            length
        )
    ) %>% 
    ungroup() %>% 
    dplyr::select(
        -c(
            method,
            alternative,
            n_genes
        )
    ) %>%
    dplyr::rename('odds.ratio'=estimate)
}

load_enrichment_data <- function(file_pattern='Enrichments-RescuedGenes-.*.tsv'){
    ENRICHMENTS_DIR %>% 
    list.files(
        pattern=file_pattern,
        full.names=TRUE
    ) %>%
    tibble(filepath=.) %>%
    mutate(
    fileinfo=
        filepath %>%
        basename() %>% 
        str_remove('.tsv')
    ) %>%
    separate_wider_delim(
        fileinfo,
        delim='-',
        names=
            c(
                NA,
                NA,
                NA,
                'Rescue.Criteria',
                'OverlapList'
            )
    ) %>% 
    rowwise() %>% 
    mutate(enrichment.results=list(read_tsv(filepath))) %>% 
    ungroup() %>% 
    unnest(enrichment.results) %>% 
    dplyr::select(-c(filepath))
}

