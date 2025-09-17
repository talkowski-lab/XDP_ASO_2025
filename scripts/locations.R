################################################
# Input Data
################################################
SAMPLE_METADATA_FILE     <- file.path(BASE_DIR, 'NSC.Sample.Metadata.tsv')
READ_COUNTS_FILE         <- file.path(BASE_DIR, 'raw_counts.txt')
UNIQ_READS_FILE          <- file.path(BASE_DIR, 'uniq_reads.txt')

################################################
# Annotation Data
################################################
ANNOTATIONS_DIR          <- file.path(BASE_DIR, 'annotations')
TAF1_TARGETS_FILE        <- file.path(ANNOTATIONS_DIR, 'TAF1.targets.txt')
NDD_GENES_FILE           <- file.path(ANNOTATIONS_DIR, 'NDG.genes.txt')
GENE_METADATA_FILE       <- file.path(ANNOTATIONS_DIR, 'gene.metadata.tsv')
GENE_ID_MAPPING_FILE     <- file.path(ANNOTATIONS_DIR, 'Symbol.EnsemblID.mapping.tsv')
MSIG_DB_DIR              <- file.path(ANNOTATIONS_DIR, "MSigDB7.4")
COMBINED_ANNOTATION_FILE <- file.path(ANNOTATIONS_DIR, 'combined_annotation_terms.tsv')

################################################
# GLM Results
################################################
GLM_RESULTS_DIR            <- file.path(BASE_DIR, 'rm2Outliers_GLMResults')
XDP_DEG_RESULTS_FILE       <- file.path(GLM_RESULTS_DIR, 'NO.CON+NO.XDP-DEGResults.tsv')
BOOTSTRAP_RESULTS_DIR      <- file.path(BASE_DIR, 'BootstrapResults')
BOOTSTRAP_RESULTS_FILE     <- file.path(BOOTSTRAP_RESULTS_DIR, 'bootstrap_results.tsv')

################################################
# Paper Results
################################################
TABLES_DIR                 <- file.path(BASE_DIR, '2025-06-22_paper.tables')
FIGURES_DIR                <- file.path(BASE_DIR, '2025-09-09_figures')
COMPUTED_DATA_DIR          <- file.path(FIGURES_DIR, 'data')
CORRECTED_COUNTS_FILE      <- file.path(COMPUTED_DATA_DIR, 'all.corrected.counts.tsv')
RESCUED_GENES_RESULTS_FILE <- file.path(COMPUTED_DATA_DIR, 'RescueData-AllTreatments-AllCriteria.tsv')
XDP_DEG_RESULTS_FILE       <- file.path(COMPUTED_DATA_DIR, 'NO.CON+NO.XDP-LFC+Var.tsv')
ENRICHMENTS_FDR_FILE       <- file.path(COMPUTED_DATA_DIR, 'Enrichments-NO.CON.vs.NO.XDP-FDR.tsv')
ENRICHMENTS_NOMINAL_FILE   <- file.path(COMPUTED_DATA_DIR, 'Enrichments-NO.CON.vs.NO.XDP-Nominal.tsv')

################################################
# Misc
################################################
# Unused results dirs
# ENRICHMENTS_DIR            <- file.path(BASE_DIR, '2025-06-18_enrichment.results')
# OUTLIER_ANALYSIS_DIR       <- file.path(BASE_DIR, 'rm2Outliers_OutlierAnalysis')
# MISC_RESULTS_DIR         <- file.path(BASE_DIR, 'MiscResults')
# L3L5_DIR                 <- file.path(MISC_RESULTS_DIR, '5vs3levels_GLMResults')
# PAIRED_TESTING_DIR       <- file.path(MISC_RESULTS_DIR, 'rm2Outliers_PairedAnalysis')
# results dirs
# SHINY_APPS_DIR           <- file.path(BASE_DIR, 'ShinyApps')
# GLM_PAIRED_RESULTS_DIR   <- file.path(BASE_DIR, 'rm2Outliers_GLMPairComparison')
# ENRICHMENTS_DIR          <- file.path(BASE_DIR, 'rm2Outliers_Enrichments')
# PAPER_FIGURES_DIR        <- file.path(BASE_DIR, 'rm2Outliers_PaperFigures')

