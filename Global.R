# ═══════════════════════════════════════════════════════════════
# IMGT/StatIGTR - Configuration Globale
# ═══════════════════════════════════════════════════════════════

# VALIDATION DU RÉPERTOIRE TEMPORAIRE
temp_dir <- tempdir()
if (!dir.exists(temp_dir)) {
	dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
	cat("✓ Temporary directory recreated:", temp_dir, "\n")
} 
# ── PACKAGES ───────────────────────────────────────────────────
# Interface
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(plotly)
library(htmlwidgets)
library(webshot2)

# Données et graphiques 
library(tidyverse)
library(scales)
library(viridis)
library(ggrepel)
library(ggridges)
library(ggseqlogo)
library(patchwork)

# Statistiques
library(broom)
library(vegan)
library(pROC)
library(igraph)

# Arbres phylogénétiques
library(ape)
library(phangorn)
library(ggtree)

# Exportoui
library(openxlsx)
library(rmarkdown)

options(shiny.maxRequestSize = Inf)  

# ── SOURCES
# Couche 0 : fichiers autonomes
source("R/app_meta.R",                     encoding = "UTF-8")
source("R/utils_core.R",                   encoding = "UTF-8")
source("R/dataset_helpers.R",              encoding = "UTF-8")
source("R/ui_helpers.R",                   encoding = "UTF-8")
source("R/cooccurrence_analysis.R",        encoding = "UTF-8")
source("R/genomic_map_data.R",             encoding = "UTF-8")

# Couche 1 : dépendent de la couche 0 uniquement
source("R/fasta_parser.R",                 encoding = "UTF-8")
source("R/subgroup_display.R",             encoding = "UTF-8")
source("R/stat_descriptives.R",            encoding = "UTF-8")
source("R/diversity_functions.R",          encoding = "UTF-8")
source("R/logistic_validation.R",          encoding = "UTF-8")
source("R/visualization_advanced.R",       encoding = "UTF-8")
source("R/cooccurrence_network.R",         encoding = "UTF-8")
source("R/partial_genes.R",                encoding = "UTF-8")

# Couche 2 : dépend de la couche 1
source("R/conservation_functions.R",       encoding = "UTF-8")

# Modules (interface + serveur)
source("modules/mod_import.R",             encoding = "UTF-8")
source("modules/mod_alluvial_ridgeline.R", encoding = "UTF-8")
source("modules/mod_cooccurrence.R",       encoding = "UTF-8")
source("modules/gmap_module.R",            encoding = "UTF-8")
source("modules/permanova_module.R",       encoding = "UTF-8")
source("modules/phylogeny_module.R",       encoding = "UTF-8")
source("modules/mod_motif_scanner.R",      encoding = "UTF-8")
# ── Interface (une fonction par onglet) ────────────────────────
source("ui/ui_home.R",            encoding = "UTF-8")
source("ui/ui_dataset.R",         encoding = "UTF-8")
source("ui/ui_partial_genes.R",   encoding = "UTF-8")
source("ui/ui_descriptive.R",     encoding = "UTF-8")
source("ui/ui_conservation.R",    encoding = "UTF-8")
source("ui/ui_similarity.R",      encoding = "UTF-8")
source("ui/ui_diversity.R",       encoding = "UTF-8")
source("ui/ui_advanced.R",        encoding = "UTF-8")
source("ui/ui_scanner.R",         encoding = "UTF-8")
source("ui/ui_resources.R",       encoding = "UTF-8")
