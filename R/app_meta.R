# ═══════════════════════════════════════════════════════════════════════════
#  Métadonnées centralisées 
# ═══════════════════════════════════════════════════════════════════════════
# Dépend de : rien
# Packages  : base R
APP_META <- list(
	name        = "IMGT/PROMISE",
	expansion   = "PROmoter Motif Identification and Sequence Explorer",
	tagline     = "Statistical analysis of regulatory regions of IG and TR V-genes",
	version     = "1.0.0",
	release_date = "2026-10-15",          
	status      = "Pre-release",
	authors     = list(
		list(name = "Boutayna Dhioui", role = "Developer",
			affiliation = "IMGT® Team, Institut de Génétique Humaine (IGH), CNRS – Université de Montpellier")
	),
	contact     = "boutayna.dhioui@igh.cnrs.fr",
	website     = "https://www.imgt.org",
	# ── je dois renseigner qd j'aurais publication ──────────────────────────────
	license     = NULL,   
	repo_url    = NULL,   
	doi         = NULL,   
	article     = NULL    
)

meta_or <- function(x, fallback = "To be defined at public release") {
	if (is.null(x) || length(x) == 0 || is.na(x[1]) || !nzchar(x[1])) fallback else x
}

APP_FOOTER <- function() {
	sprintf("%s v%s · IMGT® · %s",
		APP_META$name, APP_META$version, meta_or(APP_META$license, "License pending"))
}

#  Citation
promise_citation_text <- function() {
	a <- APP_META$authors[[1]]$name
	url <- if (!is.null(APP_META$repo_url)) paste0(" ", APP_META$repo_url) else ""
	doi <- if (!is.null(APP_META$doi)) paste0(" https://doi.org/", APP_META$doi) else ""
	sprintf("%s (%s). %s: %s (version %s) [Computer software]. IMGT®, the international ImMunoGeneTics information system®.%s%s",
		a, substr(APP_META$release_date, 1, 4), APP_META$name, APP_META$expansion,
		APP_META$version, doi, url)
}

promise_bibtex <- function() {
	fields <- c(
		"  author  = {Dhioui, Boutayna}",
		sprintf("  title   = {{%s}: %s}", APP_META$name, APP_META$expansion),
		sprintf("  year    = {%s}", substr(APP_META$release_date, 1, 4)),
		sprintf("  version = {%s}", APP_META$version),
		"  publisher = {IMGT\\textregistered{}, the international ImMunoGeneTics information system}",
		if (!is.null(APP_META$doi))      sprintf("  doi     = {%s}", APP_META$doi),
		if (!is.null(APP_META$repo_url)) sprintf("  url     = {%s}", APP_META$repo_url)
	)
	paste0("@software{dhioui_imgt_promise_", substr(APP_META$release_date, 1, 4), ",\n",
		paste(fields, collapse = ",\n"), "\n}")
}

promise_ris <- function() {
	paste(c(
		"TY  - COMP",
		"AU  - Dhioui, Boutayna",
		sprintf("TI  - %s: %s", APP_META$name, APP_META$expansion),
		sprintf("PY  - %s", substr(APP_META$release_date, 1, 4)),
		sprintf("ET  - %s", APP_META$version),
		if (!is.null(APP_META$doi))      sprintf("DO  - %s", APP_META$doi),
		if (!is.null(APP_META$repo_url)) sprintf("UR  - %s", APP_META$repo_url),
		"ER  - "
	), collapse = "\n")
}

# Sections disponibles dans le rapport 
REPORT_SECTIONS <- c(
	"Dataset overview"                = "overview",
	"Descriptive statistics"          = "descriptive",
	"Motif presence frequencies"      = "freq",
	"Motif diversity"                 = "diversity",
	"Motif distances to L-PART1"      = "distances",
	"Partial genes (no-motif)"        = "partial",
	"Motif scanner: method comparison" = "scanner"
)

#  Packages 
REQUIRED_PACKAGES <- c("shiny", "shinydashboard", "dplyr", "tidyr", "purrr", "stringr",
	"ggplot2", "plotly", "DT", "vegan", "openxlsx", "patchwork",
	"htmltools", "rmarkdown", "knitr")

# Dictionnaire de données
DATA_DICTIONARY <- data.frame(
	Column = c("ALLELE", "SPECIES", "GENE", "LOCI", "SUBGROUP", "FUNCTIONALITY",
		"LENGTH_NT", "SEQ_<MOTIF>", "POS_START_<MOTIF>", "POS_END_<MOTIF>",
		"BIN_<MOTIF>", "DIST_<MOTIF>", "SEQ_UPSTREAM"),
	Type = c("text", "text", "text", "text", "text", "text",
		"integer", "text", "integer", "integer", "0/1", "integer", "text"),
	Required = c("Yes", "Yes", "Yes", "Yes", "Derived", "Yes",
		"No", "Per motif", "Per motif", "Per motif", "Derived", "Derived", "Reference only"),
	Description = c(
		"Allele name (IMGT nomenclature), join key.",
		"Species name, join key.",
		"Gene name (allele without the *NN suffix), join key.",
		"Locus (IGHV, IGKV, IGLV, TRAV, TRBV, TRDV, TRGV).",
		"Subgroup derived from the gene name (e.g. IGHV1).",
		"F (functional), P (pseudogene) or ORF.",
		"Length of the 5′UTR sequence in nucleotides.",
		"Motif sequence; 'absent' when not detected.",
		"Motif start position in the 5′UTR (1-based).",
		"Motif end position in the 5′UTR (1-based).",
		"Motif presence (1) / absence (0), built from SEQ_<MOTIF>.",
		"POS_START_LPART_1 − POS_END_<MOTIF> (nt); negative values set to NA.",
		"Upstream region used for PWM calibration and method comparison."),
	Example = c("IGHV1-69*01", "Homo sapiens", "IGHV1-69", "IGHV", "IGHV1", "F",
		"1204", "ATGCAAAT", "812", "819", "1", "385", "...ACTG..."),
	stringsAsFactors = FALSE
)

# Glossaire
GLOSSARY <- data.frame(
	Term = c("5′UTR", "L-PART1", "Motif", "Consensus", "IUPAC code", "PWM", "Log-odds score",
		"Null distribution", "p-value", "q-value / FDR", "Sensitivity", "Specificity",
		"Precision", "F1-score", "Shannon index (H′)", "Simpson index (D)",
		"Richness (S)", "Evenness (J)", "PERMANOVA", "PCoA", "Odds ratio (OR)"),
	Definition = c(
		"5′ untranslated region upstream of the V-gene coding sequence, containing regulatory elements.",
		"First exon of the leader peptide; used as the anchor to measure motif distances.",
		"Short conserved regulatory sequence (octamer, decamer, heptanucleotide, pentadecamer, TATA-box, E-box).",
		"Fixed reference pattern for a motif, taken from the literature or IMGT curation.",
		"Ambiguity alphabet for nucleotides (e.g. N = any base).",
		"Position Weight Matrix: per-position nucleotide frequencies learned from confirmed instances.",
		"Log2 ratio of the motif model probability to the background probability; higher = better match.",
		"Distribution of the maximum PWM score on random sequences with the background composition (Monte Carlo).",
		"Probability of a result at least as extreme as observed, under the null hypothesis.",
		"p-value adjusted for multiple testing (Benjamini–Hochberg); expected proportion of false discoveries.",
		"TP / (TP + FN): proportion of confirmed motifs that are detected.",
		"TN / (TN + FP): proportion of shuffled negative controls correctly left undetected.",
		"TP / (TP + FP): proportion of detections that are correct.",
		"Harmonic mean of precision and sensitivity.",
		"Overall diversity of motif profiles: −Σ pᵢ ln(pᵢ).",
		"Probability that two draws belong to different motifs: 1 − Σ pᵢ².",
		"Number of distinct motifs detected in a group.",
		"Uniformity of the motif distribution: H′ / ln(S), between 0 and 1.",
		"Permutational multivariate analysis of variance testing differences between group centroids.",
		"Principal coordinates analysis: low-dimensional ordination of a distance matrix.",
		"exp(coefficient) of a logistic regression; >1 increases, <1 decreases the odds relative to the reference level."),
	stringsAsFactors = FALSE
)

# Changelog
CHANGELOG <- list(
	list(version = "1.0.0", date = "2026-10-15", status = "Pre-release",
		items = c("Dataset construction from 5′UTR and motif files, with Roman-numeral clan exclusion.",
			"Descriptive statistics, motif conservation, positioning, similarity, diversity and co-occurrence analyses.",
			"Partial genes analysis and literature comparison.",
			"Motif scanner: consensus baseline and PWM-based probabilistic detection, with k-fold method comparison.",
			"Dynamic analysis report (HTML / Word / PDF) with reproducibility information."))
)