# ═══════════════════════════════════════════════════════════════
# Utilitaires communs 
# ═══════════════════════════════════════════════════════════════

# Dépend de : rien 
# Packages  : ggplot2

#  Opérateur "valeur par défaut
`%||%` <- function(a, b) {
	if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a
}

# Palettes
MOTIF_COLORS <- c(
	"TATA_BOX"    = "#E63946",
	"HEPTA_PENTA" = "#F4A261",
	"OCTA_DECA"   = "#2A9D8F",
	"E_BOX"       = "#457B9D",
	"TATA_BOX2"   = "#9B2226",
	"EBOX_2"      = "#606C38",
	"LPART_1"     = "#7B2D8B"
)

FUNC_COLORS <- c("F" = "#28A745", "P" = "#DC3545", "ORF" = "#FFC107")
FUNC_LABEL_COLORS <- c("F" = "white", "P" = "white", "ORF" = "grey20")

LOCI_COLORS <- c(
	IGHV = "#1565C0", IGKV = "#2E7D32", IGLV = "#F9A825",
	TRAV = "#AD1457", TRBV = "#6A1B9A", TRDV = "#00838F", TRGV = "#E65100"
)

SP_PALETTE <- c(
	"#E63946", "#F4A261", "#2A9D8F", "#8338EC", "#FB5607",
	"#FF006E", "#06D6A0", "#FFBE0B", "#8B0000", "#457B9D",
	"#A8DADC", "#6D6875", "#B5838D", "#264653", "#C77DFF"
)

get_species_palette <- function(species_vec) {
	sp <- unique(species_vec)
	setNames(rep_len(SP_PALETTE, length(sp)), sp)
}

get_motif_palette <- function(ms) { cols <- MOTIF_COLORS[ms]; cols[is.na(cols)] <- "#999999"; cols }

# Thème ggplot commun
theme_imgt <- function(base_size = 13) {
	ggplot2::theme_minimal(base_size = base_size) +
		ggplot2::theme(
			plot.title       = ggplot2::element_text(size = base_size + 3, face = "bold",
				hjust = 0.5, margin = ggplot2::margin(b = 10)),
			plot.subtitle    = ggplot2::element_text(size = base_size - 1, hjust = 0.5,
				color = "grey50", margin = ggplot2::margin(b = 8)),
			axis.title       = ggplot2::element_text(size = base_size, face = "bold"),
			axis.text        = ggplot2::element_text(size = base_size - 2, color = "grey30"),
			legend.title     = ggplot2::element_text(size = base_size - 1, face = "bold"),
			legend.text      = ggplot2::element_text(size = base_size - 2),
			legend.position  = "bottom",
			legend.key.size  = ggplot2::unit(0.5, "cm"),
			strip.background = ggplot2::element_rect(fill = "#001f54", color = NA),
			strip.text       = ggplot2::element_text(face = "bold", color = "white",
				size = base_size - 1),
			panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.4),
			panel.grid.minor = ggplot2::element_blank(),
			panel.spacing    = ggplot2::unit(0.8, "lines"),
			plot.margin      = ggplot2::margin(10, 15, 10, 15)
		)
}

.logit_terms <- function(df_model) {
	terms <- c("LOCI", "SPECIES", "FUNCTIONALITY")
	keep <- terms[vapply(terms, function(v) nlevels(droplevels(as.factor(df_model[[v]]))) >= 2, logical(1))]
	if (length(keep) == 0) stop("No explanatory factor with at least 2 levels (single species, locus and functionality).")
	list(terms = keep, formula = stats::as.formula(paste("Presence ~", paste(keep, collapse = " + "))))
}