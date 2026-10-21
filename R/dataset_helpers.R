# ═══════════════════════════════════════════════════════════════
#Fonctions utilitaires sur le dataset
# ═══════════════════════════════════════════════════════════════

#  Nettoyage des sous-groupes 
clean_subgroup <- function(x) {
	x   <- as.character(x)
	out <- x
	pos <- regexpr("^[A-Za-z]+[0-9]+", x)
	
	matched <- regmatches(x, pos)
	idx     <- which(pos > 0 & !is.na(pos))
	out[idx] <- matched
	out
}

#Renomme la colonne LOCI en GROUP
rename_loci_display <- function(df) {
	if ("LOCI" %in% names(df)) {
		names(df)[names(df) == "LOCI"] <- "GROUP"
	}
	df
}

# Séparation IG / TR du dataset final
LOCI_IG <- c("IGHV", "IGKV", "IGLV")
LOCI_TR <- c("TRAV", "TRBV", "TRDV", "TRGV")

# Columns liées aux motifs régulateurs, non pertinentes pour les TR
get_motif_cols <- function(df) {
	names(df)[grepl("^(SEQ_|POS_START_|POS_END_|BIN_|DIST_)", names(df))]
}

split_dataset_ig_tr <- function(df) {
	if (is.null(df)) return(list(ig = NULL, tr = NULL))
	motif_cols <- get_motif_cols(df)
	df_ig <- df %>% dplyr::filter(LOCI %in% LOCI_IG)
	df_tr <- df %>% dplyr::filter(LOCI %in% LOCI_TR) %>%
		dplyr::select(-dplyr::any_of(motif_cols))
	list(ig = df_ig, tr = df_tr)
}
clean_functionality <- function(raw) {
	f <- gsub("[\\[\\]\\(\\)]", "", raw, perl = TRUE)
	f <- trimws(f)
	f <- toupper(f)
	if (!f %in% c("F", "P", "ORF")) f <- NA_character_
	return(f)
}
