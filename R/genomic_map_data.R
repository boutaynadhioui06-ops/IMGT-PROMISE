# ═══════════════════════════════════════════════════════════════
# Genomic map
# ═══════════════════════════════════════════════════════════════

# Dépend de : rien
# Packages  : dplyr, stringr

library(tidyverse)


# ══════════════════════════════════════════════════════════════
# CARTE GÉNOMIQUE DES MOTIFS 
# ══════════════════════════════════════════════════════════════

prepare_genomic_map <- function(df,
	type_data = "IG",
	species   = NULL,
	locus     = NULL,
	subgroup  = NULL,
	max_seqs  = 50) {
	
	if (type_data == "IG") {
		df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	} else {
		df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	}
	if (!is.null(species)  && species  != "" && species  != "All")
		df <- df %>% dplyr::filter(SPECIES == species)
	if (!is.null(locus)    && locus    != "" && locus    != "All")
		df <- df %>% dplyr::filter(LOCI == locus)
	if (!is.null(subgroup) && subgroup != "" && subgroup != "All")
		df <- df %>% dplyr::filter(SUBGROUP == subgroup)
	
	if (nrow(df) == 0) stop("No sequence after filtering.")
	
	pos_start_cols <- names(df)[stringr::str_starts(names(df), "POS_START_")]
	pos_end_cols   <- names(df)[stringr::str_starts(names(df), "POS_END_")]
	pos_start_cols <- pos_start_cols[!stringr::str_detect(pos_start_cols, "LPART")]
	pos_end_cols   <- pos_end_cols[!stringr::str_detect(pos_end_cols, "LPART")]
	
	if (length(pos_start_cols) == 0)
		stop("No POS_START_* column available (excluding LPART1).")
	
	motif_names <- stringr::str_remove(pos_start_cols, "^POS_START_")
	has_utr_pos <- "POS_START" %in% names(df)
	
	if (nrow(df) > max_seqs) {
		df <- df %>%
			dplyr::arrange(LOCI, SUBGROUP, ALLELE) %>%
			dplyr::slice_head(n = max_seqs)
	}
	
	rows <- list()
	
	for (i in seq_len(nrow(df))) {
		row     <- df[i, ]
		utr_ref <- if (has_utr_pos && !is.na(row$POS_START)) row$POS_START else 0L
		utr_len <- if ("LENGTH_NT" %in% names(row) && !is.na(row$LENGTH_NT))
			row$LENGTH_NT else NA_integer_
		
		rows[[paste0("UTR_", i)]] <- data.frame(
			ALLELE     = as.character(row$ALLELE),
			LOCI       = as.character(row$LOCI),
			SUBGROUP   = as.character(row$SUBGROUP),
			SPECIES    = as.character(row$SPECIES),
			Motif      = "5'UTR",
			Start_rel  = 0L,
			End_rel    = if (!is.na(utr_len)) as.integer(utr_len) else 200L,
			Length     = if (!is.na(utr_len)) as.integer(utr_len) else 200L,
			UTR_length = if (!is.na(utr_len)) as.integer(utr_len) else NA_integer_,
			Present    = TRUE,
			stringsAsFactors = FALSE
		)
		
		for (m in motif_names) {
			sc <- paste0("POS_START_", m)
			ec <- paste0("POS_END_",   m)
			bc <- paste0("BIN_", m)
			if (!sc %in% names(df) || !ec %in% names(df)) next
			ps <- row[[sc]]; pe <- row[[ec]]
			present <- if (bc %in% names(df)) !is.na(row[[bc]]) && row[[bc]] == 1L else !is.na(ps)
			ps_rel <- if (!is.na(ps) && ps > 0) as.integer(ps - utr_ref) else NA_integer_
			pe_rel <- if (!is.na(pe) && pe > 0) as.integer(pe - utr_ref) else NA_integer_
			if (!is.na(ps_rel) && ps_rel < 0) ps_rel <- 0L
			if (!is.na(ps_rel) && !is.na(pe_rel) && pe_rel > ps_rel) {
				rows[[paste0(m, "_", i)]] <- data.frame(
					ALLELE = as.character(row$ALLELE), LOCI = as.character(row$LOCI),
					SUBGROUP = as.character(row$SUBGROUP), SPECIES = as.character(row$SPECIES),
					Motif = m, Start_rel = ps_rel, End_rel = pe_rel,
					Length = pe_rel - ps_rel,
					UTR_length = if (!is.na(utr_len)) as.integer(utr_len) else NA_integer_,
					Present = present, stringsAsFactors = FALSE
				)
			}
		}
	}
	
	if (length(rows) == 0) stop("No valid position found.")
	dplyr::bind_rows(rows)
}
