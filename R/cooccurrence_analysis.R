# ═══════════════════════════════════════════════════════════════
# Coocurrence analysis
# ═══════════════════════════════════════════════════════════════

# Dépend de : rien
# Packages  : dplyr, tidyr, stringr, igraph

library(tidyverse)
library(igraph)
library(vegan)

# ══════════════════════════════════════════════════════════════
# CO-OCCURRENCE DES MOTIFS
# ══════════════════════════════════════════════════════════════

calculate_cooccurrence <- function(df,
	type_data    = "IG",
	species      = NULL,
	locus        = NULL,
	func_filter  = NULL) {
	
	if (type_data == "IG") {
		df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	} else {
		df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	}
	if (!is.null(species) && species != "" && species != "All")
		df <- df %>% dplyr::filter(SPECIES == species)
	if (!is.null(locus) && locus != "" && locus != "All")
		df <- df %>% dplyr::filter(LOCI == locus)
	if (!is.null(func_filter) && func_filter != "" && func_filter != "All")
		df <- df %>% dplyr::filter(FUNCTIONALITY == func_filter)
	
	if (nrow(df) == 0) stop("No sequences after filtering.")
	
	bin_cols <- names(df)[stringr::str_starts(names(df), "BIN_")]
	if (length(bin_cols) < 2) stop("Less than 2 BIN_* motifs available.")
	
	motif_names <- stringr::str_remove(bin_cols, "^BIN_")
	n_seq       <- nrow(df)
	n_motifs    <- length(bin_cols)
	
	bin_mat <- df %>%
		dplyr::select(dplyr::all_of(bin_cols)) %>%
		as.matrix()
	colnames(bin_mat) <- motif_names
	
	mat_pct <- matrix(NA_real_, n_motifs, n_motifs, dimnames = list(motif_names, motif_names))
	mat_phi <- matrix(NA_real_, n_motifs, n_motifs, dimnames = list(motif_names, motif_names))
	
	for (i in seq_len(n_motifs)) {
		for (j in seq_len(n_motifs)) {
			xi <- bin_mat[, i]
			xj <- bin_mat[, j]
			mat_pct[i, j] <- round(mean(xi == 1 & xj == 1, na.rm = TRUE) * 100, 1)
			if (i != j) {
				cc <- suppressWarnings(cor(xi, xj, use = "complete.obs"))
				mat_phi[i, j] <- round(ifelse(is.na(cc), 0, cc), 3)
			} else {
				mat_phi[i, j] <- 1
			}
		}
	}
	
	return(list(
		matrix_pct = mat_pct,
		matrix_phi = mat_phi,
		n_seq      = n_seq,
		motifs     = motif_names,
		bin_mat    = bin_mat
	))
}


fisher_pairwise_tests <- function(cooc_result) {
	
	bin_mat <- cooc_result$bin_mat
	motifs  <- cooc_result$motifs
	n       <- length(motifs)
	results <- list()
	
	for (i in seq_len(n - 1)) {
		for (j in (i + 1):n) {
			xi <- bin_mat[, i]
			xj <- bin_mat[, j]
			
			tbl <- table(xi, xj)
			if (nrow(tbl) < 2 || ncol(tbl) < 2) next
			
			ft <- tryCatch(
				fisher.test(tbl, simulate.p.value = TRUE, B = 2000),
				error = function(e) NULL
			)
			if (is.null(ft)) next
			
			or <- round(ft$estimate, 3)
			pv <- ft$p.value
			
			freq_i <- round(mean(xi == 1, na.rm = TRUE) * 100, 1)
			freq_j <- round(mean(xj == 1, na.rm = TRUE) * 100, 1)
			cooc   <- round(mean(xi == 1 & xj == 1, na.rm = TRUE) * 100, 1)
			results[[paste(motifs[i], motifs[j], sep = "__")]] <- data.frame(
				Motif_1     = motifs[i],
				Motif_2     = motifs[j],
				Freq_Motif1 = freq_i,
				Freq_Motif2 = freq_j,
				CoOcc_Pct   = cooc,
				OddsRatio   = or,
				p_value     = round(pv, 4),
				Significant = dplyr::case_when(
					pv < 0.001 ~ "***",
					pv < 0.01  ~ "**",
					pv < 0.05  ~ "*",
					TRUE       ~ "ns"
				),
				Association = dplyr::case_when(
					pv >= 0.05  ~ "Not significant",
					or > 1      ~ "Positive co-occurrence",
					or < 1      ~ "Mutual exclusion",
					TRUE        ~ "Independent"
				),
				stringsAsFactors = FALSE
			)
		}
	}
	
	if (length(results) == 0) return(NULL)
	dplyr::bind_rows(results) %>% dplyr::arrange(p_value)
}


build_cooccurrence_graph <- function(fisher_df,
	cooc_result,
	phi_threshold = 0.05,
	sig_only      = FALSE) {
	
	if (is.null(fisher_df) || nrow(fisher_df) == 0)
		stop("No Fisher results available.")
	
	df_edges <- fisher_df
	if (sig_only)
		df_edges <- df_edges %>% dplyr::filter(Significant != "ns")
	
	mat_phi      <- cooc_result$matrix_phi
	df_edges$phi <- mapply(function(m1, m2) {
		if (m1 %in% rownames(mat_phi) && m2 %in% colnames(mat_phi))
			mat_phi[m1, m2]
		else NA_real_
	}, df_edges$Motif_1, df_edges$Motif_2)
	
	df_edges <- df_edges %>%
		dplyr::filter(!is.na(phi), abs(phi) >= phi_threshold)
	
	motifs    <- cooc_result$motifs
	bin_mat   <- cooc_result$bin_mat
	node_freq <- sapply(motifs, function(m) {
		round(mean(bin_mat[, m] == 1, na.rm = TRUE) * 100, 1)
	})
	vertices_df <- data.frame(
		name = motifs,
		freq = node_freq,
		stringsAsFactors = FALSE
	)
	
	if (nrow(df_edges) == 0) {
		
		g <- igraph::graph_from_data_frame(
			d        = data.frame(from = character(0), to = character(0)),
			vertices = vertices_df,
			directed = FALSE
		)
		message(paste0(
			"build_cooccurrence_graph: no edge with |phi| >= ",
			phi_threshold,
			". Graph returned with nodes only."
		))
		return(g)
	}
	
	edges_df <- df_edges %>%
		dplyr::select(
			from     = Motif_1,
			to       = Motif_2,
			weight   = phi,
			cooc_pct = CoOcc_Pct,
			or       = OddsRatio,
			p_value  = p_value,
			sig      = Significant,
			assoc    = Association
		) %>%
		dplyr::mutate(
			abs_weight = abs(weight),
			edge_color = dplyr::case_when(
				sig   == "ns"                       ~ "#CBD5E1",
				assoc == "Positive co-occurrence"   ~ "#16A34A",
				assoc == "Mutual exclusion"         ~ "#DC2626",
				TRUE                                ~ "#6B7280"
			)
		)
	
	g <- igraph::graph_from_data_frame(
		d        = edges_df,
		vertices = vertices_df,
		directed = FALSE
	)
	return(g)
}

