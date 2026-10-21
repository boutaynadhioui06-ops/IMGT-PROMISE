# ═══════════════════════════════════════════════════════════════
# Module Gènes Partiels 
# ═══════════════════════════════════════════════════════════════

# Dépend de : utils_core.R, dataset_helpers.R
# Packages  : tidyverse, plotly

# ══════════════════════════════════════════════════════════════
#  STATISTIQUES DESCRIPTIVES 
# ══════════════════════════════════════════════════════════════

# Columns BIN_ correspondant à de vrais motifs régulateurs 
get_regulatory_bin_cols <- function(df) {
	bin_cols <- names(df)[startsWith(names(df), "BIN_")]
	bin_cols[!grepl("^BIN_LPART[_]?1$", bin_cols, ignore.case = TRUE)]
}

# Ajoute nb de motifs régulateurs présents 
add_motif_count <- function(df) {
	bin_cols <- get_regulatory_bin_cols(df)
	if (length(bin_cols) == 0) { df$N_MOTIFS <- 0L; return(df) }
	mat <- as.matrix(df[, bin_cols, drop = FALSE])
	mat[is.na(mat)] <- 0
	df$N_MOTIFS <- as.integer(rowSums(mat))
	df
}

# Sequences IG SANS aucun motif régulateur 
identify_no_motif_sequences <- function(df_complet) {
	if (is.null(df_complet) || nrow(df_complet) == 0) return(data.frame())
	df_ig <- df_complet %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	if (nrow(df_ig) == 0) return(data.frame())
	df_ig <- add_motif_count(df_ig)
	df_no <- df_ig %>% dplyr::filter(N_MOTIFS == 0)
	cols_keep <- intersect(
		c("ALLELE","GENE","SUBGROUP","LOCI","SPECIES","FUNCTIONALITY","LENGTH_NT"),
		names(df_no))
	df_no %>%
		dplyr::select(dplyr::all_of(cols_keep)) %>%
		dplyr::arrange(SPECIES, LOCI, SUBGROUP, ALLELE)
}

# Sequences IG AVEC au moins un motif régulateur 
identify_with_motif_sequences <- function(df_complet) {
	if (is.null(df_complet) || nrow(df_complet) == 0) return(data.frame())
	df_ig <- df_complet %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	if (nrow(df_ig) == 0) return(data.frame())
	df_ig <- add_motif_count(df_ig)
	df_with <- df_ig %>% dplyr::filter(N_MOTIFS > 0)
	
	bin_cols    <- get_regulatory_bin_cols(df_ig)
	motif_names <- gsub("^BIN_", "", bin_cols)
	dist_cols   <- intersect(paste0("DIST_", motif_names), names(df_with))
	
	cols_keep <- intersect(
		c("ALLELE","GENE","SUBGROUP","LOCI","SPECIES","FUNCTIONALITY",
			"N_MOTIFS", dist_cols, "LENGTH_NT"),
		names(df_with))
	df_with %>%
		dplyr::select(dplyr::all_of(cols_keep)) %>%
		dplyr::arrange(SPECIES, LOCI, SUBGROUP, ALLELE)
}

# Filtre Species -> Locus -> Subgroup 
apply_partial_cascade_filters <- function(df, species_sel = NULL, locus_sel = NULL, subgroup_sel = NULL) {
	if (is.null(df) || nrow(df) == 0) return(df)
	if (!is.null(species_sel)  && species_sel  != "" && species_sel  != "All species")
		df <- df %>% dplyr::filter(SPECIES == species_sel)
	if (!is.null(locus_sel)    && locus_sel    != "" && locus_sel    != "All loci")
		df <- df %>% dplyr::filter(LOCI == locus_sel)
	if (!is.null(subgroup_sel) && subgroup_sel != "" && subgroup_sel != "All subgroups")
		df <- df %>% dplyr::filter(SUBGROUP == subgroup_sel)
	df
}

# Résumé par espèce/locus — séquences SANS motif
summary_no_motif_by_species <- function(df_no_motif) {
	if (is.null(df_no_motif) || nrow(df_no_motif) == 0) return(data.frame())
	df_no_motif %>%
		dplyr::group_by(SPECIES, LOCI) %>%
		dplyr::summarise(
			N_alleles    = dplyr::n(),
			Fonctionnels = sum(FUNCTIONALITY == "F",   na.rm = TRUE),
			Pseudogenes  = sum(FUNCTIONALITY == "P",   na.rm = TRUE),
			ORF          = sum(FUNCTIONALITY == "ORF", na.rm = TRUE),
			Length_mean  = as.integer(round(mean(LENGTH_NT, na.rm = TRUE), 0)),
			Length_min   = suppressWarnings(min(LENGTH_NT, na.rm = TRUE)),
			Length_max   = suppressWarnings(max(LENGTH_NT, na.rm = TRUE)),
			.groups = "drop"
		) %>%
		dplyr::arrange(SPECIES, LOCI)
}

# Résumé par espèce/locus — séquences AVEC motif
summary_with_motif_by_species <- function(df_with_motif) {
	if (is.null(df_with_motif) || nrow(df_with_motif) == 0) return(data.frame())
	df_with_motif %>%
		dplyr::group_by(SPECIES, LOCI) %>%
		dplyr::summarise(
			N_alleles     = dplyr::n(),
			Mean_N_motifs = round(mean(N_MOTIFS, na.rm = TRUE), 2),
			Min_N_motifs  = suppressWarnings(min(N_MOTIFS, na.rm = TRUE)),
			Max_N_motifs  = suppressWarnings(max(N_MOTIFS, na.rm = TRUE)),
			Length_mean   = as.integer(round(mean(LENGTH_NT, na.rm = TRUE), 0)),
			.groups = "drop"
		) %>%
		dplyr::arrange(SPECIES, LOCI)
}

# Graphique : séquences SANS motif, par espèce et fonctionnalité
plot_no_motif_by_species <- function(df_no_motif, interactive = TRUE) {
	if (is.null(df_no_motif) || nrow(df_no_motif) == 0)
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No sequences without motif for this selection", size = 6, color = "grey50") +
				ggplot2::theme_void())
	
	data_plot <- df_no_motif %>%
		dplyr::filter(!is.na(FUNCTIONALITY)) %>%
		dplyr::group_by(SPECIES, LOCI, FUNCTIONALITY) %>%
		dplyr::summarise(N = dplyr::n(), .groups = "drop")
	if (nrow(data_plot) == 0)
		return(ggplot2::ggplot() + ggplot2::theme_void())
	
	ordre_sp <- data_plot %>% dplyr::group_by(SPECIES) %>%
		dplyr::summarise(tot = sum(N), .groups = "drop") %>%
		dplyr::arrange(tot) %>% dplyr::pull(SPECIES)
	data_plot$SPECIES       <- factor(data_plot$SPECIES, levels = ordre_sp)
	data_plot$FUNCTIONALITY <- factor(data_plot$FUNCTIONALITY, levels = c("F","P","ORF"))
	
	FUNC_COLORS <- c("F" = "#28A745", "P" = "#DC3545", "ORF" = "#FFC107")
	n_sp     <- length(ordre_sp)
	n_facets <- min(length(unique(data_plot$LOCI)), 3)
	
	p <- ggplot2::ggplot(data_plot, ggplot2::aes(x = SPECIES, y = N, fill = FUNCTIONALITY)) +
		ggplot2::geom_bar(stat = "identity", width = 0.65, position = "stack") +
		ggplot2::coord_flip() +
		ggplot2::facet_wrap(~ LOCI, ncol = n_facets, scales = "free_x") +
		ggplot2::scale_fill_manual(values = FUNC_COLORS, name = "Functionality", drop = FALSE) +
		ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.28))) +
		ggplot2::labs(
			title    = "Sequences with NO regulatory motif (LPART1 excluded)",
			subtitle = paste0(sum(data_plot$N), " sequences across ", n_sp, " species"),
			x = NULL, y = "Number of sequences"
		) +
		theme_imgt(base_size = 13) +
		ggplot2::theme(
			axis.text.y     = ggplot2::element_text(face = "bold", size = 11),
			strip.text      = ggplot2::element_text(face = "bold", size = 12),
			legend.position = "bottom",
			panel.spacing   = ggplot2::unit(1.3, "lines")
		)
	
	h <- max(420, n_sp * 42)
	if (interactive) p <- plotly::ggplotly(p, tooltip = c("x","y","fill"), height = min(h, 1400) + 60) %>%
		plotly::layout(legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.18),
			margin = list(b = 90, r = 60))
	p
}

# Distribution des longueurs de 5'UTR — séquences SANS motif
plot_no_motif_length_dist <- function(df_no_motif, interactive = TRUE) {
	if (is.null(df_no_motif) || nrow(df_no_motif) == 0)
		return(ggplot2::ggplot() + ggplot2::theme_void())
	df <- df_no_motif %>% dplyr::filter(!is.na(LENGTH_NT))
	if (nrow(df) == 0)
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No length data for this selection", size = 5, color = "grey50") +
				ggplot2::theme_void())
	
	rng <- range(df$LENGTH_NT, na.rm = TRUE)
	bw  <- max(1, ceiling(diff(rng) / 20))
	n_loci_present <- length(unique(df$LOCI))
	
	p <- ggplot2::ggplot(df, ggplot2::aes(x = LENGTH_NT, fill = LOCI)) +
		ggplot2::geom_histogram(binwidth = bw, color = "white", alpha = 0.85) +
		ggplot2::scale_fill_manual(values = LOCI_COLORS, name = "Locus") +
		ggplot2::facet_wrap(~ LOCI, scales = "free_y", ncol = min(n_loci_present, 3)) +
		ggplot2::labs(
			title    = "5'UTR length — sequences with NO regulatory motif",
			subtitle = paste0(nrow(df), " sequences"),
			x = "5'UTR length (nt)", y = "Count"
		) +
		theme_imgt(base_size = 13) +
		ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(face = "bold", size = 12))
	
	h <- if (n_loci_present <= 1) 420 else 460
	if (interactive) p <- plotly::ggplotly(p, tooltip = c("x","y"), height = h) %>%
		plotly::layout(margin = list(t = 60, r = 40))
	p
}

# Distribution du nombre de motifs par allèle — séquences AVEC motif
plot_with_motif_count <- function(df_with_motif, interactive = TRUE) {
	if (is.null(df_with_motif) || nrow(df_with_motif) == 0)
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No sequences with motif for this selection", size = 5, color = "grey50") +
				ggplot2::theme_void())
	
	data_plot <- df_with_motif %>% dplyr::count(SPECIES, LOCI, N_MOTIFS)
	
	p <- ggplot2::ggplot(data_plot, ggplot2::aes(x = factor(N_MOTIFS), y = n, fill = LOCI)) +
		ggplot2::geom_col(position = "dodge") +
		ggplot2::facet_wrap(~ SPECIES, scales = "free_y") +
		ggplot2::scale_fill_manual(values = LOCI_COLORS, name = "Locus") +
		ggplot2::labs(
			title = "Number of regulatory motifs per allele (LPART1 excluded)",
			x = "Number of motifs", y = "Number of alleles"
		) +
		theme_imgt(base_size = 12) +
		ggplot2::theme(legend.position = "bottom")
	
	if (interactive) p <- plotly::ggplotly(p, tooltip = c("x","y","fill"), height = 500)
	p
}

# Distribution des distances motif -> LPART1 — séquences AVEC motif
plot_with_motif_distances <- function(df_with_motif, interactive = TRUE) {
	if (is.null(df_with_motif) || nrow(df_with_motif) == 0)
		return(ggplot2::ggplot() + ggplot2::theme_void())
	
	dist_cols <- names(df_with_motif)[startsWith(names(df_with_motif), "DIST_")]
	if (length(dist_cols) == 0)
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No distance data available for this selection", size = 5, color = "grey50") +
				ggplot2::theme_void())
	
	df_long <- df_with_motif %>%
		dplyr::select(SPECIES, LOCI, dplyr::all_of(dist_cols)) %>%
		tidyr::pivot_longer(cols = dplyr::all_of(dist_cols), names_to = "Motif", values_to = "Distance") %>%
		dplyr::mutate(Motif = sub("^DIST_", "", Motif)) %>%
		dplyr::filter(!is.na(Distance))
	
	if (nrow(df_long) == 0)
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No distance values available (motifs may be missing POS_END)", size = 5, color = "grey50") +
				ggplot2::theme_void())
	
	p <- ggplot2::ggplot(df_long, ggplot2::aes(x = Motif, y = Distance, fill = Motif)) +
		ggplot2::geom_boxplot(alpha = 0.8, outlier.size = 0.8) +
		ggplot2::labs(
			title = "Distance between motif and LPART1",
			x = NULL, y = "Distance to LPART1 (nt)"
		) +
		theme_imgt(base_size = 12) +
		ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
	
	if (interactive) p <- plotly::ggplotly(p, tooltip = c("x","y"), height = 480)
	p
}
# ══════════════════════════════════════════════════════════════
# PARSER FASTA LITTÉRATURE
# ══════════════════════════════════════════════════════════════

parse_literature_fasta <- function(fasta_text) {
	
	lines <- strsplit(fasta_text, "\n")[[1]]
	lines <- trimws(lines)
	lines <- lines[nchar(lines) > 0] 
	
	header_idx <- which(startsWith(lines, ">"))
	if (length(header_idx) == 0)
		stop("No FASTA header found (line starting with '>')")
	
	n           <- length(header_idx)
	seq_end_idx <- c(header_idx[-1] - 1, length(lines))
	
	result <- data.frame(
		ACCESSION     = character(n),
		ALLELE        = character(n),
		GENE          = character(n),
		SUBGROUP      = character(n),
		LOCI          = character(n),
		SPECIES       = character(n),
		FUNCTIONALITY = character(n),
		SEQUENCE      = character(n),
		LENGTH_LIT    = integer(n),
		stringsAsFactors = FALSE
	)
	
	for (i in seq_len(n)) {
		h      <- sub("^>", "", lines[header_idx[i]])
		fields <- strsplit(h, "\\|")[[1]]
		
		#  ACCESSION
		result$ACCESSION[i] <- trimws(fields[1])
		
		#  ALLELE
		if (length(fields) >= 2) {
			ga <- trimws(fields[2])
			result$LOCI[i]     <- gsub("^([A-Za-z]+).*", "\\1", ga)
			result$SUBGROUP[i] <- clean_subgroup(gsub("[-*].*$", "", ga))
			result$GENE[i]     <- strsplit(ga, "\\*")[[1]][1]
			result$ALLELE[i]   <- ga
		}
		
		# SPECIES
		if (length(fields) >= 3) {
			sp <- trimws(fields[3])
			result$SPECIES[i] <- gsub("_.*$", "", sp)
		}
		
		#  FUNCTIONALITY
		if (length(fields) >= 4) {
			result$FUNCTIONALITY[i] <- clean_functionality(fields[4])
		}
		
		# Séquence 
		seq_lines <- if (header_idx[i] + 1 <= seq_end_idx[i]) {
			lines[(header_idx[i] + 1):seq_end_idx[i]]
		} else {
			character(0)   
		}
		seq_lines         <- seq_lines[!startsWith(seq_lines, ">")]
		full_seq          <- paste(trimws(seq_lines), collapse = "")
		result$SEQUENCE[i] <- toupper(full_seq)
		result$LENGTH_LIT[i] <- nchar(full_seq)
	}
	
	# Garder uniquement  IG
	result <- result[grepl("V$", result$LOCI, ignore.case = FALSE) &
			result$LOCI %in% c("IGHV", "IGKV", "IGLV"), ]
	
	return(result)
}


# ══════════════════════════════════════════════════════════════
# COMPARAISON DES SÉQUENCES (300 derniers nt)
# ══════════════════════════════════════════════════════════════

.last_n_nt <- function(seq, n = 300) {
	if (is.na(seq) || nchar(seq) == 0) return(NA_character_)
	#  ne garder que A/C/G/T 
	seq_clean <- gsub("[^ACGTacgt]", "", seq)
	seq_clean <- toupper(seq_clean)
	if (nchar(seq_clean) == 0) return(NA_character_)
	len <- nchar(seq_clean)
	if (len <= n) return(seq_clean)
	substr(seq_clean, len - n + 1, len)
}

.pairwise_identity <- function(s1, s2) {
	if (is.na(s1) || is.na(s2)) return(NA_real_)
	
	# Nettoyage de sécurité 
	s1 <- gsub("[^ACGTacgt]", "", s1)
	s2 <- gsub("[^ACGTacgt]", "", s2)
	s1 <- toupper(s1)
	s2 <- toupper(s2)
	
	if (nchar(s1) == 0 || nchar(s2) == 0) return(NA_real_)
	
	len1 <- nchar(s1)
	len2 <- nchar(s2)
	min_len <- min(len1, len2)
	
	# Seuil de fiabilité : 5 nt minimum
	if (min_len < 5) return(NA_real_)
	
	s1_tail <- substr(s1, len1 - min_len + 1, len1)
	s2_tail <- substr(s2, len2 - min_len + 1, len2)
	
	s1c <- strsplit(s1_tail, "", fixed = TRUE)[[1]]
	s2c <- strsplit(s2_tail, "", fixed = TRUE)[[1]]
	
	round(sum(s1c == s2c) / min_len * 100, 1)
}

prepare_comparison <- function(df_lit, allele_sel, n_last = 300) {
	
	df_allele <- df_lit %>%
		dplyr::filter(ALLELE == allele_sel, !is.na(SEQUENCE), nchar(SEQUENCE) > 0)
	
	if (nrow(df_allele) == 0)
		return(list(sequences = data.frame(), identity_matrix = NULL))
	
	# Extraire les N derniers nt 
	df_allele <- df_allele %>%
		dplyr::mutate(
			SEQ_LAST   = sapply(SEQUENCE, .last_n_nt, n = n_last),
			LEN_ACTUAL = ifelse(is.na(SEQ_LAST), 0L, nchar(SEQ_LAST)),
			SEQ_LABEL  = paste0(ACCESSION, " (", LENGTH_LIT, "nt)")
		)
	
	n <- nrow(df_allele)
	if (n < 2) return(list(sequences = df_allele, identity_matrix = NULL))
	
	# Matrice d'identité pairwise
	mat <- matrix(NA_real_, nrow = n, ncol = n,
		dimnames = list(df_allele$SEQ_LABEL, df_allele$SEQ_LABEL))
	for (i in seq_len(n)) {
		for (j in seq_len(n)) {
			mat[i, j] <- .pairwise_identity(df_allele$SEQ_LAST[i], df_allele$SEQ_LAST[j])
		}
	}
	valid_seq <- !is.na(df_allele$SEQ_LAST) & nchar(df_allele$SEQ_LAST) >= 5
	diag(mat) <- ifelse(valid_seq, 100, NA_real_)
	
	list(sequences = df_allele, identity_matrix = mat)
}

plot_comparison_heatmap <- function(comparison_result, allele_sel, n_last = 300,
	interactive = TRUE) {
	
	mat <- comparison_result$identity_matrix
	if (is.null(mat) || nrow(mat) < 2) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "At least 2 sequences required for comparison",
					size = 5, color = "grey50") +
				ggplot2::theme_void())
	}
	
	n_seq <- nrow(mat)
	n_na  <- sum(is.na(mat[upper.tri(mat)]))
	
	mat_df <- expand.grid(
		Seq1             = rownames(mat),
		Seq2             = colnames(mat),
		KEEP.OUT.ATTRS   = FALSE,
		stringsAsFactors = FALSE
	)
	mat_df$Identity <- as.vector(mat)
	
	mat_dist <- 100 - mat
	mat_dist[is.na(mat_dist)] <- 100
	hc    <- hclust(as.dist(mat_dist), method = "complete")
	order <- rownames(mat)[hc$order]
	mat_df$Seq1 <- factor(mat_df$Seq1, levels = order)
	mat_df$Seq2 <- factor(mat_df$Seq2, levels = rev(order))
	
	tile_size <- dplyr::case_when(
		n_seq <= 8  ~ 4.5,
		n_seq <= 15 ~ 3.5,
		n_seq <= 25 ~ 2.8,
		TRUE        ~ 2.2
	)
	
	fill_colors <- c("#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
		"#FDDBC7", "#D6604D", "#B2182B")
	fill_stops  <- c(0, 40, 60, 75, 85, 95, 100)
	
	# VERSION STATIQUE pour export  ──────────────
	p <- ggplot2::ggplot(mat_df,
		ggplot2::aes(x = Seq1, y = Seq2, fill = Identity)) +
		ggplot2::geom_tile(color = "white", linewidth = 0.4) +
		ggplot2::geom_text(ggplot2::aes(
			label = ifelse(is.na(Identity), "—", paste0(round(Identity, 1), "%")),
			color = dplyr::case_when(
				is.na(Identity)  ~ "grey40",
				Identity < 70    ~ "white",
				TRUE             ~ "grey15"
			)),
			size = tile_size, fontface = "bold") +
		ggplot2::scale_color_identity() +
		ggplot2::scale_fill_gradientn(
			colours  = fill_colors,
			values   = scales::rescale(fill_stops),
			limits   = c(0, 100),
			na.value = "grey85",
			name     = "Identity (%)",
			breaks   = c(0, 25, 50, 75, 100),
			labels   = c("0%", "25%", "50%", "75%", "100%")
		) +
		ggplot2::labs(
			title    = paste0("Pairwise identity — ", allele_sel),
			subtitle = paste0("Last ", n_last, " nt | ", n_seq, " literature sequences",
				if (n_na > 0) paste0(" | ⚠ ", n_na, " pair(s) too short/unreliable (grey)") else ""),
			x = NULL, y = NULL
		) +
		theme_imgt(base_size = 10) +
		ggplot2::theme(
			axis.text.x       = ggplot2::element_text(angle = 40, hjust = 1, size = 8),
			axis.text.y       = ggplot2::element_text(size = 8),
			panel.grid        = ggplot2::element_blank(),
			legend.position   = "right",
			legend.key.height = ggplot2::unit(1.8, "cm")
		)
	
	if (!interactive) return(p)
	
	z_mat    <- mat[order, order, drop = FALSE]
	text_mat <- matrix(ifelse(is.na(z_mat), "—", paste0(round(z_mat, 1), "%")),
		nrow = nrow(z_mat), dimnames = dimnames(z_mat))
	
	full_colorscale <- purrr::map2(
		scales::rescale(fill_stops), fill_colors,
		function(pos, col) list(pos, col)
	)
	
	h_px <- max(400, n_seq * 42 + 120)
	
	plotly::plot_ly(
		x = colnames(z_mat), y = rownames(z_mat), z = z_mat,
		type = "heatmap",
		zmin = 0, zmax = 100, zauto = FALSE,
		colorscale = full_colorscale,
		text = text_mat, texttemplate = "%{text}", textfont = list(size = 11),
		hovertemplate = "Seq1: %{y}<br>Seq2: %{x}<br>Identity: %{z}%<extra></extra>",
		colorbar = list(title = "Identity (%)")
	) %>%
		plotly::layout(
			title  = list(text = paste0("Pairwise identity — ", allele_sel,
				"<br><sub>Last ", n_last, " nt | ", n_seq, " sequences",
				if (n_na > 0) paste0(" | ", n_na, " pair(s) unreliable") else "",
				"</sub>")),
			xaxis  = list(tickangle = 40, title = ""),
			yaxis  = list(title = "", autorange = "reversed"),
			height = min(h_px, 900)
		)
}


prepare_alignment_table <- function(comparison_result) {
	df <- comparison_result$sequences
	if (is.null(df) || nrow(df) == 0) return(data.frame())
	
	df %>%
		dplyr::select(
			dplyr::any_of(c("SEQ_LABEL", "ACCESSION", "SPECIES",
				"FUNCTIONALITY", "LENGTH_LIT", "LEN_ACTUAL", "SEQ_LAST"))
		) %>%
		dplyr::rename(
			dplyr::any_of(c(
				Sequence    = "SEQ_LABEL",
				Accession   = "ACCESSION",
				Species     = "SPECIES",
				Function    = "FUNCTIONALITY",
				Total_nt    = "LENGTH_LIT",
				Compared_nt = "LEN_ACTUAL",
				Alignment   = "SEQ_LAST"
			))
		)
}


summary_comparison <- function(comparison_result) {
	mat <- comparison_result$identity_matrix
	if (is.null(mat)) return(NULL)
		vals <- mat[upper.tri(mat)]
	if (all(is.na(vals))) {
		return(data.frame(
			N_sequences   = nrow(mat),
			Mean_identity = NA_real_,
			Min_identity  = NA_real_,
			Max_identity  = NA_real_,
			SD_identity   = NA_real_,
			Note          = "All pairwise comparisons too short/unreliable"
		))
	}
	data.frame(
		N_sequences  = nrow(mat),
		Mean_identity = round(mean(vals, na.rm = TRUE), 1),
		Min_identity  = round(min(vals,  na.rm = TRUE), 1),
		Max_identity  = round(max(vals,  na.rm = TRUE), 1),
		SD_identity   = round(sd(vals,   na.rm = TRUE), 2)
	)
}

plot_identity_distribution <- function(comparison_result, allele_sel = "",
	n_last = 300, interactive = TRUE) {
	
	mat  <- comparison_result$identity_matrix
	seqs <- comparison_result$sequences
	if (is.null(mat) || nrow(mat) < 2) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "At least 2 sequences required for comparison",
					size = 5, color = "grey50") +
				ggplot2::theme_void())
	}
	
	n_seq <- nrow(mat)
	mat_no_diag <- mat
	diag(mat_no_diag) <- NA
	
	safe_min <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else min(x) }
	safe_max <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else max(x) }
	
	df_summary <- data.frame(
		Sequence = rownames(mat),
		Mean_id  = apply(mat_no_diag, 1, mean, na.rm = TRUE),
		Min_id   = apply(mat_no_diag, 1, safe_min),
		Max_id   = apply(mat_no_diag, 1, safe_max),
		stringsAsFactors = FALSE
	)
	df_summary$Mean_id[is.nan(df_summary$Mean_id)] <- NA_real_
	
	seuil_fiable <- max(50, round(n_last * 0.5))
	if (!is.null(seqs) && all(c("SEQ_LABEL", "LEN_ACTUAL") %in% names(seqs))) {
		df_summary <- df_summary %>%
			dplyr::left_join(
				seqs %>% dplyr::select(SEQ_LABEL, LEN_ACTUAL) %>%
					dplyr::rename(Sequence = SEQ_LABEL),
				by = "Sequence"
			)
	} else {
		df_summary$LEN_ACTUAL <- NA_integer_
	}
	df_summary$Fiable <- ifelse(is.na(df_summary$LEN_ACTUAL), TRUE,
		df_summary$LEN_ACTUAL >= seuil_fiable)
	
	df_summary <- df_summary %>%
		dplyr::mutate(
			Band = dplyr::case_when(
				is.na(Mean_id) ~ "No data",
				Mean_id >= 95  ~ "> 95% (near-identical)",
				Mean_id >= 85  ~ "85–95% (high)",
				Mean_id >= 70  ~ "70–85% (moderate)",
				TRUE           ~ "< 70% (divergent)"
			),
			Band = factor(Band, levels = c("> 95% (near-identical)", "85–95% (high)",
				"70–85% (moderate)", "< 70% (divergent)", "No data"))
		) %>%
		dplyr::arrange(Mean_id)
	df_summary$Sequence <- factor(df_summary$Sequence, levels = df_summary$Sequence)
	
	band_colors <- c(
		"> 95% (near-identical)" = "#1B7F3B",
		"85–95% (high)"          = "#66BB6A",
		"70–85% (moderate)"      = "#F4C300",
		"< 70% (divergent)"      = "#C0392B",
		"No data"                = "grey70"
	)
	
	overall_mean <- mean(mat_no_diag, na.rm = TRUE)
	n_unreliable <- sum(!df_summary$Fiable, na.rm = TRUE)
	
	p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = Sequence, y = Mean_id, color = Band)) +
		{ if (is.finite(overall_mean))
			ggplot2::geom_hline(yintercept = overall_mean, linetype = "dashed",
				color = "#0B2A5B", linewidth = 0.6)
		} +
		ggplot2::geom_segment(ggplot2::aes(x = Sequence, xend = Sequence,
			y = Min_id, yend = Max_id), linewidth = 1, alpha = 0.35, na.rm = TRUE) +
		ggplot2::geom_point(ggplot2::aes(
			shape = ifelse(Fiable, "Reliable", "Short fragment (unreliable)")),
			size = 3.2, na.rm = TRUE) +
		ggplot2::scale_shape_manual(
			values = c("Reliable" = 16, "Short fragment (unreliable)" = 4),
			name   = paste0("Compared length (\u2265 ", seuil_fiable, " nt)")
		) +
		ggplot2::scale_color_manual(values = band_colors, name = "Mean identity band") +
		ggplot2::coord_flip() +
		ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
		ggplot2::labs(
			title    = paste0("Per-sequence identity profile — ", allele_sel),
			subtitle = paste0("Last ", n_last, " nt | ", n_seq,
				" sequences | dot = mean identity vs others, line = min\u2013max",
				if (n_unreliable > 0)
					paste0(" | \u26a0 ", n_unreliable, " sequence(s) compared on < ",
						seuil_fiable, " nt")
				else ""),
			x = NULL, y = "Identity vs other sequences (%)"
		) +
		theme_imgt(base_size = 12) +
		ggplot2::theme(legend.position = "bottom")
	
	h_px <- max(380, n_seq * 32 + 140)
	if (interactive) p <- plotly::ggplotly(p, tooltip = c("x", "y", "colour"),
		height = min(h_px, 900))
	return(p)
}

