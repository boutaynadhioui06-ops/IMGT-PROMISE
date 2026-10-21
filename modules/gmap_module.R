# ═══════════════════════════════════════════════════════════════
#  gmap_module
# ═══════════════════════════════════════════════════════════════

# Dépend de : genomic_map_data.R, utils_core.R
# Packages  : shiny, ggplot2, ggseqlogo, patchwork, purrr, dplyr, DT

.motif_color <- function(motif_names) {
	cols <- MOTIF_COLORS[motif_names]
	cols[is.na(cols)] <- "#94A3B8"   # gris pour tout motif non listé
	unname(cols)
}


# ══════════════════════════════════════════════════════════════
# HELPERS — construction des logos
# ══════════════════════════════════════════════════════════════

# Filtre un vecteur de séquences : ne garde que ACGT
.build_logo_seqs <- function(seq_vec) {
	seq_vec <- toupper(as.character(seq_vec))
	seq_vec <- seq_vec[!is.na(seq_vec) & nchar(seq_vec) > 0]
	seq_vec <- seq_vec[grepl("^[ACGT]+$", seq_vec)]
	if (length(seq_vec) == 0) return(NULL)
	len_tab   <- table(nchar(seq_vec))
	modal_len <- as.integer(names(len_tab)[which.max(len_tab)])
	seq_vec[nchar(seq_vec) == modal_len]
}


.make_logo_grob <- function(seqs, min_n = 1) {
	if (is.null(seqs) || length(seqs) < min_n) return(NULL)
	p <- tryCatch(
		
		ggseqlogo::ggseqlogo(seqs, method = "prob", seq_type = "dna") +
			ggplot2::theme_void() +
			ggplot2::theme(
				legend.position = "none",
				plot.margin     = ggplot2::margin(0, 0, 0, 0)
			),
		error = function(e) NULL
	)
	if (is.null(p)) return(NULL)
	tryCatch(ggplot2::ggplotGrob(p), error = function(e) NULL)
}

.build_seq_lookup <- function(map_data, df_scope) {
	all_seq_cols <- names(df_scope)[stringr::str_starts(names(df_scope), "SEQ_")]
	if (length(all_seq_cols) == 0) return(list())
	
	motifs_from_df <- stringr::str_remove(all_seq_cols, "^SEQ_")
		seq_long <- df_scope %>%
		dplyr::select(ALLELE, dplyr::all_of(all_seq_cols)) %>%
		tidyr::pivot_longer(-ALLELE, names_to = "Motif_col", values_to = "Seq") %>%
		dplyr::mutate(Motif = stringr::str_remove(Motif_col, "^SEQ_")) %>%
		dplyr::select(-Motif_col) %>%
		dplyr::filter(!is.na(Seq), nchar(trimws(Seq)) > 0)
	
	# Joint avec map_data pour récupérer SUBGROUP par allèle
	allele_sg <- map_data %>%
		dplyr::select(ALLELE, SUBGROUP) %>%
		dplyr::distinct()
	
	joined <- seq_long %>%
		dplyr::left_join(allele_sg, by = "ALLELE")
	
	if (nrow(joined) == 0) return(list())
	keys_allele <- paste(trimws(joined$ALLELE), joined$Motif, sep = "__ALLELE__")
	raw_allele  <- split(joined$Seq, keys_allele)
	
	# Clés par sous-groupe
	joined_sg <- joined %>% dplyr::filter(!is.na(SUBGROUP))
	keys_sg   <- paste(trimws(joined_sg$SUBGROUP), joined_sg$Motif, sep = "__")
	raw_sg    <- split(joined_sg$Seq, keys_sg)
	
	# Clés globales (
	keys_all <- paste0("__", joined$Motif)
	raw_all  <- split(joined$Seq, keys_all)
	
	result <- c(raw_allele, raw_sg, raw_all)
	lapply(result, .build_logo_seqs)
}


# ══════════════════════════════════════════════════════════════
#  HELPER — distances inter-motifs 
# ══════════════════════════════════════════════════════════════

.compute_motif_gaps <- function(sub_long, motif_order, dist_lpart = list()) {
	if (length(motif_order) < 2) return(data.frame())
	
	wide <- sub_long %>%
		dplyr::filter(Motif %in% motif_order, Present) %>%
		dplyr::select(ALLELE, Motif, Start_rel, End_rel) %>%
		dplyr::distinct(ALLELE, Motif, .keep_all = TRUE)
	
	purrr::map_dfr(seq_len(length(motif_order) - 1), function(i) {
		m1 <- motif_order[i]; m2 <- motif_order[i + 1]
		
		if (m2 %in% c("LPART_1", "LPART1")) {
			dval <- dist_lpart[[m1]]
			if (!is.null(dval) && is.finite(dval)) {
				return(data.frame(Motif1 = m1, Motif2 = m2,
					gap_min = round(dval), gap_max = round(dval), n = NA_integer_))
			}
			return(NULL)
		}
		
		d1 <- wide %>% dplyr::filter(Motif == m1) %>% dplyr::select(ALLELE, End1 = End_rel)
		d2 <- wide %>% dplyr::filter(Motif == m2) %>% dplyr::select(ALLELE, Start2 = Start_rel)
		mg <- dplyr::inner_join(d1, d2, by = "ALLELE") %>%
			dplyr::mutate(gap = round(Start2 - End1)) %>%
			dplyr::filter(!is.na(gap))
		
		if (nrow(mg) == 0) {
			
			if (nrow(d1) > 0 && nrow(d2) > 0) {
				gap_est <- round(median(d2$Start2, na.rm = TRUE) - median(d1$End1, na.rm = TRUE))
				if (is.finite(gap_est))
					return(data.frame(Motif1 = m1, Motif2 = m2,
						gap_min = gap_est, gap_max = gap_est, n = 0L))
			}
			return(NULL)
		}
		data.frame(Motif1 = m1, Motif2 = m2,
			gap_min = min(mg$gap), gap_max = max(mg$gap), n = nrow(mg))
	})
}


# ══════════════════════════════════════════════════════════════
# Positionnement des motifs
# ══════════════════════════════════════════════════════════════

.compute_slot_positions <- function(motif_order, med_pos, scale_mode = "schematic") {
	n <- length(motif_order)
	if (n == 0) return(numeric(0))
	
	if (n == 1 || !identical(scale_mode, "real")) {
		return(stats::setNames(seq_len(n), motif_order))
	}
	
	is_lpart <- motif_order %in% c("LPART_1", "LPART1")
	
	lpart_med_start <- {
		v <- med_pos$med_start[med_pos$Motif %in% c("LPART_1", "LPART1")]
		if (length(v) > 0 && is.finite(v[1])) v[1] else NA_real_
	}
	
	raw <- vapply(motif_order, function(m) {
		if (m %in% c("LPART_1", "LPART1")) return(0)
		
		row <- med_pos[med_pos$Motif == m, , drop = FALSE]
		d   <- if (nrow(row) > 0) row$order_key[1] else NA_real_
		
		#  distance réelle nt 
		if (is.finite(d) && d > 0) return(d)
		
		# écart de position médiane vs LPART1
		ms <- if (nrow(row) > 0) row$med_start[1] else NA_real_
		if (is.finite(lpart_med_start) && is.finite(ms)) {
			d2 <- lpart_med_start - ms
			if (is.finite(d2) && d2 > 0) return(d2)
		}
		
		NA_real_
	}, numeric(1))
	
	finite_others <- raw[!is_lpart & is.finite(raw)]
	default_raw   <- if (length(finite_others) > 0) max(finite_others) + 100 else 100
	raw[!is.finite(raw)] <- default_raw
	raw[raw < 0] <- 0
	
	
	compressed <- sqrt(raw) / 2
	
	x <- max(compressed) - compressed + 1
	
	for (i in seq_len(n)[-1]) {
		if (x[i] - x[i - 1] < 1) x[i] <- x[i - 1] + 1
	}
	
	stats::setNames(x, motif_order)
}


# ══════════════════════════════════════════════════════════════
#  PANEL — un sous-groupe (ou  allèle) = un schéma
# ══════════════════════════════════════════════════════════════


NO_LOGO_MOTIFS <- c("LPART_1")

.plot_subgroup_panel <- function(sub_long,
	subgroup_name,               # libellé AFFICHÉ 
	seq_lookup,
	bin_motifs  = character(0),  # motifs détectés
	bin_pos     = list(),        
	bin_npres   = list(),        # n allèles présents par motif 
	dist_lpart  = list(),        # médiane de distance réelle à LPART1
	show_logos  = TRUE,
	show_dist   = TRUE,
	min_logo_n  = 1,
	label_size  = 3.0,
	scale_mode    = "schematic", 
	sg_raw        = NULL,        
	allele_filter = NULL) {       
	
	sg_key <- if (!is.null(sg_raw)) sg_raw else subgroup_name
	
	motifs_from_map <- sub_long %>%
		dplyr::filter(Motif != "5'UTR", Present) %>%
		dplyr::pull(Motif) %>% unique()
	
	motifs_present <- union(motifs_from_map, bin_motifs)
	
	n_tot <- dplyr::n_distinct(sub_long$ALLELE)
	
	# ── Construction de la table des slots
	med_map <- sub_long %>%
		dplyr::filter(Motif %in% motifs_from_map, Present) %>%
		dplyr::group_by(Motif) %>%
		dplyr::summarise(
			med_start = median(Start_rel, na.rm = TRUE),
			n_pres    = dplyr::n_distinct(ALLELE),
			.groups = "drop")
	
	extra_rows <- purrr::map_dfr(bin_motifs, function(em) {
		data.frame(
			Motif     = em,
			med_start = if (!is.null(bin_pos[[em]])) bin_pos[[em]] else Inf,
			n_pres    = if (!is.null(bin_npres[[em]])) bin_npres[[em]] else 0L
		)
	})
	med_pos <- dplyr::bind_rows(med_map, extra_rows)
	
	med_pos$order_key <- vapply(seq_len(nrow(med_pos)), function(i) {
		m <- med_pos$Motif[i]
		if (m %in% c("LPART_1", "LPART1")) return(-Inf)
		d <- dist_lpart[[m]]
		if (!is.null(d) && is.finite(d)) return(d)
		ms <- med_pos$med_start[i]
		if (is.finite(ms)) return(-ms)
		-Inf
	}, numeric(1))
	
	med_pos <- med_pos %>% dplyr::arrange(dplyr::desc(order_key))
	
	if (nrow(med_pos) == 0) {
		return(
			ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = paste0(subgroup_name, " — no motif detected (n=", n_tot, " alleles)"),
					size = 4, color = "grey50") +
				ggplot2::labs(title = subgroup_name) +
				ggplot2::theme_void() +
				ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12, color = "#0B2A5B"))
		)
	}
	
	motif_order <- med_pos$Motif
	
	x_pos <- .compute_slot_positions(motif_order, med_pos, scale_mode)
	
	slots <- data.frame(
		Motif = motif_order,
		x     = unname(x_pos[motif_order]),
		pct   = round(med_pos$n_pres / n_tot * 100, 0),
		n_p   = med_pos$n_pres
	)
	
	hw  <- 0.34
	pal <- setNames(.motif_color(motif_order), motif_order)
	
	p <- ggplot2::ggplot()
	
	# Backbone
	p <- p + ggplot2::geom_segment(
		ggplot2::aes(x = 0.5, xend = max(slots$x) + 0.5, y = 0, yend = 0),
		color = "#334155", linewidth = 0.9)
	
	# Blocs motifs
	p <- p + ggplot2::geom_rect(
		data = slots,
		ggplot2::aes(xmin = x - hw, xmax = x + hw, ymin = -0.16, ymax = 0.16, fill = Motif),
		color = "white", linewidth = 0.5) +
		ggplot2::scale_fill_manual(values = pal, guide = "none")
	
	# Étiquettes motif + % présence
	p <- p +
		ggplot2::geom_text(
			data = slots,
			ggplot2::aes(x = x, y = -0.32, label = Motif),
			size = label_size, fontface = "bold", color = "#1F2937") +
		ggplot2::geom_text(
			data = slots,
			ggplot2::aes(x = x, y = -0.46, label = paste0(pct, "% (n=", n_p, "/", n_tot, ")")),
			size = label_size * 0.72, color = "#6B7280")
	
	if (show_logos) {
		for (i in seq_len(nrow(slots))) {
			if (slots$Motif[i] %in% NO_LOGO_MOTIFS) next   # pas de logo pour LPART_1
			
			seqs <- NULL
			
			# priorité absolue à l'allèle sélectionné
			if (!is.null(allele_filter) && allele_filter != "") {
				key_allele <- paste(trimws(allele_filter), slots$Motif[i], sep = "__ALLELE__")
				seqs <- seq_lookup[[key_allele]]
			}
			
			#  repli : consensus du VRAI sous-groupe 
			if (is.null(seqs) || length(seqs) < min_logo_n) {
				key_sg <- paste(trimws(sg_key), slots$Motif[i], sep = "__")
				seqs   <- seq_lookup[[key_sg]]
			}
			
			# dernier repli : motif dans le périmètre déjà filtré
			if (is.null(seqs) || length(seqs) < min_logo_n) {
				key_all <- paste0("__", slots$Motif[i])
				seqs    <- seq_lookup[[key_all]]
			}
			
			g <- .make_logo_grob(seqs, min_n = min_logo_n)
			if (!is.null(g)) {
				p <- p + ggplot2::annotation_custom(
					g,
					xmin = slots$x[i] - hw, xmax = slots$x[i] + hw,
					ymin = 0.22,            ymax = 0.95)
			}
		}
	}
	
	# Distances inter-motifs 
	if (show_dist && nrow(slots) >= 2) {
		gaps <- .compute_motif_gaps(sub_long, motif_order, dist_lpart)
		if (nrow(gaps) > 0) {		gaps <- gaps %>%
			dplyr::mutate(
				x1  = slots$x[match(Motif1, slots$Motif)] + hw,
				x2  = slots$x[match(Motif2, slots$Motif)] - hw,
				xm  = (x1 + x2) / 2,
				lbl = ifelse(gap_min == gap_max,
					ifelse(!is.na(n) & n == 0,
						paste0("~", gap_min, " nt"),
						paste0(gap_min, " nt")),
					paste0("[", gap_min, "-", gap_max, "] nt"))
			)
		p <- p +
			ggplot2::geom_segment(
				data = gaps,
				ggplot2::aes(x = x1, xend = x2, y = -0.62, yend = -0.62),
				color = "#94A3B8", linewidth = 0.4,
				arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "inches"),
					ends = "both", type = "open")) +
			ggplot2::geom_text(
				data = gaps,
				ggplot2::aes(x = xm, y = -0.74, label = lbl),
				size = label_size * 0.68, color = "#475569")
		}
	}
	
	subtitle_txt <- paste0(n_tot, " allele(s) | ", length(motifs_present), " motif(s) detected")
	if (identical(scale_mode, "real"))
		subtitle_txt <- paste0(subtitle_txt, " | scale: \u221ant to LPART1 (compressed, real nt distances)")
	
	p <- p +
		ggplot2::coord_cartesian(
			xlim   = c(0.3, max(slots$x) + 0.7),
			ylim   = c(-0.95, 1.0),
			expand = FALSE) +
		ggplot2::labs(
			title    = subgroup_name,
			subtitle = subtitle_txt
		) +
		ggplot2::theme_void() +
		ggplot2::theme(
			plot.title    = ggplot2::element_text(face = "bold", size = 12.5, color = "#0B2A5B",
				hjust = 0, margin = ggplot2::margin(b = 2)),
			plot.subtitle = ggplot2::element_text(size = 9, color = "#6B7280",
				hjust = 0, margin = ggplot2::margin(b = 6)),
			plot.margin   = ggplot2::margin(8, 12, 6, 12)
		)
	
	p
}


# ══════════════════════════════════════════════════════════════
# Allèles SANS aucun motif régulateur (LPART1 exclu)
# ══════════════════════════════════════════════════════════════

.alleles_without_motif <- function(df) {
	bin_cols <- names(df)[stringr::str_starts(names(df), "BIN_")]
	bin_cols <- setdiff(bin_cols, c("BIN_LPART_1", "BIN_LPART1"))
	if (length(bin_cols) == 0 || nrow(df) == 0) return(character(0))
	has_motif <- Reduce(`|`, lapply(bin_cols, function(cc) {
		v <- df[[cc]]
		!is.na(v) & v == 1
	}))
	unique(df$ALLELE[!has_motif])
}


# ══════════════════════════════════════════════════════════════
# FONCTION PRINCIPAL
# ══════════════════════════════════════════════════════════════

plot_genomic_map_logos <- function(df_final,
	type_data   = "IG",
	species     = NULL,
	subgroup    = NULL,
	func_filter = NULL,
	allele      = NULL,     
	scale_mode  = "schematic", 
	show_logos  = TRUE,
	show_dist   = TRUE,
	max_seqs    = 300,
	max_panels  = 40) {
	
	if (is.null(species) || species == "")
		stop("Please select a species.")
	
	df_use <- df_final
	if (!is.null(func_filter) && func_filter != "")
		df_use <- df_use %>% dplyr::filter(FUNCTIONALITY == func_filter)
	
	map_data <- prepare_genomic_map(df_use,
		type_data = type_data,
		species   = species,
		locus     = NULL,
		subgroup  = if (is.null(subgroup) || subgroup == "") NULL else subgroup,
		max_seqs  = max_seqs)
	
	if (is.null(map_data) || nrow(map_data) == 0)
		stop("No data for these filters.")
	
	
	df_scope_motif <- df_final %>% dplyr::filter(SPECIES == species)
	if (!is.null(subgroup) && subgroup != "")
		df_scope_motif <- df_scope_motif %>% dplyr::filter(SUBGROUP == subgroup)
	if (!is.null(func_filter) && func_filter != "")
		df_scope_motif <- df_scope_motif %>% dplyr::filter(FUNCTIONALITY == func_filter)
	
	alleles_no_motif   <- .alleles_without_motif(df_scope_motif)
	alleles_with_motif <- setdiff(unique(df_scope_motif$ALLELE), alleles_no_motif)
	
	no_motif_table <- df_scope_motif %>%
		dplyr::filter(ALLELE %in% alleles_no_motif) %>%
		dplyr::select(dplyr::any_of(c("ALLELE", "SPECIES", "LOCI", "SUBGROUP",
			"FUNCTIONALITY", "LENGTH_NT"))) %>%
		dplyr::distinct()
	
	
	if (!is.null(allele) && allele != "") {
		map_data <- map_data %>% dplyr::filter(ALLELE == allele)
		if (is.null(map_data) || nrow(map_data) == 0)
			stop("No data available for this allele (check the species/subgroup/functionality filters).")
	} else {
	
		if (length(alleles_with_motif) == 0)
			stop("No allele with at least one regulatory motif (LPART1 excluded) for these filters.")
		map_data <- map_data %>% dplyr::filter(ALLELE %in% alleles_with_motif)
	}
	
	
	df_filtered <- df_final %>%
		dplyr::filter(SPECIES == species)
	if (!is.null(func_filter) && func_filter != "")
		df_filtered <- df_filtered %>% dplyr::filter(FUNCTIONALITY == func_filter)
	if (!is.null(allele) && allele != "") {           
		df_filtered <- df_filtered %>% dplyr::filter(ALLELE == allele)
	} else {

		df_filtered <- df_filtered %>% dplyr::filter(ALLELE %in% alleles_with_motif)
	}
	
	seq_lookup <- .build_seq_lookup(map_data, df_filtered)
	
	bin_cols <- names(df_final)[stringr::str_starts(names(df_final), "BIN_")]
	pos_cols <- names(df_final)[stringr::str_starts(names(df_final), "POS_START_")]
	motifs_in_map <- setdiff(unique(map_data$Motif), "5'UTR")
	extra_bin_motifs <- setdiff(
		stringr::str_remove(bin_cols, "^BIN_"),
		motifs_in_map
	)
	
	subgroups_to_plot <- sort(unique(map_data$SUBGROUP))
	if (length(subgroups_to_plot) > max_panels) {
		message(paste0(
			"plot_genomic_map_logos: ", length(subgroups_to_plot),
			" subgroups detected — display limited to the first ", max_panels,
			" (filter by subgroup to see all)."))
		subgroups_to_plot <- subgroups_to_plot[seq_len(max_panels)]
	}
	
	panels <- lapply(subgroups_to_plot, function(sg) {
		sub_long <- map_data %>% dplyr::filter(SUBGROUP == sg)
		
		sg_alleles <- df_filtered %>%
			dplyr::filter(SUBGROUP == sg) %>%
			dplyr::pull(ALLELE)
		
		df_sg <- df_filtered %>% dplyr::filter(ALLELE %in% sg_alleles)
		
		bin_motifs_sg <- character(0)
		bin_pos_sg    <- list()
		bin_npres_sg  <- list()
		
		for (em in extra_bin_motifs) {
			bin_col <- paste0("BIN_", em)
			pos_col <- paste0("POS_START_", em)
			
			if (!bin_col %in% names(df_sg)) next
			n_pres <- sum(df_sg[[bin_col]] == 1, na.rm = TRUE)
			if (n_pres == 0) next
			
			med_start <- if (pos_col %in% names(df_sg)) {
				median(df_sg[[pos_col]][df_sg[[bin_col]] == 1], na.rm = TRUE)
			} else {
				Inf   
			}
			
			bin_motifs_sg          <- c(bin_motifs_sg, em)
			bin_pos_sg[[em]]       <- med_start
			bin_npres_sg[[em]]     <- n_pres
		}
		
		motifs_all_sg <- union(
			sub_long %>% dplyr::filter(Motif != "5'UTR", Present) %>%
				dplyr::pull(Motif) %>% unique(),
			bin_motifs_sg
		)
		dist_lpart_sg <- list()
		for (mo in motifs_all_sg) {
			if (mo %in% c("LPART_1", "LPART1")) next
			dist_col <- paste0("DIST_", mo)
			if (!dist_col %in% names(df_sg)) next
			dval <- suppressWarnings(median(df_sg[[dist_col]], na.rm = TRUE))
			if (is.finite(dval)) dist_lpart_sg[[mo]] <- dval
		}
		
		
		panel_label <- if (!is.null(allele) && allele != "")
			paste0(sg, " \u2022 ", allele) else sg
		
		.plot_subgroup_panel(sub_long, panel_label, seq_lookup,
			bin_motifs  = bin_motifs_sg,
			bin_pos     = bin_pos_sg,
			bin_npres   = bin_npres_sg,
			dist_lpart  = dist_lpart_sg,
			show_logos  = show_logos,
			show_dist   = show_dist,
			scale_mode  = scale_mode,
			sg_raw        = sg,      
			allele_filter = allele)  
	})
	
	combined <- patchwork::wrap_plots(panels, ncol = 1)
	attr(combined, "n_panels")            <- length(panels)
	attr(combined, "no_motif_table")      <- no_motif_table   
	attr(combined, "n_alleles_no_motif")  <- length(alleles_no_motif)   
	attr(combined, "n_alleles_with_motif")<- length(alleles_with_motif) 
	combined
}


# ══════════════════════════════════════════════════════════════
#  TABLEAU récapitulatif 
# ══════════════════════════════════════════════════════════════

summarize_motif_table <- function(df_final,
	type_data   = "IG",
	species     = NULL,
	subgroup    = NULL,
	func_filter = NULL,
	allele      = NULL,     
	max_seqs    = 300) {
	
	df_use <- df_final
	if (!is.null(func_filter) && func_filter != "")
		df_use <- df_use %>% dplyr::filter(FUNCTIONALITY == func_filter)
	
	map_data <- prepare_genomic_map(df_use,
		type_data = type_data,
		species   = species,
		locus     = NULL,
		subgroup  = if (is.null(subgroup) || subgroup == "") NULL else subgroup,
		max_seqs  = max_seqs)
	
	if (!is.null(allele) && allele != "") {           
		map_data <- map_data %>% dplyr::filter(ALLELE == allele)
	} else {
	.
		df_scope_motif <- df_final %>% dplyr::filter(SPECIES == species)
		if (!is.null(subgroup) && subgroup != "")
			df_scope_motif <- df_scope_motif %>% dplyr::filter(SUBGROUP == subgroup)
		if (!is.null(func_filter) && func_filter != "")
			df_scope_motif <- df_scope_motif %>% dplyr::filter(FUNCTIONALITY == func_filter)
		alleles_with_motif <- setdiff(unique(df_scope_motif$ALLELE),
			.alleles_without_motif(df_scope_motif))
		map_data <- map_data %>% dplyr::filter(ALLELE %in% alleles_with_motif)
	}
	
	map_data %>%
		dplyr::filter(Motif != "5'UTR") %>%
		dplyr::group_by(SUBGROUP, LOCI, Motif) %>%
		dplyr::summarise(
			N_alleles_total = dplyr::n_distinct(ALLELE),
			N_present       = sum(Present, na.rm = TRUE),
			Pct_presence    = round(100 * N_present / N_alleles_total, 1),
			Median_Start    = suppressWarnings(round(median(Start_rel[Present], na.rm = TRUE))),
			Median_End      = suppressWarnings(round(median(End_rel[Present],   na.rm = TRUE))),
			Median_Length   = suppressWarnings(round(median(Length[Present],    na.rm = TRUE))),
			.groups = "drop") %>%
		dplyr::arrange(SUBGROUP, Median_Start)
}


# ══════════════════════════════════════════════════════════════
#  UI — mod_gmap_ui
# ══════════════════════════════════════════════════════════════

mod_gmap_ui <- function(id) {
	ns <- shiny::NS(id)
	shiny::tagList(
		
		shiny::div(
			style = "background: linear-gradient(90deg,#EFF6FF,#F0FDF4);
               border-radius:6px; padding:10px 14px; font-size:12px;
               color:#374151; margin-bottom:12px;"
		),
		
		shiny::fluidRow(
			
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("filter"), shiny::strong("Filters")),
					shiny::selectInput(ns("gmap_species"), "Species",
						choices = NULL, selected = NULL),
					shiny::selectInput(ns("gmap_subgroup"), "Subgroup",
						choices = c("All subgroups" = ""), selected = ""),
					
			
					shiny::conditionalPanel(
						condition = paste0("input['", ns("gmap_subgroup"), "'] != ''"),
						shiny::selectInput(ns("gmap_allele"),
							"Allele (optional — individual map)",
							choices = c("All alleles (subgroup consensus)" = ""),
							selected = ""),
						shiny::tags$small(shiny::icon("info-circle"),
							" Select one allele to display its own (real) motif positions",
							" instead of the subgroup consensus.")
					),
					
					shiny::selectInput(ns("gmap_func"), "Functionality (optional)",
						choices = c("All" = "", "F" = "F", "P" = "P", "ORF" = "ORF"),
						selected = "")
				)
			),
			
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("sliders-h"), shiny::strong("Display options")),
					shiny::checkboxInput(ns("gmap_show_logos"),
						"Show consensus sequence logos", value = TRUE),
					shiny::checkboxInput(ns("gmap_show_dist"),
						"Show inter-motif distances (nt)", value = TRUE),
										shiny::radioButtons(ns("gmap_scale_mode"), "Scale",
						choices  = c("Schematic (equal spacing)" = "schematic",
							"To scale (\u221ant distance to LPART1)" = "real"),
						selected = "schematic"),
					shiny::tags$small(shiny::icon("lightbulb"),
						" 'To scale' spaces motifs according to their real nt distance to LPART1",
						" Exact nt distances are always labeled between blocks.")
				)
			),
			
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("chart-bar"), shiny::strong("Summary")),
					shiny::uiOutput(ns("gmap_quick_stats")),
					shiny::tags$small(shiny::icon("info-circle"),
						" LPART1 is not counted as a regulatory motif.",
						" Alleles below.")
				)
			),
			
			shiny::column(3,
				shiny::div(style = "padding-top:16px;",
					shiny::actionButton(ns("btn_run_gmap"),
						"🗺️ Generate map",
						icon  = shiny::icon("play"),
						class = "btn-success btn-lg",
						style = "width:100%; height:65px; font-size:15px;")
				)
			)
		),
		
		shiny::hr(),
		shiny::uiOutput(ns("gmap_status")),
		
		shiny::tabsetPanel(
			id = ns("gmap_tabs"),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("map"), " Genomic map"),
				shiny::br(),
				shiny::uiOutput(ns("gmap_plot_ui")),
				shiny::br()
			),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("table"), " Summary table"),
				shiny::br(),
				DT::DTOutput(ns("gmap_summary_table")),
				shiny::br()
			),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("ban"), " No-motif alleles"),
				shiny::br(),
				shiny::div(
					style = "background-color:#FFFBEB; border-left:4px solid #D97706;
                     border-radius:6px; padding:10px 14px; margin-bottom:14px; font-size:12.5px;",
					shiny::icon("info-circle"),
					shiny::strong(" LPART1 is not a regulatory motif — "),
					"alleles listed below carry NONE of the regulatory motifs ",
					"(TATA_BOX, HEPTA_PENTA, OCTA_DECA, E_BOX, TATA_BOX2, EBOX_2). ",
					"They are automatically excluded from the subgroup-consensus genomic map ",
					"(percentages, median positions, logos) above, but each can still be viewed ",
					"individually via the ", shiny::strong("Allele"), " filter."
				),
				shiny::uiOutput(ns("gmap_no_motif_status")),
				shiny::br(),
				DT::DTOutput(ns("gmap_no_motif_table")),
				shiny::br(),
				shiny::downloadButton(ns("dl_gmap_no_motif_csv"),
					"📥 CSV — Alleles without motif",
					class = "btn-warning")
			),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("download"), " Exports"),
				shiny::br(),
				shiny::fluidRow(
					shiny::column(4,
						shiny::div(
							style = "text-align:center; padding:14px; background:#EFF6FF;
                       border-radius:8px; border:1px solid #93C5FD;",
							shiny::icon("image",
								style = "font-size:24px; color:#1D4ED8; display:block; margin-bottom:6px;"),
							shiny::p(style = "font-size:11px; color:#374151; margin:0 0 8px 0;",
								"High-resolution map (PNG 300 dpi)"),
							shiny::downloadButton(ns("dl_gmap_png"),
								"📥 PNG — Map",
								class = "btn-primary btn-block", style = "width:100%;")
						)
					),
					shiny::column(4,
						shiny::div(
							style = "text-align:center; padding:14px; background:#F5F3FF;
                       border-radius:8px; border:1px solid #C4B5FD;",
							shiny::icon("file-pdf",
								style = "font-size:24px; color:#7C3AED; display:block; margin-bottom:6px;"),
							shiny::p(style = "font-size:11px; color:#374151; margin:0 0 8px 0;",
								"Vector map (PDF, publication)"),
							shiny::downloadButton(ns("dl_gmap_pdf"),
								"📄 PDF — Map",
								class = "btn-info btn-block", style = "width:100%;")
						)
					),
					shiny::column(4,
						shiny::div(
							style = "text-align:center; padding:14px; background:#FFFBEB;
                       border-radius:8px; border:1px solid #FCD34D;",
							shiny::icon("file-csv",
								style = "font-size:24px; color:#D97706; display:block; margin-bottom:6px;"),
							shiny::p(style = "font-size:11px; color:#374151; margin:0 0 8px 0;",
								"Summary table (CSV)"),
							shiny::downloadButton(ns("dl_gmap_csv"),
								"📥 CSV — Summary",
								class = "btn-warning btn-block", style = "width:100%;")
						)
					)
				)
			)
		)
	)
}


# ══════════════════════════════════════════════════════════════
# SERVER — mod_gmap_server
# ══════════════════════════════════════════════════════════════

mod_gmap_server <- function(id, df_final) {
	shiny::moduleServer(id, function(input, output, session) {
		
		ns <- session$ns
		gmap_store <- shiny::reactiveValues(
			plot = NULL, summary = NULL, done = FALSE, n_panels = 1,
			allele_label = NULL,   
			no_motif_table = NULL, n_no_motif = 0, n_with_motif = 0   
		)
		
		#Sélecteur espèce
		shiny::observe({
			shiny::req(df_final())
			df      <- df_final() %>% dplyr::filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
			sp_list <- sort(unique(df$SPECIES))
			shiny::isolate({
				shiny::updateSelectInput(session, "gmap_species",
					choices  = setNames(sp_list, sp_list),
					selected = sp_list[1])
			})
		})
		
		# Sélecteur sous-groupe dépend de l'espèc
		shiny::observeEvent(input$gmap_species, {
			shiny::req(df_final(), input$gmap_species)
			df <- df_final() %>%
				dplyr::filter(
					LOCI %in% c("IGHV", "IGKV", "IGLV"),
					SPECIES == input$gmap_species)
			sg_list <- sort(unique(df$SUBGROUP))
			shiny::updateSelectInput(session, "gmap_subgroup",
				choices  = c("All subgroups" = "", setNames(sg_list, sg_list)),
				selected = "")
		}, ignoreInit = FALSE)
		
		shiny::observeEvent(list(input$gmap_species, input$gmap_subgroup), {
			shiny::req(df_final(), input$gmap_species)
			
			if (is.null(input$gmap_subgroup) || input$gmap_subgroup == "") {
				shiny::updateSelectInput(session, "gmap_allele",
					choices  = c("All alleles (subgroup consensus)" = ""),
					selected = "")
				return()
			}
			
			df <- df_final() %>%
				dplyr::filter(
					LOCI %in% c("IGHV", "IGKV", "IGLV"),
					SPECIES  == input$gmap_species,
					SUBGROUP == input$gmap_subgroup)
			al_list <- sort(unique(df$ALLELE))
			
			shiny::updateSelectInput(session, "gmap_allele",
				choices  = c("All alleles (subgroup consensus)" = "", setNames(al_list, al_list)),
				selected = "")
		}, ignoreInit = FALSE, ignoreNULL = FALSE)
		
		# Stats rapides
		output$gmap_quick_stats <- shiny::renderUI({
			shiny::req(df_final(), input$gmap_species)
			df <- df_final() %>%
				dplyr::filter(LOCI %in% c("IGHV", "IGKV", "IGLV"), SPECIES == input$gmap_species)
			if (!is.null(input$gmap_subgroup) && input$gmap_subgroup != "")
				df <- df %>% dplyr::filter(SUBGROUP == input$gmap_subgroup)
			if (!is.null(input$gmap_allele) && input$gmap_allele != "")   
				df <- df %>% dplyr::filter(ALLELE == input$gmap_allele)
			
			seq_cols <- names(df)[stringr::str_starts(names(df), "SEQ_")]
			
			n_no_motif_live <- length(.alleles_without_motif(df))
			
			shiny::tags$table(style = "width:100%; font-size:12px;",
				shiny::tags$tr(
					shiny::tags$td(shiny::icon("dna"), " Alleles:"),
					shiny::tags$td(shiny::strong(nrow(df)))),
				shiny::tags$tr(
					shiny::tags$td(shiny::icon("layer-group"), " Subgroups:"),
					shiny::tags$td(shiny::strong(length(unique(df$SUBGROUP))))),
				shiny::tags$tr(
					shiny::tags$td(shiny::icon("puzzle-piece"), " Motifs (sequence available):"),
					shiny::tags$td(shiny::strong(length(seq_cols)))),
				shiny::tags$tr(
					shiny::tags$td(shiny::icon("ban"), " Alleles with NO motif (excl. LPART1):"),
					shiny::tags$td(shiny::strong(
						style = if (n_no_motif_live > 0) "color:#D97706;" else NULL,
						n_no_motif_live)))
			)
		})
		
		# Génération
		shiny::observeEvent(input$btn_run_gmap, {
			shiny::req(df_final(), input$gmap_species)
			shiny::withProgress(message = "Generating genomic map…", value = 0.1, {
				tryCatch({
					fn <- if (is.null(input$gmap_func)     || input$gmap_func     == "") NULL else input$gmap_func
					sg <- if (is.null(input$gmap_subgroup) || input$gmap_subgroup == "") NULL else input$gmap_subgroup
					al <- if (is.null(input$gmap_allele)   || input$gmap_allele   == "") NULL else input$gmap_allele   
					sm <- if (is.null(input$gmap_scale_mode)) "schematic" else input$gmap_scale_mode                   
					
					shiny::incProgress(0.3, detail = "Building logos…")
					p <- plot_genomic_map_logos(df_final(),
						type_data   = "IG",
						species     = input$gmap_species,
						subgroup    = sg,
						func_filter = fn,
						allele      = al,          
						scale_mode  = sm,          
						show_logos  = isTRUE(input$gmap_show_logos),
						show_dist   = isTRUE(input$gmap_show_dist))
					
					shiny::incProgress(0.7, detail = "Summary table…")
					tbl <- summarize_motif_table(df_final(),
						type_data   = "IG",
						species     = input$gmap_species,
						subgroup    = sg,
						func_filter = fn,
						allele      = al)         
					
					shiny::incProgress(1)
					gmap_store$plot         <- p
					gmap_store$summary      <- tbl
					gmap_store$n_panels     <- attr(p, "n_panels")
					if (is.null(gmap_store$n_panels)) gmap_store$n_panels <- 1
					gmap_store$allele_label <- al     
					gmap_store$no_motif_table <- attr(p, "no_motif_table")   
					gmap_store$n_no_motif     <- attr(p, "n_alleles_no_motif") %||% 0   
					gmap_store$n_with_motif   <- attr(p, "n_alleles_with_motif") %||% 0 
					gmap_store$done <- TRUE
					
				}, error = function(e) {
					shiny::showNotification(paste0("❌ Genomic map: ", e$message),
						type = "error", duration = 8)
				})
			})
		})
		
		# Statut
		output$gmap_status <- shiny::renderUI({
			if (!gmap_store$done)
				shiny::div(class = "alert alert-info", shiny::icon("info-circle"),
					" Configure the filters and click ",
					shiny::strong("🗺️ Generate map"), ".")
			else {
				allele_txt <- if (!is.null(gmap_store$allele_label) && gmap_store$allele_label != "")
					paste0(" — individual allele: ", gmap_store$allele_label) else ""
				excl_txt <- if ((is.null(gmap_store$allele_label) || gmap_store$allele_label == "") &&
						isTRUE(gmap_store$n_no_motif > 0))
					paste0(" | ", gmap_store$n_no_motif,
						" allele(s) with no regulatory motif excluded from the consensus (see '🚫 No-motif alleles' tab)")
				else ""
				shiny::div(class = "alert alert-success", shiny::icon("check-circle"),
					shiny::strong(paste0(" Map generated", allele_txt, " — ")),
					gmap_store$n_panels, " panel(s) displayed.",
					shiny::tags$span(style = "color:#B45309;", excl_txt))
			}
		})
		
		# Plot dynamique
		output$gmap_plot_ui <- shiny::renderUI({
			shiny::req(gmap_store$done, gmap_store$plot)
			h <- max(260, gmap_store$n_panels * 230 + 40)
			shiny::plotOutput(ns("gmap_plot"), height = paste0(h, "px"))
		})
		
		output$gmap_plot <- shiny::renderPlot({
			shiny::req(gmap_store$done, gmap_store$plot)
			gmap_store$plot
		}, res = 110)
		
		# Tableau récapitulatif
		output$gmap_summary_table <- DT::renderDT({
			shiny::req(gmap_store$summary)
			DT::datatable(gmap_store$summary,
				options = list(pageLength = 25, scrollX = TRUE,
					dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
				extensions = "Buttons", rownames = FALSE,
				class = "cell-border stripe hover compact", filter = "top") %>%
				DT::formatStyle("Motif", fontWeight = "bold", color = "#1565C0") %>%
				DT::formatStyle("Pct_presence",
					background = DT::styleColorBar(c(0, 100), "#93C5FD"),
					backgroundSize = "100% 80%", backgroundRepeat = "no-repeat",
					backgroundPosition = "center")
		})
		
		output$gmap_no_motif_status <- shiny::renderUI({
			if (!gmap_store$done)
				shiny::div(class = "alert alert-info", shiny::icon("info-circle"),
					" Generate a map first to list the alleles excluded from the consensus.")
			else if (gmap_store$n_no_motif == 0)
				shiny::div(class = "alert alert-success", shiny::icon("check-circle"),
					" No allele without a regulatory motif in this selection.")
			else
				shiny::div(class = "alert alert-warning", shiny::icon("exclamation-triangle"),
					shiny::strong(paste0(" ", gmap_store$n_no_motif, " allele(s) ")),
					"with no regulatory motif (out of ",
					gmap_store$n_no_motif + gmap_store$n_with_motif, " total) — excluded from the map above.")
		})
		
		output$gmap_no_motif_table <- DT::renderDT({
			shiny::req(gmap_store$done)
			df <- gmap_store$no_motif_table
			if (is.null(df) || nrow(df) == 0)
				return(DT::datatable(data.frame(Message = "No allele without a regulatory motif."),
					rownames = FALSE, options = list(dom = "t")))
			DT::datatable(df,
				options = list(pageLength = 25, scrollX = TRUE,
					dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
				extensions = "Buttons", rownames = FALSE,
				class = "cell-border stripe hover compact", filter = "top")
		})
		
		output$dl_gmap_no_motif_csv <- shiny::downloadHandler(
			filename = function() {
				suf <- if (!is.null(input$gmap_subgroup) && input$gmap_subgroup != "")
					paste0("_", input$gmap_subgroup) else ""
				paste0("Alleles_NoMotif_", input$gmap_species, suf, "_", Sys.Date(), ".csv")
			},
			content = function(file) {
				shiny::req(gmap_store$no_motif_table)
				write.csv(gmap_store$no_motif_table, file, row.names = FALSE)
				shiny::showNotification("✅ CSV downloaded!", type = "message", duration = 3)
			})
		
		#Downloads
		output$dl_gmap_png <- shiny::downloadHandler(
			filename = function() {
				suf <- if (!is.null(input$gmap_allele) && input$gmap_allele != "")
					paste0("_", input$gmap_allele) else ""   
				paste0("GenomicMap_Logos_", input$gmap_species, suf, "_", Sys.Date(), ".png")
			},
			content = function(file) {
				shiny::req(gmap_store$plot)
				h_in <- max(4, gmap_store$n_panels * 3.0)
				ggplot2::ggsave(file, plot = gmap_store$plot,
					width = 11, height = h_in, dpi = 300, bg = "white")
				shiny::showNotification("✅ PNG downloaded!", type = "message", duration = 3)
			})
		
		output$dl_gmap_pdf <- shiny::downloadHandler(
			filename = function() {
				suf <- if (!is.null(input$gmap_allele) && input$gmap_allele != "")
					paste0("_", input$gmap_allele) else ""   
				paste0("GenomicMap_Logos_", input$gmap_species, suf, "_", Sys.Date(), ".pdf")
			},
			content = function(file) {
				shiny::req(gmap_store$plot)
				h_in <- max(4, gmap_store$n_panels * 3.0)
				ggplot2::ggsave(file, plot = gmap_store$plot,
					width = 11, height = h_in, device = "pdf", bg = "white")
				shiny::showNotification("✅ PDF downloaded!", type = "message", duration = 3)
			})
		
		output$dl_gmap_csv <- shiny::downloadHandler(
			filename = function() {
				suf <- if (!is.null(input$gmap_allele) && input$gmap_allele != "")
					paste0("_", input$gmap_allele) else ""   
				paste0("GenomicMap_Summary_", input$gmap_species, suf, "_", Sys.Date(), ".csv")
			},
			content = function(file) {
				shiny::req(gmap_store$summary)
				write.csv(gmap_store$summary, file, row.names = FALSE)
				shiny::showNotification("✅ CSV downloaded!", type = "message", duration = 3)
			})
	})
}
