# ═══════════════════════════════════════════════════════════════
# IMGT/StatIGTR — Module d'affichage  pour les sous-groupes
# ═══════════════════════════════════════════════════════════════

# Dépend de : utils_core.R
# Packages  : dplyr, ggplot2, plotly, stringr, tidyr, shiny

# ══════════════════════════════════════════════════════════════
#  CONSTANTES
# ══════════════════════════════════════════════════════════════

.SG_SMALL   <- 20
.SG_MEDIUM  <- 40
.SG_LARGE   <- 80

.PX_PER_ROW_SMALL  <- 32
.PX_PER_ROW_MEDIUM <- 26
.PX_PER_ROW_LARGE  <- 20
.PX_PER_ROW_XLARGE <- 16

.EXPORT_DPI <- 300


# ══════════════════════════════════════════════════════════════
#  HELPERS
# ══════════════════════════════════════════════════════════════

sg_text_size <- function(n_sg) {
	dplyr::case_when(
		n_sg <= .SG_SMALL  ~ 9.5,
		n_sg <= .SG_MEDIUM ~ 8.0,
		n_sg <= .SG_LARGE  ~ 6.5,
		TRUE               ~ 5.0
	)
}

sg_plot_height_px <- function(n_sg, n_col = 1, min_h = 400) {
	px_per_row <- dplyr::case_when(
		n_sg <= .SG_SMALL  ~ .PX_PER_ROW_SMALL,
		n_sg <= .SG_MEDIUM ~ .PX_PER_ROW_MEDIUM,
		n_sg <= .SG_LARGE  ~ .PX_PER_ROW_LARGE,
		TRUE               ~ .PX_PER_ROW_XLARGE
	)
	header <- 120
	per_sg <- px_per_row * max(1, n_col)
	max(min_h, header + n_sg * per_sg)
}

sg_export_dims <- function(n_sg, n_col = 1, n_facet = 1) {
	h_body <- n_sg * 0.30 * max(1, n_col) * max(1, n_facet)
	h      <- max(8, h_body + 3)
	label_extra <- dplyr::case_when(
		n_sg <= .SG_SMALL  ~ 3.0,
		n_sg <= .SG_MEDIUM ~ 3.5,
		n_sg <= .SG_LARGE  ~ 4.0,
		TRUE               ~ 4.5
	)
	w <- max(12, 10 + label_extra + n_col * 0.8)
	list(width = min(w, 50), height = min(h, 120))
}

sg_margin_left_px <- function(labels) {
	max_chars <- max(nchar(as.character(labels)), na.rm = TRUE)
	dplyr::case_when(
		max_chars <= 8  ~  80,
		max_chars <= 14 ~ 110,
		max_chars <= 20 ~ 150,
		max_chars <= 28 ~ 190,
		TRUE            ~ 230
	)
}

sg_use_heatmap <- function(n_sg) n_sg > .SG_LARGE

sg_tile_label_size <- function(n_sg) {
	dplyr::case_when(
		n_sg <= 15  ~ 3.5,
		n_sg <= 30  ~ 2.8,
		n_sg <= 50  ~ 2.2,
		n_sg <= 80  ~ 1.8,
		n_sg <= 120 ~ 1.4,
		TRUE        ~ 1.1   
	)
}




# ══════════════════════════════════════════════════════════════
# HEATMAP DE DISTRIBUTION DES SÉQUENCES
# ══════════════════════════════════════════════════════════════

plot_sg_heatmap_seq_count <- function(df, species = "", interactive = TRUE) {
	
	if (nrow(df) == 0 || !("SUBGROUP" %in% names(df)) || !("LOCI" %in% names(df))) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available",
					size = 6, color = "grey50") +
				ggplot2::theme_void())
	}
	
	mat <- df %>%
		dplyr::filter(!is.na(SUBGROUP), !is.na(LOCI)) %>%
		dplyr::group_by(SUBGROUP, LOCI) %>%
		dplyr::summarise(Count = dplyr::n(), .groups = "drop")
	
	if (nrow(mat) == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No data available",
					size = 6, color = "grey50") +
				ggplot2::theme_void())
	}
	
	n_sg   <- length(unique(mat$SUBGROUP))
	n_loci <- length(unique(mat$LOCI))
	
	sg_order <- mat %>%
		dplyr::group_by(SUBGROUP) %>%
		dplyr::summarise(tot = sum(Count), .groups = "drop") %>%
		dplyr::arrange(tot) %>%
		dplyr::pull(SUBGROUP)
	mat$SUBGROUP <- factor(mat$SUBGROUP, levels = sg_order)
	
	txt_size        <- sg_text_size(n_sg)
	tile_label_size <- sg_tile_label_size(n_sg) 
	
	p <- ggplot2::ggplot(mat,
		ggplot2::aes(x = LOCI, y = SUBGROUP, fill = Count)) +
		ggplot2::geom_tile(color = "white", linewidth = 0.5) +
		ggplot2::geom_text(ggplot2::aes(                        
			label = Count,
			color = ifelse(Count <= max(Count, na.rm = TRUE) * 0.6, "white", "grey20")),
			size = tile_label_size, fontface = "bold") +
		ggplot2::scale_color_identity() +
		ggplot2::scale_fill_viridis_c(
			option   = "viridis",
			name     = "Nb sequences",
			na.value = "#F9FAFB"
		) +
		ggplot2::labs(title    = paste0("Sequence distribution by subgroup",
			if (species != "") paste0(" — ", species) else ""),
			subtitle = paste0(n_sg, " subgroups | Intensity = number of sequences"),
			x = "Locus", y = NULL) +
		theme_imgt(base_size = 10) +
		ggplot2::theme(
			axis.text.y       = ggplot2::element_text(size = txt_size),
			axis.text.x       = ggplot2::element_text(size = 10, face = "bold"),
			panel.grid        = ggplot2::element_blank(),
			legend.position   = "right",
			legend.key.height = ggplot2::unit(1.5, "cm")
		)
	
	if (interactive) {
		h  <- sg_plot_height_px(n_sg, min_h = 450)
		ml <- sg_margin_left_px(unique(mat$SUBGROUP))
		p  <- plotly::ggplotly(p, tooltip = c("x", "y", "fill"), height = h) %>%
			plotly::layout(margin = list(l = ml))
	}
	return(p)
}


# Heatmap fonctionnalité × sous-groupes
plot_sg_heatmap_functionality <- function(df, species = "", interactive = TRUE) {
	
	if (nrow(df) == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No data available", size = 6, color = "grey50") +
				ggplot2::theme_void())
	}
	
	data_mat <- df %>%
		dplyr::filter(!is.na(SUBGROUP), !is.na(FUNCTIONALITY)) %>%
		dplyr::group_by(SUBGROUP, FUNCTIONALITY) %>%
		dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
		dplyr::group_by(SUBGROUP) %>%
		dplyr::mutate(Pct = round(n / sum(n) * 100, 1)) %>%
		dplyr::ungroup()
	
	if (nrow(data_mat) == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No data available", size = 6, color = "grey50") +
				ggplot2::theme_void())
	}
	
	n_sg <- length(unique(data_mat$SUBGROUP))
	
	sg_order <- data_mat %>%
		dplyr::filter(FUNCTIONALITY == "F") %>%
		dplyr::arrange(dplyr::desc(Pct)) %>%
		dplyr::pull(SUBGROUP)
	sg_order <- c(sg_order, setdiff(unique(data_mat$SUBGROUP), sg_order))
	data_mat$SUBGROUP <- factor(data_mat$SUBGROUP, levels = rev(sg_order))
	
	txt_size        <- sg_text_size(n_sg)
	tile_label_size <- sg_tile_label_size(n_sg)  
	
	p <- ggplot2::ggplot(data_mat,
		ggplot2::aes(x = FUNCTIONALITY, y = SUBGROUP, fill = Pct)) +
		ggplot2::geom_tile(color = "white", linewidth = 0.4) +
		ggplot2::geom_text(ggplot2::aes(                        
			label = ifelse(Pct > 0, paste0(Pct, "%"), ""),
			color = ifelse(Pct <= 55, "white", "grey20")),
			size = tile_label_size, fontface = "bold") +
		ggplot2::scale_color_identity() +
		ggplot2::scale_fill_viridis_c(
			option = "viridis",
			limits = c(0, 100),
			name   = "Proportion (%)",
			na.value = "#F9FAFB"
		) +
		ggplot2::facet_wrap(~ FUNCTIONALITY, nrow = 1, scales = "free_x") +
		ggplot2::labs(title    = paste0("Functionality by subgroup",
			if (species != "") paste0(" — ", species) else ""),
			subtitle = paste0(n_sg, " subgroups | Columns: F / P / ORF"),
			x = NULL, y = NULL) +
		theme_imgt(base_size = 10) +
		ggplot2::theme(
			axis.text.y    = ggplot2::element_text(size = txt_size),
			axis.text.x    = ggplot2::element_blank(),
			panel.grid     = ggplot2::element_blank(),
			legend.position = "right",
			strip.text     = ggplot2::element_text(size = 11, face = "bold")
		)
	
	if (interactive) {
		h  <- sg_plot_height_px(n_sg, min_h = 450)
		ml <- sg_margin_left_px(unique(data_mat$SUBGROUP))
		p  <- plotly::ggplotly(p, tooltip = c("x", "y", "fill"), height = h) %>%
			plotly::layout(margin = list(l = ml))
	}
	return(p)
}


# ══════════════════════════════════════════════════════════════
# WRAPPER SÉQUENCES
# ══════════════════════════════════════════════════════════════

plot_seq_by_subgroup_smart <- function(df,
	type_data      = "IG",
	species_filter = NULL,
	locus_filter   = NULL,
	interactive    = TRUE,
	top_n          = NULL,
	force_mode     = "auto") {
	
	if (type_data == "IG") df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	else                   df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	
	if (!is.null(species_filter) && species_filter != "")
		df <- df %>% dplyr::filter(SPECIES == species_filter)
	
	if (!is.null(locus_filter) && locus_filter != "" && locus_filter != "All")
		df <- df %>% dplyr::filter(LOCI == locus_filter)
	
	if (nrow(df) == 0 || length(unique(df$SUBGROUP)) == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No data available for this selection.",
					size = 5, color = "grey50") +
				ggplot2::theme_void() +
				ggplot2::labs(title = "Sequence distribution by subgroup"))
	}
	
	seq_sg <- df %>%
		dplyr::group_by(LOCI, SUBGROUP) %>%
		dplyr::summarise(Count = dplyr::n(), .groups = "drop")
	
	n_sg <- length(unique(seq_sg$SUBGROUP))
	
	if (!is.null(top_n) && is.numeric(top_n) && top_n > 0 && top_n < n_sg) {
		top_sgs <- seq_sg %>%
			dplyr::group_by(SUBGROUP) %>%
			dplyr::summarise(tot = sum(Count), .groups = "drop") %>%
			dplyr::arrange(dplyr::desc(tot)) %>%
			dplyr::slice_head(n = top_n) %>%
			dplyr::pull(SUBGROUP)
		seq_sg <- seq_sg %>% dplyr::filter(SUBGROUP %in% top_sgs)
		df     <- df     %>% dplyr::filter(SUBGROUP %in% top_sgs)
		n_sg   <- length(unique(seq_sg$SUBGROUP))
	}
	
	use_hm <- (force_mode == "heatmap") ||
		(force_mode == "auto" && sg_use_heatmap(n_sg))
	
	if (use_hm) {
		return(plot_sg_heatmap_seq_count(df,
			species     = species_filter %||% "",
			interactive = interactive))
	}
	
	txt_size <- sg_text_size(n_sg)
	h_px     <- sg_plot_height_px(n_sg, min_h = 380)
	ml_px    <- sg_margin_left_px(unique(seq_sg$SUBGROUP))
	
	sg_order <- seq_sg %>%
		dplyr::group_by(SUBGROUP) %>%
		dplyr::summarise(tot = sum(Count), .groups = "drop") %>%
		dplyr::arrange(tot) %>%
		dplyr::pull(SUBGROUP)
	seq_sg$SUBGROUP <- factor(seq_sg$SUBGROUP, levels = sg_order)
	
	subtitle_parts <- c()
	if (!is.null(species_filter) && species_filter != "")
		subtitle_parts <- c(subtitle_parts, paste("Species:", species_filter))
	if (!is.null(locus_filter) && locus_filter != "" && locus_filter != "All")
		subtitle_parts <- c(subtitle_parts, paste("Locus:", locus_filter))
	if (!is.null(top_n))
		subtitle_parts <- c(subtitle_parts, paste0("Top-", top_n))
	
	p <- ggplot2::ggplot(seq_sg,
		ggplot2::aes(x = SUBGROUP, y = Count, fill = LOCI)) +
		ggplot2::geom_bar(stat = "identity", width = 0.75) +
		ggplot2::geom_text(ggplot2::aes(label = Count),
			hjust = -0.15, size = 2.8, color = "grey20", fontface = "bold") +
		ggplot2::coord_flip() +
		ggplot2::scale_fill_manual(values = LOCI_COLORS, name = "Locus") +
		ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.22))) +
		ggplot2::labs(title    = paste("Sequence distribution by subgroup —", type_data),
			subtitle = if (length(subtitle_parts) > 0) paste(subtitle_parts, collapse = " | ") else NULL,
			x = NULL, y = "Number of sequences") +
		theme_imgt(base_size = 10) +
		ggplot2::theme(
			axis.text.y     = ggplot2::element_text(size = txt_size),
			legend.position = "right",
			panel.spacing   = ggplot2::unit(0.5, "lines")
		)
	
	if (interactive) {
		p <- plotly::ggplotly(p, tooltip = c("x","y","fill"), height = h_px) %>%
			plotly::layout(margin = list(l = ml_px, r = 90))
	}
	return(p)
}


# ══════════════════════════════════════════════════════════════
#  WRAPPER FONCTIONNALITÉ
# ══════════════════════════════════════════════════════════════

plot_func_by_subgroup_smart <- function(df,
	type_data      = "IG",
	species_filter = NULL,
	locus_filter   = NULL,
	interactive    = TRUE,
	top_n          = NULL,
	force_mode     = "auto") {
	
	if (type_data == "IG") df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	else                   df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	
	if (!is.null(species_filter) && species_filter != "")
		df <- df %>% dplyr::filter(SPECIES == species_filter)
	
	if (!is.null(locus_filter) && locus_filter != "" && locus_filter != "All")
		df <- df %>% dplyr::filter(LOCI == locus_filter)
	
	if (nrow(df) == 0 || length(unique(df$SUBGROUP)) == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No data available for this selection.",
					size = 5, color = "grey50") +
				ggplot2::theme_void() +
				ggplot2::labs(title = "Functionality distribution by subgroup"))
	}
	
	df_func <- df %>% dplyr::filter(!is.na(FUNCTIONALITY))
	if (nrow(df_func) == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No functionality data for this selection.",
					size = 5, color = "grey50") +
				ggplot2::theme_void())
	}
	
	n_sg_all <- length(unique(df$SUBGROUP))
	
	if (!is.null(top_n) && is.numeric(top_n) && top_n > 0 && top_n < n_sg_all) {
		top_sgs <- df %>%
			dplyr::group_by(SUBGROUP) %>%
			dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
			dplyr::arrange(dplyr::desc(n)) %>%
			dplyr::slice_head(n = top_n) %>%
			dplyr::pull(SUBGROUP)
		df <- df %>% dplyr::filter(SUBGROUP %in% top_sgs)
	}
	
	n_sg <- length(unique(df$SUBGROUP))
	
	use_hm <- (force_mode == "heatmap") ||
		(force_mode == "auto" && sg_use_heatmap(n_sg))
	
	if (use_hm) {
		return(plot_sg_heatmap_functionality(df,
			species     = species_filter %||% "",
			interactive = interactive))
	}
	
	txt_size <- sg_text_size(n_sg)
	h_px     <- sg_plot_height_px(n_sg, min_h = 380)
	ml_px    <- sg_margin_left_px(unique(df$SUBGROUP))
	
	data_plot <- df %>%
		dplyr::group_by(SUBGROUP, LOCI, FUNCTIONALITY) %>%
		dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
		dplyr::group_by(SUBGROUP) %>%
		dplyr::mutate(Total = sum(n), Proportion = n / Total) %>%
		dplyr::ungroup() %>%
		dplyr::filter(!is.na(FUNCTIONALITY)) %>%
		dplyr::arrange(LOCI, SUBGROUP)
	
	if (nrow(data_plot) == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No data available.", size = 5, color = "grey50") +
				ggplot2::theme_void())
	}
	
	sg_order <- data_plot %>%
		dplyr::filter(FUNCTIONALITY == "F") %>%
		dplyr::arrange(Proportion) %>%
		dplyr::pull(SUBGROUP) %>%
		unique()
	sg_order <- c(sg_order, setdiff(unique(data_plot$SUBGROUP), sg_order))
	data_plot$SUBGROUP <- factor(data_plot$SUBGROUP, levels = sg_order)
	
	subtitle_parts <- c()
	if (!is.null(species_filter) && species_filter != "")
		subtitle_parts <- c(subtitle_parts, paste("Species:", species_filter))
	if (!is.null(locus_filter) && locus_filter != "" && locus_filter != "All")
		subtitle_parts <- c(subtitle_parts, paste("Locus:", locus_filter))
	if (!is.null(top_n))
		subtitle_parts <- c(subtitle_parts, paste0("Top-", top_n))
	
	p <- ggplot2::ggplot(data_plot,
		ggplot2::aes(x = SUBGROUP, y = Proportion, fill = FUNCTIONALITY)) +
		ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.75) +
		ggplot2::geom_text(
			ggplot2::aes(
				label = ifelse(n >= 2, n, ""),
				color = FUNCTIONALITY),
			position = ggplot2::position_stack(vjust = 0.5),
			size = 2.5, fontface = "bold", show.legend = FALSE) +
		ggplot2::scale_color_manual(
			values = c("F" = "white", "P" = "white", "ORF" = "white"),
			guide  = "none") +
		ggplot2::coord_flip() +
		ggplot2::scale_y_continuous(
			labels = scales::percent_format(accuracy = 1),
			expand = ggplot2::expansion(mult = c(0, 0.04))) +
		ggplot2::scale_fill_manual(
			values = FUNC_COLORS, breaks = c("F","P","ORF"),
			labels = c("Functional(F)", "Pseudogene (P)", "ORF")) +
		ggplot2::labs(
			title    = paste("Functionality by subgroup  —", type_data),
			subtitle = if (length(subtitle_parts) > 0) paste(subtitle_parts, collapse = " | ") else NULL,
			x = NULL, y = "Proportion", fill = "Functionality"
		) +
		theme_imgt(base_size = 10) +
		ggplot2::theme(
			axis.text.y     = ggplot2::element_text(size = txt_size),
			legend.position = "bottom",
			panel.spacing   = ggplot2::unit(0.5, "lines")
		)
	
	if (interactive) {
		p <- plotly::ggplotly(p, tooltip = c("x","y","fill"), height = h_px + 60) %>%
			plotly::layout(
				legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.1),
				margin = list(l = ml_px, b = 80)
			)
	}
	return(p)
}


# ══════════════════════════════════════════════════════════════
# WRAPPER FRÉQUENCES MOTIFS sous-groupe
# ══════════════════════════════════════════════════════════════


plot_motif_freq_subgroup_smart <- function(data,
	species_filter = NULL,
	interactive    = TRUE,
	top_n          = NULL,
	force_mode     = "auto") {
	
	if (!is.null(species_filter) && species_filter != "")
		data <- data %>% dplyr::filter(SPECIES == species_filter)
	
	pct_cols  <- names(data)[stringr::str_starts(names(data), "Pct_")]
	n_sg_all  <- length(unique(data$SUBGROUP))
	n_motifs  <- length(pct_cols)
	motif_pal <- get_motif_palette(stringr::str_remove(pct_cols, "^Pct_"))
	
	if (n_sg_all == 0 || n_motifs == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No data available.", size = 5, color = "grey50") +
				ggplot2::theme_void())
	}
	
	if (!is.null(top_n) && is.numeric(top_n) && top_n > 0 && top_n < n_sg_all) {
		data_long_tmp <- data %>%
			tidyr::pivot_longer(dplyr::all_of(pct_cols), names_to = "Motif", values_to = "Percentage")
		top_sgs <- data_long_tmp %>%
			dplyr::group_by(SUBGROUP) %>%
			dplyr::summarise(m = mean(Percentage, na.rm = TRUE), .groups = "drop") %>%
			dplyr::arrange(dplyr::desc(m)) %>%
			dplyr::slice_head(n = top_n) %>%
			dplyr::pull(SUBGROUP)
		data <- data %>% dplyr::filter(SUBGROUP %in% top_sgs)
	}
	
	n_sg <- length(unique(data$SUBGROUP))
	
	use_hm <- (force_mode == "heatmap") ||
		(force_mode == "auto" && sg_use_heatmap(n_sg))
	
	#  HEATMAP
	if (use_hm) {
		data_long <- data %>%
			tidyr::pivot_longer(dplyr::all_of(pct_cols), names_to = "BIN_col", values_to = "Pct") %>%
			dplyr::mutate(Motif = stringr::str_remove(BIN_col, "^Pct_"))
		
		sg_order <- data_long %>%
			dplyr::group_by(SUBGROUP) %>%
			dplyr::summarise(m = mean(Pct, na.rm = TRUE), .groups = "drop") %>%
			dplyr::arrange(dplyr::desc(m)) %>%
			dplyr::pull(SUBGROUP)
		data_long$SUBGROUP <- factor(data_long$SUBGROUP, levels = rev(sg_order))
		
		txt_size        <- sg_text_size(n_sg)
		tile_label_size <- sg_tile_label_size(n_sg)  
		
		p <- ggplot2::ggplot(data_long,
			ggplot2::aes(x = Motif, y = SUBGROUP, fill = Pct)) +
			ggplot2::geom_tile(color = "white", linewidth = 0.4) +
			ggplot2::geom_text(ggplot2::aes(                      
				label = ifelse(Pct > 0, paste0(Pct, "%"), ""),
				color = ifelse(Pct <= 55, "white", "grey20")),
				size = tile_label_size, fontface = "bold") +
			ggplot2::scale_color_identity() +
			ggplot2::scale_fill_viridis_c(
				option = "viridis",
				limits = c(0, 100),
				name   = "Presence (%)",
				breaks = c(0, 25, 50, 75, 100),
				na.value = "#F9FAFB"
			) +
			ggplot2::labs(title    = paste0("Motif presence frequencies by subgroup",
				if (!is.null(species_filter) && species_filter != "")
					paste0(" — ", species_filter) else ""),
				subtitle = paste0(n_sg, " subgroups × ", n_motifs, " motifs"),
				x = "Motif", y = NULL) +
			theme_imgt(base_size = 10) +
			ggplot2::theme(
				axis.text.y       = ggplot2::element_text(size = txt_size),
				axis.text.x       = ggplot2::element_text(angle = 40, hjust = 1),
				panel.grid        = ggplot2::element_blank(),
				legend.position   = "right",
				legend.key.height = ggplot2::unit(1.5, "cm")
			)
		
		if (interactive) {
			h  <- sg_plot_height_px(n_sg, min_h = 450)
			ml <- sg_margin_left_px(unique(data_long$SUBGROUP))
			p  <- plotly::ggplotly(p, tooltip = c("x","y","fill"), height = h) %>%
				plotly::layout(margin = list(l = ml))
		}
		return(p)
	}
	
	data_long <- data %>%
		tidyr::pivot_longer(dplyr::all_of(pct_cols), names_to = "Motif", values_to = "Percentage") %>%
		dplyr::mutate(Motif = stringr::str_remove(Motif, "^Pct_"))
	
	txt_size <- sg_text_size(n_sg)
	h_px     <- max(400, n_sg * 22)
	ml_px    <- sg_margin_left_px(unique(data_long$SUBGROUP))
	
	sg_order <- data_long %>%
		dplyr::group_by(SUBGROUP) %>%
		dplyr::summarise(total = sum(Percentage, na.rm = TRUE), .groups = "drop") %>%
		dplyr::arrange(total) %>%
		dplyr::pull(SUBGROUP)
	data_long$SUBGROUP <- factor(data_long$SUBGROUP, levels = sg_order)
	
	sp_label <- if (!is.null(species_filter) && species_filter != "") species_filter else ""
	
	p <- ggplot2::ggplot(data_long,
		ggplot2::aes(x = SUBGROUP, y = Percentage, fill = Motif)) +
		ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.80) +
		ggplot2::geom_text(
			ggplot2::aes(label = ifelse(Percentage > 8, paste0(round(Percentage, 0), "%"), "")),
			position = ggplot2::position_stack(vjust = 0.5),
			size = 2.5, fontface = "bold", color = "white") +
		ggplot2::coord_flip() +
		ggplot2::scale_fill_manual(values = motif_pal, name = "Motif") +
		ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.04))) +
		ggplot2::labs(
			title    = "Motif presence frequencies by subgroup",
			subtitle = paste0(
				if (sp_label != "") paste0("Species: ", sp_label) else "Each color = one motif",
				if (!is.null(top_n)) paste0(" | Top-", top_n) else "",
				"  |  Stacked = cumulative % per subgroup"
			),
			x = NULL, y = "Cumulative presence (%)") +
		theme_imgt(base_size = 10) +
		ggplot2::theme(
			axis.text.y     = ggplot2::element_text(size = txt_size),
			legend.position = "right"
		)
	
	if (interactive) {
		p <- plotly::ggplotly(p, tooltip = c("x","y","fill"), height = h_px) %>%
			plotly::layout(margin = list(l = ml_px))
	}
	return(p)
}


# ══════════════════════════════════════════════════════════════
#  EXPORTS
# ══════════════════════════════════════════════════════════════

save_sg_plot_png <- function(p, file, n_sg = 0, n_col = 1, n_facet = 1,
	grouping = "subgroup", n_sp = 1) {
	if (grouping == "subgroup") {
		dims <- sg_export_dims(n_sg, n_col = n_col, n_facet = n_facet)
	} else if (grouping == "locus") {
		dims <- list(width = 13, height = max(6, n_sp * 3.8))
	} else {
		dims <- list(width = 10, height = 7)
	}
	ggplot2::ggsave(filename = file, plot = p,
		width = dims$width, height = dims$height, units = "in", dpi = .EXPORT_DPI, bg = "white")
}

save_sg_plot_pdf <- function(p, file, n_sg = 0, n_col = 1, n_facet = 1,
	grouping = "subgroup", n_sp = 1) {
	if (grouping == "subgroup") {
		dims <- sg_export_dims(n_sg, n_col = n_col, n_facet = n_facet)
	} else if (grouping == "locus") {
		dims <- list(width = 13, height = max(6, n_sp * 3.8))
	} else {
		dims <- list(width = 10, height = 7)
	}
	ggplot2::ggsave(filename = file, plot = p,
		device = grDevices::cairo_pdf,
		width = dims$width, height = dims$height, units = "in")
}


# ══════════════════════════════════════════════════════════════
#  UI HELPERS
# ══════════════════════════════════════════════════════════════

sg_alert_banner <- function(n_sg, top_n = NULL) {
	if (n_sg <= .SG_SMALL) return(NULL)
	top_active <- !is.null(top_n) && is.numeric(top_n) && top_n > 0 && top_n < n_sg
	use_hm     <- sg_use_heatmap(n_sg)
	
	if (n_sg <= .SG_MEDIUM) {
		cls <- "alert alert-info"
		ico <- shiny::icon("info-circle")
		msg <- paste0(n_sg, " subgroups detected — compact display activated.")
	} else if (n_sg <= .SG_LARGE) {
		cls <- "alert alert-warning"
		ico <- shiny::icon("exclamation-triangle")
		msg <- paste0(n_sg, " subgroups — reduced readability. ",
			"Use the Top-N filter or switch to Heatmap view.")
	} else {
		cls <- "alert alert-danger"
		ico <- shiny::icon("exclamation-circle")
		msg <- paste0(n_sg, " subgroups — Heatmap view activated automatically. ",
			"Use the Top-N filter to go back to barplot.")
	}
	
	shiny::div(class = cls, ico, shiny::strong(" "), msg,
		if (top_active)
			shiny::span(style = "margin-left:10px; color:#155724; font-weight:bold;",
				shiny::icon("filter"),
				paste0(" Active filter: Top-", top_n, " subgroups"))
	)
}


