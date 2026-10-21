# ═══════════════════════════════════════════════════════════════
#  Advanced Visualisations
# ═══════════════════════════════════════════════════════════════

# Dépend de : utils_core.R
# Packages  : tidyverse, plotly, ggridges, viridis

library(tidyverse)
library(visNetwork)
library(ggplot2)
library(plotly)

# ══════════════════════════════════════════════════════════════
# MODULE A — VISUALISATIONS CO-OCCURRENCE
# ══════════════════════════════════════════════════════════════

#  Heatmap de corrélation phi
plot_cooccurrence_heatmap <- function(cooc_result,
	mode        = "phi",
	interactive = TRUE) {
	
	mat <- if (mode == "phi") cooc_result$matrix_phi else cooc_result$matrix_pct
	motifs <- cooc_result$motifs
	n      <- length(motifs)
	
	mat_lower <- mat
	mat_lower[upper.tri(mat_lower, diag = TRUE)] <- NA
	
	data_long <- as.data.frame(mat_lower) %>%
		tibble::rownames_to_column("Motif_1") %>%
		tidyr::pivot_longer(-Motif_1, names_to = "Motif_2", values_to = "Value") %>%
		dplyr::filter(!is.na(Value)) %>%
		dplyr::mutate(
			Motif_1 = factor(Motif_1, levels = motifs),
			Motif_2 = factor(Motif_2, levels = motifs),
			Label   = if (mode == "phi") sprintf("%.2f", Value) else paste0(Value, "%")
		)
	
	if (mode == "phi") {
		fill_scale <- ggplot2::scale_fill_gradient2(
			low = "#DC2626", mid = "#F9FAFB", high = "#16A34A",
			midpoint = 0, limits = c(-1, 1), name = "φ Correlation",
			breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1","-0.5","0","0.5","1")
		)
		title_str    <- paste0("Phi correlations between motifs (N = ", cooc_result$n_seq, " seq.)")
		subtitle_str <- "Green = positive co-occurrence | Red = mutual exclusion"
	} else {
		fill_scale <- ggplot2::scale_fill_viridis_c(
			option = "viridis", limits = c(0, 100), name = "Co-occurrence (%)",
			breaks = c(0, 25, 50, 75, 100)
		)
		title_str    <- paste0("Motif co-occurrence — % sequences (N = ", cooc_result$n_seq, ")")
		subtitle_str <- "Percentage of sequences simultaneously carrying both motifs"
	}
	
	tile_size <- dplyr::case_when(n <= 4 ~ 4.5, n <= 6 ~ 3.8, n <= 8 ~ 3.0, TRUE ~ 2.5)
	
	p <- ggplot2::ggplot(data_long,
		ggplot2::aes(x = Motif_2, y = Motif_1, fill = Value)) +
		ggplot2::geom_tile(color = "white", linewidth = 0.8) +
		ggplot2::geom_text(
			ggplot2::aes(
				label = Label,
				color = ifelse(
					if (mode == "phi") abs(Value) >= 0.4 else Value >= 50,
					"white", "grey25"
				)
			),
			size = tile_size, fontface = "bold"
		) +
		ggplot2::scale_color_identity() +
		fill_scale +
		ggplot2::labs(title = title_str, subtitle = subtitle_str, x = NULL, y = NULL) +
		theme_imgt(base_size = 12) +
		ggplot2::theme(
			axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, size = 11),
			axis.text.y = ggplot2::element_text(size = 11),
			panel.grid  = ggplot2::element_blank(),
			legend.position = "right",
			legend.key.height = ggplot2::unit(1.5, "cm")
		)
	
	if (interactive)
		p <- plotly::ggplotly(p, tooltip = c("x","y","fill"), height = 480)
	
	return(p)
}




# Barplot Fisher — associations significatives
plot_fisher_associations <- function(fisher_df, interactive = TRUE) {
	
	if (is.null(fisher_df) || nrow(fisher_df) == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No Fisher results available.", size = 5, color = "grey50") +
				ggplot2::theme_void())
	}
	
	df_plot <- fisher_df %>%
		dplyr::mutate(
			Pair    = paste0(Motif_1, " × ", Motif_2),
			log10_p = -log10(pmax(p_value, 1e-6)),
			sig_label = dplyr::case_when(
				Significant == "***" ~ "p < 0.001",
				Significant == "**"  ~ "p < 0.01",
				Significant == "*"   ~ "p < 0.05",
				TRUE                 ~ "ns"
			)
		) %>%
		dplyr::arrange(dplyr::desc(log10_p)) %>%
		dplyr::slice_head(n = 30)
	
	df_plot$Pair <- factor(df_plot$Pair, levels = rev(df_plot$Pair))
	
	#  valeurs présentes uniquement
	assoc_vals_present <- unique(df_plot$Association)
	assoc_pal_full <- c(
		"Positive co-occurrence" = "#16A34A",
		"Mutual exclusion"       = "#DC2626",
		"Not significant"        = "#CBD5E1",
		"Independent"            = "#94A3B8"
	)
	assoc_pal <- assoc_pal_full[names(assoc_pal_full) %in% assoc_vals_present]
	
	p <- ggplot2::ggplot(df_plot,
		ggplot2::aes(x = Pair, y = log10_p, fill = Association)) +
		ggplot2::geom_bar(stat = "identity", width = 0.72) +
		ggplot2::geom_hline(yintercept = -log10(0.05),
			linetype = "dashed", color = "#F59E0B", linewidth = 0.8) +
		ggplot2::geom_hline(yintercept = -log10(0.001),
			linetype = "dotted", color = "#EF4444", linewidth = 0.7) +
		ggplot2::geom_text(
			ggplot2::aes(label = ifelse(Significant != "ns", Significant, "")),
			hjust = -0.2, size = 3.5, fontface = "bold", color = "grey20"
		) +
		ggplot2::coord_flip() +
		ggplot2::scale_fill_manual(values = assoc_pal, name = "Association") +
		ggplot2::scale_y_continuous(
			expand = ggplot2::expansion(mult = c(0, 0.20)),
			name   = "-log₁₀(p-value)"
		) +
		ggplot2::annotate("text",
			x = 1, y = -log10(0.05) + 0.05,
			label = "p = 0.05", size = 3, color = "#F59E0B", hjust = -0.1) +
		ggplot2::labs(
			title    = "Fisher Tests — Associations between motif pairs",
			subtitle = "Top 30 pairs | Yellow = p<0.05 | Red = p<0.001",
			x = NULL
		) +
		theme_imgt(base_size = 11) +
		ggplot2::theme(legend.position = "right")
	
	if (interactive) {
		h <- max(400, nrow(df_plot) * 22)
		p <- plotly::ggplotly(p, tooltip = c("x","y","fill"), height = min(h, 800))
	}
	return(p)
}


# ══════════════════════════════════════════════════════════════
# MODULE B — ALLUVIAL 
# ══════════════════════════════════════════════════════════════

plot_alluvial_parcats <- function(df,
	type_data = "IG",
	species   = NULL) {
	
	if (type_data == "IG") {
		df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	} else {
		df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	}
	if (!is.null(species) && species != "" && species != "All")
		df <- df %>% dplyr::filter(SPECIES == species)
	df <- df %>% dplyr::filter(!is.na(FUNCTIONALITY))
	if (nrow(df) == 0) stop("No data after filtering.")
	
	df_agg <- df %>%
		dplyr::group_by(SPECIES, LOCI, FUNCTIONALITY) %>%
		dplyr::summarise(N = dplyr::n(), .groups = "drop") %>%
		dplyr::mutate(FUNCTIONALITY = factor(FUNCTIONALITY, levels = c("F","P","ORF")))
	
	fig <- plotly::plot_ly(
		type = "parcats",
		dimensions = list(
			list(label = "Species",       values = df_agg$SPECIES,                  categoryorder = "category ascending"),
			list(label = "Locus",         values = df_agg$LOCI,                     categoryorder = "category ascending"),
			list(label = "Functionality", values = as.character(df_agg$FUNCTIONALITY), categoryarray = c("F","P","ORF"))
		),
		counts = df_agg$N,
		line = list(
			color = as.integer(factor(df_agg$FUNCTIONALITY, levels = c("F","P","ORF"))),
			colorscale = list(c(0,"#16A34A"), c(0.5,"#DC2626"), c(1,"#F59E0B")),
			shape = "hspline"
		),
		labelfont   = list(size = 13, family = "Arial", color = "#1E293B"),
		arrangement = "freeform",
		bundlecolors = TRUE
	) %>%
		plotly::layout(
			title = list(
				text = paste0("Flow SPECIES → LOCI → FUNCTIONALITY",
					if (!is.null(species) && species != "") paste0(" — ", species) else ""),
				font = list(size = 15, color = "#1E293B")
			),
			margin = list(t = 60, b = 30, l = 50, r = 50),
			height = 520
		)
	return(fig)
}


# ══════════════════════════════════════════════════════════════
# MODULE D — RIDGELINE DES DISTANCES
# ══════════════════════════════════════════════════════════════

plot_distance_ridgeline <- function(df,
	motif       = NULL,
	group_var   = "SPECIES",
	species     = NULL,
	type_data   = "IG",
	interactive = TRUE) {
	
	if (type_data == "IG") {
		df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	} else {
		df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	}
	if (!is.null(species) && species != "" && species != "All")
		df <- df %>% dplyr::filter(SPECIES == species)
	
	dist_cols <- names(df)[stringr::str_starts(names(df), "DIST_")]
	if (length(dist_cols) == 0) stop("No DIST_* columns available.")
	
	if (!is.null(motif) && motif != "" && motif != "All") {
		dist_col_name <- paste0("DIST_", motif)
		if (dist_col_name %in% dist_cols) dist_cols <- dist_col_name
	}
	
	df_long <- df %>%
		dplyr::select(dplyr::all_of(c(group_var, dist_cols))) %>%
		tidyr::pivot_longer(
			cols = dplyr::all_of(dist_cols),
			names_to  = "Motif",
			values_to = "Distance"
		) %>%
		dplyr::mutate(Motif = stringr::str_remove(Motif, "^DIST_")) %>%
		dplyr::filter(!is.na(Distance), Distance >= 0)
	
	if (nrow(df_long) == 0) stop("No valid distances after filtering.")
	
	n_groups <- length(unique(df_long[[group_var]]))
	n_motifs <- length(unique(df_long$Motif))
	txt_size <- dplyr::case_when(n_groups <= 5 ~ 10, n_groups <= 10 ~ 9, TRUE ~ 8)
	
	# Resolved palette — supports SPECIES, LOCI, SUBGROUP
	group_vals <- unique(df_long[[group_var]])
	pal <- if (group_var == "SPECIES") {
		get_species_palette(group_vals)
	} else if (group_var == "LOCI") {
		loci_pal_full <- tryCatch(get("LOCI_COLORS", envir = .GlobalEnv), error = function(e) NULL)
		if (!is.null(loci_pal_full) && all(group_vals %in% names(loci_pal_full))) {
			loci_pal_full[group_vals]
		} else {
			setNames(viridis::viridis(length(group_vals), option = "D"), group_vals)
		}
	} else {
		setNames(viridis::viridis(length(group_vals), option = "turbo"), group_vals)
	}
	pal[is.na(pal)] <- "#94A3B8"
	
	group_label <- switch(group_var,
		SPECIES  = "Species",
		LOCI     = "Locus",
		SUBGROUP = "Subgroup",
		group_var)
	
	use_ridges <- requireNamespace("ggridges", quietly = TRUE)
	
	if (use_ridges && n_motifs == 1) {
		p <- ggplot2::ggplot(df_long,
			ggplot2::aes(x = Distance, y = .data[[group_var]],
				fill = .data[[group_var]])) +
			ggridges::geom_density_ridges(
				alpha = 0.75, scale = 1.2, rel_min_height = 0.01, color = "white"
			) +
			ggplot2::scale_fill_manual(values = pal, guide = "none") +
			ggplot2::labs(
				title    = paste0("Distance to LPART1 \u2014 ", unique(df_long$Motif)),
				subtitle = paste0("Grouping: ", group_label,
					" | N = ", nrow(df_long), " observations"),
				x = "Distance (nt)", y = NULL
			) +
			theme_imgt(base_size = 11) +
			ggplot2::theme(axis.text.y = ggplot2::element_text(size = txt_size))
	} else {
		med_df <- df_long %>%
			dplyr::group_by(.data[[group_var]], Motif) %>%
			dplyr::summarise(MedianDist = median(Distance, na.rm = TRUE), .groups = "drop")
		df_long <- df_long %>%
			dplyr::left_join(med_df, by = c(group_var, "Motif")) %>%
			dplyr::mutate(!!group_var := forcats::fct_reorder(.data[[group_var]], MedianDist, .desc = TRUE))
		
		p <- ggplot2::ggplot(df_long,
			ggplot2::aes(x = .data[[group_var]], y = Distance,
				fill = MedianDist)) +
			ggplot2::geom_violin(trim = FALSE, alpha = 0.82, color = "white", linewidth = 0.4) +
			ggplot2::geom_jitter(width = 0.12, alpha = 0.25, size = 1, color = "#374151") +
			ggplot2::stat_summary(fun = median, geom = "crossbar",
				width = 0.5, color = "#1F2937", linewidth = 0.8) +
			ggplot2::scale_fill_viridis_c(option = "plasma", direction = -1,
				name = "Median distance (nt)") +
			ggplot2::facet_wrap(~ Motif, scales = "free_y") +
			ggplot2::labs(
				title    = "Distance distribution to LPART1",
				subtitle = paste0("Grouping: ", group_label, " | Bar = median"),
				x = NULL, y = "Distance (nt)"
			) +
			theme_imgt(base_size = 11) +
			ggplot2::theme(
				axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = txt_size)
			)
	}
	
	if (interactive) {
		h <- max(480, n_groups * 40)
		p <- plotly::ggplotly(p, tooltip = c("x","y","fill"), height = min(h, 800))
	}
	return(p)
}

