# ═══════════════════════════════════════════════════════════════
#PERMANOVA MODULE 
# ═══════════════════════════════════════════════════════════════

# Dépend de : utils_core.R
# Packages  : shiny, vegan, plotly, ggrepel, purrr, tibble, DT

# ══════════════════════════════════════════════════════════════
# PALETTES
# ══════════════════════════════════════════════════════════════

PERM_MOTIF_COLORS <- c(
	"TATA_BOX"    = "#E63946",
	"HEPTA_PENTA" = "#F4A261",
	"OCTA_DECA"   = "#2A9D8F",
	"E_BOX"       = "#457B9D",
	"TATA_BOX2"   = "#9B2226",
	"EBOX_2"      = "#606C38",
	"LPART_1"     = "#7B2D8B",
	"LPART1"      = "#7B2D8B"
)

.perm_motif_color <- function(motif_names) {
	base   <- PERM_MOTIF_COLORS
	unk    <- setdiff(motif_names, names(base))
	extras <- c("#1565C0","#2E7D32","#6A1B9A","#00838F","#D84315",
		"#558B2F","#4527A0","#00695C","#AD1457")
	if (length(unk) > 0) {
		base <- c(base, setNames(
			extras[seq_len(min(length(unk), length(extras)))],
			unk[seq_len(min(length(unk), length(extras)))]))
	}
	cols <- base[motif_names]
	cols[is.na(cols)] <- "#94A3B8"
	unname(cols)
}

.perm_group_palette <- function(group_vals, group_var) {
	if (group_var == "LOCI") {
		pal <- tryCatch({
			lc  <- get("LOCI_COLORS", envir = .GlobalEnv)
			out <- lc[group_vals]
			out[is.na(out)] <- "#64748B"
			out
		}, error = function(e)
			setNames(viridis::viridis(length(group_vals)), group_vals))
	} else if (group_var == "FUNCTIONALITY") {
		pal_f <- c("F" = "#16A34A", "P" = "#DC2626", "ORF" = "#F59E0B")
		pal   <- pal_f[group_vals]
		pal[is.na(pal)] <- "#94A3B8"
		names(pal) <- group_vals
	} else {
		pal <- tryCatch(
			get_species_palette(group_vals),
			error = function(e)
				setNames(viridis::viridis(length(group_vals)), group_vals))
		pal[is.na(pal)] <- "#64748B"
	}
	pal <- pal[group_vals]
	names(pal) <- group_vals
	pal[is.na(pal)] <- "#64748B"
	pal
}

.perm_empty_plot <- function(msg) {
	ggplot2::ggplot() +
		ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
			size = 4.5, color = "#DC2626", hjust = 0.5) +
		ggplot2::theme_void()
}


# ══════════════════════════════════════════════════════════════
#  run_permanova 
# ══════════════════════════════════════════════════════════════

run_permanova <- function(df,
	type_data = "IG",
	group_var = "SPECIES",
	species   = NULL,
	n_perm    = 999) {
	
	if (type_data == "IG") {
		df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	} else {
		df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	}
	if (!is.null(species) && species != "" && species != "All")
		df <- df %>% dplyr::filter(SPECIES == species)
	
	bin_cols <- names(df)[stringr::str_starts(names(df), "BIN_")]
	if (length(bin_cols) < 2) stop("Less than 2 BIN_* motifs available.")
	if (!group_var %in% names(df))
		stop(paste("Variable", group_var, "not found in dataset."))
	
	df <- df %>% dplyr::filter(!is.na(.data[[group_var]]))
	
	n_groups <- length(unique(df[[group_var]]))
	if (n_groups < 2) stop("Less than 2 distinct groups for PERMANOVA.")
	if (nrow(df) < n_groups + 2) stop("Not enough sequences for PERMANOVA.")
	
	bin_mat <- df %>%
		dplyr::select(dplyr::all_of(bin_cols)) %>%
		dplyr::mutate(dplyr::across(dplyr::everything(),
			~ as.integer(tidyr::replace_na(., 0L)))) %>%
		as.matrix()
	colnames(bin_mat) <- stringr::str_remove(bin_cols, "^BIN_")
	
	col_var <- apply(bin_mat, 2, function(x) stats::var(x, na.rm = TRUE))
	bin_mat <- bin_mat[, !is.na(col_var) & col_var > 0, drop = FALSE]
	if (ncol(bin_mat) < 2)
		stop("Less than 2 variable BIN_* columns — PERMANOVA impossible.")
	bin_mat[is.na(bin_mat)] <- 0L
	
	
	row_sums  <- rowSums(bin_mat)
	keep_rows <- row_sums > 0
	n_removed <- sum(!keep_rows)
	if (n_removed > 0) {
		message(paste0("PERMANOVA: removing ", n_removed,
			" sequences with all-zero motif profiles (undefined Jaccard distance)."))
	}
	bin_mat <- bin_mat[keep_rows, , drop = FALSE]
	df      <- df[keep_rows, ]
	
	group_vec <- as.character(df[[group_var]])
	n_groups  <- length(unique(group_vec))
	if (n_groups < 2) stop("Less than 2 distinct groups after removing empty rows.")
	if (nrow(bin_mat) < n_groups + 2) stop("Not enough sequences after removing empty rows.")
	
	stopifnot(nrow(bin_mat) == length(group_vec), nrow(bin_mat) == nrow(df))
	
	dist_mat <- vegan::vegdist(bin_mat, method = "jaccard", binary = TRUE)
	
	adonis_res <- tryCatch(
		vegan::adonis2(dist_mat ~ group_vec, permutations = n_perm, method = "jaccard"),
		error = function(e) stop(paste("PERMANOVA failed:", e$message))
	)
	
	summary_df <- as.data.frame(adonis_res) %>%
		tibble::rownames_to_column("Term") %>%
		dplyr::mutate(
			Term = dplyr::case_when(
				!Term %in% c("Residual", "Total") ~ group_var,
				TRUE ~ Term
			),
			Significant = dplyr::case_when(
				is.na(`Pr(>F)`)  ~ "\u2014",
				`Pr(>F)` < 0.001 ~ "***",
				`Pr(>F)` < 0.01  ~ "**",
				`Pr(>F)` < 0.05  ~ "*",
				TRUE             ~ "ns"
			)
		) %>%
		dplyr::rename(SS = SumOfSqs, F_stat = `F`, p_value = `Pr(>F)`)
	
	betadisp <- tryCatch({
		bd           <- vegan::betadisper(dist_mat, group_vec)
		permutest_bd <- vegan::permutest(bd, permutations = 499)
		pairwise_bd  <- tryCatch(
			vegan::permutest(bd, pairwise = TRUE, permutations = 499),
			error = function(e) NULL)
		list(
			distances   = bd$distances,
			group_means = tapply(bd$distances, group_vec, mean),
			p_value     = permutest_bd$tab$`Pr(>F)`[1],
			pairwise    = if (!is.null(pairwise_bd)) pairwise_bd$pairwise else NULL
		)
	}, error = function(e) NULL)
	
	df_filtered <- df %>%
		dplyr::select(dplyr::any_of(c("ALLELE","SPECIES","LOCI","SUBGROUP",
			"FUNCTIONALITY", group_var)))
	
	list(
		adonis_result   = adonis_res,
		summary_df      = summary_df,
		distance_matrix = dist_mat,
		bin_mat         = bin_mat,
		betadisper      = betadisp,
		group_var       = group_var,
		group_vec       = group_vec,
		df_filtered     = df_filtered,
		n_seq           = nrow(bin_mat),
		n_groups        = n_groups
	)
}


# ══════════════════════════════════════════════════════════════
#  plot_permanova_pcoa  
# ══════════════════════════════════════════════════════════════

plot_permanova_pcoa <- function(permanova_result, interactive = TRUE) {
	
	dist_mat  <- permanova_result$distance_matrix
	group_var <- permanova_result$group_var
	groups    <- permanova_result$group_vec
	
	pcoa_res <- tryCatch(stats::cmdscale(dist_mat, k = 2, eig = TRUE), error = function(e) NULL)
	if (is.null(pcoa_res) || nrow(pcoa_res$points) == 0)
		return(.perm_empty_plot("PCoA not available."))
	
	n_pts   <- nrow(pcoa_res$points)
	eig_pos <- pcoa_res$eig[pcoa_res$eig > 0]
	var_exp <- if (length(eig_pos) >= 2)
		round(eig_pos[1:2] / sum(eig_pos) * 100, 1) else c(0, 0)
	
	coords        <- as.data.frame(pcoa_res$points)
	names(coords) <- c("PCoA1", "PCoA2")
	
	if (length(groups) != n_pts) {
		return(.perm_empty_plot(paste0(
			"Group size mismatch (dist=", n_pts, " vs groups=", length(groups),
			").\nPlease re-run PERMANOVA.")))
	}
	coords$Group <- as.character(groups)
	
	if (!is.null(permanova_result$df_filtered) &&
			nrow(permanova_result$df_filtered) == n_pts) {
		meta <- permanova_result$df_filtered %>%
			dplyr::select(dplyr::any_of(c("ALLELE","SPECIES","LOCI","SUBGROUP")))
		coords <- cbind(coords, meta)
	}
	
	tt <- paste0("<b>", coords$Group, "</b>")
	if ("ALLELE"   %in% names(coords)) tt <- paste0(tt, "<br>Allele: ",   coords$ALLELE)
	if ("SPECIES"  %in% names(coords)) tt <- paste0(tt, "<br>Species: ",  coords$SPECIES)
	if ("LOCI"     %in% names(coords)) tt <- paste0(tt, "<br>Locus: ",    coords$LOCI)
	if ("SUBGROUP" %in% names(coords)) tt <- paste0(tt, "<br>Subgroup: ", coords$SUBGROUP)
	tt <- paste0(tt,
		"<br>PCoA1: ", round(coords$PCoA1, 3),
		"<br>PCoA2: ", round(coords$PCoA2, 3))
	coords$tooltip <- tt
	
	group_vals   <- sort(unique(coords$Group))
	pal          <- .perm_group_palette(group_vals, group_var)
	group_counts <- coords %>%
		dplyr::group_by(Group) %>%
		dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
		dplyr::filter(n >= 4)
	
	if (!interactive) {
		p <- ggplot2::ggplot(coords,
			ggplot2::aes(x = PCoA1, y = PCoA2, color = Group, text = tooltip)) +
			ggplot2::geom_point(size = 3.0, alpha = 0.80) +
			ggplot2::scale_color_manual(values = pal, name = group_var) +
			ggplot2::labs(
				title    = paste0("PCoA \u2014 Motif profiles by ", group_var),
				subtitle = paste0("Jaccard | PCoA1=", var_exp[1], "% | PCoA2=",
					var_exp[2], "% | N=", n_pts),
				x = paste0("PCoA1 (", var_exp[1], "%)"),
				y = paste0("PCoA2 (", var_exp[2], "%)")) +
			theme_imgt(base_size = 12) +
			ggplot2::theme(legend.position = "right")
		if (nrow(group_counts) > 0)
			p <- p + ggplot2::stat_ellipse(
				data = coords %>% dplyr::filter(Group %in% group_counts$Group),
				ggplot2::aes(fill = Group), geom = "polygon",
				level = 0.95, alpha = 0.08, show.legend = FALSE) +
				ggplot2::scale_fill_manual(values = pal, guide = "none")
		return(p)
	}
	
	fig <- plotly::plot_ly(height = 540)
	
	for (grp in group_counts$Group) {
		df_grp  <- coords %>% dplyr::filter(Group == grp)
		col_grp <- if (grp %in% names(pal)) pal[[grp]] else "#64748B"
		ell <- tryCatch({
			theta <- seq(0, 2 * pi, length.out = 100)
			mu    <- colMeans(df_grp[, c("PCoA1","PCoA2")])
			sigma <- stats::cov(df_grp[, c("PCoA1","PCoA2")])
			chi2  <- stats::qchisq(0.95, df = 2)
			eg    <- eigen(sigma)
			sdev  <- sqrt(chi2 * eg$values)
			pts   <- t(mu + eg$vectors %*% (sdev * rbind(cos(theta), sin(theta))))
			data.frame(x = pts[,1], y = pts[,2])
		}, error = function(e) NULL)
		if (!is.null(ell))
			fig <- plotly::add_trace(fig,
				type      = "scatter", mode = "lines",
				x         = c(ell$x, ell$x[1]), y = c(ell$y, ell$y[1]),
				fill      = "toself",
				fillcolor = paste0(col_grp, "1A"),
				line      = list(color = col_grp, width = 1.0, dash = "dot"),
				name      = grp, showlegend = FALSE, hoverinfo = "skip")
	}
	
	for (grp in group_vals) {
		df_grp  <- coords %>% dplyr::filter(Group == grp)
		col_grp <- if (grp %in% names(pal)) pal[[grp]] else "#64748B"
		fig <- plotly::add_trace(fig,
			type      = "scatter", mode = "markers",
			x         = df_grp$PCoA1, y = df_grp$PCoA2,
			text      = df_grp$tooltip, hoverinfo = "text",
			name      = grp, showlegend = TRUE,
			marker    = list(color = col_grp, size = 9, opacity = 0.85,
				line = list(color = "white", width = 0.5)))
	}
	
	fig %>%
		plotly::layout(
			title  = list(
				text = paste0("<b>PCoA \u2014 Motif profiles by ", group_var, "</b>"),
				font = list(size = 14, color = "#1E293B")),
			xaxis  = list(title = paste0("PCoA1 (", var_exp[1], "%)"),
				zeroline = FALSE, showgrid = TRUE, gridcolor = "#E5E7EB"),
			yaxis  = list(title = paste0("PCoA2 (", var_exp[2], "%)"),
				zeroline = FALSE, showgrid = TRUE, gridcolor = "#E5E7EB"),
			legend = list(title = list(text = paste0("<b>", group_var, "</b>")),
				bgcolor = "rgba(255,255,255,0.85)",
				bordercolor = "#CBD5E1", borderwidth = 1),
			paper_bgcolor = "white", plot_bgcolor = "#F8FAFF",
			annotations = list(list(
				text = paste0("Jaccard distance | N = ", n_pts, " sequences"),
				xref = "paper", yref = "paper",
				x = 0, y = -0.08, xanchor = "left", yanchor = "top",
				showarrow = FALSE,
				font = list(size = 10, color = "#6B7280")))
		) %>%
		plotly::config(displayModeBar = TRUE,
			modeBarButtonsToRemove = c("lasso2d","select2d","autoScale2d"))
}


# ══════════════════════════════════════════════════════════════
#  plot_permanova_pcoa_subgroup 
# ══════════════════════════════════════════════════════════════

plot_permanova_pcoa_subgroup <- function(permanova_result, color_by = "SPECIES", interactive = TRUE) {
	
	dist_mat <- permanova_result$distance_matrix
	meta     <- permanova_result$df_filtered
	
	if (is.null(meta) || !color_by %in% names(meta)) color_by <- "SPECIES"
	
	if (is.null(meta) || !all(c("SUBGROUP", color_by) %in% names(meta)))
		return(.perm_empty_plot(
			paste0("SUBGROUP and/or ", color_by, " not available in this result.\nPlease re-run PERMANOVA.")))
	
	pcoa_res <- tryCatch(stats::cmdscale(dist_mat, k = 2, eig = TRUE), error = function(e) NULL)
	if (is.null(pcoa_res)) return(.perm_empty_plot("PCoA not available."))
	
	n_pts <- nrow(pcoa_res$points)
	if (nrow(meta) != n_pts)
		return(.perm_empty_plot(paste0(
			"Inconsistent size (dist=", n_pts, " vs meta=", nrow(meta),
			").\nPlease re-run PERMANOVA.")))
	
	eig_pos <- pcoa_res$eig[pcoa_res$eig > 0]
	var_exp <- if (length(eig_pos) >= 2)
		round(eig_pos[1:2] / sum(eig_pos) * 100, 1) else c(0, 0)
	
	coords <- as.data.frame(pcoa_res$points)
	names(coords) <- c("PCoA1", "PCoA2")
	coords$SUBGROUP <- meta$SUBGROUP
	coords$ColorVar <- meta[[color_by]]
	
	centroids <- coords %>%
		dplyr::filter(!is.na(SUBGROUP), !is.na(ColorVar)) %>%
		dplyr::group_by(SUBGROUP, ColorVar) %>%
		dplyr::summarise(
			PCoA1 = mean(PCoA1), PCoA2 = mean(PCoA2),
			N = dplyr::n(), .groups = "drop") %>%
		dplyr::mutate(tooltip = paste0(
			"<b>", SUBGROUP, "</b><br>", color_by, ": ", ColorVar,
			"<br>N alleles: ", N,
			"<br>PCoA1: ", round(PCoA1, 3), "<br>PCoA2: ", round(PCoA2, 3)))
	
	if (nrow(centroids) == 0)
		return(.perm_empty_plot(paste0("No centroid could be computed (missing SUBGROUP/", color_by, ").")))
	
	group_vals <- sort(unique(centroids$ColorVar))
	pal        <- .perm_group_palette(group_vals, color_by)
	
	p <- ggplot2::ggplot(centroids,
		ggplot2::aes(x = PCoA1, y = PCoA2, color = ColorVar, size = N, text = tooltip)) +
		ggplot2::geom_point(alpha = 0.85) +
		ggrepel::geom_text_repel(
			ggplot2::aes(label = SUBGROUP), size = 2.6, max.overlaps = 30,
			show.legend = FALSE, color = "#374151") +
		ggplot2::scale_color_manual(values = pal, name = color_by) +
		ggplot2::scale_size_continuous(range = c(2, 9), name = "N alleles") +
		ggplot2::labs(
			title    = paste0("PCoA \u2014 Subgroups (centroids) colored by ", tolower(color_by)),
			subtitle = paste0(
				"Each point = centroid of one subgroup | ",
				nrow(centroids), " subgroup(s) | size = number of alleles"),
			x = paste0("PCoA1 (", var_exp[1], "%)"),
			y = paste0("PCoA2 (", var_exp[2], "%)")) +
		theme_imgt(base_size = 12) +
		ggplot2::theme(legend.position = "right")
	
	if (interactive)
		p <- plotly::ggplotly(p, tooltip = "text", height = 580) %>%
		plotly::layout(legend = list(title = list(text = color_by)))
	
	p
}

# ══════════════════════════════════════════════════════════════
# plot_betadisper_boxplot
# ══════════════════════════════════════════════════════════════

plot_betadisper_boxplot <- function(permanova_result, interactive = TRUE) {
	
	dist_mat  <- permanova_result$distance_matrix
	group_var <- permanova_result$group_var
	groups    <- permanova_result$group_vec
	
	n_dist <- attr(dist_mat, "Size")
	if (length(groups) != n_dist) {
		return(.perm_empty_plot(paste0(
			"Size mismatch: dist=", n_dist, " vs groups=", length(groups),
			".\nPlease re-run PERMANOVA.")))
	}
	groups <- as.character(groups)
	
	bd_res <- tryCatch(vegan::betadisper(dist_mat, groups), error = function(e) NULL)
	if (is.null(bd_res)) return(.perm_empty_plot("Beta-dispersion computation failed."))
	
	df_bd <- data.frame(
		Group    = groups,
		Distance = bd_res$distances,
		tooltip  = paste0("<b>", groups, "</b><br>Dist. to centroid: ",
			round(bd_res$distances, 3))
	)
	
	group_vals <- sort(unique(df_bd$Group))
	pal        <- .perm_group_palette(group_vals, group_var)
	
	bd_store <- permanova_result$betadisper
	bd_pval  <- if (!is.null(bd_store) && !is.null(bd_store$p_value))
		bd_store$p_value else NA
	
	pval_label <- if (!is.na(bd_pval)) {
		sig <- dplyr::case_when(
			bd_pval < 0.001 ~ "***",
			bd_pval < 0.01  ~ "**",
			bd_pval < 0.05  ~ "*",
			TRUE            ~ "ns")
		paste0("Homogeneity test: p = ", format(round(bd_pval, 4), nsmall = 4), " ", sig)
	} else "Homogeneity test: p = N/A"
	
	pairwise_annot <- NULL
	if (length(group_vals) >= 3) {
		pw <- tryCatch(
			stats::pairwise.wilcox.test(df_bd$Distance, df_bd$Group, p.adjust.method = "BH"),
			error = function(e) NULL)
		if (!is.null(pw)) {
			pw_df <- as.data.frame(as.table(pw$p.value)) %>%
				dplyr::rename(G1 = Var1, G2 = Var2, p = Freq) %>%
				dplyr::filter(!is.na(p), p < 0.05) %>%
				dplyr::mutate(
					G1  = as.character(G1), G2 = as.character(G2),
					sig = dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*"))
			pairwise_annot <- pw_df
		}
	}
	
	n_groups  <- length(group_vals)
	lbl_angle <- if (n_groups > 5) 35 else 0
	
	p <- ggplot2::ggplot(df_bd,
		ggplot2::aes(x = Group, y = Distance, text = tooltip)) +
		ggplot2::geom_jitter(ggplot2::aes(color = Group),
			width = 0.14, alpha = 0.35, size = 1.6) +
		ggplot2::geom_boxplot(ggplot2::aes(color = Group),
			fill = NA, outlier.shape = NA,
			linewidth = 0.9, width = 0.52) +
		ggplot2::scale_color_manual(values = pal, guide = "none") +
		ggplot2::labs(
			title    = paste0("Beta-dispersion by ", group_var),
			subtitle = paste0(
				"Jaccard distance of each sequence to its group centroid "),
			x = NULL, y = "Distance to centroid (Jaccard)") +
		ggplot2::annotate("text", x = Inf, y = Inf, label = pval_label,
			hjust = 1.05, vjust = 1.5, size = 3.2,
			color = if (!is.na(bd_pval) && bd_pval < 0.05) "#DC2626" else "#374151",
			fontface = "bold") +
		theme_imgt(base_size = 12) +
		ggplot2::theme(axis.text.x = ggplot2::element_text(
			angle = lbl_angle, hjust = if (lbl_angle > 0) 1 else 0.5))
	
	if (!is.null(pairwise_annot) && nrow(pairwise_annot) > 0) {
		y_max    <- max(df_bd$Distance, na.rm = TRUE)
		y_range  <- diff(range(df_bd$Distance, na.rm = TRUE))
		y_step   <- max(y_range * 0.10, 0.01)
		x_levels <- levels(factor(group_vals))
		
		for (k in seq_len(nrow(pairwise_annot))) {
			g1 <- pairwise_annot$G1[k]; g2 <- pairwise_annot$G2[k]
			x1 <- which(x_levels == g1); x2 <- which(x_levels == g2)
			if (length(x1) == 0 || length(x2) == 0) next
			y_br <- y_max + k * y_step
			p <- p +
				ggplot2::annotate("segment", x = x1, xend = x2, y = y_br, yend = y_br,
					color = "#374151", linewidth = 0.5) +
				ggplot2::annotate("text", x = (x1 + x2) / 2, y = y_br + y_step * 0.3,
					label = pairwise_annot$sig[k], size = 3.5, color = "#DC2626", fontface = "bold")
		}
		y_lim <- y_max + (nrow(pairwise_annot) + 1) * y_step
		p <- p + ggplot2::coord_cartesian(ylim = c(NA, y_lim))
	}
	
	if (interactive)
		p <- plotly::ggplotly(p, tooltip = "text", height = 480) %>%
		plotly::layout(xaxis = list(tickangle = if (n_groups > 5) -35 else 0))
	
	p
}


# ══════════════════════════════════════════════════════════════
# plot_motif_contribution
# ══════════════════════════════════════════════════════════════

plot_motif_contribution <- function(permanova_result, interactive = TRUE) {
	
	bin_mat   <- permanova_result$bin_mat
	group_vec <- permanova_result$group_vec
	group_var <- permanova_result$group_var
	
	if (is.null(bin_mat) || is.null(group_vec) || nrow(bin_mat) != length(group_vec))
		return(.perm_empty_plot(
			"bin_mat / group_vec missing or desynchronized.\nPlease re-run PERMANOVA."))
	
	groups <- as.character(group_vec)
	
	simp <- tryCatch(
		vegan::simper(bin_mat, group = groups, permutations = 0),
		error = function(e) { message("SIMPER failed: ", e$message); NULL })
	
	simp_summary <- NULL
	if (!is.null(simp)) {
		simp_summary <- tryCatch(summary(simp), error = function(e) {
			message("SIMPER summary failed: ", e$message); NULL })
	}
	
	if (!is.null(simp_summary)) {
		all_contribs <- lapply(names(simp_summary), function(pair) {
			s <- simp_summary[[pair]]
			data.frame(Motif = rownames(s), Contrib = s$average, Cumul = s$cumsum,
				Pair = pair, stringsAsFactors = FALSE)
		})
		simper_df <- dplyr::bind_rows(all_contribs) %>%
			dplyr::group_by(Motif) %>%
			dplyr::summarise(
				Mean_contrib = mean(Contrib, na.rm = TRUE),
				SD_contrib   = sd(Contrib,   na.rm = TRUE),
				N_pairs      = dplyr::n(), .groups = "drop") %>%
			dplyr::arrange(dplyr::desc(Mean_contrib)) %>%
			dplyr::mutate(Motif = factor(Motif, levels = rev(Motif)))
		method_label <- "SIMPER (vegan) \u2014 average contribution to Jaccard dissimilarity"
	} else {
		motif_names <- colnames(bin_mat)
		simper_df <- purrr::map_dfr(motif_names, function(m) {
			vals <- bin_mat[, m]
			g    <- factor(groups)
			if (length(unique(vals)) < 2) return(NULL)
			grp_means  <- tapply(vals, g, mean, na.rm = TRUE)
			grand_mean <- mean(vals, na.rm = TRUE)
			ss_between <- sum(table(g) * (grp_means - grand_mean)^2)
			ss_total   <- sum((vals - grand_mean)^2)
			eta2       <- if (ss_total > 0) ss_between / ss_total else 0
			data.frame(Motif = m, Mean_contrib = eta2, SD_contrib = 0, N_pairs = 1)
		})
		if (is.null(simper_df) || nrow(simper_df) == 0)
			return(.perm_empty_plot(
				"Could not compute any contribution (both SIMPER and \u03b7\u00b2 failed)."))
		simper_df <- simper_df %>%
			dplyr::arrange(dplyr::desc(Mean_contrib)) %>%
			dplyr::mutate(Motif = factor(Motif, levels = rev(Motif)))
		method_label <- "Eta\u00b2 (ANOVA) \u2014 proportion of between-group variance per motif"
	}
	
	if (nrow(simper_df) == 0)
		return(.perm_empty_plot("No motif contribution data available."))
	
	motifs_top <- as.character(levels(simper_df$Motif))
	motif_pal  <- setNames(.perm_motif_color(motifs_top), motifs_top)
	
	p1 <- ggplot2::ggplot(simper_df,
		ggplot2::aes(x = Motif, y = Mean_contrib, fill = Motif,
			text = paste0("<b>", Motif, "</b><br>Contribution: ",
				round(Mean_contrib * 100, 2), "%"))) +
		ggplot2::geom_col(width = 0.6, alpha = 0.92) +
		ggplot2::geom_errorbar(
			ggplot2::aes(ymin = pmax(0, Mean_contrib - SD_contrib), ymax = Mean_contrib + SD_contrib),
			width = 0.22, color = "grey40", linewidth = 0.5) +
		ggplot2::scale_fill_manual(values = motif_pal, guide = "none") +
		ggplot2::scale_y_continuous(
			labels = scales::percent_format(accuracy = 0.1),
			expand = ggplot2::expansion(mult = c(0, 0.10))) +
		ggplot2::coord_flip() +
		ggplot2::labs(
			title    = paste0("Motif contribution to group differences \u2014 ", group_var),
			subtitle = method_label, x = NULL, y = "Contribution (%)") +
		theme_imgt(base_size = 13) +
		ggplot2::theme(
			axis.text.y  = ggplot2::element_text(face = "bold", size = 11),
			axis.text.x  = ggplot2::element_text(size = 10),
			panel.grid.major.y = ggplot2::element_blank(),
			panel.grid.major.x = ggplot2::element_line(
				color = "#E5E7EB", linewidth = 0.4, linetype = "dashed"))
	
	if (interactive) {
		plot_h <- max(420, nrow(simper_df) * 55 + 120)
		cols   <- motif_pal[as.character(simper_df$Motif)]
		
		fig <- plotly::plot_ly(
			data        = simper_df,
			x           = ~Mean_contrib,
			y           = ~Motif,
			type        = "bar",
			orientation = "h",
			marker      = list(color = cols, line = list(color = "white", width = 0.5)),
			error_x     = list(
				type       = "data",
				symmetric  = FALSE,
				array      = simper_df$SD_contrib,
				arrayminus = pmin(simper_df$SD_contrib, simper_df$Mean_contrib),
				color      = "#6B7280",
				thickness  = 1.1,
				width      = 6
			),
			text      = paste0("<b>", simper_df$Motif, "</b><br>Contribution: ",
				round(simper_df$Mean_contrib * 100, 2), "%"),
			hoverinfo = "text",
			height    = plot_h
		) %>%
			plotly::layout(
				title  = list(
					text = paste0("<b>Motif contribution to group differences \u2014 ",
						group_var, "</b>"),
					font = list(size = 14)),
				xaxis  = list(title = "Contribution (%)", tickformat = ".1%",
					zeroline = FALSE, showgrid = TRUE, gridcolor = "#E5E7EB"),
				yaxis  = list(title = "", autorange = "reversed"),
				margin = list(l = 10, r = 30, t = 60, b = 40),
				paper_bgcolor = "white", plot_bgcolor = "#F8FAFF",
				annotations = list(list(
					text = method_label,
					xref = "paper", yref = "paper",
					x = 0, y = 1.04, xanchor = "left", yanchor = "bottom",
					showarrow = FALSE,
					font = list(size = 10, color = "#6B7280")))
			)
		
		return(fig)
	}
	
	p1
}
# ══════════════════════════════════════════════════════════════
# UI 
# ══════════════════════════════════════════════════════════════

mod_permanova_ui <- function(id) {
	ns <- shiny::NS(id)
	shiny::tagList(
		
		shiny::div(
			style = "background:linear-gradient(90deg,#EFF6FF,#F0FDF4);
               border-radius:6px;padding:10px 14px;font-size:12px;
               color:#374151;margin-bottom:12px;",
			shiny::icon("info-circle"),
			shiny::strong(" PERMANOVA (Permutational ANOVA) \u2014 "),
			"Tests whether motif profiles differ significantly between groups. ",
			"Jaccard distance on binary motif presence is used. ",
			shiny::strong("IG receptors only"), " (TR loci have no annotated regulatory motifs)."
		),
		
		shiny::fluidRow(
			shiny::column(4,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("sliders-h"), shiny::strong("Parameters")),
					shiny::selectInput(ns("perm_species"), "Species (optional)",
						choices = c("All" = ""), selected = ""),
					shiny::selectInput(ns("perm_group_var"), "Grouping variable",
						choices = c(
							"Species"       = "SPECIES",
							"Locus"         = "LOCI",
							"Functionality" = "FUNCTIONALITY",
							"Subgroup"      = "SUBGROUP"
						), selected = "SPECIES"),
					shiny::numericInput(ns("perm_nperm"), "Number of permutations",
						value = 99, min = 99, max = 999, step = 100),
					shiny::tags$small(style = "color:#6B7280;",
						shiny::icon("clock"),
						" 99 is fastest. Above 99, computation can take noticeably longer.")
				)
			),
			shiny::column(4,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("eye"), shiny::strong("PCoA view")),
					shiny::radioButtons(ns("perm_pcoa_view"), NULL,
						choices = c("Per subgroup (centroids)" = "SUBGROUP_SPECIES"),
						selected = "SUBGROUP_SPECIES"),
					shiny::selectInput(ns("perm_color_by"), "Color by",
						choices = c("Species" = "SPECIES", "Locus" = "LOCI", "Subgroup" = "SUBGROUP"),
						selected = "SPECIES")
				)
			),
			shiny::column(4,
				shiny::div(style = "padding-top:8px;",
					shiny::actionButton(ns("btn_run_perm"),
						"\U0001F52C Run PERMANOVA",
						icon  = shiny::icon("play"),
						class = "btn-success btn-lg",
						style = "width:100%;height:60px;font-size:14px;"),
					shiny::br(), shiny::br(),
					shiny::div(
						style = "background:#FEF9C3;border-radius:6px;padding:8px;
                     font-size:11px;color:#78350F;",
						shiny::icon("exclamation-triangle"),
						shiny::strong(" Note: "),
						"Sequences with no detected motif (all-zero profile) are ",
						"automatically removed (Jaccard distance requires this)."
					)
				)
			)
		),
		
		shiny::hr(),
		shiny::uiOutput(ns("perm_status")),
		
		shiny::tabsetPanel(
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("table"), " PERMANOVA Results"),
				shiny::br(),
				DT::DTOutput(ns("perm_table")),
				shiny::br()
			),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("project-diagram"), " PCoA"),
				shiny::br(),
				shiny::uiOutput(ns("perm_pcoa_hint")),
				plotly::plotlyOutput(ns("perm_pcoa"), height = "580px"),
				shiny::br()
			),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("chart-bar"), " Beta-dispersion"),
				shiny::br(),
				plotly::plotlyOutput(ns("perm_betadisp"), height = "500px"),
				shiny::br()
			),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("puzzle-piece"), " Motif Contributions"),
				shiny::br(),
				shiny::uiOutput(ns("perm_contrib_status")),
				shiny::h6(shiny::strong("\U0001F4CA Contribution per motif (SIMPER / \u03b7\u00b2)")),
				plotly::plotlyOutput(ns("perm_motif_contrib"), height = "520px"),
				shiny::br()
			),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("download"), " Exports"),
				shiny::br(),
				shiny::fluidRow(
					shiny::column(3,
						shiny::div(
							style = "text-align:center;padding:14px;background:#EFF6FF;
                       border-radius:8px;border:1px solid #93C5FD;",
							shiny::icon("file-csv",
								style = "font-size:24px;color:#1D4ED8;display:block;margin-bottom:6px;"),
							shiny::p(style = "font-size:11px;color:#374151;margin:0 0 8px 0;",
								"PERMANOVA table (CSV)"),
							shiny::downloadButton(ns("dl_perm_csv"), "\U0001F4E5 CSV \u2014 Results",
								class = "btn-primary btn-block", style = "width:100%;")
						)
					),
					shiny::column(3,
						shiny::div(
							style = "text-align:center;padding:14px;background:#EFF6FF;
                       border-radius:8px;border:1px solid #93C5FD;",
							shiny::icon("image",
								style = "font-size:24px;color:#1D4ED8;display:block;margin-bottom:6px;"),
							shiny::p(style = "font-size:11px;color:#374151;margin:0 0 8px 0;",
								"High-resolution PCoA (PNG)"),
							shiny::downloadButton(ns("dl_pcoa_png"), "\U0001F4E5 PNG \u2014 PCoA",
								class = "btn-info btn-block", style = "width:100%;")
						)
					),
					shiny::column(3,
						shiny::div(
							style = "text-align:center;padding:14px;background:#F0FDF4;
                       border-radius:8px;border:1px solid #86EFAC;",
							shiny::icon("image",
								style = "font-size:24px;color:#16A34A;display:block;margin-bottom:6px;"),
							shiny::p(style = "font-size:11px;color:#374151;margin:0 0 8px 0;",
								"Beta-dispersion (PNG)"),
							shiny::downloadButton(ns("dl_betadisp_png"), "\U0001F4E5 PNG \u2014 Beta-dispersion",
								class = "btn-success btn-block", style = "width:100%;")
						)
					),
					shiny::column(3,
						shiny::div(
							style = "text-align:center;padding:14px;background:#FFFBEB;
                       border-radius:8px;border:1px solid #FCD34D;",
							shiny::icon("image",
								style = "font-size:24px;color:#D97706;display:block;margin-bottom:6px;"),
							shiny::p(style = "font-size:11px;color:#374151;margin:0 0 8px 0;",
								"Motif contributions (PNG)"),
							shiny::downloadButton(ns("dl_contrib_png"), "\U0001F4E5 PNG \u2014 Contributions",
								class = "btn-warning btn-block", style = "width:100%;")
						)
					)
				)
			)
		)
	)
}


# ══════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════

mod_permanova_server <- function(id, df_final) {
	shiny::moduleServer(id, function(input, output, session) {
		
		ns         <- session$ns
		perm_store <- shiny::reactiveValues(result = NULL, done = FALSE)
		
		shiny::observe({
			shiny::req(df_final())
			df      <- df_final() %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
			sp_list <- sort(unique(df$SPECIES))
			shiny::isolate(
				shiny::updateSelectInput(session, "perm_species",
					choices  = c("All" = "", setNames(sp_list, sp_list)),
					selected = ""))
		})
		
		shiny::observeEvent(input$btn_run_perm, {
			shiny::req(df_final())
			
			np_check <- if (!is.null(input$perm_nperm) && !is.na(input$perm_nperm))
				as.integer(input$perm_nperm) else 99L
			if (np_check > 99) {
				shiny::showNotification(
					paste0("\u23F3 ", np_check, " permutations selected \u2014 this can take ",
						"noticeably longer. Please wait until the computation finishes."),
					type = "warning", duration = 8)
			}
			
			
			perm_store$result <- NULL
			perm_store$done   <- FALSE
			
			shiny::withProgress(message = "Running PERMANOVA\u2026", value = 0.1, {
				tryCatch({
					sp  <- if (is.null(input$perm_species) || input$perm_species == "") NULL else input$perm_species
					np  <- if (!is.null(input$perm_nperm)  && !is.na(input$perm_nperm))
						as.integer(input$perm_nperm) else 99L
					
					df_use <- df_final()
					
					shiny::incProgress(0.3, detail = "PERMANOVA\u2026")
					result <- run_permanova(df_use,
						type_data = "IG",
						group_var = input$perm_group_var %||% "SPECIES",
						species   = sp,
						n_perm    = np)
					
					shiny::incProgress(1)
					perm_store$result <- result
					perm_store$done   <- TRUE
					
				}, error = function(e) {
					shiny::showNotification(paste0("\u274c PERMANOVA: ", e$message),
						type = "error", duration = 10)
				})
			})
		})
		
		output$perm_status <- shiny::renderUI({
			if (!perm_store$done)
				shiny::div(class = "alert alert-info", shiny::icon("info-circle"),
					" Configure the parameters and click ",
					shiny::strong("\U0001F52C Run PERMANOVA"), ".")
			else {
				res <- perm_store$result
				
				r2_vals <- res$summary_df$R2[res$summary_df$Term == res$group_var]
				r2 <- if (length(r2_vals) == 1 && !is.na(r2_vals)) round(r2_vals, 3) else NA
				
				pv_vals <- res$summary_df$p_value[res$summary_df$Term == res$group_var]
				pv <- if (length(pv_vals) == 1 && !is.na(pv_vals)) round(pv_vals, 4) else NA
				
				sig <- if (!is.na(pv) && pv < 0.05) " \u2705 Significant"
				else if (is.na(pv)) " \u2753 Not available"
				else " \u26a0\ufe0f Not significant"
				
				shiny::div(class = "alert alert-success", shiny::icon("check-circle"),
					shiny::strong(" PERMANOVA done \u2014 "),
					res$n_seq, " sequences | ", res$n_groups, " groups | ",
					"R\u00b2 = ", ifelse(is.na(r2), "N/A", r2), " | p = ",
					ifelse(is.na(pv), "N/A", pv), sig)
			}
		})
		
		output$perm_table <- DT::renderDT({
			shiny::req(perm_store$result)
			df <- perm_store$result$summary_df
			DT::datatable(df,
				options = list(pageLength = 10, scrollX = TRUE, dom = "Bfrtip",
					buttons = c("copy","csv")),
				extensions = "Buttons", rownames = FALSE,
				class = "cell-border stripe hover compact") %>%
				DT::formatRound(columns = intersect(c("SS","R2","F_stat","p_value"), names(df)),
					digits = 4) %>%
				DT::formatStyle("Significant",
					backgroundColor = DT::styleEqual(
						c("***","**","*","ns","\u2014"),
						c("#D4EDDA","#D4EDDA","#FFF3CD","#F8D7DA","#F8F9FA")),
					fontWeight = "bold")
		})
		
		output$perm_pcoa_hint <- shiny::renderUI({
			if (identical(input$perm_pcoa_view, "SUBGROUP_SPECIES"))
				shiny::div(style = "background:#F0FDF4;border-left:4px solid #16A34A;
                   border-radius:6px;padding:8px 14px;font-size:12px;margin-bottom:10px;",
					shiny::icon("info-circle"),
					" Aggregated view: one point = one subgroup (centroid of its alleles), ",
					"color = species, size = number of alleles in the subgroup.")
			else
				shiny::div(style = "background:#F0FDF4;border-left:4px solid #16A34A;
                   border-radius:6px;padding:8px 14px;font-size:12px;margin-bottom:10px;",
					shiny::icon("info-circle"),
					" Each point = one sequence. Ellipses = 95% confidence per group (if \u2265 4).")
		})
		
		output$perm_pcoa <- plotly::renderPlotly({
			shiny::req(perm_store$result)
			tryCatch({
				if (identical(input$perm_pcoa_view, "SUBGROUP_SPECIES"))
					plot_permanova_pcoa_subgroup(perm_store$result, color_by = input$perm_color_by %||% "SPECIES", interactive = TRUE)
				else
					plot_permanova_pcoa(perm_store$result, interactive = TRUE)
			}, error = function(e) plotly::ggplotly(.perm_empty_plot(paste("PCoA error:", e$message))))
		})
		
		output$perm_betadisp <- plotly::renderPlotly({
			shiny::req(perm_store$result)
			tryCatch(
				plot_betadisper_boxplot(perm_store$result, interactive = TRUE),
				error = function(e) plotly::ggplotly(.perm_empty_plot(paste("Beta-disp error:", e$message))))
		})
		
		contrib_cache <- shiny::reactive({
			shiny::req(perm_store$result)
			tryCatch(
				list(res = plot_motif_contribution(perm_store$result, interactive = TRUE), err = NULL),
				error = function(e) list(res = NULL, err = e$message))
		})
		
		output$perm_contrib_status <- shiny::renderUI({
			cc <- contrib_cache()
			if (!is.null(cc$err))
				shiny::div(class = "alert alert-danger", shiny::icon("exclamation-circle"),
					shiny::strong(" Motif Contributions error: "), cc$err)
			else shiny::tagList()
		})
		
		output$perm_motif_contrib <- plotly::renderPlotly({
			cc <- contrib_cache()
			if (!is.null(cc$res)) cc$res
			else plotly::ggplotly(.perm_empty_plot("No contribution data available"))
		})
		
		output$dl_perm_csv <- shiny::downloadHandler(
			filename = function() paste0("PERMANOVA_", Sys.Date(), ".csv"),
			content = function(file) {
				shiny::req(perm_store$result)
				write.csv(perm_store$result$summary_df, file, row.names = FALSE)
			})
		
		output$dl_pcoa_png <- shiny::downloadHandler(
			filename = function() paste0("PCoA_", perm_store$result$group_var, "_", Sys.Date(), ".png"),
			content = function(file) {
				shiny::req(perm_store$result)
				p <- if (identical(input$perm_pcoa_view, "SUBGROUP_SPECIES"))
					plot_permanova_pcoa_subgroup(perm_store$result, color_by = input$perm_color_by %||% "SPECIES", interactive = TRUE)
				else
					plot_permanova_pcoa(perm_store$result, interactive = FALSE)
				ggplot2::ggsave(file, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
			})
		output$dl_betadisp_png <- shiny::downloadHandler(
			filename = function() paste0("BetaDispersion_",
				if (!is.null(perm_store$result)) perm_store$result$group_var else "na",
				"_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
			content = function(file) {
				shiny::req(perm_store$result)
				tryCatch({
					p <- plot_betadisper_boxplot(perm_store$result, interactive = FALSE)
					n_groups <- perm_store$result$n_groups
					w <- max(10, n_groups * 1.6)   
					h <- 7                          
					ggplot2::ggsave(file, plot = p, width = w, height = h, dpi = 300, bg = "white")
					shiny::showNotification("✅ Beta-dispersion PNG downloaded!", type = "message", duration = 3)
				}, error = function(e) {
					shiny::showNotification(paste0("❌ Beta-dispersion export: ", e$message),
						type = "error", duration = 8)
				})
			})
		output$dl_contrib_png <- shiny::downloadHandler(
			filename = function() paste0("MotifContrib_", perm_store$result$group_var, "_", Sys.Date(), ".png"),
			content = function(file) {
				shiny::req(perm_store$result)
				p <- plot_motif_contribution(perm_store$result, interactive = FALSE)
				ggplot2::ggsave(file, plot = p, width = 10, height = 7, dpi = 300, bg = "white")
			})
	})
}