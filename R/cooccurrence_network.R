# ═══════════════════════════════════════════════════════════════
#  Réseau Co-occurrence PLOTLY 
# ═══════════════════════════════════════════════════════════════

# Dépend de : utils_core.R, cooccurrence_analysis.R
# Packages  : plotly, igraph, dplyr

plot_cooccurrence_network_plotly <- function(graph,
	fisher_df     = NULL,
	cooc_result   = NULL,
	layout_seed   = 42) {
	
	if (is.null(graph) || igraph::vcount(graph) == 0)
		stop("Empty or NULL graph.")
	
	# Layout graphe 
	set.seed(layout_seed)
	w_layout <- if (igraph::ecount(graph) > 0) pmax(abs(igraph::E(graph)$weight), 0.01) else NULL
	layout_mat <- igraph::layout_with_fr(graph, weights = w_layout)
	
	motifs    <- igraph::V(graph)$name
	node_freq <- igraph::V(graph)$freq
	n_nodes   <- length(motifs)
	
	# Colorsnœuds
	default_pal <- c("#1565C0","#16A34A","#DC2626","#F59E0B",
		"#7C3AED","#0891B2","#EA580C","#65A30D","#DB2777")
	node_colors <- sapply(seq_along(motifs), function(i) {
		tryCatch({
			col_map <- get("MOTIF_COLORS", envir = .GlobalEnv)
			m <- motifs[i]
			if (m %in% names(col_map)) as.character(col_map[[m]])
			else default_pal[((i - 1) %% length(default_pal)) + 1]
		}, error = function(e) {
			default_pal[((i - 1) %% length(default_pal)) + 1]
		})
	})
	
	# Données nœuds
	node_x    <- layout_mat[, 1]
	node_y    <- layout_mat[, 2]
	node_size <- pmax(as.numeric(node_freq), 5) / 100 * 40 + 10  
	
	node_text <- paste0(
		"<b>", motifs, "</b><br>",
		"Frequency: ", node_freq, "%<br>",
		"N sequences: ", if (!is.null(cooc_result)) cooc_result$n_seq else "?"
	)
	
	# Arêtes
	edge_list <- igraph::as_data_frame(graph, what = "edges")
	
	# Initialiser la figure
	fig <- plotly::plot_ly(type = "scatter", mode = "markers")
	
	if (nrow(edge_list) > 0) {
		motif_idx <- setNames(seq_along(motifs), motifs)
		
		for (k in seq_len(nrow(edge_list))) {
			e   <- edge_list[k, ]
			i   <- motif_idx[e$from]
			j   <- motif_idx[e$to]
			
			# Couleur arête
			edge_col <- dplyr::case_when(
				e$assoc == "Positive co-occurrence" ~ "rgba(22,163,74,0.7)",
				e$assoc == "Mutual exclusion"       ~ "rgba(220,38,38,0.7)",
				e$sig   == "ns"                     ~ "rgba(203,213,225,0.5)",
				TRUE                                ~ "rgba(107,114,128,0.6)"
			)
			
			# Épaisseur ~ |phi|
			lw <- max(1, abs(e$weight) * 10)
			
			# Tooltip 
			mid_x <- (node_x[i] + node_x[j]) / 2
			mid_y <- (node_y[i] + node_y[j]) / 2
			
			edge_hover <- paste0(
				"<b>", e$from, " × ", e$to, "</b><br>",
				"φ = ", round(e$weight, 3), "<br>",
				"Co-occurrence: ", e$cooc_pct, "%<br>",
				"OR = ", round(e$or, 2), "<br>",
				"p = ", round(e$p_value, 4), " [", e$sig, "]<br>",
				"<i>", e$assoc, "</i>"
			)
			
			# Ligne arête
			fig <- fig %>% plotly::add_trace(
				x    = c(node_x[i], node_x[j], NA),
				y    = c(node_y[i], node_y[j], NA),
				type = "scatter", mode = "lines",
				line = list(color = edge_col, width = lw),
				hoverinfo = "none",
				showlegend = FALSE,
				inherit = FALSE
			)
			fig <- fig %>% plotly::add_annotations(
				x         = mid_x,
				y         = mid_y,
				text      = paste0("\u03c6=", round(e$weight, 2)),
				showarrow = FALSE,
				font      = list(size = 9, color = "#1E293B", family = "Arial"),
				bgcolor   = "rgba(255,255,255,0.75)",
				bordercolor = "rgba(0,0,0,0.08)",
				borderwidth = 1,
				borderpad   = 1
			)
			# Point fantôme au milieu pour afficher le tooltip
			fig <- fig %>% plotly::add_trace(
				x         = mid_x,
				y         = mid_y,
				type      = "scatter",
				mode      = "markers",
				marker    = list(size = 6, color = edge_col, opacity = 0),
				text      = edge_hover,
				hoverinfo = "text",
				showlegend = FALSE,
				inherit = FALSE
			)
		}
	}
	
	# Nœuds
	fig <- fig %>% plotly::add_trace(
		x          = node_x,
		y          = node_y,
		type       = "scatter",
		mode       = "markers+text",
		marker     = list(
			size    = node_size,
			color   = node_colors,
			line    = list(color = "#1E293B", width = 2),
			opacity = 0.92
		),
		text      = motifs,
		textposition = "top center",
		textfont  = list(size = 12, color = "#1E293B", family = "Arial"),
		customdata = node_text,
		hovertemplate = "%{customdata}<extra></extra>",
		showlegend = FALSE,
		inherit    = FALSE
	)
	
	# Legendmanuelle 
	legend_items <- list(
		list(x = 0.01, y = 0.98, text = "● Positive co-occurrence", color = "#16A34A"),
		list(x = 0.01, y = 0.93, text = "● Mutual exclusion",       color = "#DC2626"),
		list(x = 0.01, y = 0.88, text = "● Not significant",        color = "#CBD5E1")
	)
	
	annotations <- lapply(legend_items, function(li) {
		list(
			x          = li$x, y = li$y,
			xref       = "paper", yref = "paper",
			text       = li$text,
			showarrow  = FALSE,
			font       = list(size = 11, color = li$color, family = "Arial"),
			align      = "left"
		)
	})
	
	# Layout final
	n_edges <- nrow(edge_list)
	subtitle_text <- if (n_edges == 0) {
		"No edges — lower the |φ| threshold in the settings"
	} else {
		paste0(
			n_nodes, " motifs | ", n_edges, " associations | ",
			"Green = co-occurrence | Red = exclusion | Size \u221d frequency | Width \u221d |φ|"
		)
	}
	
	fig <- fig %>% plotly::layout(
		title = list(
			text = paste0("<b> motif co-occurrence network</b><br>",
				"<span style='font-size:11px;color:#64748B;'>",
				subtitle_text, "</span>"),
			font = list(size = 14, color = "#1E293B", family = "Arial"),
			x = 0.5, xanchor = "center"
		),
		xaxis = list(visible = FALSE, zeroline = FALSE, showgrid = FALSE),
		yaxis = list(visible = FALSE, zeroline = FALSE, showgrid = FALSE),
		plot_bgcolor  = "#F8FAFF",
		paper_bgcolor = "#FFFFFF",
		hovermode     = "closest",
		annotations   = annotations,
		margin        = list(t = 80, b = 20, l = 20, r = 20),
		height        = 580
	)
	
	return(fig)
}

#  rendu pour le server patch
render_cooc_network_plotly <- function(cooc_result, fisher_df,
	phi_threshold, sig_only) {
	
	if (is.null(cooc_result))
		return(plotly::plot_ly() %>%
				plotly::layout(title = "No data available"))
	
	# Reconstruire le graphe
	graph_res <- tryCatch(
		build_cooccurrence_graph(
			fisher_df,
			cooc_result,
			phi_threshold = phi_threshold,
			sig_only      = sig_only
		),
		error = function(e) NULL
	)
	
	if (is.null(graph_res) || igraph::vcount(graph_res) == 0) {
		return(plotly::plot_ly() %>%
				plotly::layout(
					title = paste0("Empty network — |φ| threshold = ", phi_threshold,
						" too high or no significant association"),
					plot_bgcolor = "#F8FAFF"))
	}
	
	plot_cooccurrence_network_plotly(
		graph       = graph_res,
		fisher_df   = fisher_df,
		cooc_result = cooc_result
	)
}
