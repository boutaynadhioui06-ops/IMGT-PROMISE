# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  DIVERSITÉ DES MOTIFS  ██████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════

div_results <- reactiveValues(
	data      = NULL,
	tests     = NULL,
	grouping  = NULL,
	df_snap   = NULL,
	done      = FALSE
)

observe({
	req(store$df_final)
	isolate({
		df_ig <- store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
		sp_list <- sort(unique(df_ig$SPECIES))
		updateSelectInput(session, "div_species_filter",
			choices  = c("All species" = "All species", setNames(sp_list, sp_list)),
			selected = "All species")
		updateSelectInput(session, "div_locus_filter",
			choices  = c("All loci", "IGHV", "IGKV", "IGLV"),
			selected = "All loci")
	})
})

observeEvent(store$df_final, {
	req(store$df_final)
	isolate({
		grouping     <- input$div_grouping %||% "species"
		#  force pour afficher toutes les espèces
		sp_filter    <- "All species"
		locus_filter <- "All loci"
		.run_diversity(grouping, sp_filter, locus_filter)
	})
}, ignoreNULL = TRUE)

observeEvent(input$btn_launch_div, {
	req(store$df_final)
	isolate({
		grouping     <- input$div_grouping %||% "species"
		sp_filter    <- input$div_species_filter %||% "All species"
		locus_filter <- input$div_locus_filter %||% "All loci"
		.run_diversity(grouping, sp_filter, locus_filter)
	})
})

.run_diversity <- function(grouping, sp_filter, locus_filter) {
	withProgress(message = "Computing diversity indices…", value = 0, {
		tryCatch({
			df <- store$df_final
			sp_val    <- if (is.null(sp_filter)    || sp_filter    == "All species") NULL else sp_filter
			locus_val <- if (is.null(locus_filter) || locus_filter == "All loci")   NULL else locus_filter
			
			incProgress(0.3, detail = "Calculating indices…")
			diversity_data <- calculate_diversity_indices(df, "IG", grouping,
				species_filter = sp_val, locus_filter = locus_val)
			
			if (is.null(diversity_data) || nrow(diversity_data) == 0) {
				showNotification("⚠️ No data after filtering.", type = "warning", duration = 5)
				return()
			}
			
			incProgress(0.7, detail = "Statistical tests…")
			tests_data <- tryCatch(
				diversity_statistical_tests(df, "IG", grouping,
					species_filter = sp_val, locus_filter = locus_val),
				error = function(e) NULL
			)
			
			div_results$data     <- diversity_data
			div_results$tests    <- tests_data
			div_results$grouping <- grouping
			div_results$df_snap  <- df
			div_results$done     <- TRUE
			incProgress(1, detail = "Done!")
			
		}, error = function(e) {
			showNotification(paste0("❌ Diversity: ", e$message), type = "error", duration = 8)
		})
	})
}

# Status 
output$div_status <- renderUI({
	if (!isTRUE(div_results$done))
		div(class = "alert alert-info", icon("info-circle"),
			" Build the dataset then click ", strong("🚀 Launch analysis"), ".")
	else {
		n_groups <- nrow(div_results$data)
		grp      <- div_results$grouping %||% "species"
		sp_label    <- if (grp != "species") input$div_species_filter else NULL
		locus_label <- if (grp == "subgroup") input$div_locus_filter else NULL
		info_parts <- c()
		if (!is.null(sp_label) && sp_label != "All species")
			info_parts <- c(info_parts, paste0("Species: ", sp_label))
		if (!is.null(locus_label) && locus_label != "All loci")
			info_parts <- c(info_parts, paste0("Locus: ", locus_label))
		div(class = "alert alert-success", icon("check-circle"),
			paste0(" Analysis complete — ", n_groups, " group(s)",
				if (length(info_parts) > 0) paste0(" | ", paste(info_parts, collapse = " | ")) else ""))
	}
})

# Heights dynamiques
output$div_barplot_ui <- renderUI({
	req(isTRUE(div_results$done))
	grp <- div_results$grouping %||% "species"
	h <- if (grp == "subgroup") {
		n_sg <- length(unique(div_results$data$SUBGROUP))
		min(max(500, n_sg * 20), 1200)
	} else 560
	tagList(
		# Sélecteur d'indice
		div(
			style = "background-color:#EFF6FF; padding:10px 14px; border-radius:6px;
                 border-left:4px solid #3B82F6; margin-bottom:12px;
                 display:flex; align-items:center; gap:16px;",
			tags$span(icon("chart-bar"), strong(" Index to display:"),
				style = "color:#1E40AF; white-space:nowrap;"),
			div(style = "width:220px;",
				selectInput("div_barplot_index", label = NULL,
					choices  = c("All 4 indices" = "all",
						"Shannon (H')"  = "Shannon",
						"Simpson (D)"   = "Simpson",
						"Richness (S)"  = "Richness",
						"Evenness (J)"  = "Evenness"),
					selected = "all", width = "100%")
			),
			tags$small(style = "color:#6B7280; font-size:11px;",
				icon("info-circle"),
				" 'All 4 indices' shows faceted panels — single index for a larger, cleaner view.")
		),
		plotlyOutput("div_barplot", height = paste0(h, "px"))
	)
})

#  Rendus des plots ─
output$div_barplot <- renderPlotly({
	req(isTRUE(div_results$done), div_results$data)
	idx <- input$div_barplot_index %||% "all"
	tryCatch(
		plot_diversity_barplot(div_results$data, div_results$grouping,
			interactive = TRUE, index_filter = idx),
		error = function(e) plotly::plot_ly() %>% plotly::layout(title = paste("Error:", e$message))
	)
})

output$div_radar <- renderPlotly({
	req(isTRUE(div_results$done), div_results$data)
	tryCatch(
		plot_diversity_radar(div_results$data, div_results$grouping),
		error = function(e) plotly::plot_ly() %>% plotly::layout(title = paste("Error:", e$message))
	)
})

output$btn_download_div_radar <- downloadHandler(
	filename = function() paste0("Diversity_Radar_", input$div_grouping, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(div_results$data)
		df   <- div_results$data
		grp  <- div_results$grouping
		n_gr <- nrow(df)
		
		#  radar 
		idx_cols <- intersect(c("Shannon","Simpson","Richness","Evenness"), names(df))
		if (length(idx_cols) == 0) {
			showNotification("❌ No index columns found for radar.", type = "error", duration = 5)
			return()
		}
		
		# Variable de groupe
		grp_col <- switch(grp,
			species  = "SPECIES",
			locus    = "LOCI",
			subgroup = "SUBGROUP",
			"SPECIES"
		)
		if (!grp_col %in% names(df)) grp_col <- names(df)[1]
		
		# Normalisation 0-1 par indice
		df_norm <- df
		for (col in idx_cols) {
			mn <- min(df[[col]], na.rm = TRUE)
			mx <- max(df[[col]], na.rm = TRUE)
			df_norm[[col]] <- if (mx > mn) (df[[col]] - mn) / (mx - mn) else 0.5
		}
		
		# Format long pour ggplot
		df_long <- df_norm %>%
			dplyr::select(dplyr::all_of(c(grp_col, idx_cols))) %>%
			tidyr::pivot_longer(cols = dplyr::all_of(idx_cols),
				names_to = "Index", values_to = "Value") %>%
			dplyr::rename(Group = dplyr::all_of(grp_col))
		
		# Fermer le polygone 
		df_closed <- df_long %>%
			dplyr::group_by(Group) %>%
			dplyr::bind_rows(df_long %>%
					dplyr::group_by(Group) %>%
					dplyr::slice_head(n = 1)) %>%
			dplyr::ungroup()
		
		n_idx <- length(idx_cols)
		angles <- seq(0, 2 * pi, length.out = n_idx + 1)[1:n_idx]
		
		# Coordonnées cartésiennes 
		df_cart <- df_long %>%
			dplyr::left_join(
				data.frame(Index = idx_cols, angle = angles),
				by = "Index"
			) %>%
			dplyr::mutate(
				x = Value * cos(angle),
				y = Value * sin(angle)
			)
		
		# Grille de fond 
		grid_df <- expand.grid(
			r     = c(0.25, 0.5, 0.75, 1.0),
			theta = seq(0, 2 * pi, length.out = 100)
		) %>%
			dplyr::mutate(x = r * cos(theta), y = r * sin(theta))
		
		# Lignes des axes
		axis_df <- data.frame(
			Index = idx_cols, angle = angles
		) %>%
			dplyr::mutate(xend = cos(angle), yend = sin(angle))
		
		# Labels des axes
		label_df <- data.frame(
			Index = idx_cols, angle = angles
		) %>%
			dplyr::mutate(
				x = 1.18 * cos(angle),
				y = 1.18 * sin(angle)
			)
		
		# Palette
		n_groups <- length(unique(df_long$Group))
		pal <- if (n_groups <= 12) {
			setNames(scales::hue_pal()(n_groups), unique(df_long$Group))
		} else {
			setNames(viridis::viridis(n_groups), unique(df_long$Group))
		}
		
		p <- ggplot2::ggplot() +
			# Grille
			ggplot2::geom_path(data = grid_df,
				ggplot2::aes(x = x, y = y, group = r),
				color = "grey80", linewidth = 0.3, linetype = "dashed") +
			# Axes
			ggplot2::geom_segment(data = axis_df,
				ggplot2::aes(x = 0, y = 0, xend = xend, yend = yend),
				color = "grey60", linewidth = 0.4) +
			# Labels axes
			ggplot2::geom_text(data = label_df,
				ggplot2::aes(x = x, y = y, label = Index),
				size = 4, fontface = "bold", color = "#374151") +
			# Polygones groupes
			ggplot2::geom_polygon(data = df_cart,
				ggplot2::aes(x = x, y = y, group = Group, fill = Group, color = Group),
				alpha = 0.15, linewidth = 0.9) +
			ggplot2::geom_point(data = df_cart,
				ggplot2::aes(x = x, y = y, color = Group),
				size = 3) +
			ggplot2::scale_fill_manual(values  = pal, name = grp_col) +
			ggplot2::scale_color_manual(values = pal, name = grp_col) +
			ggplot2::coord_equal(xlim = c(-1.35, 1.35), ylim = c(-1.35, 1.35)) +
			ggplot2::labs(
				title    = paste0("Diversity Radar — ", grp_col),
				subtitle = paste0("Normalized indices (0=min, 1=max) | ",
					n_gr, " groups | Shannon, Simpson, Richness, Evenness")
			) +
			theme_imgt(base_size = 12) +
			ggplot2::theme(
				axis.text     = ggplot2::element_blank(),
				axis.ticks    = ggplot2::element_blank(),
				axis.title    = ggplot2::element_blank(),
				panel.grid    = ggplot2::element_blank(),
				legend.position = "right"
			)
		
		w <- if (n_gr <= 6) 10 else 13
		h <- if (n_gr <= 6) 8  else 10
		ggplot2::ggsave(file, plot = p, width = w, height = h, dpi = 300, bg = "white")
		showNotification("✅ Radar PNG downloaded!", type = "message", duration = 3)
	}
)

output$div_tests_table <- renderDT({
	req(isTRUE(div_results$done), div_results$tests)
	df <- div_results$tests
	dt <- datatable(df,
		options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv")),
		extensions = "Buttons", rownames = FALSE,
		class = "cell-border stripe hover compact") %>%
		formatRound(columns = intersect(c("Statistic","p.value"), names(df)), digits = 4)
	if ("Significant" %in% names(df)) {
		dt <- dt %>% formatStyle("Significant",
			backgroundColor = styleEqual(
				c("***","**","*","ns","—"),
				c("#D4EDDA","#D4EDDA","#FFF3CD","#F8D7DA","#F9FAFB")),
			fontWeight = "bold")
	}
	dt
})

# Exports 
output$btn_download_div_csv <- downloadHandler(
	filename = function() paste0("Diversity_", input$div_grouping, "_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(div_results$data)
		write.csv(div_results$data, file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	}
)

output$btn_download_div_barplot <- downloadHandler(
	filename = function() {
		idx <- input$div_barplot_index %||% "all"
		sfx <- if (idx == "all") "All" else idx
		paste0("Diversity_Barplot_", sfx, "_", input$div_grouping, "_", Sys.Date(), ".png")
	},
	content  = function(file) {
		req(div_results$data)
		idx  <- input$div_barplot_index %||% "all"
		p    <- plot_diversity_barplot(div_results$data, div_results$grouping,
			interactive = FALSE, index_filter = idx)
		grp  <- div_results$grouping
		n_sg <- if (grp == "subgroup") length(unique(div_results$data$SUBGROUP)) else 1
		if (idx == "all") {
			h <- if (grp == "subgroup") max(10, n_sg * 0.28 + 2) else 7
			w <- 16
		} else {
			h <- if (grp == "subgroup") max(8, n_sg * 0.22 + 2) else 6
			w <- 12
		}
		ggsave(file, plot = p, width = min(w, 40), height = min(h, 60),
			dpi = 300, bg = "white")
		showNotification("✅ PNG downloaded!", type = "message", duration = 3)
	}
)
