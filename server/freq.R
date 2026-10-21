# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  CONSERVATION DES MOTIFS  ██████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════
# 1. FRÉQUENCES DE PRÉSENCE 
# ══════════════════════════════════════════════════════════════

freq_results <- reactiveValues(data=NULL, plot=NULL, done=FALSE, grouping=NULL, df_snap=NULL,
	prev_grouping=NULL, stat_tests=NULL, pairwise_tests=NULL)

observe({
	req(store$df_final)
	df_ig   <- store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	sp_list <- sort(unique(df_ig$SPECIES))
	# Subgroup mode species selector
	updateSelectInput(session, "freq_species_filter",
		choices  = setNames(sp_list, sp_list),
		selected = if (length(sp_list) > 0) sp_list[1] else NULL)
	# Locus mode species selector
	updateSelectInput(session, "freq_locus_species_filter",
		choices  = c("All species" = "", setNames(sp_list, sp_list)),
		selected = "")
	# Locus filter for subgroup mode
	loci_list <- sort(unique(df_ig$LOCI))
	updateSelectInput(session, "freq_locus_filter",
		choices  = c("All loci", loci_list),
		selected = "All loci")
})

.run_freq <- function(grouping, sp_sel, locus_sel, sp_locus_filter, sg_mode, topn_val) {
	tryCatch({
		df_ig <- store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
		if (nrow(df_ig) == 0) return()
		
		if (grouping == "subgroup") {
			sp_list <- sort(unique(df_ig$SPECIES))
			
			if (is.null(sp_sel) || sp_sel == "" || !(sp_sel %in% sp_list)) {
				sp_sel <- sp_list[1]
				updateSelectInput(session, "freq_species_filter",
					choices = setNames(sp_list, sp_list), selected = sp_sel)
			}
			df_filt <- df_ig %>% dplyr::filter(SPECIES == sp_sel)
			if (!is.null(locus_sel) && locus_sel != "" && locus_sel != "All loci")
				df_filt <- df_filt %>% dplyr::filter(LOCI == locus_sel)
			if (nrow(df_filt) == 0) return()
			top_n_val   <- if (!is.null(topn_val) && !is.na(topn_val) && topn_val > 0) topn_val else NULL
			data_filt   <- calculate_motif_frequencies(df_filt, grouping="subgroup", type_data="IG")
			result_data <- calculate_motif_frequencies(df_ig,   grouping="subgroup", type_data="IG")
			result_plot <- plot_motif_freq_subgroup_smart(data_filt, species_filter=sp_sel,
				interactive=TRUE, top_n=top_n_val, force_mode=sg_mode %||% "auto")
			stat_tests_res <- tryCatch(
				calculate_frequency_stat_tests(df_filt, grouping="subgroup", group_filter=NULL),
				error = function(e) NULL)
			
		} else if (grouping == "locus") {
			df_l <- if (!is.null(sp_locus_filter) && sp_locus_filter != "")
				df_ig %>% dplyr::filter(SPECIES == sp_locus_filter) else df_ig
			result_data <- calculate_motif_frequencies(df_l, grouping="locus", type_data="IG")
			result_plot <- tryCatch(
				plot_motif_frequencies(result_data, grouping="locus", interactive=TRUE),
				error = function(e) plotly::plot_ly() %>% plotly::layout(title=paste("Error:", e$message)))
			stat_tests_res <- tryCatch(
				calculate_frequency_stat_tests(df_l, grouping="locus", group_filter=NULL),
				error = function(e) NULL)
			# know exactly WHICH loci differ significantly 
			pairwise_tests_res <- tryCatch(
				calculate_locus_pairwise_tests(df_l),
				error = function(e) NULL)
			
		} else {  # species
			result_data <- calculate_motif_frequencies(df_ig, grouping="species", type_data="IG")
			result_plot <- plot_motif_frequencies(result_data, grouping="species", interactive=TRUE)
			stat_tests_res <- tryCatch(
				calculate_frequency_stat_tests(df_ig, grouping="species", group_filter=NULL),
				error = function(e) NULL)
			pairwise_tests_res <- NULL
		}
		
		if (grouping == "subgroup") pairwise_tests_res <- NULL
		
		freq_results$data           <- result_data
		freq_results$plot           <- result_plot
		freq_results$grouping       <- grouping
		freq_results$prev_grouping  <- grouping
		freq_results$df_snap        <- df_ig
		freq_results$stat_tests     <- stat_tests_res
		freq_results$pairwise_tests <- pairwise_tests_res
		freq_results$done           <- TRUE
	}, error = function(e) {
		showNotification(paste0("❌ Freq: ", e$message), type="error", duration=6)
	})
}

observe({
	req(store$df_final)
	grp          <- input$freq_grouping
	sp_sel       <- input$freq_species_filter
	locus_sel    <- input$freq_locus_filter
	sp_loc_filt  <- input$freq_locus_species_filter
	sg_mode      <- input$freq_sg_mode
	topn_val     <- input$freq_topn
	.run_freq(grp, sp_sel, locus_sel, sp_loc_filt, sg_mode, topn_val)
})

observeEvent(input$btn_launch_freq, {
	req(store$df_final)
	.run_freq(
		grouping     = input$freq_grouping,
		sp_sel       = input$freq_species_filter,
		locus_sel    = input$freq_locus_filter,
		sp_locus_filter = input$freq_locus_species_filter,
		sg_mode      = input$freq_sg_mode,
		topn_val     = input$freq_topn
	)
})

output$freq_sg_alert <- renderUI({
	req(freq_results$done, freq_results$grouping == "subgroup")
	n_sg <- if (!is.null(freq_results$data) && "SUBGROUP" %in% names(freq_results$data))
		length(unique(freq_results$data$SUBGROUP)) else 0
	top_n_val <- if (!is.na(input$freq_topn) && input$freq_topn > 0) input$freq_topn else NULL
	sg_alert_banner(n_sg, top_n_val)
})

output$freq_plot_ui <- renderUI({
	req(freq_results$done)
	grp <- freq_results$grouping
	h <- if (grp == "subgroup") {
		n_sg <- if (!is.null(freq_results$data) && "SUBGROUP" %in% names(freq_results$data))
			length(unique(freq_results$data$SUBGROUP)) else 20
		sg_plot_height_px(n_sg, min_h = 400)
	} else if (grp == "locus") {
		n_panels <- if (!is.null(freq_results$data) && "SPECIES" %in% names(freq_results$data) &&
				"LOCI" %in% names(freq_results$data))
			length(unique(paste(freq_results$data$SPECIES, freq_results$data$LOCI)))
		else 6
		min(max(450, n_panels * 220), 900)
	} else {
		550
	}
	plotly::plotlyOutput("freq_plot", height = paste0(h, "px"))
})

output$freq_plot <- renderPlotly({
	req(freq_results$plot, freq_results$done)
	p <- freq_results$plot
	if (inherits(p, "plotly")) return(p)
	plotly::ggplotly(p)
})

output$btn_download_freq_plot <- downloadHandler(
	filename = function() paste0("Freq_Motifs_", freq_results$grouping, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(freq_results$df_snap, freq_results$grouping)
		grp       <- freq_results$grouping
		top_n_val <- if (!is.na(input$freq_topn) && input$freq_topn > 0) input$freq_topn else NULL
		
		withProgress(message = "Generating PNG…", value = 0.5, {
			tryCatch({
				if (grp == "subgroup") {
					sp        <- input$freq_species_filter
					locus_sel <- input$freq_locus_filter
					df_filt   <- freq_results$df_snap
					if (!is.null(sp) && sp != "") df_filt <- df_filt %>% dplyr::filter(SPECIES == sp)
					if (!is.null(locus_sel) && locus_sel != "All loci") df_filt <- df_filt %>% dplyr::filter(LOCI == locus_sel)
					data <- calculate_motif_frequencies(df_filt, grouping="subgroup", type_data="IG")
					p    <- plot_motif_freq_subgroup_smart(data, species_filter=sp,
						interactive=FALSE, top_n=top_n_val, force_mode=input$freq_sg_mode %||% "auto")
					n_sg <- length(unique(data$SUBGROUP))
					h_in <- max(6, n_sg * 0.22 + 2)
					ggplot2::ggsave(file, plot=p, width=14, height=min(h_in,60),
						units="in", dpi=300, bg="white")
					
				} else if (grp == "locus") {
					data <- freq_results$data
					p    <- plot_motif_frequencies(data, grouping="locus", interactive=FALSE)
					n_panels <- if ("SPECIES" %in% names(data) && "LOCI" %in% names(data))
						length(unique(paste(data$SPECIES, data$LOCI))) else 6
					ggplot2::ggsave(file, plot=p, width=14, height=max(5, n_panels * 2.2 + 2),
						units="in", dpi=300, bg="white")
					
				} else {
					data <- freq_results$data
					p    <- plot_motif_frequencies(data, grouping="species", interactive=FALSE)
					ggplot2::ggsave(file, plot=p, width=12, height=7,
						units="in", dpi=300, bg="white")
				}
				incProgress(0.5)
			}, error = function(e) {
				showNotification(paste0("❌ PNG error: ", e$message), type="error", duration=8)
			})
		})
		showNotification("✅ PNG downloaded!", type="message", duration=3)
	}
)

output$btn_download_freq_species_png <- downloadHandler(
	filename = function() {
		sp        <- input$freq_species_filter
		locus_sel <- input$freq_locus_filter
		suffix    <- if (!is.null(locus_sel) && locus_sel != "All loci") paste0("_", locus_sel) else ""
		paste0("Freq_Motifs_subgroup_", sp, suffix, "_", Sys.Date(), ".png")
	},
	content = function(file) {
		req(freq_results$df_snap, input$freq_species_filter)
		sp        <- input$freq_species_filter
		locus_sel <- input$freq_locus_filter
		top_n_val <- if (!is.na(input$freq_topn) && input$freq_topn > 0) input$freq_topn else NULL
		
		withProgress(message = "Generating PNG…", value = 0.5, {
			tryCatch({
				df_filt <- freq_results$df_snap %>% dplyr::filter(SPECIES == sp)
				if (!is.null(locus_sel) && locus_sel != "All loci") df_filt <- df_filt %>% dplyr::filter(LOCI == locus_sel)
				data <- calculate_motif_frequencies(df_filt, grouping="subgroup", type_data="IG")
				p    <- plot_motif_freq_subgroup_smart(data, species_filter=sp,
					interactive=FALSE, top_n=top_n_val, force_mode=input$freq_sg_mode %||% "auto")
				n_sg <- length(unique(data$SUBGROUP))
				h_in <- max(6, n_sg * 0.22 + 2)
				ggplot2::ggsave(file, plot=p, width=14, height=min(h_in,60),
					units="in", dpi=300, bg="white")
				incProgress(0.5)
			}, error = function(e) {
				showNotification(paste0("❌ PNG error: ", e$message), type="error", duration=8)
			})
		})
		showNotification(paste0("✅ PNG for ", sp, " downloaded!"), type="message", duration=3)
	}
)

output$btn_download_freq_csv <- downloadHandler(
	filename = function() paste0("Freq_Motifs_", input$freq_grouping, "_", Sys.Date(), ".csv"),
	content  = function(file) { req(freq_results$data); write.csv(freq_results$data, file, row.names=FALSE) }
)

# ── Statistical tests output for frequencies ──────────────
output$freq_stat_tests_table <- renderDT({
	if (is.null(freq_results$stat_tests)) {
		return(datatable(
			data.frame(Message = "Run the analysis first to compute statistical tests."),
			rownames = FALSE, options = list(dom = "t")
		))
	}
	df <- freq_results$stat_tests
	numeric_cols <- intersect(c("Statistic","p.value","p.adjusted"), names(df))
	dt <- datatable(df,
		options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv")),
		extensions = "Buttons", rownames = FALSE,
		class = "cell-border stripe hover compact")
	if (length(numeric_cols) > 0)
		dt <- dt %>% formatRound(columns = numeric_cols, digits = 5)
	if ("Significant" %in% names(df)) {
		dt <- dt %>% formatStyle("Significant",
			backgroundColor = styleEqual(
				c("***","**","*","ns","—"),
				c("#D4EDDA","#D4EDDA","#FFF3CD","#F8D7DA","#F9FAFB")),
			fontWeight = "bold", color = styleEqual(
				c("***","**","*","ns","—"),
				c("#155724","#155724","#856404","#721c24","#6B7280")))
	}
	dt
})

output$btn_download_freq_stat_csv <- downloadHandler(
	filename = function() paste0("Freq_StatTests_", input$freq_grouping, "_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(freq_results$stat_tests)
		write.csv(freq_results$stat_tests, file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	}
)

# ── Pairwise post-hoc tests between loci  ──────
output$freq_pairwise_tests_table <- renderDT({
	if (is.null(freq_results$pairwise_tests) || freq_results$grouping != "locus") {
		return(datatable(
			data.frame(Message = "Switch the grouping level to 'By locus' and launch the analysis to see pairwise comparisons between loci."),
			rownames = FALSE, options = list(dom = "t")
		))
	}
	df <- freq_results$pairwise_tests
	numeric_cols <- intersect(c("Statistic","p.value","p.adjusted"), names(df))
	dt <- datatable(df,
		options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv")),
		extensions = "Buttons", rownames = FALSE,
		class = "cell-border stripe hover compact")
	if (length(numeric_cols) > 0)
		dt <- dt %>% formatRound(columns = numeric_cols, digits = 5)
	if ("Significant" %in% names(df)) {
		dt <- dt %>% formatStyle("Significant",
			backgroundColor = styleEqual(
				c("***","**","*","ns"),
				c("#D4EDDA","#D4EDDA","#FFF3CD","#F8D7DA")),
			fontWeight = "bold", color = styleEqual(
				c("***","**","*","ns"),
				c("#155724","#155724","#856404","#721c24")))
	}
	dt
})

output$btn_download_freq_pairwise_csv <- downloadHandler(
	filename = function() paste0("Freq_PairwiseLocusTests_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(freq_results$pairwise_tests)
		write.csv(freq_results$pairwise_tests, file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	}
)
