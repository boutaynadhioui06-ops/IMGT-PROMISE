# ══════════════════════════════════════════════════════════════
# 4. CONSERVATION / FONCTIONNALITÉ
# ══════════════════════════════════════════════════════════════

cons_func_results <- reactiveValues(data=NULL, plot=NULL, chi=NULL, stat_tests=NULL, done=FALSE, grouping=NULL, motif_filter=NULL, df_snap=NULL)

observe({
	req(store$df_final)
	isolate({
		bin_cols    <- names(store$df_final)[str_starts(names(store$df_final),"BIN_")]
		motif_names <- str_remove(bin_cols, "^BIN_")
		updateSelectInput(session, "cons_func_motif", choices=motif_names, selected=motif_names[1])
		sp_list <- sort(unique(store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV")) %>% pull(SPECIES)))
		updateSelectInput(session, "cons_func_species", choices=sp_list, selected=sp_list[1])
	})
})

observeEvent(input$btn_launch_cons_func, {
	req(store$df_final)
	isolate({
		withProgress(message = "Computing conservation…", value = 0.3, {
			tryCatch({
				df            <- store$df_final
				grouping      <- input$cons_func_grouping
				analysis_type <- input$cons_func_analysis_type
				motif_filter  <- if (analysis_type == "single") input$cons_func_motif else NULL
				
				if (grouping %in% c("locus","subgroup")) {
					req(input$cons_func_species)
					df_filtered <- df %>% dplyr::filter(SPECIES==input$cons_func_species)
				} else {
					df_filtered <- df
				}
				
				locus_sel_cf <- input$cons_func_locus_filter
				if (grouping == "subgroup" && !is.null(locus_sel_cf) && locus_sel_cf != "" && locus_sel_cf != "All loci") {
					df_filtered <- df_filtered %>% dplyr::filter(LOCI==locus_sel_cf)
					if (nrow(df_filtered) == 0) { showNotification("⚠️ No data for this locus.", type="warning"); return() }
				}
				if (grouping == "subgroup" && length(unique(na.omit(df_filtered$SUBGROUP))) == 0) return()
				
				incProgress(0.5, detail = "Running analysis…")
				results <- conservation_by_functionality(df_filtered, grouping=grouping, type_data="IG", motif_filter=motif_filter)
				
				if (is.null(results$data) || nrow(results$data) == 0) {
					showNotification("⚠️ No data available for this selection.", type="warning", duration=5)
					return()
				}
				
				result_plot <- plot_conservation_functionality(results$data, grouping=grouping, interactive=TRUE, motif_filter=motif_filter)
				
				cons_func_results$data         <- results$data
				cons_func_results$plot         <- result_plot
				cons_func_results$chi          <- results$chi_tests
				cons_func_results$stat_tests   <- results$stat_tests
				cons_func_results$grouping     <- grouping
				cons_func_results$motif_filter <- motif_filter
				cons_func_results$df_snap      <- df_filtered
				cons_func_results$done         <- TRUE
				incProgress(1)
				
			}, error=function(e) {
				showNotification(paste0("❌ Conservation: ", e$message), type="error", duration=8)
			})
		})
	})
})

output$cons_func_status <- renderUI({
	if (!cons_func_results$done)
		div(class="alert alert-info", icon("info-circle"),
			" Configure parameters and click ", strong("🚀 Launch analysis"), ".")
	else
		div(class="alert alert-success", icon("check-circle"), " Analysis complete!")
})

output$cons_func_plot_ui <- renderUI({
	req(cons_func_results$done)
	grp  <- cons_func_results$grouping
	data <- cons_func_results$data
	
	if (grp == "subgroup") {
		df_snap   <- cons_func_results$df_snap
		locus_sel <- input$cons_func_locus_filter
		if (!is.null(df_snap)) {
			df_for_count <- df_snap
			if (!is.null(locus_sel) && locus_sel != "" && locus_sel != "All loci")
				df_for_count <- df_for_count %>% dplyr::filter(LOCI == locus_sel)
			n_sg <- length(unique(df_for_count$SUBGROUP))
		} else {
			n_sg <- length(unique(data$SUBGROUP))
		}
		top_n_val <- if (!is.na(input$cons_func_topn) && input$cons_func_topn > 0) input$cons_func_topn else NULL
		if (!is.null(top_n_val) && top_n_val < n_sg) n_sg <- top_n_val
		h <- sg_plot_height_px(n_sg, min_h = 400)
		tagList(
			sg_alert_banner(n_sg, top_n_val),
			plotlyOutput("cons_func_plot", height=paste0(h,"px"))
		)
	} else if (grp == "locus") {
		n_sp_v <- length(unique(data$SPECIES))
		n_col  <- if (n_sp_v <= 3) 1 else 2
		n_rows <- ceiling(n_sp_v / n_col)
		h      <- if (n_sp_v == 1) 500 else max(500, n_rows * 500)
		plotlyOutput("cons_func_plot", height=paste0(min(h, 900),"px"))
	} else {
		n_sp_v <- length(unique(data$SPECIES))
		plotlyOutput("cons_func_plot", height=paste0(min(max(450, n_sp_v*120), 650),"px"))
	}
})


output$cons_func_plot <- renderPlotly({
	req(cons_func_results$data, cons_func_results$grouping)
	grp <- cons_func_results$grouping
	motif_filter <- cons_func_results$motif_filter
	data <- cons_func_results$data
	
	if (grp == "subgroup") {
		top_n_val  <- if (!is.na(input$cons_func_topn) && input$cons_func_topn > 0) input$cons_func_topn else NULL
		force_mode <- input$cons_func_sg_mode %||% "auto"
		n_sg       <- length(unique(data$SUBGROUP))
		if (!is.null(top_n_val) && is.numeric(top_n_val) && top_n_val > 0 && top_n_val < n_sg) {
			pct_cols <- names(data)[str_starts(names(data),"Pct_")]
			top_sgs  <- data %>% dplyr::group_by(SUBGROUP) %>%
				dplyr::summarise(mean_pct=mean(dplyr::across(dplyr::all_of(pct_cols),~.x),na.rm=TRUE),.groups="drop") %>%
				dplyr::arrange(dplyr::desc(mean_pct)) %>% dplyr::slice_head(n=top_n_val) %>% dplyr::pull(SUBGROUP)
			data <- data %>% dplyr::filter(SUBGROUP %in% top_sgs)
		}
		p <- plot_conservation_functionality(data, grouping=grp, interactive=TRUE,
			motif_filter=motif_filter, force_mode=force_mode)
	} else {
		p <- plot_conservation_functionality(data, grouping=grp, interactive=TRUE, motif_filter=motif_filter)
	}
	
	if ("plotly" %in% class(p)) return(p)
	if ("ggplot" %in% class(p)) return(ggplotly(p))
})


output$cons_func_chi_table <- renderDT({
	req(cons_func_results$chi)
	df <- cons_func_results$chi
	numeric_cols <- intersect(c('Chi2','p.value'),names(df))
	dt <- datatable(df, options=list(pageLength=25,scrollX=TRUE), rownames=FALSE, class='cell-border stripe hover')
	if (length(numeric_cols) > 0) dt <- dt %>% formatRound(columns=numeric_cols, digits=4)
	dt %>% formatStyle('Significant', backgroundColor=styleEqual(c('***','ns'),c('#D4EDDA','#F8D7DA')))
})

# ── Tests statistiques avancés 
output$cons_func_stat_tests_table <- renderDT({
	req(cons_func_results$stat_tests)
	df <- cons_func_results$stat_tests
	if (is.null(df) || nrow(df) == 0) return(DT::datatable(data.frame(Message="No tests available"), rownames=FALSE))
	numeric_cols <- intersect(c('Statistic','p.value','p.adjusted'), names(df))
	dt <- datatable(df,
		options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv")),
		extensions = "Buttons", rownames = FALSE,
		class = "cell-border stripe hover compact")
	if (length(numeric_cols) > 0) dt <- dt %>% formatRound(columns = numeric_cols, digits = 4)
	if ("Significant" %in% names(df)) {
		dt <- dt %>% formatStyle("Significant",
			backgroundColor = styleEqual(
				c("***","**","*","ns","—"),
				c("#D4EDDA","#D4EDDA","#FFF3CD","#F8D7DA","#F9FAFB")),
			fontWeight = "bold")
	}
	if ("Test" %in% names(df)) {
		dt <- dt %>% formatStyle("Test", color = "#1565C0", fontWeight = "bold")
	}
	dt
})

output$btn_download_cons_stat_csv <- downloadHandler(
	filename = function() paste0("ConservStat_Tests_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(cons_func_results$stat_tests)
		write.csv(cons_func_results$stat_tests, file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	}
)

output$btn_download_cons_func_plot <- downloadHandler(
	filename = function() paste0("ConservFunc_", cons_func_results$grouping, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(cons_func_results$data)
		grp          <- cons_func_results$grouping
		data         <- cons_func_results$data
		motif_filter <- cons_func_results$motif_filter
		
		if (grp == "subgroup") {
			top_n_val <- if (!is.na(input$cons_func_topn) && input$cons_func_topn > 0) input$cons_func_topn else NULL
			if (!is.null(top_n_val)) {
				pct_cols <- names(data)[str_starts(names(data),"Pct_")]
				top_sgs  <- data %>% dplyr::group_by(SUBGROUP) %>%
					dplyr::summarise(m=mean(dplyr::across(dplyr::all_of(pct_cols),~.x),na.rm=TRUE),.groups="drop") %>%
					dplyr::arrange(dplyr::desc(m)) %>% dplyr::slice_head(n=top_n_val) %>% dplyr::pull(SUBGROUP)
				data <- data %>% dplyr::filter(SUBGROUP %in% top_sgs)
			}
		}
		n_sg     <- if ("SUBGROUP" %in% names(data)) length(unique(data$SUBGROUP)) else 0
		n_sp_v   <- length(unique(data$SPECIES))
		n_loci   <- if ("LOCI" %in% names(data)) length(unique(data$LOCI)) else 1
		p <- plot_conservation_functionality(data, grouping=grp, interactive=FALSE, motif_filter=motif_filter)
		if (grp == "subgroup") {
			dims  <- sg_export_dims(n_sg, n_col=1)
			h_in  <- dims$height; w_in <- dims$width
		} else if (grp == "locus") {
			h_in <- max(8, n_sp_v*n_loci*2.8+2); w_in <- 16
		} else {
			h_in <- max(7, n_sp_v*2.5+2); w_in <- 14
		}
		ggsave(file, plot=p, width=min(w_in,40), height=min(h_in,80), units="in", dpi=300, bg="white")
		showNotification("✅ PNG downloaded!", type="message", duration=3)
	}
)

output$btn_download_cons_func_csv <- downloadHandler(
	filename = function() paste0("ConservFunc_", input$cons_func_grouping, "_", Sys.Date(), ".csv"),
	content  = function(file) { req(cons_func_results$data); write.csv(cons_func_results$data, file, row.names=FALSE) }
)

