# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  HEATMAP D'IDENTITÉ  ██████████████████ 
# ═══════════════════════════════════════════════════════════════════════════

identity_results <- reactiveValues(data=NULL, plot=NULL, mode="species", done=FALSE)

observe({
	req(store$df_final)
	isolate({
		df_ig <- store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
		seq_cols    <- names(df_ig)[str_starts(names(df_ig),"SEQ_")]
		motif_names <- str_remove(seq_cols, "^SEQ_")
		motif_names <- motif_names[sapply(motif_names, function(m) {
			col <- paste0("SEQ_",m); sum(!is.na(df_ig[[col]]) & df_ig[[col]] != "" & df_ig[[col]] != "absent") > 0
		})]
		updateSelectInput(session, "identity_motif", choices=motif_names, selected=motif_names[1])
		sp_list <- sort(unique(df_ig$SPECIES))
		updateSelectInput(session, "identity_species", choices=sp_list, selected=sp_list[1])
	})
})

observeEvent(input$btn_launch_identity, {
	req(store$df_final, input$identity_motif, nchar(input$identity_motif) > 0)
	isolate({
		withProgress(message = "Computing identity heatmap…", value = 0.4, {
			tryCatch({
				df <- store$df_final
				motif_name <- input$identity_motif
				mode       <- input$identity_mode
				species    <- input$identity_species
				if (mode == "species") {
					result_plot <- plot_sequence_identity_heatmap(df, motif_name, interactive=TRUE)
					result_data <- calculate_sequence_identity(df, motif_name)
				} else if (mode == "locus") {
					result_plot <- plot_locus_identity_heatmap(df, motif_name, species, interactive=TRUE)
					result_data <- calculate_locus_identity(df, motif_name, species)
				} else {
					result_plot <- plot_subgroup_identity_heatmap(df, motif_name, species, interactive=TRUE)
					result_data <- calculate_subgroup_identity(df, motif_name, species)
				}
				identity_results$data <- result_data
				identity_results$plot <- result_plot
				identity_results$mode <- mode
				identity_results$done <- TRUE
				incProgress(1)
			}, error=function(e) {
				showNotification(paste0("❌ Identity: ", e$message), type="error", duration=8)
			})
		})
	})
})

output$identity_status <- renderUI({
	if (!identity_results$done)
		div(class="alert alert-info", icon("info-circle"),
			" Select a motif and mode, then click ", strong("🚀 Generate heatmap"), ".")
	else {
		mode_label <- switch(identity_results$mode,
			species  = "between species",
			locus    = paste0("between loci — ", input$identity_species),
			subgroup = paste0("between subgroups — ", input$identity_species))
		div(class="alert alert-success", icon("check-circle"), " Heatmap generated — motif: ",
			strong(input$identity_motif), " | ", mode_label)
	}
})

output$identity_plot_ui <- renderUI({
	req(identity_results$done)
	if (identity_results$mode == "subgroup") {
		n_sg <- length(unique(identity_results$data[[1]]))
		plotlyOutput("identity_plot", height=paste0(max(450,n_sg*28),"px"))
	} else {
		plotlyOutput("identity_plot", height="480px")
	}
})

output$identity_plot <- renderPlotly({ req(identity_results$plot); identity_results$plot })

output$identity_table <- renderDT({
	req(identity_results$data)
	mode <- identity_results$mode
	col2 <- switch(mode, species="Species_2", locus="Locus_2", subgroup="Subgroup_2")
	datatable(
		identity_results$data %>% mutate(Identity=paste0(Identity,"%")) %>%
			pivot_wider(names_from=all_of(col2), values_from=Identity),
		options=list(pageLength=25,scrollX=TRUE), rownames=FALSE, class="cell-border stripe hover"
	)
})

output$btn_download_identity_png <- downloadHandler(
	filename = function() paste0("Identity_",input$identity_motif,"_",input$identity_mode,"_",Sys.Date(),".png"),
	content  = function(file) {
		req(store$df_final, input$identity_motif)
		df   <- store$df_final; mode <- input$identity_mode; sp <- input$identity_species
		if (mode == "species") {
			p <- plot_sequence_identity_heatmap(df, input$identity_motif, interactive=FALSE)
			n <- length(unique(df$SPECIES[df$LOCI %in% c("IGHV","IGKV","IGLV")]))
			w <- max(8,n*1.2); h <- max(6,n*1.0)
		} else if (mode == "locus") {
			p <- plot_locus_identity_heatmap(df, input$identity_motif, sp, interactive=FALSE); w <- 7; h <- 5
		} else {
			p <- plot_subgroup_identity_heatmap(df, input$identity_motif, sp, interactive=FALSE)
			n_sg <- tryCatch(length(unique(calculate_subgroup_identity(df,input$identity_motif,sp)$Subgroup_1)),error=function(e)10)
			w <- max(10,n_sg*0.5); h <- max(8,n_sg*0.4)
		}
		ggsave(file, plot=p, width=min(w,40), height=min(h,60), dpi=300, bg="white")
		showNotification("✅ PNG downloaded!", type="message", duration=3)
	}
)

output$btn_download_identity_csv <- downloadHandler(
	filename = function() paste0("Identity_",input$identity_motif,"_",input$identity_mode,"_",Sys.Date(),".csv"),
	content  = function(file) { req(identity_results$data); write.csv(identity_results$data, file, row.names=FALSE) }
)
