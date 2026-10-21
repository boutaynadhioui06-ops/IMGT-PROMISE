# ══════════════════════════════════════════════════════════════
# 3. LOGOS DE SÉQUENCES 
# ══════════════════════════════════════════════════════════════

logo_results <- reactiveValues(data=NULL, plot=NULL, done=FALSE)

observe({
	req(store$df_final)
	isolate({
		df_ig <- store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
		seq_cols    <- names(df_ig)[str_starts(names(df_ig),"SEQ_")]
		motif_names <- str_remove(seq_cols, "^SEQ_")
		motif_names <- motif_names[sapply(motif_names, function(m) {
			col <- paste0("SEQ_",m); sum(!is.na(df_ig[[col]]) & df_ig[[col]] != "" & df_ig[[col]] != "absent") > 0
		})]
		updateSelectInput(session, "logo_motif", choices=motif_names, selected=motif_names[1])
		sp_list <- sort(unique(df_ig$SPECIES))
		updateSelectInput(session, "logo_species", choices=c("All species"="",setNames(sp_list,sp_list)), selected="")
		loci_list <- sort(unique(df_ig$LOCI))
		updateSelectInput(session, "logo_locus", choices=c("All loci"="",setNames(loci_list,loci_list)), selected="")
	})
})

observeEvent(input$btn_launch_logo, {
	req(store$df_final, input$logo_motif, nchar(input$logo_motif) > 0)
	isolate({
		withProgress(message = "Generating sequence logos…", value = 0.2, {
			tryCatch({
				df <- store$df_final; motif_name <- input$logo_motif
				seqs <- prepare_sequences_for_logo(df, motif_name, "IG")
				if (!is.null(input$logo_species) && input$logo_species != "") seqs <- seqs %>% dplyr::filter(SPECIES==input$logo_species)
				if (!is.null(input$logo_locus)   && input$logo_locus   != "") seqs <- seqs %>% dplyr::filter(LOCI==input$logo_locus)
				if (nrow(seqs) == 0) { showNotification("⚠️ No sequences found.", type="warning"); return() }
				
				incProgress(0.4, detail = "Building logos…")
				library(ggseqlogo)
				seq_counts <- seqs %>% group_by(SPECIES,LOCI) %>% summarise(n=n(),.groups="drop")
				valid_groups <- seq_counts %>% dplyr::filter(n >= 3)
				if (nrow(valid_groups) == 0) return()
				
				seq_list <- seqs %>%
					inner_join(valid_groups %>% select(SPECIES,LOCI), by=c("SPECIES","LOCI")) %>%
					group_by(SPECIES,LOCI) %>% summarise(Sequences=list(Sequence),.groups="drop")
				
				plots <- list()
				for (i in 1:nrow(seq_list)) {
					species   <- seq_list$SPECIES[i]; locus <- seq_list$LOCI[i]
					sequences <- toupper(trimws(seq_list$Sequences[[i]]))
					
					if (grepl("TATA", motif_name, ignore.case = TRUE)) {
						sequences <- substr(sequences, 1, 7)
					}
					sequences <- sequences[grepl("^[ATGC]+$",sequences)]
					if (length(sequences) < 3) next
					seq_lengths <- nchar(sequences)
					if (length(unique(seq_lengths)) > 1) {
						most_common <- as.numeric(names(sort(table(seq_lengths),decreasing=TRUE))[1])
						sequences   <- sequences[seq_lengths == most_common]
						if (length(sequences) < 3) next
					}
					tryCatch({
						p <- ggseqlogo(sequences, method='prob', seq_type='dna') +
							ggtitle(paste(species,"-",locus)) + theme_minimal() +
							theme(plot.title=element_text(face="bold",hjust=0.5))
						plots[[paste(species,locus,sep="_")]] <- p
					}, error=function(e) NULL)
				}
				
				if (length(plots) == 0) return()
				library(patchwork)
				logo_results$data <- seqs
				n_plots <- length(plots)
				ncols   <- if (n_plots == 1) 1 else 2
				logo_results$plot <- wrap_plots(plots, ncol = ncols)
				logo_results$done <- TRUE
				incProgress(1)
				
			}, error=function(e) {
				showNotification(paste0("❌ Logos: ", e$message), type="error", duration=8)
			})
		})
	})
})

output$logo_status <- renderUI({
	if (!logo_results$done)
		div(class="alert alert-info", icon("info-circle"),
			" Select a motif and click ", strong("🚀 Generate logos"), ".")
	else
		div(class="alert alert-success", icon("check-circle"), " Logos generated for ", strong(input$logo_motif))
})

output$logo_plot <- renderPlot({ req(logo_results$plot); print(logo_results$plot) })
output$logo_plot_ui <- renderUI({
	req(logo_results$data)
	n_combos <- tryCatch(
		length(unique(paste(logo_results$data$SPECIES, logo_results$data$LOCI))),
		error = function(e) 2
	)
	n_cols  <- if (n_combos == 1) 1 else 2
	h_px    <- max(300, ceiling(n_combos / n_cols) * 300)
	w_pct   <- if (n_combos == 1) "55%" else "100%"
	div(
		style = "display:flex; justify-content:center;",
		plotOutput("logo_plot", height = paste0(h_px, "px"), width = w_pct)
	)
})
output$logo_sequences_table <- renderDT({
	req(logo_results$data)
	datatable(logo_results$data, options=list(pageLength=25,scrollX=TRUE), rownames=FALSE, class='cell-border stripe hover')
})

output$btn_download_logo_plot <- downloadHandler(
	filename = function() paste0("Logos_",input$logo_motif,"_",Sys.Date(),".png"),
	content  = function(file) { req(logo_results$plot); ggsave(file,plot=logo_results$plot,width=14,height=10,dpi=300) }
)
output$btn_download_logo_csv <- downloadHandler(
	filename = function() paste0("Sequences_",input$logo_motif,"_",Sys.Date(),".csv"),
	content  = function(file) { req(logo_results$data); write.csv(logo_results$data,file,row.names=FALSE) }
)