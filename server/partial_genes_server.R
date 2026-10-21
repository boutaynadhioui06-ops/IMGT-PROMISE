
# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  GÈNES PARTIELS  ████████████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════

partial_lit     <- reactiveVal(NULL)   
partial_compare <- reactiveVal(NULL)   

# ── Launch analysis ────────────────────────────────────────────────
partial_desc_results <- reactiveValues(done = FALSE)

observeEvent(store$df_no_motif, {
	partial_desc_results$done <- FALSE
})

observeEvent(input$btn_launch_partial, {
	req(store$df_final, nrow(store$df_final) > 0)
	withProgress(message = "Analyzing motif presence…", value = 0.3, {
		incProgress(0.6, detail = "Computing summaries…")
		partial_desc_results$done <- TRUE
		incProgress(1)
	})
})

output$partial_desc_status <- renderUI({
	if (is.null(store$df_final) || nrow(store$df_final) == 0) return(NULL)
	if (!partial_desc_results$done)
		div(class = "alert alert-info", icon("info-circle"),
			" Optionally filter by species/locus/subgroup below, then click ",
			strong("🚀 Launch analysis"), " to display results.")
	else
		div(class = "alert alert-success", icon("check-circle"),
			" Analysis complete — ", nrow(store$df_no_motif), " sequences without motif, ",
			nrow(store$df_with_motif), " sequences with motif(s).")
})

output$warn_no_partiels <- renderUI({
	if (is.null(store$df_final) || nrow(store$df_final) == 0)
		div(class = "alert alert-warning",
			icon("exclamation-triangle"), strong(" Warning : "),
			"Build the dataset first (📊 Final Dataset tab).")
})

# ══════════════════════════════════════════════════════════════
# SOUS-ONGLET A : PRÉSENCE / ABSENCE DE MOTIF RÉGULATEUR
# ══════════════════════════════════════════════════════════════

# ── Filtres en cascade : Species Locus Subgroup
observe({
	req(store$df_final)
	df_ig <- store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	sp_list <- sort(unique(df_ig$SPECIES))
	updateSelectInput(session, "partial_species_filter",
		choices = c("All species" = "All species", setNames(sp_list, sp_list)),
		selected = "All species")
})

observeEvent(input$partial_species_filter, {
	req(store$df_final)
	df_ig <- store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	if (!is.null(input$partial_species_filter) &&
			input$partial_species_filter != "All species")
		df_ig <- df_ig %>% dplyr::filter(SPECIES == input$partial_species_filter)
	loci_list <- sort(unique(df_ig$LOCI))
	updateSelectInput(session, "partial_locus_filter",
		choices = c("All loci" = "All loci", setNames(loci_list, loci_list)),
		selected = "All loci")
}, ignoreNULL = FALSE)

observeEvent(list(input$partial_species_filter, input$partial_locus_filter), {
	req(store$df_final)
	df_ig <- store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	if (!is.null(input$partial_species_filter) &&
			input$partial_species_filter != "All species")
		df_ig <- df_ig %>% dplyr::filter(SPECIES == input$partial_species_filter)
	if (!is.null(input$partial_locus_filter) &&
			input$partial_locus_filter != "All loci")
		df_ig <- df_ig %>% dplyr::filter(LOCI == input$partial_locus_filter)
	sg_list <- sort(unique(na.omit(df_ig$SUBGROUP)))
	updateSelectInput(session, "partial_subgroup_filter",
		choices = c("All subgroups" = "All subgroups", setNames(sg_list, sg_list)),
		selected = "All subgroups")
}, ignoreNULL = FALSE)

# ── Données filtrées  
filtered_no_motif <- reactive({
	req(partial_desc_results$done, store$df_no_motif)
	apply_partial_cascade_filters(store$df_no_motif,
		input$partial_species_filter, input$partial_locus_filter, input$partial_subgroup_filter)
})
filtered_with_motif <- reactive({
	req(partial_desc_results$done, store$df_with_motif)
	apply_partial_cascade_filters(store$df_with_motif,
		input$partial_species_filter, input$partial_locus_filter, input$partial_subgroup_filter)
})

# ── Value boxes
output$vbox_partial_total <- renderValueBox({
	n <- if (isTRUE(partial_desc_results$done)) nrow(store$df_no_motif %||% data.frame()) else 0
	valueBox(n, "Sequences with NO motif", icon = icon("ban"), color = "red")
})
output$vbox_partial_with_motif <- renderValueBox({
	n <- if (isTRUE(partial_desc_results$done)) nrow(store$df_with_motif %||% data.frame()) else 0
	valueBox(n, "Sequences WITH motif(s)", icon = icon("check"), color = "green")
})
output$vbox_partial_species <- renderValueBox({
	n <- if (isTRUE(partial_desc_results$done) && !is.null(store$df_no_motif) && nrow(store$df_no_motif) > 0)
		length(unique(store$df_no_motif$SPECIES)) else 0
	valueBox(n, "Species concerned (no motif)", icon = icon("globe"), color = "orange")
})

# ── Tableau 1 : résumé — sans motif
output$partial_summary_nomotif_table <- renderDT({
	req(partial_desc_results$done)
	df_sum <- summary_no_motif_by_species(store$df_no_motif)
	if (nrow(df_sum) == 0) return(datatable(data.frame(Message = "No data"), rownames = FALSE))
	df_sum <- rename_loci_display(df_sum)
	datatable(df_sum,
		options = list(pageLength = 30, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv","excel")),
		extensions = "Buttons", rownames = FALSE,
		class = "cell-border stripe hover compact",
		colnames = c("Species","Group","N alleles","F","P","ORF","5'UTR Mean length (nt)","5'UTR Min length (nt)","5'UTR Max length")
	) %>%
		formatStyle("GROUP",
			backgroundColor = styleEqual(c("IGHV","IGKV","IGLV"), c("#DBEAFE","#DCFCE7","#FEF9C3")),
			fontWeight = "bold")
})

# ── Tableau 2 : liste détaillée — sans motif
output$partial_alleles_nomotif_table <- renderDT({
	df_list <- filtered_no_motif()
	if (is.null(df_list) || nrow(df_list) == 0) return(datatable(data.frame(Message = "No data"), rownames = FALSE))
	df_list <- rename_loci_display(df_list)
	df_list <- rename_loci_display(df_list)
	if ("LENGTH_NT" %in% names(df_list))
		names(df_list)[names(df_list) == "LENGTH_NT"] <- "5'UTR_LENGTH_NT"
	datatable(df_list,
		options = list(pageLength = 50, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv","excel")),
		extensions = "Buttons", rownames = FALSE, filter = "top",
		class = "cell-border stripe hover compact"
	) %>%
		formatStyle("FUNCTIONALITY",
			backgroundColor = styleEqual(c("F","P","ORF"), c("#D4EDDA","#F8D7DA","#FFF3CD"))) %>%
		formatStyle("GROUP", fontWeight = "bold")
})

# ── Tableau 3 : résumé avec  motif
output$partial_summary_withmotif_table <- renderDT({
	req(partial_desc_results$done)
	df_sum <- summary_with_motif_by_species(store$df_with_motif)
	if (nrow(df_sum) == 0) return(datatable(data.frame(Message = "No data"), rownames = FALSE))
	df_sum <- rename_loci_display(df_sum)
	datatable(df_sum,
		options = list(pageLength = 30, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv","excel")),
		extensions = "Buttons", rownames = FALSE,
		class = "cell-border stripe hover compact",
		colnames = c("Species","GROUP","N alleles","Mean N motifs","Min","Max","5'UTR Mean length (nt)")
	) %>%
		formatStyle("GROUP",
			backgroundColor = styleEqual(c("IGHV","IGKV","IGLV"), c("#DBEAFE","#DCFCE7","#FEF9C3")),
			fontWeight = "bold")
})

# ── Tableau 4 : liste détaillée avec motif 
output$partial_alleles_withmotif_table <- renderDT({
	df_list <- filtered_with_motif()
	if (is.null(df_list) || nrow(df_list) == 0) return(datatable(data.frame(Message = "No data"), rownames = FALSE))
	df_list <- rename_loci_display(df_list)
	df_list <- rename_loci_display(df_list)
	if ("LENGTH_NT" %in% names(df_list))
		names(df_list)[names(df_list) == "LENGTH_NT"] <- "5'UTR_LENGTH_NT"
	dt <- datatable(df_list,
		options = list(pageLength = 50, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv","excel")),
		extensions = "Buttons", rownames = FALSE, filter = "top",
		class = "cell-border stripe hover compact"
	) %>%
		formatStyle("GROUP", fontWeight = "bold") %>%
		formatStyle("N_MOTIFS",
			background = styleColorBar(c(0, max(df_list$N_MOTIFS, na.rm = TRUE)), "#93C5FD"),
			backgroundSize = "100% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")
	if ("FUNCTIONALITY" %in% names(df_list))
		dt <- dt %>% formatStyle("FUNCTIONALITY",
			backgroundColor = styleEqual(c("F","P","ORF"), c("#D4EDDA","#F8D7DA","#FFF3CD")))
	dt
})

# ── Graphiques : sans motif
output$partial_plot_nomotif_species_ui <- renderUI({
	df <- filtered_no_motif(); req(df)
	n_sp <- length(unique(df$SPECIES))
	plotlyOutput("partial_plot_nomotif_species",
		height = paste0(min(max(420, n_sp*42), 1400) + 60, "px"))  
})
output$partial_plot_nomotif_species <- renderPlotly({
	plot_no_motif_by_species(filtered_no_motif(), interactive = TRUE)
})
output$partial_plot_nomotif_length_ui <- renderUI({
	df <- filtered_no_motif(); req(df)
	n_loci <- length(unique(df$LOCI))
	plotlyOutput("partial_plot_nomotif_length", height = paste0(if (n_loci<=1) 420 else 460, "px"))
})
output$partial_plot_nomotif_length <- renderPlotly({
	plot_no_motif_length_dist(filtered_no_motif(), interactive = TRUE)
})

# ── Graphiques : avec motif
output$partial_plot_withmotif_count_ui <- renderUI({
	req(filtered_with_motif())
	plotlyOutput("partial_plot_withmotif_count", height = "520px")
})
output$partial_plot_withmotif_count <- renderPlotly({
	plot_with_motif_count(filtered_with_motif(), interactive = TRUE)
})
output$partial_plot_withmotif_dist_ui <- renderUI({
	req(filtered_with_motif())
	plotlyOutput("partial_plot_withmotif_dist", height = "500px")
})
output$partial_plot_withmotif_dist <- renderPlotly({
	plot_with_motif_distances(filtered_with_motif(), interactive = TRUE)
})

# Downloads
output$btn_dl_partial_summary_nomotif_csv <- downloadHandler(
	filename = function() paste0("NoMotif_Summary_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(store$df_no_motif)
		write.csv(summary_no_motif_by_species(store$df_no_motif), file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	})
output$btn_dl_partial_alleles_nomotif_csv <- downloadHandler(
	filename = function() paste0("NoMotif_Alleles_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(filtered_no_motif())
		write.csv(filtered_no_motif(), file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	})
output$btn_dl_partial_summary_withmotif_csv <- downloadHandler(
	filename = function() paste0("WithMotif_Summary_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(store$df_with_motif)
		write.csv(summary_with_motif_by_species(store$df_with_motif), file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	})
output$btn_dl_partial_alleles_withmotif_csv <- downloadHandler(
	filename = function() paste0("WithMotif_Alleles_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(filtered_with_motif())
		write.csv(filtered_with_motif(), file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	})
output$btn_dl_partial_plot_nomotif_species_png <- downloadHandler(
	filename = function() paste0("NoMotif_BySpecies_", Sys.Date(), ".png"),
	content  = function(file) {
		df <- filtered_no_motif(); req(df)
		p <- plot_no_motif_by_species(df, interactive = FALSE)
		n_sp <- length(unique(df$SPECIES))
		ggplot2::ggsave(file, plot = p, width = 13, height = min(max(6, n_sp*0.5+2), 30), dpi = 300, bg = "white")
		showNotification("✅ PNG downloaded!", type = "message", duration = 3)
	})
output$btn_dl_partial_plot_nomotif_length_png <- downloadHandler(
	filename = function() paste0("NoMotif_LengthDist_", Sys.Date(), ".png"),
	content  = function(file) {
		p <- plot_no_motif_length_dist(filtered_no_motif(), interactive = FALSE)
		ggplot2::ggsave(file, plot = p, width = 12, height = 7, dpi = 300, bg = "white")
		showNotification("✅ PNG downloaded!", type = "message", duration = 3)
	})
output$btn_dl_partial_plot_withmotif_count_png <- downloadHandler(
	filename = function() paste0("WithMotif_Count_", Sys.Date(), ".png"),
	content  = function(file) {
		p <- plot_with_motif_count(filtered_with_motif(), interactive = FALSE)
		ggplot2::ggsave(file, plot = p, width = 12, height = 7, dpi = 300, bg = "white")
		showNotification("✅ PNG downloaded!", type = "message", duration = 3)
	})
output$btn_dl_partial_plot_withmotif_dist_png <- downloadHandler(
	filename = function() paste0("WithMotif_Distances_", Sys.Date(), ".png"),
	content  = function(file) {
		p <- plot_with_motif_distances(filtered_with_motif(), interactive = FALSE)
		ggplot2::ggsave(file, plot = p, width = 10, height = 7, dpi = 300, bg = "white")
		showNotification("✅ PNG downloaded!", type = "message", duration = 3)
	})


# ══════════════════════════════════════════════════════════════
# SOUS-ONGLET B : COMPARAISON AVEC SÉQUENCES LITTÉRATURE
# ══════════════════════════════════════════════════════════════

# ── Import FASTA littérature ──────────────────────────────────
# Source 1 : copier-coller texte
observeEvent(input$btn_parse_lit_text, {
	req(nchar(trimws(input$lit_fasta_text)) > 10)
	tryCatch({
		df_lit <- parse_literature_fasta(input$lit_fasta_text)
		if (nrow(df_lit) == 0) {
			showNotification("⚠️ No valid IG V-gene sequences found in the pasted FASTA.",
				type = "warning", duration = 6)
			return()
		}
		partial_lit(df_lit)
		showNotification(
			HTML(paste0("<strong>✅ Literature FASTA parsed</strong><br>",
				"• ", nrow(df_lit), " sequences<br>",
				"• ", length(unique(df_lit$ALLELE)), " alleles<br>",
				"• ", length(unique(df_lit$SPECIES)), " species")),
			type = "message", duration = 6)
	}, error = function(e) {
		showNotification(paste0("❌ Parsing error: ", e$message),
			type = "error", duration = 8)
	})
})

# Source 2 : upload fichier )
observeEvent(input$btn_parse_lit_file, {
	req(input$lit_fasta_file)
	tryCatch({
		txt    <- paste(readLines(input$lit_fasta_file$datapath, warn = FALSE),
			collapse = "\n")
		df_lit <- parse_literature_fasta(txt)
		if (nrow(df_lit) == 0) {
			showNotification("⚠️ No valid IG V-gene sequences found in the file.",
				type = "warning", duration = 6)
			return()
		}
		partial_lit(df_lit)
		showNotification(
			HTML(paste0("<strong>✅ Literature FASTA loaded</strong><br>",
				"• ", nrow(df_lit), " sequences<br>",
				"• ", length(unique(df_lit$ALLELE)), " alleles")),
			type = "message", duration = 5)
	}, error = function(e) {
		showNotification(paste0("❌ File error: ", e$message),
			type = "error", duration = 8)
	})
})

# ── Statut import littérature
output$lit_import_status <- renderUI({
	df_lit <- partial_lit()
	if (is.null(df_lit) || nrow(df_lit) == 0) {
		div(class = "alert alert-info",
			icon("info-circle"),
			" Paste your FASTA sequences or upload a file above, then click ",
			strong("Parse FASTA"), ".")
	} else {
		n_alleles  <- length(unique(df_lit$ALLELE))
		n_seqs     <- nrow(df_lit)
		n_species  <- length(unique(df_lit$SPECIES))
		div(class = "alert alert-success",
			icon("check-circle"),
			strong(paste0(" ", n_seqs, " sequences loaded")),
			paste0(" | ", n_alleles, " alleles | ", n_species, " species"),
			br(),
			tags$small(
				icon("info-circle"),
				" Select an allele below to compare its sequences."
			))
	}
})

#  Sélecteur allèle (intersection littérature × séquences retirées)
observe({
	df_lit  <- partial_lit()
	df_part <- store$df_partiels
	req(!is.null(df_lit) && nrow(df_lit) > 0)
	
	alleles_lit <- sort(unique(df_lit$ALLELE))
	
	# Mettre en évidence les allèles qui sont aussi dans df_partiels
	if (!is.null(df_part) && nrow(df_part) > 0) {
		norm <- function(x) toupper(trimws(x))
		
		key_part <- unique(paste(norm(df_part$ALLELE), norm(df_part$SPECIES)))
		
		allele_species_lit <- df_lit %>%
			dplyr::distinct(ALLELE, SPECIES) %>%
			dplyr::mutate(key = paste(norm(ALLELE), norm(SPECIES)))
		
		matched_alleles <- unique(allele_species_lit$ALLELE[allele_species_lit$key %in% key_part])
		
		alleles_matched <- alleles_lit[alleles_lit %in% matched_alleles]
		alleles_other   <- alleles_lit[!(alleles_lit %in% matched_alleles)]
		choices <- c(
			if (length(alleles_matched) > 0)
				setNames(alleles_matched, paste0("★ ", alleles_matched, " [also in removed]"))
			else NULL,
			if (length(alleles_other) > 0)
				setNames(alleles_other, alleles_other)
			else NULL
		)
	} else {
		choices <- setNames(alleles_lit, alleles_lit)
	}
	
	updateSelectInput(session, "lit_allele_select",
		choices  = choices,
		selected = if (length(choices) > 0) choices[1] else NULL)
})

output$lit_removed_allele_info <- renderUI({
	req(input$lit_allele_select)
	df_part <- store$df_partiels
	if (is.null(df_part) || nrow(df_part) == 0) return(NULL)
	
	rows <- df_part %>% dplyr::filter(ALLELE == input$lit_allele_select)
	if (nrow(rows) == 0) return(NULL)   # cet allèle n'est pas dans les séquences retirées
	
	div(
		style = "background:#FFFBEB; border-left:4px solid #D97706;
                 border-radius:8px; padding:12px 16px; margin-top:10px; font-size:13px;",
		icon("star", style = "color:#D97706;"),
		strong(paste0(" ", input$lit_allele_select,
			" was removed from the final dataset (5'UTR < 40 nt).")),
		" Details of the removed sequence(s), by species:",
		tags$ul(style = "margin:8px 0 0 0;",
			lapply(seq_len(nrow(rows)), function(i) {
				r <- rows[i, ]
				tags$li(
					strong(r$SPECIES), " | Length: ", r$LENGTH_NT, " nt",
					" | Functionality: ", r$FUNCTIONALITY %||% "—"
				)
			})
		),
		tags$small(style = "color:#6B7280; display:block; margin-top:6px;",
			icon("info-circle"),
			" Full list: tab '🧩 Partial Genes → 🔍 Allele list'.")
	)
})
# ── Lancer la comparaison
observeEvent(input$btn_run_comparison, {
	req(partial_lit(), input$lit_allele_select)
	tryCatch({
		n_last <- as.integer(input$lit_n_last %||% 300)
		result <- prepare_comparison(partial_lit(), input$lit_allele_select, n_last)
		if (nrow(result$sequences) == 0) {
			showNotification("⚠️ No sequences found for this allele.",
				type = "warning", duration = 5)
			return()
		}
		partial_compare(result)
		showNotification(
			paste0("✅ Comparison done — ", nrow(result$sequences),
				" sequences × last ", n_last, " nt"),
			type = "message", duration = 4)
	}, error = function(e) {
		showNotification(paste0("❌ Comparison error: ", e$message),
			type = "error", duration = 8)
	})
})

# ── Résultats comparaison

output$compare_status <- renderUI({
	res <- partial_compare()
	if (is.null(res)) {
		div(class = "alert alert-info", icon("info-circle"),
			" Select an allele and click ",
			strong("🚀 Run comparison"), ".")
	} else {
		n   <- nrow(res$sequences)
		smr <- summary_comparison(res)
		div(class = "alert alert-success", icon("check-circle"),
			strong(paste0(" ", n, " sequences compared")),
			if (!is.null(smr))
				paste0(" | Mean identity: ", smr$Mean_identity, "% ",
					"(min ", smr$Min_identity, "% — max ", smr$Max_identity, "%)")
		)
	}
})

# Heatmap identité pairwise
output$compare_heatmap_ui <- renderUI({
	res <- partial_compare()
	if (is.null(res) || is.null(res$identity_matrix)) return(NULL)
	n   <- nrow(res$identity_matrix)
	h   <- min(max(420, n * 45 + 120), 900)
	plotlyOutput("compare_heatmap", height = paste0(h, "px"))
})

output$compare_heatmap <- renderPlotly({
	req(partial_compare(), input$lit_allele_select)
	plot_comparison_heatmap(partial_compare(), input$lit_allele_select,
		as.integer(input$lit_n_last %||% 300),
		interactive = TRUE)
})

# Tableau dalignement 
output$compare_alignment_table <- renderDT({
	req(partial_compare())
	df_aln <- prepare_alignment_table(partial_compare())
	if (nrow(df_aln) == 0) return(datatable(data.frame(Message="No data"),
		rownames = FALSE))
	# On masque la colonne Alignment trop longue 
	dt <- datatable(df_aln,
		options = list(
			pageLength = 20, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv","excel"),
			columnDefs = list(
				list(targets = which(names(df_aln) == "Alignment") - 1,
					visible = FALSE)  # masqué visuellement, présent en export
			)
		),
		extensions = "Buttons", rownames = FALSE,
		class = "cell-border stripe hover compact"
	)
	if ("Function" %in% names(df_aln)) {
		dt <- dt %>% formatStyle("Function",
			backgroundColor = styleEqual(c("F","P","ORF"),
				c("#D4EDDA","#F8D7DA","#FFF3CD")))
	}
	dt
})

# Table résumé statistique
output$compare_summary_table <- renderDT({
	req(partial_compare())
	smr <- summary_comparison(partial_compare())
	if (is.null(smr)) return(datatable(data.frame(Message="Need ≥ 2 sequences"),
		rownames = FALSE))
	datatable(smr,
		options = list(dom = "t", scrollX = TRUE),
		rownames = FALSE, class = "cell-border stripe hover compact")
})

# Distribution des identités pairwise vue par séquence 
output$compare_dist_plot_ui <- renderUI({
	req(partial_compare())
	mat <- partial_compare()$identity_matrix
	n   <- if (!is.null(mat)) nrow(mat) else 0
	h   <- max(380, n * 32 + 140)
	plotlyOutput("compare_dist_plot", height = paste0(min(h, 900), "px"))
})

output$compare_dist_plot <- renderPlotly({
	req(partial_compare())
	allele_lbl <- input$lit_allele_select %||% ""
	plot_identity_distribution(
		partial_compare(),
		allele_sel  = allele_lbl,
		n_last      = input$lit_n_last %||% 300,
		interactive = TRUE
	)
})

# Downloads sous-onglet B
output$btn_dl_compare_heatmap_png <- downloadHandler(
	filename = function() paste0("Comparison_Heatmap_",
		input$lit_allele_select, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(partial_compare(), input$lit_allele_select)
		p <- plot_comparison_heatmap(partial_compare(), input$lit_allele_select,
			as.integer(input$lit_n_last %||% 300),
			interactive = FALSE)
		n <- nrow(partial_compare()$identity_matrix)
		ggplot2::ggsave(file, plot = p,
			width  = max(8,  n * 0.6 + 3),
			height = max(6, n * 0.55 + 3),
			dpi = 300, bg = "white")
		showNotification("✅ PNG downloaded!", type = "message", duration = 3)
	}
)
output$btn_dl_compare_csv <- downloadHandler(
	filename = function() paste0("Comparison_Sequences_",
		input$lit_allele_select, "_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(partial_compare())
		write.csv(prepare_alignment_table(partial_compare()), file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	}
)
output$btn_dl_compare_dist_png <- downloadHandler(
	filename = function() paste0("Comparison_IdentityProfile_",
		input$lit_allele_select, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(partial_compare(), input$lit_allele_select)
		p <- plot_identity_distribution(
			partial_compare(),
			allele_sel  = input$lit_allele_select,
			n_last      = input$lit_n_last %||% 300,
			interactive = FALSE
		)
		n <- nrow(partial_compare()$identity_matrix)
		ggplot2::ggsave(file, plot = p, width = 10,
			height = max(6, n * 0.35 + 2), dpi = 300, bg = "white")
		showNotification("✅ PNG downloaded!", type = "message", duration = 3)
	}
)
output$btn_dl_lit_parsed_csv <- downloadHandler(
	filename = function() paste0("Literature_Sequences_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(partial_lit())
		df_export <- partial_lit() %>% dplyr::select(-SEQUENCE)
		write.csv(df_export, file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	}
)
