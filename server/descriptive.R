# ══ REACTIVE AUTO-TRIGGERS ═══════════════════════════════════
# Sequence distribution auto update when type or grouping changes
observe({
	req(store$df_final)
	input$seq_type; input$seq_grouping
	input$seq_species_filter; input$seq_locus_filter
	input$seq_locus_species_filter; input$seq_sg_mode; input$seq_topn
	isolate({
		tryCatch({
			df        <- store$df_final
			type_data <- isolate(input$seq_type)
			grouping  <- isolate(input$seq_grouping)
			if (type_data == "IG") df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
			else                   df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
			if (nrow(df) == 0) return()
			if (grouping == "species") {
				result_data <- df %>% group_by(SPECIES) %>% summarise(Count=n(), .groups="drop")
				result_plot <- plot_seq_by_species(df, type_data, interactive=TRUE)
			} else if (grouping == "locus") {
				sp_filter <- isolate(input$seq_locus_species_filter)
				df_l <- if (!is.null(sp_filter) && sp_filter != "") df %>% dplyr::filter(SPECIES==sp_filter) else df
				result_data <- df_l %>% group_by(LOCI,SPECIES) %>% summarise(Count=n(),.groups="drop")
				result_plot <- plot_seq_by_locus(df_l, type_data, interactive=TRUE)
			} else {
				result_data <- df %>% group_by(SUBGROUP,SPECIES) %>% summarise(Count=n(),.groups="drop")
				result_plot <- NULL
			}
			seq_results$data     <- result_data
			seq_results$plot     <- result_plot
			seq_results$grouping <- grouping
			seq_results$df_snap  <- df
			seq_results$done     <- TRUE
		}, error = function(e) NULL)
	})
})

# Functionality distribution auto-update 
observe({
	req(store$df_final)
	input$func_type; input$func_grouping
	input$func_species_filter; input$func_locus_filter
	input$func_locus_species_filter; input$func_sg_mode; input$func_topn
	isolate({
		tryCatch({
			df        <- store$df_final
			type_data <- isolate(input$func_type)
			grouping  <- isolate(input$func_grouping)
			if (type_data == "IG") df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
			else                   df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
			if (nrow(df) == 0 || all(is.na(df$FUNCTIONALITY))) return()
			if (grouping == "species") {
				result_data <- df %>% group_by(SPECIES,FUNCTIONALITY) %>% summarise(n=n(),.groups="drop") %>%
					group_by(SPECIES) %>% mutate(Proportion=n/sum(n)) %>% ungroup()
				result_plot <- plot_func_by_species(df, type_data, interactive=TRUE)
			} else if (grouping == "locus") {
				sp_filter <- isolate(input$func_locus_species_filter)
				df_l <- if (!is.null(sp_filter) && sp_filter != "") df %>% dplyr::filter(SPECIES==sp_filter) else df
				result_data <- df_l %>% group_by(LOCI,SPECIES,FUNCTIONALITY) %>% summarise(n=n(),.groups="drop") %>%
					group_by(LOCI,SPECIES) %>% mutate(Proportion=n/sum(n)) %>% ungroup()
				result_plot <- plot_func_by_locus(df_l, type_data, interactive=TRUE)
			} else {
				result_data <- df %>% group_by(SUBGROUP,SPECIES,FUNCTIONALITY) %>% summarise(n=n(),.groups="drop") %>%
					group_by(SUBGROUP,SPECIES) %>% mutate(Proportion=n/sum(n)) %>% ungroup()
				result_plot <- NULL
			}
			func_results$data     <- result_data
			func_results$plot     <- result_plot
			func_results$grouping <- grouping
			func_results$df_snap  <- df
			func_results$done     <- TRUE
		}, error = function(e) NULL)
	})
})

# Descriptive table — auto-update
observe({
	req(store$df_final)
	input$tableau_type; input$tableau_species
	isolate({
		tryCatch({
			species_sel <- isolate(input$tableau_species)
			type_sel    <- isolate(input$tableau_type)
			req(!is.null(species_sel), nchar(species_sel) > 0)
			result_data <- generate_descriptive_table_species(store$df_final, type_sel, species_sel)
			tableau_results$data <- result_data
			tableau_results$done <- TRUE
		}, error = function(e) NULL)
	})
})

# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  STATISTIQUES DESCRIPTIVES  ██████████████████ 
# ═══════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════
# A. TABLEAU DESCRIPTIF
# ══════════════════════════════════════════════════════════════

tableau_results <- reactiveValues(data = NULL, done = FALSE)

observe({
	req(store$df_final)
	isolate({
		df_filtered <- store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
		species_list <- unique(df_filtered$SPECIES)
		updateSelectInput(session, "tableau_species",
			choices = species_list, selected = species_list[1])
	})
})

observeEvent(input$btn_launch_tableau, {
	req(store$df_final)
	isolate({
		species_sel <- input$tableau_species
		type_sel    <- input$tableau_type
		req(species_sel, nchar(species_sel) > 0)
		withProgress(message = "Generating table…", value = 0.5, {
			tryCatch({
				result_data <- generate_descriptive_table_species(
					store$df_final, type_sel, species_sel
				)
				tableau_results$data <- result_data
				tableau_results$done <- TRUE
				incProgress(1)
			}, error = function(e) {
				showNotification(paste0("❌ Tableau error: ", e$message), type="error", duration=8)
			})
		})
	})
})

output$tableau_status <- renderUI({
	if (!tableau_results$done) {
		div(class="alert alert-info", icon("info-circle"),
			" Select a species above and click ", strong("🚀 Generate table"), ".")
	} else {
		div(class="alert alert-success", icon("check-circle"),
			" Table generated for ", strong(input$tableau_species))
	}
})

output$tableau_table <- renderDT({
	req(tableau_results$data)
	df <- tableau_results$data
	grp       <- get_tableau_col_groups(df)
	taux_cols <- grp$taux
	dist_cols <- grp$dist
	
	df_display <- df
	for (col in c(taux_cols, dist_cols)) {
		if (col %in% names(df_display)) {
			df_display[[col]] <- ifelse(is.na(df_display[[col]]), "—", as.character(df_display[[col]]))
		}
	}
	for (col in dist_cols) {
		if (col %in% names(df_display) && is.numeric(df[[col]])) {
			df_display[[col]] <- ifelse(is.na(df[[col]]), "—",
				as.character(as.integer(round(df[[col]], 0))))
		}
	}
	
	n_id    <- length(grp$id);   n_eff  <- length(grp$eff)
	n_taux  <- length(taux_cols); n_dist <- length(dist_cols); n_taille <- length(grp$taille)
	
	th_group <- function(label, colspan, bg, border_right = TRUE) {
		br_style <- if (border_right) "border-right:2px solid #ffd700;" else ""
		sprintf('<th colspan="%d" style="background:%s;color:white;font-weight:bold;text-align:center;vertical-align:middle;%s">%s</th>',
			colspan, bg, br_style, label)
	}
	th_col <- function(name) sprintf('<th style="text-align:center;vertical-align:bottom;">%s</th>', name)
	
	# En-tête ligne 1 : groupes de colonnes (motifs masqués pour TR)
	row1_cells <- c(th_group("Identification", n_id, "#001f54"),
		th_group("Counts", n_eff, "#1565C0"))
	if (n_taux   > 0) row1_cells <- c(row1_cells, th_group("Motif presence rate (%)",       n_taux,  "#166534"))
	if (n_dist   > 0) row1_cells <- c(row1_cells, th_group("Mean distance to LPART1 (nt)", n_dist,  "#854D0E"))
	if (n_taille > 0) row1_cells <- c(row1_cells, th_group("5'UTR length (nt)",            n_taille, "#581C87", border_right = FALSE))
	
	# En-tête ligne 2 : labels courts
	all_cols   <- c(grp$id, grp$eff, taux_cols, dist_cols, grp$taille)
	col_labels <- all_cols
	col_labels <- sub("^SPECIES$",        "Species",       col_labels)
	col_labels <- sub("^LOCI$",           "Locus",         col_labels)
	col_labels <- sub("^SUBGROUP$",       "Subgroup",      col_labels)
	col_labels <- sub("^Taux_",           "",              col_labels)
	col_labels <- sub("^Dist_",           "",              col_labels)
	col_labels <- sub("Taille_5UTR_moy",  "Mean",          col_labels)
	col_labels <- sub("Taille_5UTR_med",  "Median",        col_labels)
	col_labels <- sub("Nb_alleles",       "Alleles",       col_labels)
	col_labels <- sub("Fonctionnels",     "F",             col_labels)
	col_labels <- sub("Pseudogenes",      "P",             col_labels)
	row2_cells <- sapply(col_labels, th_col)
	
	container_html <- htmltools::withTags(table(class = "display",
		thead(
			tr(HTML(paste(row1_cells, collapse = ""))),
			tr(HTML(paste(row2_cells, collapse = "")))
		)))
	
	dt <- datatable(df_display, container = container_html,
		options = list(pageLength = 50, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","excel","csv"),
			columnDefs = list(list(className = "dt-center", targets = "_all"))),
		extensions = "Buttons", rownames = FALSE, filter = "top",
		class = "cell-border stripe hover compact")
	
	# Style colonne LOCI
	if ("LOCI" %in% names(df_display)) {
		dt <- dt %>%
			formatStyle("LOCI",
				backgroundColor = styleEqual(
					c("IGHV","IGKV","IGLV","TRAV","TRBV","TRDV","TRGV"),
					c("#DBEAFE","#DCFCE7","#FEF9C3","#FCE7F3","#FFE4E6","#E0F2FE","#F3E8FF")),
				fontWeight = "bold",
				color = styleEqual(
					c("IGHV","IGKV","IGLV","TRAV","TRBV","TRDV","TRGV"),
					c("#1E3A8A","#166534","#854D0E","#9D174D","#991B1B","#0C4A6E","#581C87")))
	}
	if ("SUBGROUP" %in% names(df_display))
		dt <- dt %>% formatStyle("SUBGROUP", fontWeight = "bold", textAlign = "left")
	if ("SPECIES" %in% names(df_display))
		dt <- dt %>% formatStyle("SPECIES",  fontStyle  = "italic", color = "#374151")
	
	# Barres de couleur pour les effectifs
	for (col_name in c("Nb_alleles")) {
		if (!col_name %in% names(df)) next
		rng <- range(df[[col_name]], na.rm = TRUE)
		if (rng[1] < rng[2])
			dt <- dt %>% formatStyle(col_name,
				background = styleColorBar(rng, "#93C5FD"),
				backgroundSize = "100% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")
	}
	for (col_name in c("Fonctionnels","Pseudogenes","ORF")) {
		if (!col_name %in% names(df)) next
		mx <- max(df[[col_name]], na.rm = TRUE)
		if (is.na(mx) || mx == 0) next
		col_bar <- switch(col_name, Fonctionnels = "#16A34A", Pseudogenes = "#EF4444", "#F59E0B")
		dt <- dt %>% formatStyle(col_name,
			background = styleColorBar(c(0, mx), col_bar),
			backgroundSize = "100% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "center")
	}
	# Style colonnes Taux 
	for (col in taux_cols) {
		if (!col %in% names(df_display)) next
		dt <- dt %>% formatStyle(col, fontWeight = "bold",
			color           = styleEqual("—", "grey", default = "#166534"),
			backgroundColor = styleEqual("—", "#F9FAFB", default = "#DCFCE7"))
	}
	# Style colonnes Dist 
	for (col in dist_cols) {
		if (!col %in% names(df_display)) next
		dt <- dt %>% formatStyle(col,
			backgroundColor = styleEqual("—", "#F9FAFB", default = "#FEF9C3"),
			color           = styleEqual("—", "grey", default = "#854D0E"))
	}
	for (col in grp$taille) {
		if (!col %in% names(df_display)) next
		dt <- dt %>% formatStyle(col, backgroundColor = "#F3E8FF", color = "#581C87")
	}
	return(dt)
})

output$btn_download_tableau_csv <- downloadHandler(
	filename = function() paste0("Tableau_", input$tableau_species, "_", Sys.Date(), ".csv"),
	content  = function(file) { req(tableau_results$data); write.csv(tableau_results$data, file, row.names=FALSE) }
)

output$btn_download_tableau_excel <- downloadHandler(
	filename = function() paste0("Tableau_", input$tableau_species, "_", Sys.Date(), ".xlsx"),
	content  = function(file) {
		req(tableau_results$data)
		tryCatch({
			library(openxlsx)
			wb  <- createWorkbook()
			sn  <- substr(input$tableau_species, 1, 31)
			addWorksheet(wb, sn, gridLines=FALSE)
			df  <- tableau_results$data
			grp <- get_tableau_col_groups(df)
			writeData(wb, sn, df, startRow=2, rowNames=FALSE)
			
			group_defs <- list(
				list(label="Identification", cols=which(names(df) %in% grp$id), fill="#001f54", font="#FFFFFF"),
				list(label="Counts", cols=which(names(df) %in% grp$eff), fill="#1565C0", font="#FFFFFF"),
				list(label="Motif occurrence rate (%)", cols=which(names(df) %in% grp$taux), fill="#166534", font="#FFFFFF"),
				list(label="Mean distance to LPART1 (nt)", cols=which(names(df) %in% grp$dist), fill="#854D0E", font="#FFFFFF"),
				list(label="5′ UTR length", cols=which(names(df) %in% grp$taille), fill="#581C87", font="#FFFFFF")
			)
			for (g in group_defs) {
				if (length(g$cols) == 0) next
				style_g <- createStyle(fontName="Calibri", fontSize=11, fontColour=g$font, fgFill=g$fill,
					halign="CENTER", valign="CENTER", textDecoration="bold",
					border="LeftRight", borderColour="#ffd700", borderStyle="medium")
				writeData(wb, sn, g$label, startRow=1, startCol=min(g$cols), rowNames=FALSE, colNames=FALSE)
				addStyle(wb, sn, style_g, rows=1, cols=g$cols, gridExpand=TRUE)
				if (length(g$cols) > 1) mergeCells(wb, sn, rows=1, cols=g$cols)
			}
			
			style_header <- createStyle(fontName="Calibri", fontSize=10, fontColour="#FFFFFF", fgFill="#1F4E79",
				halign="CENTER", valign="CENTER", textDecoration="bold", wrapText=TRUE,
				border="Bottom", borderColour="#ffd700", borderStyle="thin")
			addStyle(wb, sn, style_header, rows=2, cols=1:ncol(df), gridExpand=TRUE)
			
			style_id   <- createStyle(fontName="Calibri", halign="LEFT", fontColour="#1E3A8A", textDecoration="bold")
			style_sp   <- createStyle(fontName="Calibri", halign="LEFT", fontColour="#374151", textDecoration="italic")
			style_eff  <- createStyle(fontName="Calibri", halign="CENTER", numFmt="0")
			style_F    <- createStyle(fgFill="#D4EDDA", fontColour="#155724", fontName="Calibri", halign="CENTER")
			style_P    <- createStyle(fgFill="#F8D7DA", fontColour="#721C24", fontName="Calibri", halign="CENTER")
			style_O    <- createStyle(fgFill="#FFF3CD", fontColour="#856404", fontName="Calibri", halign="CENTER")
			style_taux <- createStyle(fgFill="#DCFCE7", fontColour="#166534", fontName="Calibri", halign="CENTER", numFmt="0.0")
			style_dist <- createStyle(fgFill="#FEF9C3", fontColour="#854D0E", fontName="Calibri", halign="CENTER", numFmt="0")
			style_taille <- createStyle(fgFill="#F3E8FF", fontColour="#581C87", fontName="Calibri", halign="CENTER", numFmt="0")
			style_alt  <- createStyle(fgFill="#F8FAFF", fontName="Calibri")
			rows_data  <- 3:(nrow(df)+2)
			
			for (r in seq(3, nrow(df)+2, by=2)) addStyle(wb, sn, style_alt, rows=r, cols=1:ncol(df), gridExpand=TRUE)
			sp_idx <- which(names(df) == "SPECIES")
			if (length(sp_idx) > 0) addStyle(wb, sn, style_sp, rows=rows_data, cols=sp_idx, gridExpand=TRUE)
			id_no_sp <- grp$id[grp$id != "SPECIES"]
			if (length(id_no_sp) > 0) addStyle(wb, sn, style_id, rows=rows_data, cols=which(names(df) %in% id_no_sp), gridExpand=TRUE)
			if (length(grp$eff) > 0) addStyle(wb, sn, style_eff, rows=rows_data, cols=which(names(df) %in% grp$eff), gridExpand=TRUE)
			for (col_name in c("Fonctionnels","Pseudogenes","ORF")) {
				if (col_name %in% names(df)) {
					sty <- switch(col_name, Fonctionnels=style_F, Pseudogenes=style_P, ORF=style_O)
					addStyle(wb, sn, sty, rows=rows_data, cols=which(names(df)==col_name), gridExpand=TRUE)
				}
			}
			taux_idx <- which(names(df) %in% grp$taux)
			if (length(taux_idx) > 0) addStyle(wb, sn, style_taux, rows=rows_data, cols=taux_idx, gridExpand=TRUE)
			dist_idx <- which(names(df) %in% grp$dist)
			if (length(dist_idx) > 0) addStyle(wb, sn, style_dist, rows=rows_data, cols=dist_idx, gridExpand=TRUE)
			taille_idx <- which(names(df) %in% grp$taille)
			if (length(taille_idx) > 0) addStyle(wb, sn, style_taille, rows=rows_data, cols=taille_idx, gridExpand=TRUE)
			
			setColWidths(wb, sn, cols=1:ncol(df), widths="auto")
			freezePane(wb, sn, firstActiveRow=3, firstActiveCol=3)
			saveWorkbook(wb, file, overwrite=TRUE)
			showNotification("✅ Excel downloaded!", type="message", duration=4)
		}, error = function(e) {
			showNotification(paste0("❌ Excel error: ", e$message), type="error", duration=8)
		})
	}
)

# ══════════════════════════════════════════════════════════════
# B. DISTRIBUTION DES SÉQUENCES — RÉACTIVE AUTOMATIQUE
# ══════════════════════════════════════════════════════════════

seq_results <- reactiveValues(data=NULL, plot=NULL, done=FALSE, grouping=NULL, df_snap=NULL)

observe({
	req(store$df_final)
	input$seq_type  
	isolate({
		df <- store$df_final
		if (input$seq_type == "IG") df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
		else                        df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
		sp_list <- sort(unique(df$SPECIES))
		updateSelectInput(session, "seq_species_filter", choices=setNames(sp_list,sp_list), selected=sp_list[1])
		updateSelectInput(session, "seq_locus_filter", choices=c("All loci", sort(unique(df$LOCI))), selected="All loci")
		updateSelectInput(session, "seq_locus_species_filter",
			choices  = c("All species" = "", setNames(sp_list, sp_list)),
			selected = "")
	})
})

observeEvent(input$btn_launch_seq, {
	req(store$df_final)
	isolate({
		tryCatch({
			df        <- store$df_final
			type_data <- input$seq_type
			grouping  <- input$seq_grouping
			
			if (type_data == "IG") df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
			else                   df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
			if (nrow(df) == 0) return()
			
			if (grouping == "species") {
				result_data <- df %>% group_by(SPECIES) %>% summarise(Count=n(), .groups="drop")
				result_plot <- plot_seq_by_species(df, type_data, interactive=TRUE)
				
			} else if (grouping == "locus") {
				sp_locus_filter <- isolate(input$seq_locus_species_filter)
				df_for_locus <- if (!is.null(sp_locus_filter) && sp_locus_filter != "")
					df %>% dplyr::filter(SPECIES == sp_locus_filter) else df
				result_data <- df_for_locus %>% group_by(LOCI, SPECIES) %>% summarise(Count=n(), .groups="drop")
				result_plot <- plot_seq_by_locus(df_for_locus, type_data, interactive=TRUE)
				
			} else {
				result_data <- df %>% group_by(SUBGROUP, SPECIES) %>% summarise(Count=n(), .groups="drop")
				result_plot <- NULL
			}
			
			seq_results$data     <- result_data
			seq_results$plot     <- result_plot
			seq_results$grouping <- grouping
			seq_results$df_snap  <- df
			seq_results$done     <- TRUE
			
		}, error = function(e) {
			showNotification(paste0("❌ Sequences: ", e$message), type="error", duration=8)
		})
	})
})

output$seq_sg_alert <- renderUI({
	req(seq_results$done, seq_results$grouping == "subgroup")
	df   <- seq_results$df_snap
	sp   <- input$seq_species_filter
	if (!is.null(sp) && sp != "") df <- df %>% dplyr::filter(SPECIES == sp)
	n_sg <- length(unique(df$SUBGROUP))
	top_n_val <- if (!is.na(input$seq_topn) && input$seq_topn > 0) input$seq_topn else NULL
	sg_alert_banner(n_sg, top_n_val)
})

output$seq_plot_ui <- renderUI({
	req(seq_results$done)
	grp <- seq_results$grouping
	if (grp == "subgroup") {
		df   <- seq_results$df_snap
		sp   <- input$seq_species_filter
		if (!is.null(sp) && sp != "") df <- df %>% dplyr::filter(SPECIES == sp)
		n_sg <- length(unique(df$SUBGROUP))
		h    <- sg_plot_height_px(n_sg, min_h=400)
		plotly::plotlyOutput("seq_plot", height=paste0(h,"px"))
	} else if (grp == "locus") {
		sp_filter <- input$seq_locus_species_filter
		df_count  <- seq_results$df_snap
		if (!is.null(sp_filter) && sp_filter != "")
			df_count <- df_count %>% dplyr::filter(SPECIES == sp_filter)
		n_sp_v <- length(unique(df_count$SPECIES))
		n_col  <- if (n_sp_v <= 3) 1 else 2
		n_rows <- ceiling(n_sp_v / n_col)
		h      <- if (n_sp_v == 1) 420 else max(420, n_rows * 340)
		plotly::plotlyOutput("seq_plot", height=paste0(h,"px"))
	} else {
		
		n_sp_v <- length(unique(seq_results$df_snap$SPECIES))
		h      <- min(max(400, n_sp_v * 55), 700)
		plotly::plotlyOutput("seq_plot", height=paste0(h,"px"))
	}
})

output$seq_plot <- renderPlotly({
	req(seq_results$df_snap, seq_results$grouping)
	grp <- seq_results$grouping
	df  <- seq_results$df_snap
	
	if (grp == "subgroup") {
		sp_sel    <- input$seq_species_filter; req(sp_sel)
		locus_sel <- input$seq_locus_filter
		top_n_val <- if (!is.na(input$seq_topn) && input$seq_topn > 0) input$seq_topn else NULL
		p <- plot_seq_by_subgroup_smart(df, input$seq_type,
			species_filter = sp_sel,
			locus_filter   = if (!is.null(locus_sel) && locus_sel != "All loci") locus_sel else NULL,
			interactive    = TRUE, top_n=top_n_val, force_mode=input$seq_sg_mode %||% "auto")
	} else {
		p <- seq_results$plot
	}
	if (inherits(p,"plotly")) return(p)
	plotly::ggplotly(p)
})

output$btn_download_seq_plot_png <- downloadHandler(
	filename = function() paste0("DistSeq_", seq_results$grouping, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(seq_results$df_snap, seq_results$grouping)
		grp <- seq_results$grouping; df <- seq_results$df_snap
		top_n_val <- if (!is.na(input$seq_topn) && input$seq_topn > 0) input$seq_topn else NULL
		if (grp == "subgroup") {
			sp <- input$seq_species_filter
			if (!is.null(sp) && sp != "") df <- df %>% dplyr::filter(SPECIES==sp)
			p <- plot_seq_by_subgroup_smart(df, input$seq_type, species_filter=sp, interactive=FALSE, top_n=top_n_val, force_mode=input$seq_sg_mode %||% "auto")
			save_sg_plot_png(p, file, n_sg=length(unique(df$SUBGROUP)), grouping="subgroup")
		} else if (grp == "locus") {
			sp_filter_dl <- input$seq_locus_species_filter
			df_dl <- if (!is.null(sp_filter_dl) && sp_filter_dl != "")
				df %>% dplyr::filter(SPECIES == sp_filter_dl) else df
			p <- plot_seq_by_locus(df_dl, input$seq_type, interactive=FALSE)
			save_sg_plot_png(p, file, grouping="locus", n_sp=length(unique(df_dl$SPECIES)))
		} else {
			p <- plot_seq_by_species(df, input$seq_type, interactive=FALSE)
			save_sg_plot_png(p, file, grouping="species")
		}
		showNotification("✅ PNG downloaded!", type="message", duration=3)
	}
)

output$btn_download_seq_plot_pdf <- downloadHandler(
	filename = function() paste0("DistSeq_", seq_results$grouping, "_", Sys.Date(), ".pdf"),
	content  = function(file) {
		req(seq_results$df_snap, seq_results$grouping)
		grp <- seq_results$grouping; df <- seq_results$df_snap
		top_n_val <- if (!is.na(input$seq_topn) && input$seq_topn > 0) input$seq_topn else NULL
		if (grp == "subgroup") {
			sp <- input$seq_species_filter
			if (!is.null(sp) && sp != "") df <- df %>% dplyr::filter(SPECIES==sp)
			p <- plot_seq_by_subgroup_smart(df, input$seq_type, species_filter=sp, interactive=FALSE, top_n=top_n_val, force_mode=input$seq_sg_mode %||% "auto")
			save_sg_plot_pdf(p, file, n_sg=length(unique(df$SUBGROUP)), grouping="subgroup")
		} else if (grp == "locus") {
			sp_filter_dl <- input$seq_locus_species_filter
			df_dl <- if (!is.null(sp_filter_dl) && sp_filter_dl != "")
				df %>% dplyr::filter(SPECIES == sp_filter_dl) else df
			p <- plot_seq_by_locus(df_dl, input$seq_type, interactive=FALSE)
			save_sg_plot_pdf(p, file, grouping="locus", n_sp=length(unique(df_dl$SPECIES)))
		} else {
			p <- plot_seq_by_species(df, input$seq_type, interactive=FALSE)
			save_sg_plot_pdf(p, file, grouping="species")
		}
		showNotification("✅ PDF downloaded!", type="message", duration=3)
	}
)

output$btn_download_seq_csv <- downloadHandler(
	filename = function() paste0("DistSeq_", input$seq_grouping, "_", Sys.Date(), ".csv"),
	content  = function(file) { req(seq_results$data); write.csv(seq_results$data, file, row.names=FALSE) }
)

observeEvent(input$seq_locus_species_filter, {
	req(seq_results$done, seq_results$grouping == "locus", store$df_final)
	isolate({
		tryCatch({
			df <- seq_results$df_snap
			sp_filter <- input$seq_locus_species_filter
			df_for_locus <- if (!is.null(sp_filter) && sp_filter != "")
				df %>% dplyr::filter(SPECIES == sp_filter) else df
			result_plot <- plot_seq_by_locus(df_for_locus, input$seq_type, interactive=TRUE)
			seq_results$plot <- result_plot
		}, error = function(e) {
			showNotification(paste0("❌ Locus species filter: ", e$message), type="error", duration=6)
		})
	})
})

# ══════════════════════════════════════════════════════════════
# C. DISTRIBUTION FONCTIONNALITÉ — RÉACTIVE AUTOMATIQUE
# ══════════════════════════════════════════════════════════════

func_results <- reactiveValues(data=NULL, plot=NULL, done=FALSE, grouping=NULL, df_snap=NULL)

observe({
	req(store$df_final)
	input$func_type  
	isolate({
		df <- store$df_final
		if (input$func_type == "IG") df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
		else                         df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
		sp_list <- sort(unique(df$SPECIES))
		updateSelectInput(session, "func_species_filter", choices=setNames(sp_list,sp_list), selected=sp_list[1])
		updateSelectInput(session, "func_locus_filter", choices=c("All loci",sort(unique(df$LOCI))), selected="All loci")
		updateSelectInput(session, "func_locus_species_filter",
			choices  = c("All species" = "", setNames(sp_list, sp_list)),
			selected = "")
	})
})

observeEvent(input$btn_launch_func, {
	req(store$df_final)
	isolate({
		tryCatch({
			df        <- store$df_final
			type_data <- input$func_type
			grouping  <- input$func_grouping
			
			if (type_data == "IG") df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
			else                   df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
			if (nrow(df) == 0 || all(is.na(df$FUNCTIONALITY))) return()
			
			if (grouping == "species") {
				result_data <- df %>% group_by(SPECIES,FUNCTIONALITY) %>% summarise(n=n(),.groups="drop") %>%
					group_by(SPECIES) %>% mutate(Proportion=n/sum(n)) %>% ungroup()
				result_plot <- plot_func_by_species(df, type_data, interactive=TRUE)
				
			} else if (grouping == "locus") {
				sp_locus_filter <- isolate(input$func_locus_species_filter)
				df_for_locus <- if (!is.null(sp_locus_filter) && sp_locus_filter != "")
					df %>% dplyr::filter(SPECIES == sp_locus_filter) else df
				result_data <- df_for_locus %>% group_by(LOCI,SPECIES,FUNCTIONALITY) %>% summarise(n=n(),.groups="drop") %>%
					group_by(LOCI,SPECIES) %>% mutate(Proportion=n/sum(n)) %>% ungroup()
				result_plot <- plot_func_by_locus(df_for_locus, type_data, interactive=TRUE)
				
			} else {
				if (length(unique(na.omit(df$SUBGROUP))) == 0) return()
				result_data <- df %>% group_by(SUBGROUP,SPECIES,FUNCTIONALITY) %>% summarise(n=n(),.groups="drop") %>%
					group_by(SUBGROUP,SPECIES) %>% mutate(Proportion=n/sum(n)) %>% ungroup()
				result_plot <- NULL
			}
			
			func_results$data     <- result_data
			func_results$plot     <- result_plot
			func_results$grouping <- grouping
			func_results$df_snap  <- df
			func_results$done     <- TRUE
			
		}, error = function(e) {
			showNotification(paste0("❌ Functionality: ", e$message), type="error", duration=8)
		})
	})
})

output$func_sg_alert <- renderUI({
	req(func_results$done, func_results$grouping=="subgroup")
	df <- func_results$df_snap; sp <- input$func_species_filter
	if (!is.null(sp) && sp != "") df <- df %>% dplyr::filter(SPECIES==sp)
	n_sg <- length(unique(df$SUBGROUP))
	top_n_val <- if (!is.na(input$func_topn) && input$func_topn > 0) input$func_topn else NULL
	sg_alert_banner(n_sg, top_n_val)
})

output$func_plot_ui <- renderUI({
	req(func_results$done)
	grp <- func_results$grouping
	if (grp == "subgroup") {
		df <- func_results$df_snap; sp <- input$func_species_filter
		if (!is.null(sp) && sp != "") df <- df %>% dplyr::filter(SPECIES==sp)
		n_sg <- length(unique(df$SUBGROUP))
		
		plotly::plotlyOutput("func_plot", height=paste0(sg_plot_height_px(n_sg, min_h=380) + 60,"px"))
	} else if (grp == "locus") {
		sp_filter <- input$func_locus_species_filter
		df_count  <- func_results$df_snap
		if (!is.null(sp_filter) && sp_filter != "")
			df_count <- df_count %>% dplyr::filter(SPECIES == sp_filter)
		n_sp_v <- length(unique(df_count$SPECIES))
		n_col  <- if (n_sp_v <= 3) 1 else 2
		n_rows <- ceiling(n_sp_v / n_col)
		h      <- if (n_sp_v == 1) 420 else max(420, n_rows * 340)
		plotly::plotlyOutput("func_plot", height=paste0(h + 60,"px"))
	} else {
		
		n_sp_v <- length(unique(func_results$df_snap$SPECIES))
		h      <- min(max(450, n_sp_v * 55), 700)
		plotly::plotlyOutput("func_plot", height=paste0(h,"px"))
	}
})

output$func_plot <- renderPlotly({
	req(func_results$df_snap, func_results$grouping)
	grp <- func_results$grouping; df <- func_results$df_snap
	if (grp == "subgroup") {
		sp_sel    <- input$func_species_filter; req(sp_sel)
		locus_sel <- input$func_locus_filter
		top_n_val <- if (!is.na(input$func_topn) && input$func_topn > 0) input$func_topn else NULL
		p <- plot_func_by_subgroup_smart(df, input$func_type,
			species_filter = sp_sel,
			locus_filter   = if (!is.null(locus_sel) && locus_sel != "All loci") locus_sel else NULL,
			interactive=TRUE, top_n=top_n_val, force_mode=input$func_sg_mode %||% "auto")
	} else { p <- func_results$plot }
	if (inherits(p,"plotly")) return(p)
	plotly::ggplotly(p)
})

output$btn_download_func_plot <- downloadHandler(
	filename = function() paste0("DistFunc_", func_results$grouping, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(func_results$df_snap, func_results$grouping)
		grp <- func_results$grouping; df <- func_results$df_snap
		top_n_val <- if (!is.na(input$func_topn) && input$func_topn > 0) input$func_topn else NULL
		if (grp == "subgroup") {
			sp <- input$func_species_filter
			if (!is.null(sp) && sp != "") df <- df %>% dplyr::filter(SPECIES==sp)
			p <- plot_func_by_subgroup_smart(df, input$func_type, species_filter=sp, interactive=FALSE, top_n=top_n_val, force_mode=input$func_sg_mode %||% "auto")
			save_sg_plot_png(p, file, n_sg=length(unique(df$SUBGROUP)), grouping="subgroup")
		} else if (grp == "locus") {
			sp_filter_dl <- input$func_locus_species_filter
			df_dl <- if (!is.null(sp_filter_dl) && sp_filter_dl != "")
				df %>% dplyr::filter(SPECIES == sp_filter_dl) else df
			p <- plot_func_by_locus(df_dl, input$func_type, interactive=FALSE)
			save_sg_plot_png(p, file, grouping="locus", n_sp=length(unique(df_dl$SPECIES)))
		} else {
			p <- plot_func_by_species(df, input$func_type, interactive=FALSE)
			save_sg_plot_png(p, file, grouping="species")
		}
		showNotification("✅ PNG downloaded!", type="message", duration=3)
	}
)



output$btn_download_func_csv <- downloadHandler(
	filename = function() paste0("DistFunc_", input$func_grouping, "_", Sys.Date(), ".csv"),
	content  = function(file) { req(func_results$data); write.csv(func_results$data, file, row.names=FALSE) }
)

observeEvent(input$func_locus_species_filter, {
	req(func_results$done, func_results$grouping == "locus", store$df_final)
	isolate({
		tryCatch({
			df <- func_results$df_snap
			sp_filter <- input$func_locus_species_filter
			df_for_locus <- if (!is.null(sp_filter) && sp_filter != "")
				df %>% dplyr::filter(SPECIES == sp_filter) else df
			result_plot <- plot_func_by_locus(df_for_locus, input$func_type, interactive=TRUE)
			func_results$plot <- result_plot
		}, error = function(e) {
			showNotification(paste0("❌ Locus species filter: ", e$message), type="error", duration=6)
		})
	})
})
