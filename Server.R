# ═══════════════════════════════════════════════════════════════
# IMGT/PROMISE - Serveur
# ═══════════════════════════════════════════════════════════════
filter <- dplyr::filter

server <- function(input, output, session) {
	
	# ══ STOCKAGE GLOBAL ══════════════════════════════════════════
	store <- reactiveValues(
		df_5utr  = NULL,
		df_final = NULL
	)
	
	# ══  IMPORT 5'UTR ════════════════════════════════════
	imported_utr <- mod_import_server("import_module")
	observeEvent(imported_utr(), {
		req(imported_utr())
		store$df_5utr <- imported_utr()
		store$df_5utr <- store$df_5utr %>%
			mutate(
				FUNCTIONALITY = str_remove_all(FUNCTIONALITY, "[\\(\\)\\[\\]]"),
				FUNCTIONALITY = str_trim(FUNCTIONALITY),
				FUNCTIONALITY = ifelse(FUNCTIONALITY == "", NA, FUNCTIONALITY)
			)
		showNotification(
			HTML(paste0(
				"<strong>✅ 5' UTR file successfully loaded!</strong><br>",
				"• ", nrow(imported_utr()), " headers extracted<br>",
				"• ", length(unique(imported_utr()$LOCI)), " LOCI detected<br>",
				"• ", length(unique(imported_utr()$SPECIES)), " species<br>",
				"<em style='color: #666;'>DNA sequences ignored.</em>"
			)),
			type = "message", duration = 8
		)
	})
	
	# ══  IMPORT MOTIFS ════════════════════════════════════
	motifs_result <- mod_motifs_server(
		"motifs_module",
		data_5utr = reactive(store$df_5utr)
	)
	
	# ══ AVERTISSEMENTS ══════════════════════════════════════════
	output$warn_no_utr <- renderUI({
		if (is.null(store$df_5utr))
			div(class = "alert alert-warning",
				icon("exclamation-triangle"), strong(" Warning : "),
				"Please import the 5′ UTR file first (onglet 📂 Import 5'UTR).")
	})
	
	# ══ CONSTRUCTION DU DATASET  ══════════════════════════
	output$dataset_build_status <- renderUI({
		n_utr       <- if (!is.null(store$df_5utr)) nrow(store$df_5utr) else 0
		n_motifs    <- length(motifs_result$motif_data())
		motif_names <- names(motifs_result$motif_data())
		sp_names    <- if (!is.null(store$df_5utr)) sort(unique(store$df_5utr$SPECIES)) else character(0)
		loci_names  <- if (!is.null(store$df_5utr)) sort(unique(store$df_5utr$LOCI))    else character(0)
		
		tagList(
			p(icon("file-code"), " 5'UTR : ",
				if (n_utr > 0) span(style="color:green;", strong(n_utr), " sequences imported ✅")
				else           span(style="color:orange;", "Not imported yet ⏳")),
			if (n_utr > 0) tagList(
				p(style = "font-size:12px; color:#374151; margin:2px 0 2px 22px;",
					strong("Species: "), paste(sp_names, collapse = ", ")),
				p(style = "font-size:12px; color:#374151; margin:2px 0 8px 22px;",
					strong("Loci: "), paste(loci_names, collapse = ", "))
			),
			p(icon("puzzle-piece"), " Motifs : ",
				if (n_motifs > 0) span(style="color:green;", strong(n_motifs), " motifs imported ✅")
				else              span(style="color:orange;", "No motif imported yet ⏳")),
			if (n_motifs > 0)
				p(style = "font-size:12px; color:#374151; margin:2px 0 0 22px;",
					strong("Motifs: "), paste(motif_names, collapse = ", "))
		)
	})
	
	observeEvent(input$btn_build_dataset, {
		if (is.null(store$df_5utr)) {
			showNotification("❌ Please import the 5'UTR file first.", type = "error", duration = 5)
			return()
		}
		withProgress(message = "Building dataset...", value = 0, {
			df <- store$df_5utr
			md <- motifs_result$motif_data()
			
			if ("SUBGROUP" %in% colnames(df)) {
				df$SUBGROUP <- clean_subgroup(df$SUBGROUP)
			}
			incProgress(0.2, detail = "5'UTR base loaded")
			for (motif_name in names(md)) {
				df_m <- md[[motif_name]]
				cols_join <- c("ALLELE","SPECIES", "GENE",
					paste0("SEQ_",       motif_name),
					paste0("POS_START_", motif_name),
					paste0("POS_END_",   motif_name))
				cols_join <- intersect(cols_join, colnames(df_m))
				df <- merge(df, df_m[, cols_join], by = c("ALLELE","SPECIES", "GENE"), all.x = TRUE)
			}
			incProgress(0.5, detail = "Join complete")
			lpart_col <- if ("POS_START_LPART_1" %in% colnames(df)) "POS_START_LPART_1"
			else if ("POS_START_LPART1" %in% colnames(df)) "POS_START_LPART1"
			else NULL
			if (!is.null(lpart_col)) {
				for (motif_name in names(md)) {
					if (motif_name %in% c("LPART_1", "LPART1")) next
					pos_end_col <- paste0("POS_END_",   motif_name)
					dist_col    <- paste0("DIST_",      motif_name)
					if (pos_end_col %in% colnames(df)) {
						df[[dist_col]] <- df[[lpart_col]] - df[[pos_end_col]]
						df[[dist_col]] <- ifelse(df[[dist_col]] < 0, NA_real_, df[[dist_col]])
					}
				}
			}
			incProgress(0.7, detail = "Distances calculated")
			for (motif_name in names(md)) {
				seq_col <- paste0("SEQ_", motif_name)
				bin_col <- paste0("BIN_", motif_name)
				if (seq_col %in% colnames(df)) {
					df[[bin_col]] <- ifelse(is.na(df[[seq_col]]) | df[[seq_col]] == "", 0L, 1L)
					df[[seq_col]] <- ifelse(is.na(df[[seq_col]]) | df[[seq_col]] == "", "absent", df[[seq_col]])
				}
			}
			incProgress(0.85, detail = "Analyzing motif presence")
			
			store$df_final <- df
			store$df_final   <- df
			
			store$df_no_motif   <- identify_no_motif_sequences(df)
			store$df_with_motif <- identify_with_motif_sequences(df)
			store$df_partiels   <- store$df_no_motif
			incProgress(1, detail = "Done!")
		})
		
		
		showNotification(
			HTML(paste0(
				"<strong>✅ Dataset built</strong><br>",
				"• ", nrow(store$df_final), " sequences retained<br>",
				"• ", ncol(store$df_final), " columns<br>",
				"• <em>", nrow(store$df_no_motif), " IG sequences with NO regulatory motif detected</em>"
			)),
			type = "message", duration = 10
		)
	})
	

	
	# ══ VALUE BOXES ══════════════════════════════════════════════
	output$vbox_total <- renderValueBox({
		n <- if (!is.null(store$df_final)) nrow(store$df_final) else 0
		valueBox(n, "Sequences", icon = icon("dna"), color = "blue")
	})
	output$vbox_loci <- renderValueBox({
		n <- if (!is.null(store$df_final)) length(unique(store$df_final$LOCI)) else 0
		valueBox(n, "LOCI", icon = icon("layer-group"), color = "purple")
	})
	output$vbox_species <- renderValueBox({
		n <- if (!is.null(store$df_final)) length(unique(store$df_final$SPECIES)) else 0
		valueBox(n, "Species", icon = icon("globe"), color = "green")
	})
	output$vbox_func <- renderValueBox({
		n <- if (!is.null(store$df_final)) sum(store$df_final$FUNCTIONALITY == "F", na.rm = TRUE) else 0
		valueBox(n, "Functional Genes (F)", icon = icon("check-circle"), color = "teal")
	})
	
	
	
	source("server/dataset_tables.R", local = TRUE, encoding = "UTF-8")
	# ══ MESSAGE DE BIENVENUE ════════════════════════════════════
	showModal(modalDialog(
		title = "🧬 Welcome to IMGT-PROMISE",
		HTML("
      <p>IMGT-PROMISE is an integrated environment for exploring and analyzing
      <strong>regulatory features of 5\u2032UTR regions</strong> in IG and TR V-genes.
      It combines <strong>statistical analyses</strong>, <strong>visualization</strong>,
      <strong>comparative approaches</strong>, and <strong>automatic motif discovery</strong>
      for species lacking curated motif annotations.</p>
      <hr>
      <p><strong>Getting started</strong></p>
      <ol>
        <li>Build your dataset (📊 Dataset Construction)</li>
        <li>Explore the analyses in the left menu, organized by research question</li>
        <li><strong>Export</strong> results, tables and figures at every step</li>
      </ol>
    "),
		easyClose = TRUE, footer = modalButton("Launch →")
	))
	
	# ══ SECTIONS DU SERVEUR 
	source("server/descriptive.R",           local = TRUE, encoding = "UTF-8")
	source("server/freq.R",                  local = TRUE, encoding = "UTF-8")
	source("server/regression.R",            local = TRUE, encoding = "UTF-8")
	source("server/logos.R",                 local = TRUE, encoding = "UTF-8")
	source("server/cons_func.R",             local = TRUE, encoding = "UTF-8")
	source("server/clustering.R",            local = TRUE, encoding = "UTF-8")
	source("server/diversity.R",             local = TRUE, encoding = "UTF-8")
	source("server/identity.R",              local = TRUE, encoding = "UTF-8")
	source("server/module_calls.R",          local = TRUE, encoding = "UTF-8")
	source("server/partial_genes_server.R",  local = TRUE, encoding = "UTF-8")
	source("server/resources.R",             local=TRUE,   encoding = "UTF-8")
	
	} # end server