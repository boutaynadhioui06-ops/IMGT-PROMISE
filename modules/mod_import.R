# ═══════════════════════════════════════════════════════════════
#  Import 5'UTR + Import Motifs
# ═══════════════════════════════════════════════════════════════
# Dépend de : dataset_helpers.R, fasta_parser.R
# Packages  : shiny, shinydashboard, DT 
# ══════════════════════════════════════════════════════════════
# UI DU MODULE IMPORT 5'UTR
# ══════════════════════════════════════════════════════════════
mod_import_ui <- function(id) {
	ns <- NS(id)
	tagList(
		
		fluidRow(
			#  Upload
			box(
				width = 12, title = "📁 Select the 5′ UTR file",
				status = "primary", solidHeader = TRUE,
				
				fluidRow(
					column(6,
						fileInput(ns("file_5utr"),
							label = NULL,
							accept = c(".txt", ".fasta", ".fa"),
							buttonLabel = "Browse...",
							placeholder = "No file selected"),
						
						tags$small(icon("info-circle"),
							" Accepted formats: .txt, .fasta, .fa — Unlimited file size")
					),
					column(6,
						div(
							style = "background-color: #FFF3CD; padding: 12px; border-radius: 5px; 
                       border-left: 4px solid #FFA726;",
							h5(icon("filter"), strong("Header extraction only")),
							p(style = "margin-bottom: 5px; font-size: 13px;",
								icon("info-circle"), 
								" The application will extract", 
								strong("only the headers"), 
								"from the 5′ UTR file."),
							p(style = "margin-bottom: 0; font-size: 13px;",
								"The", strong("DNA sequences will be ignored."), 
								" Only metadata (GENE, ALLELE, SPECIES, etc.) will be retained.")
						)
					)
				),
				
				hr(),
				
				actionButton(ns("btn_parse"), "▶  Analyze file",
					icon  = icon("play"),
					class = "btn-success btn-lg",
					width = "30%")
			)
		),
		
		# Rapport + aperçu
		fluidRow(
			box(
				width = 5, title = "📋 Parsing report",
				status = "success", solidHeader = TRUE,
				verbatimTextOutput(ns("parse_report"))
			),
			box(
				width = 7, title = "👁️Preview of extracted headers (first 10 lines))",
				status = "info", solidHeader = TRUE,
				collapsible = TRUE, collapsed = FALSE,
				
				div(
					style = "background-color: #E8F5E9; padding: 10px; border-radius: 5px; 
                   margin-bottom: 10px; border-left: 3px solid #4CAF50;",
					p(style = "margin: 0; font-size: 13px;",
						icon("check-circle", style = "color: #4CAF50;"), 
						strong(" DNA sequences were ignored. "),
						" Only metadata were retained.")
				),
				
				DTOutput(ns("preview_table"))
			)
		)
	)
}

# ══════════════════════════════════════════════════════════════
# UI DU MODULE IMPORT MOTIFS 
# ══════════════════════════════════════════════════════════════
mod_motifs_ui <- function(id) {
	ns <- NS(id)
	tagList(
		
		#  Liste des motifs disponibles
		fluidRow(
			box(
				width = 12, title = "📂 Import motif files",
				status = "warning", solidHeader = TRUE,
				
				div(
					style = "background-color: #E3F2FD; padding: 12px; border-radius: 6px;
                   border-left: 4px solid #2196F3; margin-bottom: 15px;",
					p(style = "margin: 0; font-size: 13px;",
						icon("info-circle"),
						"Upload the motif files required for your analysis.",
						" Only imported motifs will be included in the final dataset.")
				),
				
				uiOutput(ns("motif_inputs"))
			)
		),
		
		#Statut global 
		fluidRow(
			box(
				width = 12, title = "✅ import status"
				,        status = "success", solidHeader = TRUE,
				uiOutput(ns("import_status"))
			)
		)
	)
}


# ══════════════════════════════════════════════════════════════
# SERVER 
# ══════════════════════════════════════════════════════════════
mod_import_server <- function(id) {
	moduleServer(id, function(input, output, session) {
		
		parsed_data <- reactiveVal(NULL)
		
		# Parsing au clic
		observeEvent(input$btn_parse, {
			req(input$file_5utr)
			
			
			withProgress(message = "Parsing in progress...", value = 0, {
				incProgress(0.3, detail = "Extracting headers only...")
				
				tryCatch({
					df <- parse_5utr_fasta(input$file_5utr$datapath)
					incProgress(0.6, detail = "Headers extracted...")
					parsed_data(df)
					incProgress(1,   detail = "Completed !")
					
					showNotification(
						HTML(paste0(
							"<strong>✅ extraction successful!</strong><br>",
							nrow(df), " headers extracted<br>",
							"<em>DNA sequences were ignored as expected.</em>"
						)),
						type = "message", 
						duration = 7
					)
					
				}, error = function(e) {
					showNotification(paste0("❌ error : ", e$message),
						type = "error", duration = 10)
				})
			})
		})
		
		# Rapport de parsing
		output$parse_report <- renderPrint({
			req(parsed_data())
			df <- parsed_data()
			
			cat("══════════════════════════════════\n")
			cat("  PARSING REPORT\n")
			cat("══════════════════════════════════\n\n")
			cat("• Total sequences :", nrow(df), "\n")
			cat("• Columns         :", ncol(df), "\n\n")
			
			cat("─── GROUP ───────────────────────\n")
			df  <- rename_loci_display(df)               
			
			tl <- table(df$GROUP)
			for (l in names(tl)) cat(" ", l, ":", tl[l], "seq.\n")
			
			cat("\n─── SPECIES (", length(unique(df$SPECIES)), " total) ─────\n", sep = "")
			sp_tbl <- sort(table(df$SPECIES), decreasing = TRUE)
			for (sp in names(sp_tbl)) cat("  -", sp, ":", sp_tbl[sp], "seq.\n")
			
			cat("\n─── FUNCTIONALITY ─────────────\n")
			tf <- table(df$FUNCTIONALITY)
			for (f in names(tf)) cat(" ", f, ":", tf[f], "\n")
			
			na_c <- colSums(is.na(df))
			na_c <- na_c[na_c > 0]
			cat("\n─── Missing values ─────────\n")
			if (length(na_c) == 0) {
				cat("  None\n")
			} else {
				for (col in names(na_c)) cat(" ", col, ":", na_c[col], "\n")
			}
			cat("══════════════════════════════════\n")
		})
		
		# Aperçu du tableau
		output$preview_table <- renderDT({
			req(parsed_data())
			df_display <- rename_loci_display(head(parsed_data(), 10))
			datatable(
				df_display,
				options  = list(scrollX = TRUE, dom = "t", pageLength = 10),
				rownames = FALSE,
				class    = "cell-border stripe hover"
			) %>%
				formatStyle("FUNCTIONALITY",
					backgroundColor = styleEqual(
						c("F",       "P",       "ORF"),
						c("#D4EDDA", "#F8D7DA", "#FFF3CD")
					),
					color = styleEqual(
						c("F",       "P",       "ORF"),
						c("#155724", "#721c24", "#856404")
					)
				) %>%
				formatStyle("GROUP", fontWeight = "bold")
		})
		
		return(parsed_data)
	})
}

# ══════════════════════════════════════════════════════════════
# SERVER 
# ══════════════════════════════════════════════════════════════
# Info-bulles pour les motifs ambigus selon le locus
.motif_info_text <- function(motif) {
	switch(motif,
		"HEPTA_PENTA" = "HEPTAMER for IGH · PENTADECAMER for IGK/IGL",
		"OCTA_DECA"   = "OCTAMER for IGH · DECAMER for IGK/IGL",
		NULL
	)
}

mod_motifs_server <- function(id, data_5utr) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns
		
		# Liste  les motifs disponibles
		ALL_MOTIFS <- c("TATA_BOX", "HEPTA_PENTA", "OCTA_DECA",
			"E_BOX", "LPART_1")
		
		# Stockage de tous les motifs importés
		motif_data <- reactiveVal(list())
		
		# ── Inputs 
		output$motif_inputs <- renderUI({
			
			boxes <- lapply(ALL_MOTIFS, function(motif) {
				info_txt <- .motif_info_text(motif)
				column(
					width = 4,
					div(
						style = "background:#FFFDF0; border:1px solid #FFC107;
                 border-radius:8px; padding:15px; margin-bottom:10px;",
						h5(icon("file-alt"), strong(motif),
							if (!is.null(info_txt))
								tags$span(
									icon("info-circle"),
									title = info_txt,
									style = "margin-left:6px; color:#2196F3; cursor:help; font-size:13px;"
								)
						),
						fileInput(
							ns(paste0("file_", motif)),
							label       = NULL,
							accept      = c(".txt", ".fasta", ".fa"),
							buttonLabel = "Import",
							placeholder = "No file"
						),
						uiOutput(ns(paste0("status_", motif)))
					)
				)
			})
			
			do.call(fluidRow, boxes)
		})
		
		#  Observer : parser chaque motif à l'upload
		lapply(ALL_MOTIFS, function(motif) {
			local({
				m         <- motif
				input_id  <- paste0("file_", m)
				status_id <- paste0("status_", m)
				
				observeEvent(input[[input_id]], {
					req(input[[input_id]])
					
					tryCatch({
						df_m    <- parse_motif_fasta(input[[input_id]]$datapath, m)
						current <- motif_data()
						current[[m]] <- df_m
						motif_data(current)
						
						output[[status_id]] <- renderUI({
							div(class = "alert alert-success",
								style = "padding:5px; margin:0;",
								icon("check-circle"),
								strong(nrow(df_m)), " Imported sequences") 
						})
						
						showNotification(
							paste0("✅ ", m, " : ", nrow(df_m), " Imported sequences)" ) ,  
							type = "message", duration = 3
						)
						
					}, error = function(e) {
						output[[status_id]] <- renderUI({
							div(class = "alert alert-danger",
								style = "padding:5px; margin:0;",
								icon("times-circle"), " Error : ", e$message)
						})
					})
				})
			})
		})
		
		# Statut global des imports
		output$import_status <- renderUI({
			md      <- motif_data()
			n_done  <- length(names(md))
			n_total <- length(ALL_MOTIFS)
			
			tagList(
				# Barre de progression
				div(style = "margin-bottom:10px;",
					tags$b(n_done, "/", n_total, " Imported motifs"),
					div(class = "progress",
						div(class = "progress-bar progress-bar-success",
							style = paste0("width:", round(n_done / n_total * 100), "%;"),
							paste0(round(n_done / n_total * 100), "%"))
					)
				),
				
				# Tableau statut
				tags$table(
					class = "table table-bordered table-sm",
					tags$thead(
						tags$tr(
							tags$th("Motif"),
							tags$th("Statut"),
							tags$th("Sequences")
						)
					),
					tags$tbody(
						lapply(ALL_MOTIFS, function(m) {
							if (m %in% names(md)) {
								tags$tr(
									tags$td(strong(m)),
									tags$td(span(class = "label label-success", "✅ Imported" )) ,  
									tags$td(nrow(md[[m]]))
								)
							} else {
								tags$tr(
									tags$td(m),
									tags$td(span(class = "label label-warning", "⏳ Pending ")),
									tags$td("—")
								)
							}
						})
					)
				),
				
				# Message si au moins 1 motif importé
				if (n_done > 0)
					div(class = "alert alert-success",
						icon("check-double"),
						strong(paste0(" ", n_done, " motif(s) imported !")),
						" You can proceed to dataset construction.")
			)
		})
		
		return(list(
			motif_data = motif_data
		))
	})
}
