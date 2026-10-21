# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  tabItem (descriptive)  ████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════
ui_tab_stat_desc_tableau <- function()
tabItem(tabName = "stat_desc_tableau",
	
	tab_header("table", "Descriptive Table — Full Dataset Characterization",
		"How are sequences and gene functionality distributed across species, loci and subgroups?",
		"#001f54", "#1565C0"),
	
	hr(),
	
	# ── Description des sections ──────────────────────────────
	fluidRow(
		column(12,
			div(
				style = "display:flex; gap:8px; flex-wrap:wrap; margin-bottom:16px;",
				div(style = "background:#001f54; color:white; padding:6px 14px;
                     border-radius:6px; font-size:12px; font-weight:bold;",
					"Identification"),
				div(style = "background:#1565C0; color:white; padding:6px 14px;
                     border-radius:6px; font-size:12px; font-weight:bold;",
					"Counts"),
				div(style = "background:#166534; color:white; padding:6px 14px;
                     border-radius:6px; font-size:12px; font-weight:bold;",
					"Motif presence rate (%)"),
				div(style = "background:#854D0E; color:white; padding:6px 14px;
                     border-radius:6px; font-size:12px; font-weight:bold;",
					"Mean distance to LPART1 (nt)"),
				div(style = "background:#581C87; color:white; padding:6px 14px;
                     border-radius:6px; font-size:12px; font-weight:bold;",
					"5\\'UTR length")
			)
		)
	),
	
	fluidRow(
		box(
			width = 12, title = "⚙️ Parameters",
			status = "primary", solidHeader = TRUE,
			
			fluidRow(
				column(3,
					div(style = "background-color: #E3F2FD; padding: 15px; border-radius: 8px;",
						h5(icon("dna"), strong("Data type")),
						radioButtons("tableau_type", label = NULL,
							choices = c("IG" = "IG", "TR" = "TR"), selected = "IG")
					)
				),
				column(4,
					div(style = "background-color: #E8F5E9; padding: 15px; border-radius: 8px;",
						h5(icon("globe"), strong("Species selection")),
						selectInput("tableau_species", label = NULL,
							choices = NULL, selected = NULL),
						tags$small(
							icon("info-circle"),
							" The table shows ALL loci and subgroups for this species.",
							br(),
							icon("info-circle"),
							" \"—\" = motif absent or not imported for this subgroup."
						)
					)
				),
				column(3, offset = 2,
					div(style = "padding-top: 30px;",
						actionButton("btn_launch_tableau", "🚀 Generate table",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height: 70px;")
					)
				)
			)
		)
	),
	
	fluidRow(
		box(
			width = 12, title = "📊 Generated Table",
			status = "info", solidHeader = TRUE,
			
			uiOutput("tableau_status"),
			hr(),
			
			# Notes de tableau
			div(
				style = "background:#F8FAFF; border-left:3px solid #1565C0;
                 padding:8px 14px; margin-bottom:10px; font-size:12px; color:#374151;",
				tags$sup("1"), " Distances calculated: ",
				strong("POS_START_LPART1 − POS_END_motif"),
				" (gap between motif end and L-PART1 start, in nucleotides).",
				br(),
				tags$sup("2"), " Rate = proportion of alleles in the subgroup carrying the motif.",
				br(),
				tags$sup("3"), " ",
				tags$span(style="color:#155724; font-weight:bold;", "F"),
				" = Functional, ",
				tags$span(style="color:#721C24; font-weight:bold;", "P"),
				" = Pseudogene, ",
				tags$span(style="color:#856404; font-weight:bold;", "ORF"),
				" = Open Reading Frame."
			),
			
			DTOutput("tableau_table"),
			br(),
			hr(),
			
			div(style = "margin-top:18px;",
				fluidRow(
					column(6,
						div(
							style = "text-align:center; padding:10px; background:#FFFBEB;
                   border-radius:8px; border:1px solid #FCD34D;",
							icon("file-csv", style="font-size:20px; color:#D97706; display:block; margin-bottom:6px;"),
							p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
								"Raw data table"),
							downloadButton("btn_download_tableau_csv",
								"📥 CSV — Data",
								class = "btn-warning btn-block",
								style = "width:100%;")
						)
					),
					column(6,
						div(
							style = "text-align:center; padding:10px; background:#EFF6FF;
                   border-radius:8px; border:1px solid #93C5FD;",
							icon("file-excel", style="font-size:20px; color:#1D4ED8; display:block; margin-bottom:6px;"),
							p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
								"Coloured sections + auto widths"),
							downloadButton("btn_download_tableau_excel",
								"📊 Structured Excel",
								class = "btn-primary btn-block",
								style = "width:100%;")
						)
					)
				)
			)
		)
	)
)
# ═══════════════════════════════════════════════════════════════
# ONGLET : b. Distribution des Sequences
# ═══════════════════════════════════════════════════════════════
ui_tab_stat_desc_sequences <- function()
tabItem(tabName = "stat_desc_sequences",
	
	tab_header("chart-bar", "Sequence Distribution",
		"How are sequences and gene functionality distributed across species, loci and subgroups?",
		"#1565C0", "#16A34A"),
	
	hr(),
	
	fluidRow(
		box(
			width = 12, title = "⚙️ Parameters",
			status = "primary", solidHeader = TRUE,
			
			fluidRow(
				column(3,
					div(style = "background-color: #E3F2FD; padding: 15px; border-radius: 8px;",
						h5(icon("dna"), strong("Data type")),
						radioButtons("seq_type", label = NULL,
							choices = c("IG" = "IG", "TR" = "TR"),
							selected = "IG")
					)
				),
				column(4,
					div(style = "background-color: #FFF9C4; padding: 15px; border-radius: 8px;",
						h5(icon("layer-group"), strong("Grouping level")),
						radioButtons("seq_grouping", label = NULL,
							choices = c("By species"     = "species",
								"By locus"       = "locus",
								"By subgroup"    = "subgroup"),
							selected = "species")
					)
				),
				column(4,
					div(style = "padding-top: 40px;",
						actionButton("btn_launch_seq", "🚀 Launch analysis",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height: 70px;")
					)
				)
			),
			
			# Filtre espèce  en mode LOCUS
			shiny::conditionalPanel(
				condition = "input.seq_grouping == 'locus'",
				shiny::hr(),
				shiny::div(
					style = "background-color: #F0FDF4; padding: 14px; border-radius: 8px;
             border-left: 4px solid #22C55E;",
					shiny::fluidRow(
						shiny::column(4,
							shiny::div(
								style = "display: flex; align-items: center; gap: 10px;",
								shiny::tags$div(shiny::icon("globe"), style = "font-size:18px; color:#166534;"),
								shiny::div(style = "flex: 1;",
									shiny::tags$label(shiny::strong("Filter by species (optional)"),
										style = "color:#166534; margin-bottom:4px; display:block;"),
									shiny::selectInput("seq_locus_species_filter", label = NULL,
										choices  = c("All species" = ""),
										selected = "", width = "100%")
								)
							)
						),
						shiny::column(8,
							shiny::tags$small(
								style = "color:#374151; font-size:12px; padding-top:8px; display:block;",
								shiny::icon("info-circle"),
								" Empty = all species displayed together.",
								shiny::br(),
								shiny::icon("lightbulb"),
								" Each panel = one species × one locus."
							)
						)
					)
				)
			),
			
			shiny::conditionalPanel(
				condition = "input.seq_grouping == 'subgroup'",
				shiny::hr(),
				shiny::div(
					style = "background-color: #EFF6FF; padding: 14px; border-radius: 8px;
             border-left: 4px solid #3B82F6;",
					shiny::fluidRow(
						# Espèce
						shiny::column(3,
							shiny::div(
								style = "display: flex; align-items: center; gap: 10px;",
								shiny::tags$div(shiny::icon("globe"), style = "font-size:18px; color:#1E40AF;"),
								shiny::div(style = "flex: 1;",
									shiny::tags$label(shiny::strong("Species"),
										style = "color:#1E40AF; margin-bottom:4px; display:block;"),
									shiny::selectInput("seq_species_filter", label = NULL,
										choices = NULL, selected = NULL, width = "100%")
								)
							)
						),
						shiny::column(3,
							shiny::div(
								style = "display: flex; align-items: center; gap: 10px;",
								shiny::tags$div(shiny::icon("layer-group"), style = "font-size:18px; color:#1E40AF;"),
								shiny::div(style = "flex: 1;",
									shiny::tags$label(shiny::strong("Locus (optional)"),
										style = "color:#1E40AF; margin-bottom:4px; display:block;"),
									shiny::selectInput("seq_locus_filter", label = NULL,
										choices = c("All loci"), selected = "All loci", width = "100%")
								)
							)
						),
						# Mode d'affichage
						shiny::column(3,
							shiny::tags$label(
								shiny::strong("Display mode"),
								shiny::tags$span(
									title = "Auto: switches to heatmap above 80 subgroups. Barplot: horizontal bars, best for < 80 subgroups. Heatmap: color matrix, ideal for many subgroups — shows sequence count per cell.",
									style = "cursor:help; margin-left:5px; color:#60A5FA; font-size:13px;",
									"ℹ️"
								),
								style = "color:#1E40AF; margin-bottom:4px; display:block;"
							),
							shiny::radioButtons("seq_sg_mode", label = NULL,
								choices = c(
									"🔄 Auto"    = "auto",
									"📊 Barplot" = "barplot",
									"🔥 Heatmap" = "heatmap"
								),
								selected = "auto", inline = TRUE
							)
						),
						# Top-N
						shiny::column(3,
							shiny::tags$label(shiny::strong("Top-N subgroups"),
								style = "color:#1E40AF; margin-bottom:4px; display:block;"),
							shiny::fluidRow(
								shiny::column(7,
									shiny::numericInput("seq_topn", label = NULL,
										value = NA, min = 5, max = 500, step = 5, width = "100%")
								),
								shiny::column(5,
									shiny::tags$small(style = "color:#6B7280; font-size:11px; padding-top:8px; display:block;",
										"Empty = all")
								)
							)
						)
					),
					shiny::uiOutput("seq_sg_alert"),
					shiny::tags$small(
						style = "color:#4B5563; font-size:11px; margin-top:4px; display:block;",
						shiny::icon("lightbulb"),
						shiny::strong(" Heatmap (> 80 SG):"),
						" shows sequence count per subgroup × locus (blue intensity).",
						shiny::br(),
						shiny::icon("filter"),
						shiny::strong(" Locus filter:"),
						" restricts barplot display to subgroups of a single locus."
					)
				)
			)
		)
	),
	
	fluidRow(
		box(
			width = 12, title = "📊 Results",
			status = "info", solidHeader = TRUE,
			
			hr(),
			
			br(),
			uiOutput("seq_plot_ui")
		)
	),
	
	# Boutons download  
	fluidRow(
		box(
			width = 12, title = "⬇️ Downloads",
			status = "success", solidHeader = FALSE,
			collapsible = TRUE,
			fluidRow(
				column(4,
					div(
						style = "text-align:center; padding:10px; background:#EFF6FF;
                   border-radius:8px; border:1px solid #93C5FD;",
						icon("image", style="font-size:20px; color:#1D4ED8; display:block; margin-bottom:6px;"),
						p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
							"Full figure"),
						downloadButton("btn_download_seq_plot_png",
							"📥 High-resolution PNG",
							class = "btn-primary btn-block",
							style = "width:100%;")
					)
				),
				column(4,
					div(
						style = "text-align:center; padding:10px; background:#F5F3FF;
                   border-radius:8px; border:1px solid #C4B5FD;",
						icon("file-pdf", style="font-size:20px; color:#7C3AED; display:block; margin-bottom:6px;"),
						p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
							"Vector format"),
						downloadButton("btn_download_seq_plot_pdf",
							"📄 Vector PDF ",
							class = "btn-info btn-block",
							style = "width:100%;")
					)
				),
				column(4,
					div(
						style = "text-align:center; padding:10px; background:#FFFBEB;
                   border-radius:8px; border:1px solid #FCD34D;",
						icon("file-csv", style="font-size:20px; color:#D97706; display:block; margin-bottom:6px;"),
						p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
							"Raw data table"),
						downloadButton("btn_download_seq_csv",
							"📥 CSV — Data",
							class = "btn-warning btn-block",
							style = "width:100%;")
					)
				)
			)
		)
	)
)
# ═══════════════════════════════════════════════════════════════
# ONGLET : c. Distribution de la Functionality
# ═══════════════════════════════════════════════════════════════
ui_tab_stat_desc_functionality <- function()
tabItem(tabName = "stat_desc_functionality",
	
	tab_header("dna", "Functionality Distribution",
		"How are sequences and gene functionality distributed across species, loci and subgroups?",
		"#16A34A", "#0B2A5B"),
	
	hr(),
	
	fluidRow(
		box(
			width = 12, title = "⚙️ Parameters",
			status = "primary", solidHeader = TRUE,
			
			fluidRow(
				column(3,
					div(style = "background-color: #E3F2FD; padding: 15px; border-radius: 8px;",
						h5(icon("dna"), strong("Data type")),
						radioButtons("func_type", label = NULL,
							choices = c("IG" = "IG", "TR" = "TR"),
							selected = "IG")
					)
				),
				column(4,
					div(style = "background-color: #FFF9C4; padding: 15px; border-radius: 8px;",
						h5(icon("layer-group"), strong("Grouping level")),
						radioButtons("func_grouping", label = NULL,
							choices = c("By species"     = "species",
								"By locus"       = "locus",
								"By subgroup"    = "subgroup"),
							selected = "species")
					)
				),
				column(4,
					div(style = "padding-top: 40px;",
						actionButton("btn_launch_func", "🚀 Launch analysis",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height: 70px;")
					)
				)
			),
			
			# ── séléct espèce  en mode locus 
			shiny::conditionalPanel(
				condition = "input.func_grouping == 'locus'",
				shiny::hr(),
				shiny::div(
					style = "background-color: #F0FDF4; padding: 14px; border-radius: 8px;
             border-left: 4px solid #22C55E;",
					shiny::fluidRow(
						shiny::column(4,
							shiny::div(
								style = "display: flex; align-items: center; gap: 10px;",
								shiny::tags$div(shiny::icon("globe"), style = "font-size:18px; color:#166534;"),
								shiny::div(style = "flex: 1;",
									shiny::tags$label(shiny::strong("Filter by species (optional)"),
										style = "color:#166534; margin-bottom:4px; display:block;"),
									shiny::selectInput("func_locus_species_filter", label = NULL,
										choices  = c("All species" = ""),
										selected = "", width = "100%")
								)
							)
						),
						shiny::column(8,
							shiny::tags$small(
								style = "color:#374151; font-size:12px; padding-top:8px; display:block;",
								shiny::icon("info-circle"),
								" Empty = all species displayed together.",
								shiny::br(),
								shiny::icon("lightbulb"),
								" Each panel = one species × one locus."
							)
						)
					)
				)
			),
			
			# ── Sélecteur espèce  en mode sous-groupe ──
			
			shiny::conditionalPanel(
				condition = "input.func_grouping == 'subgroup'",
				shiny::hr(),
				shiny::div(
					style = "background-color: #EFF6FF; padding: 14px; border-radius: 8px;
             border-left: 4px solid #3B82F6;",
					shiny::fluidRow(
						shiny::column(3,
							shiny::div(
								style = "display: flex; align-items: center; gap: 10px;",
								shiny::tags$div(shiny::icon("globe"), style = "font-size:18px; color:#1E40AF;"),
								shiny::div(style = "flex: 1;",
									shiny::tags$label(shiny::strong("Species"),
										style = "color:#1E40AF; margin-bottom:4px; display:block;"),
									shiny::selectInput("func_species_filter", label = NULL,
										choices = NULL, selected = NULL, width = "100%")
								)
							)
						),
						#  filtre locus pour fonctionnalité
						shiny::column(3,
							shiny::div(
								style = "display: flex; align-items: center; gap: 10px;",
								shiny::tags$div(shiny::icon("layer-group"), style = "font-size:18px; color:#1E40AF;"),
								shiny::div(style = "flex: 1;",
									shiny::tags$label(shiny::strong("Locus (optional)"),
										style = "color:#1E40AF; margin-bottom:4px; display:block;"),
									shiny::selectInput("func_locus_filter", label = NULL,
										choices = c("All loci"), selected = "All loci", width = "100%")
								)
							)
						),
						shiny::column(3,
							shiny::tags$label(shiny::strong("Display mode"),
								style = "color:#1E40AF; margin-bottom:4px; display:block;"),
							shiny::radioButtons("func_sg_mode", label = NULL,
								choices = c(
									"Auto"           = "auto",
									"Force barplot"  = "barplot",
									"Force heatmap"  = "heatmap"
								),
								selected = "auto", inline = TRUE
							)
						),
						shiny::column(3,
							shiny::tags$label(shiny::strong("Top-N subgroups"),
								style = "color:#1E40AF; margin-bottom:4px; display:block;"),
							shiny::fluidRow(
								shiny::column(7,
									shiny::numericInput("func_topn", label = NULL,
										value = NA, min = 5, max = 500, step = 5, width = "100%")
								),
								shiny::column(5,
									shiny::tags$small(style = "color:#6B7280; font-size:11px; padding-top:8px; display:block;",
										"Empty = all")
								)
							)
						)
					),
					shiny::uiOutput("func_sg_alert"),
					shiny::tags$small(
						style = "color:#4B5563; font-size:11px; margin-top:4px; display:block;",
						shiny::icon("lightbulb"),
						shiny::strong(" Auto :"),
						" Automatic switch to an F×P×ORF heatmap beyond 80 subgroups.",
						shiny::br(),
						shiny::icon("palette"),
						" Colors: ",
						shiny::tags$span(style = "color:#28A745; font-weight:bold;", "■ F"),
						" — ",
						shiny::tags$span(style = "color:#DC3545; font-weight:bold;", "■ P"),
						" — ",
						shiny::tags$span(style = "color:#FFC107; font-weight:bold;", "■ ORF")
					)
				)
			)
		)
	),
	
	fluidRow(
		box(
			width = 12, title = "📊 Results",
			status = "info", solidHeader = TRUE,
			
			hr(),
			
			# Legendcouleurs 
			div(
				style = "background-color: #FFF3E0; padding: 10px; border-radius: 5px;
                 margin-bottom: 15px; display: inline-block; width: 100%;",
				icon("palette"), strong(" Legend:"),
				tags$span(style = "color:#28A745; font-weight:bold; margin-left:10px;",
					"■ F = Functional"),
				tags$span(style = "color:#DC3545; font-weight:bold; margin-left:15px;",
					"■ P = Pseudogene"),
				tags$span(style = "color:#FFC107; font-weight:bold; margin-left:15px;",
					"■ ORF")
			),
			
			br(),
			uiOutput("func_plot_ui")
		)
	),
	
	# Boutons download 
	fluidRow(
		box(
			width = 12, title = "⬇️ Downloads",
			status = "success", solidHeader = FALSE,
			collapsible = TRUE,
			fluidRow(
				column(4,
					div(
						style = "text-align:center; padding:10px; background:#EFF6FF;
                   border-radius:8px; border:1px solid #93C5FD;",
						icon("image", style="font-size:20px; color:#1D4ED8; display:block; margin-bottom:6px;"),
						p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
							"Full figure"),
						downloadButton("btn_download_func_plot",
							"📥 High-resolution PNG",
							class = "btn-primary btn-block",
							style = "width:100%;")
					)
				),
				
				column(4,
					div(
						style = "text-align:center; padding:10px; background:#FFFBEB;
                   border-radius:8px; border:1px solid #FCD34D;",
						icon("file-csv", style="font-size:20px; color:#D97706; display:block; margin-bottom:6px;"),
						p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
							"Raw data table"),
						downloadButton("btn_download_func_csv",
							"📥 CSV — Data",
							class = "btn-warning btn-block",
							style = "width:100%;")
					)
				)
			)
		)
	)
)
