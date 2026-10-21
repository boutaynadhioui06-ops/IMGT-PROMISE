# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  tabItem (clustring)  ████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════
# TAB: 5. Hierarchical Clustering
# ═══════════════════════════════════════════════════════════════
ui_tab_conserv_clustering <- function()
tabItem(tabName = "conserv_clustering",
	
	tab_header("project-diagram", "Hierarchical Clustering",
		"Which genes or subgroups share similar sequence or motif profiles?",
		"#4338CA", "#0B2A5B"),
	
	hr(),
	
	imgt_ligm_note(),
	
	# ── PARAMÈTRES ──
	fluidRow(
		box(
			width = 12,
			title = "⚙️ Parameters",
			status = "primary",
			solidHeader = TRUE,
			
			fluidRow(
				column(3,
					div(
						style = "background-color: #E8F5E9; padding: 15px; border-radius: 8px;",
						h5(icon("globe"), strong("Species")),
						selectInput("clust_species", label = NULL,
							choices = NULL,
							selected = NULL)
					)
				),
				
				column(4,
					div(
						style = "background-color: #FFF9C4; padding: 15px; border-radius: 8px;",
						h5(icon("layer-group"), strong("Clustering method")),
						radioButtons("clust_method", label = NULL,
							choices = c("Motif presence/absence" = "presence",
								"Sequence similarity" = "sequence"),
							selected = "presence")
					)
				),
				
				# Conditionnel : sélection du motif si similarité séquence
				column(3,
					conditionalPanel(
						condition = "input.clust_method == 'sequence'",
						div(
							style = "background-color: #F3E5F5; padding: 15px; border-radius: 8px;",
							h5(icon("dna"), strong("Motif")),
							selectInput("clust_motif", label = NULL,
								choices = NULL,
								selected = NULL)
						)
					)
				)
			),
			
			hr(),
			
			# ── Nombre de clusters k ─────────────────────────────────
			fluidRow(
				column(6,
					div(
						style = "background-color: #E3F2FD; padding: 15px; border-radius: 8px;",
						h5(icon("object-group"), strong("Number of clusters (k)")),
						sliderInput("clust_k",
							label = NULL,
							min = 2, max = 10, value = 3, step = 1,
							width = "100%"),
						tags$small(
							icon("info-circle"),
							" Colored boxes on the dendrogram show cluster assignments.",
							br(),
							icon("lightbulb"),
							" Adjust k to explore different groupings."
						)
					)
				),
				column(6,
					div(
						style = "background-color: #FFF3E0; padding: 15px; border-radius: 8px;",
						h5(icon("table"), strong("Cluster assignments")),
						p(style = "font-size: 12px; color: #555;",
							"After generating the dendrogram, the table below shows ",
							"which cluster each subgroup belongs to for the selected k."),
						uiOutput("clust_k_summary")
					)
				)
			),
			
			hr(),
			
			fluidRow(
				column(12, style = "text-align: center;",
					actionButton("btn_launch_clust", "🚀 Generate dendrogram",
						icon = icon("play"), class = "btn-success btn-lg",
						style = "width: 300px; height: 60px; font-size: 16px;")
				)
			)
		)
	),
	
	# ── RÉSULTATS ──
	fluidRow(
		box(
			width = 12,
			title = "📊 Dendrogram",
			status = "info",
			solidHeader = TRUE,
			
			uiOutput("clust_status"),
			hr(),
			
			div(
				style = "background-color: #E3F2FD; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
				icon("info-circle"), strong(" Interpretation: "),
				"Subgroups that are close in the dendrogram have similar profiles.",
				"The height reflects the level of dissimilarity between groups."
			),
			
			tabsetPanel(
				tabPanel("🌳 Dendrogram",
					br(),
					plotOutput("clust_plot", height = "700px"),
					br(),
					div(style = "margin-top:18px;",
						fluidRow(
							column(12,
								div(
									style = "text-align:center; padding:10px; background:#EFF6FF;
                               border-radius:8px; border:1px solid #93C5FD;",
									icon("image", style="font-size:20px; color:#1D4ED8; display:block; margin-bottom:6px;"),
									p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
										"Hierarchical clustering dendrogram"),
									downloadButton("btn_download_clust_plot",
										"📥 PNG — Dendrogram",
										class = "btn-primary btn-block",
										style = "width:100%;")
								)
							)
						)
					)
				),
				tabPanel("📋 Distance matrix",
					br(),
					DTOutput("clust_dist_table")
				),
				tabPanel("🏷️ Cluster Assignments",
					br(),
					div(
						style = "background-color: #E3F2FD; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
						icon("info-circle"),
						strong(" Cluster assignments for k = "),
						textOutput("clust_k_display", inline = TRUE),
						strong(" clusters"),
						" — adjust the k slider above and regenerate to update."
					),
					DTOutput("clust_assignments_table"),
					br(),
					div(style = "text-align:left; padding:10px; background:#FFFBEB;
                     border-radius:8px; border:1px solid #FCD34D; display:inline-block;",
						icon("file-csv", style="font-size:18px; color:#D97706; margin-right:6px;"),
						downloadButton("btn_download_clust_assignments",
							"📥 CSV — Cluster Assignments",
							class = "btn-warning", style = "margin-left:8px;")
					)
				)
			)
		)
	)
)
# ══════════════════════════════════════════════════════════════
# 4. NOUVEL ONGLET : Heatmap d'identité des séquences
# ══════════════════════════════════════════════════════════════
ui_tab_conserv_identity <- function()
tabItem(tabName = "conserv_identity",
	
	tab_header("th", "Sequence Identity Heatmap",
		"Which genes or subgroups share similar sequence or motif profiles?",
		"#0EA5E9", "#1E3A5F"),
	
	hr(),
	
	imgt_ligm_note(),
	
	fluidRow(
		box(
			width = 12, title = "⚙️ Parameters",
			status = "primary", solidHeader = TRUE,
			
			div(
				style = "background-color: #DBEAFE; padding: 10px; border-radius: 6px;
                 border-left: 4px solid #3B82F6; margin-bottom: 15px;",
				icon("info-circle"),
				strong(" Principle:"),
				" The consensus sequence (most frequent) of each entity is extracted",
				" and compared pairwise to compute % identity.",
				strong(" Threshold: ≥ 75% = high conservation (dark color).")
			),
			
			fluidRow(
				column(4,
					div(style = "background-color: #FFF9C4; padding: 15px; border-radius: 8px;",
						h5(icon("dna"), strong("Motif")),
						selectInput("identity_motif", label = NULL,
							choices = NULL, selected = NULL),
						tags$small(icon("info-circle"), "Only imported motifs are available.")
					)
				),
				column(4,
					div(style = "background-color: #E8F5E9; padding: 15px; border-radius: 8px;",
						h5(icon("layer-group"), strong("Comparison mode")),
						radioButtons("identity_mode", label = NULL,
							choices = c(
								"Between species"                  = "species",
								"Between loci of one species"      = "locus",
								"Between subgroups of one species" = "subgroup"
							),
							selected = "species")
					)
				),
				column(4,
					div(style = "padding-top: 30px;",
						actionButton("btn_launch_identity", "🚀 Generate heatmap",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height: 70px;")
					)
				)
			),
			
			#  Sélecteur espèce 
			conditionalPanel(
				condition = "input.identity_mode == 'locus' || input.identity_mode == 'subgroup'",
				hr(),
				div(
					style = "background-color: #DBEAFE; padding: 12px; border-radius: 8px;",
					fluidRow(
						column(5,
							div(
								style = "display: flex; align-items: center; gap: 12px;",
								tags$div(icon("globe"), style = "font-size:18px; color:#1E40AF;"),
								div(style = "flex: 1;",
									tags$label(strong("Species"),
										style = "color:#1E40AF; margin-bottom:4px; display:block;"),
									selectInput("identity_species", label = NULL,
										choices = NULL, selected = NULL, width = "100%")
								)
							)
						),
						column(7,
							tags$small(
								style = "color:#555; font-size:12px; padding-top:6px; display:block;",
								icon("info-circle"), strong(" Locus :"),
								"Compare IGH, IGK, and IGL for the selected species.", br(),
								icon("info-circle"), strong(" Subgroups:"),
								" compares all subgroups (all loci) for the selected species."
							)
						)
					)
				)
			)
		)
	),
	
	fluidRow(
		box(
			width = 12, title = "🔥 Identity Heatmap",
			status = "info", solidHeader = TRUE,
			
			uiOutput("identity_status"),
			hr(),
			
			div(
				style = "background-color: #EFF6FF; padding: 10px; border-radius: 6px; margin-bottom: 15px;",
				icon("info-circle"),
				strong(" Reading: "),
				"Dark blue cells (≥ 75%) = high conservation. Light cells = high divergence."
			),
			
			tabsetPanel(
				tabPanel("📈 Heatmap",
					br(),
					uiOutput("identity_plot_ui"),
					br(),
					div(style = "margin-top:18px;",
						fluidRow(
							column(6,
								div(
									style = "text-align:center; padding:10px; background:#EFF6FF;
                               border-radius:8px; border:1px solid #93C5FD;",
									icon("image", style="font-size:20px; color:#1D4ED8; display:block; margin-bottom:6px;"),
									p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
										"Identity heatmap (PNG, 300 dpi)"),
									downloadButton("btn_download_identity_png",
										"📥 PNG — Identity Heatmap",
										class = "btn-primary btn-block",
										style = "width:100%;")
								)
							),
							column(6,
								div(
									style = "text-align:center; padding:10px; background:#FFFBEB;
                               border-radius:8px; border:1px solid #FCD34D;",
									icon("file-csv", style="font-size:20px; color:#D97706; display:block; margin-bottom:6px;"),
									p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
										"Pairwise identity matrix"),
									downloadButton("btn_download_identity_csv",
										"📥 CSV — Identity Matrix",
										class = "btn-warning btn-block",
										style = "width:100%;")
								)
							)
						)
					)
				),
				tabPanel("📋 Identity_table",
					br(),
					DTOutput("identity_table")
				)
			)
		)
	)
)
# ══════════════════════════════════════════════════════════════
# Tab E — Motif Phylogeny
# ══════════════════════════════════════════════════════════════
ui_tab_adv_tsne <- function()
tabItem(tabName = "adv_tsne",
	tab_header("braille", "Motif Phylogeny",
		"Which genes or subgroups share similar sequence or motif profiles?",
		"#7C3AED", "#0B2A5B"),
	uiOutput("warn_no_dataset_tsne"),
	hr(),
	imgt_ligm_note(),
	mod_tsne_ui("adv_tsne_mod")
)