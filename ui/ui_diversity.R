# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  tabItem (diversity)  ████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════
# ONGLET : Diversity des Motifs
# ═══════════════════════════════════════════════════════════════
ui_tab_diversity <- function()
tabItem(tabName = "diversity",
	
	tab_header("cubes", "Motif Diversity",
		"Do motif profiles differ between biological groups (species, loci, functionality)?",
		"#9C27B0", "#0B2A5B"),
	
	hr(),
	imgt_ligm_note(),
	
	# ── PARAMÈTRES ──────────────────────────────────────────────
	fluidRow(
		box(
			width = 12, title = "⚙️ Analysis parameters",
			status = "primary", solidHeader = TRUE,
			
			div(
				style = "background-color:#DBEAFE; padding:10px; border-radius:6px;
                 border-left:4px solid #3B82F6; margin-bottom:15px;",
				icon("info-circle"),
				strong("Sequence diversity quantified using Shannon entropy, Simpson’s index, richness, and evenness measures")
			),
			
			fluidRow(
				# Niveau de groupement
				column(3,
					div(style = "background-color:#FFF9C4; padding:15px; border-radius:8px;",
						h5(icon("layer-group"), strong("Analysis level")),
						radioButtons("div_grouping", label = NULL,
							choices = c("By species" = "species",
								"By locus"   = "locus",
								"By subgroup" = "subgroup"),
							selected = "species")
					)
				),
				# Filtre espèce en mode locus et subgroup
				column(4,
					shiny::conditionalPanel(
						condition = "input.div_grouping == 'locus' || input.div_grouping == 'subgroup'",
						div(style = "background-color:#E8F5E9; padding:15px; border-radius:8px;",
							h5(icon("globe"), strong("Species (optional)")),
							selectInput("div_species_filter", label = NULL,
								choices  = c("All species" = "All species"),
								selected = "All species",
								width    = "100%"),
							tags$small(icon("info-circle"),
								" 'All species' = all species together")
						)
					)
				),
				# Filtre locus en mode subgroup
				column(2,
					shiny::conditionalPanel(
						condition = "input.div_grouping == 'subgroup'",
						div(style = "background-color:#FFF3E0; padding:15px; border-radius:8px;",
							h5(icon("filter"), strong("Locus")),
							selectInput("div_locus_filter", label = NULL,
								choices  = c("All loci", "IGHV", "IGKV", "IGLV"),
								selected = "All loci",
								width    = "100%"),
							tags$small(icon("info-circle"),
								" Optional")
						)
					)
				),
				# Bouton
				column(3,
					div(style = "padding-top:30px;",
						actionButton("btn_launch_div", "🚀 Launch analysis",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height:65px; font-size:15px;")
					)
				)
			)
		)
	),
	
	# ── INFO INDICES ─────────────────────────────────────────────
	fluidRow(
		box(
			width = 12, title = "ℹ️ Diversity indices explained",
			status = "info", solidHeader = FALSE,
			collapsible = TRUE, collapsed = TRUE,
			fluidRow(
				column(3,
					div(style = "background-color:#E8F5E9; padding:10px; border-radius:5px;",
						h5(icon("chart-line"), strong("Shannon (H')")),
						p("Overall diversity. Higher = greater diversity."),
						tags$small("H' = -Σ(pi × ln(pi))"))),
				column(3,
					div(style = "background-color:#FFF9C4; padding:10px; border-radius:5px;",
						h5(icon("percent"), strong("Simpson (D)")),
						p("Probability 2 sequences have different motifs."),
						tags$small("D = 1 - Σ(pi²)"))),
				column(3,
					div(style = "background-color:#E3F2FD; padding:10px; border-radius:5px;",
						h5(icon("hashtag"), strong("Richness (S)")),
						p("Total number of distinct motifs detected."),
						tags$small("S = number of motifs present"))),
				column(3,
					div(style = "background-color:#F3E5F5; padding:10px; border-radius:5px;",
						h5(icon("balance-scale"), strong("Evenness (J)")),
						p("Uniformity of motif distribution (0–1)."),
						tags$small("J = H' / ln(S)")))
			)
		)
	),
	
	# ── RÉSULTATS ────────────────────────────────────────────────
	fluidRow(
		box(
			width = 12, title = "📊 Analysis results",
			status = "info", solidHeader = TRUE,
			
			uiOutput("div_status"),
			hr(),
			
			tabsetPanel(
				
				# ── Tab 1 : Barplot groupé ──────────────────────────────
				tabPanel("📊 Barplot",
					br(),
					uiOutput("div_barplot_ui"),
					br()
				),
				
				
				# ── Tab 6 : Radar chart ────────────────────────────────
				tabPanel("🎯 Radar Chart",
					br(),
					div(
						style = "background-color:#F3E5F5; padding:10px; border-radius:6px;
                     border-left:4px solid #9C27B0; margin-bottom:12px;",
						icon("info-circle"),
						strong(" Normalized radar chart"),
						" comparing all 4 indices simultaneously (values normalized 0–1)."
					),
					plotlyOutput("div_radar", height = "650px")
				),
				
				# ── Tab 7 : Tests statistiques ─────────────────────────
				tabPanel("📊 Statistical tests",
					br(),
					div(
						style = "background-color:#E8F5E9; padding:12px; border-radius:6px; margin-bottom:15px;",
						icon("info-circle"), strong(" Kruskal-Wallis tests — computed on individual sequences"),
						tags$ul(
							tags$li(strong("N motifs per sequence:"), " tests whether the total number of motifs per sequence differs between groups."),
							tags$li(strong("% motifs per sequence:"), " tests the proportion of motifs present per sequence."),
							tags$li(strong("Per motif:"), " tests whether presence of each motif varies between groups."),
							tags$li("p < 0.05 ★ | p < 0.01 ★★ | p < 0.001 ★★★ | ns = not significant")
						)
					),
					DTOutput("div_tests_table")
				)
			)
		)
	),
	
	# ── TÉLÉCHARGEMENTS  des figures 
	fluidRow(
		box(
			width = 12, title = "⬇️ Downloads",
			status = "success", solidHeader = FALSE,
			collapsible = TRUE,
			fluidRow(
				column(4,
					div(style = "text-align:center; padding:10px; background:#EFF6FF;
                   border-radius:8px; border:1px solid #93C5FD;",
						icon("image", style="font-size:20px; color:#1D4ED8; display:block; margin-bottom:6px;"),
						p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
							"Diversity barplot (selected index)"),
						downloadButton("btn_download_div_barplot", "📥 PNG — Barplot",
							class = "btn-primary btn-block", style = "width:100%;"))
				),
				column(4,
					div(style = "text-align:center; padding:10px; background:#F3E5F5;
                   border-radius:8px; border:1px solid #CE93D8;",
						icon("image", style="font-size:20px; color:#7B1FA2; display:block; margin-bottom:6px;"),
						p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
							"Radar chart (normalized profiles)"),
						downloadButton("btn_download_div_radar", "📥 PNG — Radar",
							class = "btn-info btn-block", style = "width:100%;"))
				),
				column(4,
					div(style = "text-align:center; padding:10px; background:#FFFBEB;
                   border-radius:8px; border:1px solid #FCD34D;",
						icon("file-csv", style="font-size:20px; color:#D97706; display:block; margin-bottom:6px;"),
						p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
							"Diversity indices table (CSV)"),
						downloadButton("btn_download_div_csv", "📥 CSV — Diversity data",
							class = "btn-warning btn-block", style = "width:100%;"))
				)
			)
		)
	)
) # fin tabItem diversity

# Tab F — PERMANOVA
ui_tab_adv_permanova <- function()
tabItem(tabName = "adv_permanova",
	tab_header("flask", "PERMANOVA",
		"Do motif profiles differ between biological groups (species, loci, functionality)?",
		"#7F1D1D", "#0B2A5B"),
	uiOutput("warn_no_dataset_permanova"),
	hr(),
	imgt_ligm_note(),
	mod_permanova_ui("adv_permanova_mod")
)
