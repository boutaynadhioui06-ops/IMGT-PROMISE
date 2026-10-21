
# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  tabItem(conservation)  ████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════
# ONGLET : 1. Fréquences de Présence des Motifs
# ═══════════════════════════════════════════════════════════════

ui_tab_conserv_frequencies <- function()
tabItem(tabName = "conserv_frequencies",
	
	tab_header("percent", "Motif Presence Frequencies",
		"How conserved are regulatory motifs across species, loci and subgroups?",
		"#166534", "#0B2A5B"),
	
	hr(),
	
	imgt_ligm_note(),
	
	# ── PARAMÈTRES ────────────────────────────────────────────────
	fluidRow(
		box(
			width = 12, title = "⚙️ Parameters",
			status = "primary", solidHeader = TRUE,
			
			div(
				style = "background-color: #DBEAFE; padding: 10px; border-radius: 6px;
                 border-left: 4px solid #3B82F6; margin-bottom: 15px;",
				icon("info-circle"),
				strong(" Analysis restricted to Immunoglobulins (IGH, IGK, IGL)"),
				" — TR genes currently lack sufficient 5'UTR motif data for robust analysis."
			),
			
			fluidRow(
				column(5,
					div(style = "background-color: #FFF9C4; padding: 15px; border-radius: 8px;",
						h5(icon("layer-group"), strong("Grouping level")),
						radioButtons("freq_grouping", label = NULL,
							choices = c("By species"  = "species",
								"By locus"    = "locus",
								"By subgroup" = "subgroup"),
							selected = "species")
					)
				),
				column(5,
					div(style = "padding-top: 30px;",
						actionButton("btn_launch_freq", "🚀 Launch analysis",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height: 70px;")
					)
				)
			),
			
			# Filtres  en mode LOCUS
			shiny::conditionalPanel(
				condition = "input.freq_grouping == 'locus'",
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
									shiny::selectInput("freq_locus_species_filter", label = NULL,
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
								" Each panel = one species × one locus (IGH, IGK, IGL)."
							)
						)
					)
				)
			),
			
			#  Filtres  en mode SUBGROUP ────────────────────────
			shiny::conditionalPanel(
				condition = "input.freq_grouping == 'subgroup'",
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
									shiny::selectInput("freq_species_filter", label = NULL,
										choices = NULL, selected = NULL, width = "100%")
								)
							)
						),
						# Locus
						shiny::column(3,
							shiny::div(
								style = "display: flex; align-items: center; gap: 10px;",
								shiny::tags$div(shiny::icon("layer-group"), style = "font-size:18px; color:#1E40AF;"),
								shiny::div(style = "flex: 1;",
									shiny::tags$label(shiny::strong("Locus (optional)"),
										style = "color:#1E40AF; margin-bottom:4px; display:block;"),
									selectInput("freq_locus_filter", label = NULL,
										choices  = c("All loci" = "All loci", "IGHV" = "IGHV", "IGKV" = "IGKV", "IGLV" = "IGLV"),
										selected = "All loci", width = "100%")
								)
							)
						),
						# Mode affichage
						shiny::column(3,
							shiny::tags$label(shiny::strong("Display mode"),
								shiny::tags$span(
									title = "Auto: switches to heatmap above 80 subgroups.",
									style = "cursor:help; margin-left:5px; color:#60A5FA; font-size:13px;",
									"ℹ️"),
								style = "color:#1E40AF; margin-bottom:4px; display:block;"),
							shiny::radioButtons("freq_sg_mode", label = NULL,
								choices  = c("🔄 Auto" = "auto", "📊 Barplot" = "barplot", "🔥 Heatmap" = "heatmap"),
								selected = "auto", inline = TRUE)
						),
						# Top-N
						shiny::column(3,
							shiny::tags$label(shiny::strong("Top-N subgroups"),
								style = "color:#1E40AF; margin-bottom:4px; display:block;"),
							shiny::fluidRow(
								shiny::column(7,
									shiny::numericInput("freq_topn", label = NULL,
										value = NA, min = 5, max = 500, step = 5, width = "100%")
								),
								shiny::column(5,
									shiny::tags$small(style = "color:#6B7280; font-size:11px;
                                     padding-top:8px; display:block;", "Empty = all")
								)
							)
						)
					),
					shiny::uiOutput("freq_sg_alert"),
					shiny::tags$small(
						style = "color:#4B5563; font-size:11px; margin-top:4px; display:block;",
						shiny::icon("lightbulb"), shiny::strong(" Heatmap:"),
						" Y = subgroups, X = motifs, intensity = % presence.",
						shiny::br(),
						shiny::icon("filter"), shiny::strong(" Locus:"),
						" restricts display to one locus."
					)
				)
			)
		)
	),
	
	# ── RÉSULTATS ─────────────────────────────────────────────────
	fluidRow(
		box(
			width = 12, title = "📊 Results",
			status = "info", solidHeader = TRUE,
			
			hr(),
			
			tabsetPanel(
				tabPanel("📈 Frequency Plot",
					br(),
					uiOutput("freq_plot_ui")
				),
				tabPanel("🔬 Statistical Tests",
					br(),
					div(
						style = "background-color: #E8F5E9; padding: 12px; border-radius: 6px; margin-bottom: 15px;",
						icon("info-circle"), strong(" Statistical tests on motif presence frequencies"),
						tags$ul(
							tags$li(strong("Chi-square test:"), " used when all expected cell counts ≥ 5."),
							tags$li(strong("Fisher's exact test:"), " used when any expected count < 5 (simulated p-value for tables > 2×2)."),
							tags$li(strong("Bonferroni correction:"), " applied across all motifs tested."),
							tags$li("p < 0.05 ★ | p < 0.01 ★★ | p < 0.001 ★★★ | ns = not significant")
						)
					),
					DTOutput("freq_stat_tests_table"),
					br(),
					div(style = "text-align:left; padding:10px; background:#FFFBEB;
                         border-radius:8px; border:1px solid #FCD34D; display:inline-block;",
						icon("file-csv", style="font-size:18px; color:#D97706; margin-right:6px;"),
						downloadButton("btn_download_freq_stat_csv",
							"📥 CSV — Statistical Tests",
							class = "btn-warning", style = "margin-left:8px;")
					)
				),
				tabPanel("🔬 Locus Pairwise Tests",
					br(),
					div(
						style = "background-color: #DBEAFE; padding: 12px; border-radius: 6px; margin-bottom: 15px;",
						icon("info-circle"), strong(" Which loci differ significantly?"),
						tags$ul(
							tags$li("Available only when ", strong("Grouping level = By locus"), "."),
							tags$li("Each row compares ", strong("two loci at a time"), " (e.g. IGH vs IGK) for one motif."),
							tags$li(strong("Chi-square test:"), " used when all expected cell counts \u2265 5; ",
								strong("Fisher's exact test:"), " used otherwise."),
							tags$li(strong("Bonferroni correction:"), " applied across all motif \u00d7 locus-pair comparisons."),
							tags$li("p < 0.05 \u2605 | p < 0.01 \u2605\u2605 | p < 0.001 \u2605\u2605\u2605 | ns = not significant")
						)
					),
					DTOutput("freq_pairwise_tests_table"),
					br(),
					div(style = "text-align:left; padding:10px; background:#FFFBEB;
                         border-radius:8px; border:1px solid #FCD34D; display:inline-block;",
						icon("file-csv", style="font-size:18px; color:#D97706; margin-right:6px;"),
						downloadButton("btn_download_freq_pairwise_csv",
							"📥 CSV — Pairwise Locus Tests",
							class = "btn-warning", style = "margin-left:8px;")
					)
				)
			)
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
							"Full figure (all species/loci)"),
						downloadButton("btn_download_freq_plot",
							"📥 PNG — Full figure",
							class = "btn-primary btn-block",
							style = "width:100%;")
					)
				),
				column(4,
					div(
						style = "text-align:center; padding:10px; background:#F0FDF4;
                   border-radius:8px; border:1px solid #86EFAC;",
						icon("filter", style="font-size:20px; color:#16A34A; display:block; margin-bottom:6px;"),
						p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
							"Selected species (subgroup mode)"),
						downloadButton("btn_download_freq_species_png",
							"📥 PNG — Selected species",
							class = "btn-success btn-block",
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
						downloadButton("btn_download_freq_csv",
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
# TAB: 2. Logistic Regression
# ═══════════════════════════════════════════════════════════════
ui_tab_conserv_regression <- function()

tabItem(tabName = "conserv_regression",
	
	tab_header("chart-line", "Logistic Regression",
		"How conserved are regulatory motifs across species, loci and subgroups?",
		"#7C3AED", "#0B2A5B"),
	
	hr(),
	
	imgt_ligm_note(),
	
	fluidRow(
		box(
			width = 12, title = "📖 About this method — how to read the results",
			status = "info", solidHeader = TRUE,
			collapsible = TRUE, collapsed = FALSE,
			
			HTML("
						<p><strong>What this model estimates:</strong> logistic regression models the
						probability that a given motif is <em>present</em> (coded HAS_MOTIF = 1) as a
						function of three explanatory factors: <b>LOCI</b>, <b>SPECIES</b>, and
						<b>FUNCTIONALITY</b>.</p>
 
						<p><strong>The reference level (baseline) — this is essential to interpret
						the table correctly:</strong> for each categorical factor, R automatically
						picks one level as the <em>reference</em> (by default, the first level in
						alphabetical order) and does not show a coefficient for it. Every other level
						shown in the coefficient table is a comparison <u>against this omitted
						reference</u>, not against the overall average.</p>
 
						<p style=\'background:#FFFBEB;border-left:3px solid #F59E0B;padding:8px 12px;
						border-radius:4px;\'>
						<strong>Example:</strong> if SPECIES has levels \'Homo sapiens\' and
						\'Mus musculus\', and \'Homo sapiens\' is the reference (shown nowhere in the
						table), then the row <code>SPECIESMus musculus</code> means: <em>\"compared to
						Homo sapiens, how does being Mus musculus change the odds of the motif being
						present?\"</em></p>
 
						<p><strong>Odds Ratio (OR = exp(coefficient)):</strong></p>
						<ul>
							<li><b>OR &gt; 1</b> → this level increases the odds of motif presence
							relative to the reference</li>
							<li><b>OR &lt; 1</b> → this level decreases the odds relative to the
							reference</li>
							<li><b>OR ≈ 1</b> → no meaningful difference from the reference</li>
						</ul>
 
						<p><strong>Where to find the reference level used:</strong> it is now shown
						directly in the <em>Reference</em> column of the coefficient table below, and
						also summarized in the box under \'Model summary\'.</p>
						")
		)
	),
	
	# ── PARAMÈTRES ──
	fluidRow(
		box(
			width = 12,
			title = "⚙️ Parameters",
			status = "primary",
			solidHeader = TRUE,
			
			fluidRow(
				column(4,
					div(
						style = "background-color: #FFF9C4; padding: 15px; border-radius: 8px;",
						h5(icon("dna"), strong("Motif selection")),
						selectInput("regr_motif", label = NULL,
							choices = NULL,
							selected = NULL),
						tags$small(icon("info-circle"), 
							" Dependent variable: HAS_[MOTIF]")
					)
				),
				
				column(3, offset = 2,
					div(style = "padding-top: 40px;",
						actionButton("btn_launch_regr", "🚀 Run model",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height: 70px;")
					)
				)
			),
			
			hr(),
			
			div(
				style = "background-color: #E8F5E9; padding: 10px; border-radius: 5px;",
				icon("info-circle"), strong(" Model: "),
				"HAS_[MOTIF] ~ LOCI + SPECIES + FUNCTIONALITY (factors with a single level are dropped automatically)"
			)
		)
	),
	
	# ── RÉSULTATS ──
	fluidRow(
		box(
			width = 12,
			title = "📊 Model Results",
			status = "info",
			solidHeader = TRUE,
			
			uiOutput("regr_status"),
			hr(),
			
			tabsetPanel(
				tabPanel("📋 Coefficients & P-values",
					br(),
					uiOutput("regr_reference_box"),
					br(),
					DTOutput("regr_table"),
					br(),
					div(
						style = "background-color: #FFF3E0; padding: 10px; border-radius: 5px;",
						icon("info-circle"), strong(" Interpretation:"),
						tags$ul(
							tags$li("OR > 1: factor increases the probability of motif presence"),
							tags$li("OR < 1: factor decreases the probability of motif presence"),
							tags$li("p < 0.05: statistically significant effect")
						)
					),
					br(),
					div(style = "margin-top:18px;",
						fluidRow(
							column(12,
								div(
									style = "text-align:center; padding:10px; background:#FFFBEB;
                           border-radius:8px; border:1px solid #FCD34D;",
									icon("file-csv", style="font-size:20px; color:#D97706; display:block; margin-bottom:6px;"),
									p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
										"Regression coefficients & statistics"),
									downloadButton("btn_download_regr_csv",
										"📥 CSV — Regression results",
										class = "btn-warning btn-block",
										style = "width:100%;")
								)
							)
						)
					)
				),
				tabPanel("📊 Model summary",
					br(),
					verbatimTextOutput("regr_summary")
				),
				tabPanel("📈 ROC Curve",
					br(),
					div(style="display:flex; justify-content:center;",
						plotly::plotlyOutput("regr_roc_plot", height="620px", width="620px")								),
					br(),
					div(style="text-align:center; margin-top:10px;",
						downloadButton("btn_dl_regr_roc_png",
							"📥 PNG — ROC Curve",
							class = "btn-primary", style = "margin-right:8px;")
					)
				),
				tabPanel("📊 CV Metrics",
					br(), plotly::plotlyOutput("regr_cv_plot", height="440px"),
					br(),
					div(style="text-align:left; margin-top:10px; display:flex; gap:10px;",
						downloadButton("btn_dl_regr_cv_png",
							"📥 PNG — CV Metrics",
							class = "btn-primary"),
						downloadButton("btn_dl_regr_cv_csv",
							"📥 CSV — CV Metrics (long format)",
							class = "btn-warning")
					),
					hr(),
					h4(icon("table"), " All folds combined — one table"),
					div(
						style = "background-color: #E8F5E9; padding: 10px; border-radius: 6px; margin-bottom: 12px; font-size:12px;",
						icon("info-circle"),
						" Rows = metrics, columns = folds. The ", strong("Mean"),
						" column is the average across all folds."
					),
					DTOutput("regr_cv_wide_table"),
					br(),
					div(style="text-align:left; margin-top:10px;",
						downloadButton("btn_dl_regr_cv_wide_csv",
							"📥 CSV — All folds (wide format)",
							class = "btn-warning")
					)
				),
				tabPanel("🎚️ Calibration",
					br(), plotly::plotlyOutput("regr_calibration", height="440px"),
					br(),
					div(style="text-align:left; margin-top:10px;",
						downloadButton("btn_dl_regr_calib_png",
							"📥 PNG — Calibration Curve",
							class = "btn-primary", style = "margin-right:8px;")
					)
				),
				tabPanel("🔲 Confusion Matrix",
					br(), plotly::plotlyOutput("regr_confusion", height="400px"),
					br(),
					div(style="text-align:left; margin-top:10px; display:flex; gap:10px;",
						downloadButton("btn_dl_regr_conf_png",
							"📥 PNG — Confusion Matrix",
							class = "btn-primary"),
						downloadButton("btn_dl_regr_interp_csv",
							"📥 CSV — Performance Summary",
							class = "btn-success")
					)
				),
				tabPanel("🎯 Interpretation",
					br(), uiOutput("regr_interpretation")
				)
				
			)
		)
	)
)
# ═══════════════════════════════════════════════════════════════
# TAB: 3. Sequence Logos
# ═══════════════════════════════════════════════════════════════
ui_tab_conserv_logos <- function()
tabItem(tabName = "conserv_logos",
	
	tab_header("dna", "Sequence Logos",
		"How conserved are regulatory motifs across species, loci and subgroups?",
		"#0F766E", "#1E3A5F"),
	
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
				# Motif
				column(3,
					div(
						style = "background-color: #FFF9C4; padding: 15px; border-radius: 8px;",
						h5(icon("dna"), strong("Motif")),
						selectInput("logo_motif", label = NULL,
							choices = NULL, selected = NULL),
						tags$small(icon("info-circle"), " Sequence from SEQ_[MOTIF]")
					)
				),
				# Filtre espèce 
				column(3,
					div(
						style = "background-color: #E8F5E9; padding: 15px; border-radius: 8px;",
						h5(icon("globe"), strong("Species (optional)")),
						selectInput("logo_species", label = NULL,
							choices  = c("All species" = ""),
							selected = ""),
						tags$small(
							icon("info-circle"),
							" Empty = all available species"
						)
					)
				),
				# Filtre locus — 
				column(3,
					div(
						style = "background-color: #E3F2FD; padding: 15px; border-radius: 8px;",
						h5(icon("layer-group"), strong("Locus (optional)")),
						selectInput("logo_locus", label = NULL,
							choices  = c("All loci" = ""),
							selected = ""),
						tags$small(
							icon("info-circle"),
							" Empty = IGH + IGK + IGL"
						)
					)
				),
				# Bouton
				column(3,
					div(style = "padding-top: 30px;",
						actionButton("btn_launch_logo", "🚀 Generate logos",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height: 70px;")
					)
				)
			)
		)
	),
	
	# ── RÉSULTATS ──
	fluidRow(
		box(
			width = 12,
			title = "📊 Sequence Logos by Species and Locus",
			status = "info",
			solidHeader = TRUE,
			
			uiOutput("logo_status"),
			hr(),
			
			div(
				style = "background-color: #E3F2FD; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
				icon("info-circle"), strong(" Information : "),
				"Sequence logos represent the frequency of each nucleotide at each position. ",
				"The larger the letter, the more conserved the nucleotide is at that position."
			),
			
			tabsetPanel(
				tabPanel("📈 Logos",
					br(),
					uiOutput("logo_plot_ui"),   
					br(),
					div(style = "margin-top:18px;",
						fluidRow(
							column(12,
								div(
									style = "text-align:center; padding:10px; background:#EFF6FF;
                       border-radius:8px; border:1px solid #93C5FD;",
									icon("image", style="font-size:20px; color:#1D4ED8; display:block; margin-bottom:6px;"),
									p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
										"Sequence logo figure (PNG, 300 dpi)"),
									downloadButton("btn_download_logo_plot",
										"📥 PNG — Sequence Logos",
										class = "btn-primary btn-block",
										style = "width:100%;")
								)
							)
						)
					)
				),
				tabPanel("📋 Sequences",
					br(),
					p("Sequences used to generate the logos:"),
					DTOutput("logo_sequences_table"),
					br(),
					div(style = "margin-top:18px;",
						fluidRow(
							column(12,
								div(
									style = "text-align:center; padding:10px; background:#FFFBEB;
                               border-radius:8px; border:1px solid #FCD34D;",
									icon("file-csv", style="font-size:20px; color:#D97706; display:block; margin-bottom:6px;"),
									p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
										"Sequences used for logo generation"),
									downloadButton("btn_download_logo_csv",
										"📥 CSV — Sequences",
										class = "btn-warning btn-block",
										style = "width:100%;")
								)
							)
						)
					)
				)
			)
		)
	)
)
# ═══════════════════════════════════════════════════════════════
# TAB: 4. Conservation by Functionality
# ═══════════════════════════════════════════════════════════════
ui_tab_conserv_functionality <- function()
tabItem(tabName = "conserv_functionality",
	
	tab_header("balance-scale", "Conservation by Functionality",
		"How conserved are regulatory motifs across species, loci and subgroups?",
		"#854D0E", "#1E3A5F"),
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
				strong(" Analysis restricted to Immunoglobulins (IGH, IGK, IGL)")
			),
			
			# Type d'analyse + Motif + Groupement 
			fluidRow(
				column(4,
					div(style = "background-color: #EDE9FE; padding: 15px; border-radius: 8px;",
						h5(icon("sliders-h"), strong("Analysis type")),
						radioButtons("cons_func_analysis_type", label = NULL,
							choices = c(
								"Global (all motifs)" = "global",
								"Single motif"        = "single"
							),
							selected = "global"
						)
					)
				),
				column(4,
					div(style = "background-color: #FFF9C4; padding: 15px; border-radius: 8px;",
						h5(icon("dna"), strong("Motif (if Single motif)")),
						selectInput("cons_func_motif", label = NULL,
							choices = NULL, selected = NULL),
						conditionalPanel(
							condition = "input.cons_func_analysis_type == 'global'",
							div(
								style = "background:#F3F4F6; border-radius:4px; padding:6px 10px;
                         font-size:11px; color:#9CA3AF; margin-top:4px;",
								icon("lock"), " Disabled in Global mode"
							)
						)
					)
				),
				column(4,
					div(style = "background-color: #FFF9C4; padding: 15px; border-radius: 8px;",
						h5(icon("layer-group"), strong("Grouping level")),
						radioButtons("cons_func_grouping", label = NULL,
							choices = c(
								"By species"  = "species",
								"By locus"    = "locus",
								"By subgroup" = "subgroup"
							),
							selected = "species"
						)
					)
				)
			),
			
			hr(),
			
			
			fluidRow(
				
				# Sélecteur espèce 
				column(3,
					conditionalPanel(
						condition = "input.cons_func_grouping == 'locus' ||
                         input.cons_func_grouping == 'subgroup'",
						div(
							style = "background-color: #DCFCE7; padding: 12px; border-radius: 8px;",
							div(
								style = "display:flex; align-items:center; gap:12px;",
								tags$div(icon("globe"), style = "font-size:20px; color:#166534;"),
								div(style = "flex:1;",
									tags$label(strong("Species"),
										style = "color:#166534; margin-bottom:4px; display:block;"),
									selectInput("cons_func_species", label = NULL,
										choices = NULL, selected = NULL, width = "100%")
								)
							),
							tags$small(
								style = "color:#555; font-size:11px; margin-top:6px; display:block;",
								icon("info-circle"),
								" Locus/Subgroup: one species at a time."
							)
						)
					),
					conditionalPanel(
						condition = "input.cons_func_grouping == 'species'",
						div(
							style = "background:#F0FDF4; padding:12px; border-radius:8px;
                       border:1px dashed #86EFAC;",
							icon("info-circle"),
							tags$span(style = "color:#166634; font-size:12px;",
								" Species mode: all IG species analyzed together.")
						)
					)
				),
				
				# Sélecteur locus 
				column(3,
					conditionalPanel(
						condition = "input.cons_func_grouping == 'subgroup'",
						div(
							style = "background-color: #EFF6FF; padding: 12px; border-radius: 8px;
                       border-left: 4px solid #3B82F6;",
							div(
								style = "display:flex; align-items:center; gap:10px; margin-bottom:8px;",
								tags$div(icon("layer-group"), style = "font-size:18px; color:#1E40AF;"),
								div(style = "flex:1;",
									tags$label(strong("Locus (optional)"),
										style = "color:#1E40AF; margin-bottom:4px; display:block;"),
									selectInput("cons_func_locus_filter", label = NULL,
										choices = c("All loci" = "All loci", "IGHV" = "IGHV", "IGKV" = "IGKV", "IGLV" = "IGLV"),
										selected = "All loci", width = "100%")
								)
							),
							tags$label(strong("Display mode"),
								style = "color:#1E40AF; margin-bottom:4px; display:block;"),
							radioButtons("cons_func_sg_mode", label = NULL,
								choices = c("Auto" = "auto", "Force barplot" = "barplot",
									"Force heatmap" = "heatmap"),
								selected = "auto", inline = TRUE
							),
							tags$label(strong("Top-N subgroups"),
								style = "color:#1E40AF; margin-top:8px; margin-bottom:4px; display:block;"),
							fluidRow(
								column(7,
									numericInput("cons_func_topn", label = NULL,
										value = NA, min = 5, max = 500, step = 5, width = "100%")
								),
								column(5,
									tags$small(style = "color:#6B7280; font-size:11px;
                               padding-top:8px; display:block;",
										"Empty = all")
								)
							),
							tags$small(
								style = "color:#4B5563; font-size:11px; margin-top:4px; display:block;",
								icon("lightbulb"), strong(" Auto:"),
								" switches to heatmap above 80 subgroups.", br(),
								icon("filter"), strong(" Locus:"),
								" restricts barplot to one locus."
							)
						)
					)
				),
				
				# Bouton lancement
				column(4, offset = 2,
					div(style = "padding-top: 10px;",
						actionButton("btn_launch_cons_func", "🚀 Launch analysis",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height: 70px; font-size:16px;")
					)
				)
			)
		)
	),
	
	# Résultats 
	fluidRow(
		box(
			width = 12, title = "📊 Results",
			status = "info", solidHeader = TRUE,
			
			uiOutput("cons_func_status"),
			hr(),
			
			div(
				style = "background-color: #FFF3E0; padding: 10px; border-radius: 5px;
                 margin-bottom: 15px;",
				icon("palette"), strong(" Legend:"),
				conditionalPanel(
					condition = "input.cons_func_analysis_type == 'single'",
					tagList(
						tags$span(style = "color:#28A745; font-weight:bold;", " ■ F = Functional"),
						tags$span(style = "color:#DC3545; font-weight:bold;", "  ■ P = Pseudogene"),
						tags$span(style = "color:#FFC107; font-weight:bold;", "  ■ ORF"),
						tags$em(style = "color:#6B7280; font-size:11px;",
							"  — bars colored by functionality")
					)
				),
				conditionalPanel(
					condition = "input.cons_func_analysis_type == 'global'",
					tags$span(
						style = "color:#374151; font-size:12px;",
						" Each color = one  motif. Panels separate F / P / ORF."
					)
				)
			),
			
			tabsetPanel(
				tabPanel("📈 Graphique",
					br(),
					uiOutput("cons_func_plot_ui")
				),
				tabPanel("📊 Chi² Tests",
					br(),
					div(
						style = "background-color: #E8F5E9; padding: 10px; border-radius: 5px;
                     margin-bottom: 15px;",
						icon("info-circle"), strong(" Statistical tests:"),
						p("Chi²: is motif presence independent of functionality?"),
						tags$ul(tags$li("p < 0.05: significant dependence (***)"))
					),
					DTOutput("cons_func_chi_table")
				),
				tabPanel(
					tagList(icon("flask"), " Advanced Statistical Tests"),
					br(),
					div(
						style = "background:#DBEAFE; padding:12px; border-radius:6px;
                     border-left:4px solid #1565C0; margin-bottom:14px; font-size:12px;",
						icon("info-circle"),
						strong(" Methods applied:"),
						tags$ul(
							tags$li(
								strong("Kruskal-Wallis (global): "),
								"non-parametric test of independence between motif presence and functionality group.",
								" H₀: the distribution of presence is identical across all groups (F/P/ORF)."
							),
							tags$li(
								strong("Mann-Whitney (pairwise): "),
								"applied for each significant Kruskal-Wallis (p.adjusted < 0.05).",
								" Compares each pair of functionality groups (F vs P, F vs ORF, P vs ORF)."
							),
							tags$li(
								strong("Bonferroni correction: "),
								"p-values adjusted by the number of motifs × number of pairs tested."
							),
							tags$li(
								strong("Effect size (r = |Z|/√N): "),
								"< 0.1 = negligible | 0.1–0.3 = small | 0.3–0.5 = medium | > 0.5 = large."
							)
						)
					),
					DT::DTOutput("cons_func_stat_tests_table"),
					br(),
					div(
						style = "text-align:left; padding:10px; background:#FFFBEB;
                     border-radius:8px; border:1px solid #FCD34D; display:inline-block;",
						icon("file-csv", style="font-size:18px; color:#D97706; margin-right:6px;"),
						downloadButton("btn_download_cons_stat_csv",
							"📥 CSV — Statistical Tests",
							class = "btn-warning", style = "margin-left:8px;")
					),
					br()
				)
			)
		)
	),
	
	# Boutons  
	fluidRow(
		box(
			width = 12, title = "⬇️ Downloads",
			status = "success", solidHeader = FALSE,
			collapsible = TRUE,
			fluidRow(
				column(6,
					div(
						style = "text-align:center; padding:10px; background:#EFF6FF;
                   border-radius:8px; border:1px solid #93C5FD;",
						icon("image", style="font-size:20px; color:#1D4ED8; display:block; margin-bottom:6px;"),
						p(style="font-size:11px; color:#374151; margin:0 0 6px 0;",
							"Full figure (all subgroups/species)"),
						downloadButton("btn_download_cons_func_plot",
							"📥 High-resolution PNG",
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
							"Raw data table"),
						downloadButton("btn_download_cons_func_csv",
							"📥 CSV",
							class = "btn-warning btn-block",
							style = "width:100%;")
					)
				)
			)
		)
	)
)