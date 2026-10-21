# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  tabItem (partial genes)  ████████████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════
# ════════════════════════════════════════════════════════
# PARTIAL GENES
# ════════════════════════════════════════════════════════
ui_tab_partial_desc <- function()
tabItem(tabName = "partial_desc",
	
	tab_header("cut", "Partial Genes — Motif Presence Analysis",
		"Which sequences carry no regulatory motif, and how do motif count/position relate to 5\u2032UTR length?",
		"#7F1D1D", "#B45309", "#1E3A5F"),
	
	uiOutput("warn_no_partiels"),
	
	fluidRow(
		valueBoxOutput("vbox_partial_total",       width = 4),
		valueBoxOutput("vbox_partial_with_motif",  width = 4),
		valueBoxOutput("vbox_partial_species",     width = 4)
	),
	
	# Filtres en cascade partagés par les deux analyses
	fluidRow(
		box(
			width = 12, title = "⚙️ Filters (cascading)",
			status = "primary", solidHeader = FALSE,
			fluidRow(
				column(3,
					div(style = "background:#EFF6FF; padding:12px; border-radius:8px;",
						h5(icon("globe"), strong("Species")),
						selectInput("partial_species_filter", label = NULL,
							choices = c("All species" = "All species"), selected = "All species")
					)
				),
				column(3,
					div(style = "background:#F0FDF4; padding:12px; border-radius:8px;",
						h5(icon("layer-group"), strong("Locus")),
						selectInput("partial_locus_filter", label = NULL,
							choices = c("All loci" = "All loci"), selected = "All loci"),
						tags$small(icon("info-circle"), " Updates based on Species.")
					)
				),
				column(3,
					div(style = "background:#FFF3E0; padding:12px; border-radius:8px;",
						h5(icon("sitemap"), strong("Subgroup")),
						selectInput("partial_subgroup_filter", label = NULL,
							choices = c("All subgroups" = "All subgroups"), selected = "All subgroups"),
						tags$small(icon("info-circle"), " Updates based on Species + Locus.")
					)
				),
				column(3,
					div(style = "padding-top: 8px;",
						actionButton("btn_launch_partial", "🚀 Launch analysis",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height: 90px;")
					)
				)
			),
			div(
				style = "background:#FFFBEB; padding:10px 14px; border-radius:6px;
                 border-left:4px solid #D97706; margin-top:10px; font-size:12px;",
				icon("info-circle", style = "color:#D97706;"),
				strong(" Note: "),
				"LPART1 is not counted as a regulatory motif. A sequence is 'without motif' if none of ",
				"TATA_BOX, HEPTA_PENTA, OCTA_DECA, E_BOX, TATA_BOX2, EBOX_2 are detected."
			)
		)
	),
	
	uiOutput("partial_desc_status"),
	
	# ── Tableau A : séquences SANS motif ─────────────────────────
	fluidRow(
		box(
			width = 12, title = "🚫 Sequences WITHOUT any regulatory motif",
			status = "danger", solidHeader = TRUE,
			tabsetPanel(
				tabPanel("📋 Summary by species",
					br(), DTOutput("partial_summary_nomotif_table"), br(),
					downloadButton("btn_dl_partial_summary_nomotif_csv", "📥 CSV — Summary", class = "btn-warning")
				),
				tabPanel("🔍 Allele list",
					br(), DTOutput("partial_alleles_nomotif_table"), br(),
					downloadButton("btn_dl_partial_alleles_nomotif_csv", "📥 CSV — Allele list", class = "btn-warning")
				),
				tabPanel("📈 Plots",
					br(),
					h4("Sequences without motif — by species and functionality"),
					uiOutput("partial_plot_nomotif_species_ui"), br(),
					downloadButton("btn_dl_partial_plot_nomotif_species_png", "📥 PNG", class = "btn-primary"),
					hr(),
					h4("5'UTR length distribution — sequences without motif"),
					uiOutput("partial_plot_nomotif_length_ui"), br(),
					downloadButton("btn_dl_partial_plot_nomotif_length_png", "📥 PNG", class = "btn-primary")
				)
			)
		)
	),
	
	# ── Tableau B : séquences AVEC motif ─────────────────────────
	fluidRow(
		box(
			width = 12, title = "✅ Sequences WITH at least one regulatory motif",
			status = "success", solidHeader = TRUE,
			div(
				style = "background:#F0FDF4; padding:10px; border-radius:6px;
                 border-left:4px solid #16A34A; margin-bottom:12px; font-size:12px;",
				icon("info-circle"),
				strong(" Columns: "), "N_MOTIFS = number of regulatory motifs detected | ",
				"DIST_xxx = distance (nt) between that motif and LPART1 | LENGTH_NT = 5'UTR length."
			),
			tabsetPanel(
				tabPanel("📋 Summary by species",
					br(), DTOutput("partial_summary_withmotif_table"), br(),
					downloadButton("btn_dl_partial_summary_withmotif_csv", "📥 CSV — Summary", class = "btn-warning")
				),
				tabPanel("🔍 Allele list",
					br(), DTOutput("partial_alleles_withmotif_table"), br(),
					downloadButton("btn_dl_partial_alleles_withmotif_csv", "📥 CSV — Allele list", class = "btn-warning")
				),
				tabPanel("📈 Plots",
					br(),
					h4("Number of motifs per allele"),
					uiOutput("partial_plot_withmotif_count_ui"), br(),
					downloadButton("btn_dl_partial_plot_withmotif_count_png", "📥 PNG", class = "btn-primary"),
					hr(),
					h4("Distance to LPART1, by motif"),
					uiOutput("partial_plot_withmotif_dist_ui"), br(),
					downloadButton("btn_dl_partial_plot_withmotif_dist_png", "📥 PNG", class = "btn-primary")
				)
			)
		)
	)
)

# ══════════════════════════════════════════════════════════════
# TAB B — Comparison avec séquences de la littérature
# ══════════════════════════════════════════════════════════════
ui_tab_partial_compare<- function()
tabItem(tabName = "partial_compare",
	
	tab_header("microscope", "Literature Comparison — 5\u2032UTR of Partial Genes",
		"How similar are the excluded partial sequences to homologous sequences reported in the literature?",
		"#1E3A5F", "#065F46", "#7F1D1D"),
	
	# ── Import FASTA littérature ──────────────────────────────────
	fluidRow(
		box(
			width = 12,
			title = tagList(
				tags$span(
					style = "background:#1E3A5F; color:white; border-radius:50%; width:26px; height:26px;
                   display:inline-flex; align-items:center; justify-content:center;
                   font-weight:800; font-size:13px; margin-right:8px;",
					"1"),
				"Import Literature FASTA Sequences"
			),
			status = "primary", solidHeader = TRUE,
			
			div(
				style = "background:#EFF6FF; padding:12px; border-radius:6px;
                 border-left:4px solid #3B82F6; margin-bottom:14px; font-size:12px;",
				icon("info-circle"),
				strong(" Expected FASTA header format (IMGT):"),
				br(),
				tags$code(style = "background:#F1F5F9; padding:4px 8px; border-radius:4px;",
					">AB019437.1|IGHV7-81*01|Homo sapiens|ORF|5'UTR|3501..6314|2814 nt||..."),
				br(), br(),
				"You can import multiple alleles in a single file. ",
				"Each allele can have several literature sequences — they will be grouped automatically."
			),
			
			tabsetPanel(
				tabPanel("📋 Copy-paste FASTA",
					br(),
					textAreaInput("lit_fasta_text",
						label = NULL,
						placeholder = paste0(
							">AB019437.1|IGHV7-81*01|Homo sapiens|ORF|5'UTR|3501..6314|2814 nt||||||\n",
							"ATGAAACACCTGTGGTTCTTCCTCCTGCTGGTGGCAGCTCCCAGATGGGTCCTGTCCCAG...\n",
							">AC245166.1|IGHV7-81*01|Homo sapiens|F|5'UTR|1200..3100|1900 nt||||||\n",
							"ATGAAACACCTGTGGTTCTTCCTCCTGCTGGTGGCAGCTCC..."
						),
						rows = 10, width = "100%"),
					actionButton("btn_parse_lit_text", "🚀 Launch analysis",
						icon = icon("play"), class = "btn-success btn-lg",
						width = "100%", style = "height: 60px;")
				),
				tabPanel("📁 Upload file",
					br(),
					fileInput("lit_fasta_file",
						label = "Select FASTA file (.fasta, .fa, .txt)",
						accept = c(".fasta", ".fa", ".txt"),
						buttonLabel = "Browse...",
						placeholder = "No file selected"),
					actionButton("btn_parse_lit_file", "🚀 Launch analysis",
						icon = icon("play"), class = "btn-success btn-lg",
						width = "100%", style = "height: 60px;")
				)
			),
			
			br(),
			uiOutput("lit_import_status"),
			
			# Bouton export CSV de l'import
			conditionalPanel(
				condition = "output.lit_import_status",
				div(style = "margin-top:8px;",
					downloadButton("btn_dl_lit_parsed_csv",
						"📥 CSV — Parsed literature sequences", class = "btn-warning btn-sm")
				)
			)
		)
	),
	
	# ──  Sélection allèle and paramètres ────────────────────────────
	fluidRow(
		box(
			width = 12,
			title = tagList(
				tags$span(
					style = "background:#065F46; color:white; border-radius:50%; width:26px; height:26px;
                   display:inline-flex; align-items:center; justify-content:center;
                   font-weight:800; font-size:13px; margin-right:8px;",
					"2"),
				"Select Allele & Run Comparison"
			),
			status = "success", solidHeader = TRUE,
			
			fluidRow(
				# Sélecteur allèle
				column(5,
					div(style = "background:#F0FDF4; padding:14px; border-radius:8px;",
						h5(icon("dna"), strong("Allele")),
						selectInput("lit_allele_select", label = NULL,
							choices  = NULL, selected = NULL, width = "100%"),
						div(
							style = "font-size:11px; color:#555; margin-top:6px;",
							icon("star", style = "color:#D97706;"),
							" ★ = allele also present in removed sequences"
						)
					)
				),
				# Paramètre : nombre de nt terminaux
				column(3,
					div(style = "background:#FFF9C4; padding:14px; border-radius:8px;",
						h5(icon("ruler"), strong("Last N nucleotides")),
						numericInput("lit_n_last",
							label = NULL,
							value = 300, min = 50, max = 1000, step = 50,
							width = "100%"),
						tags$small(icon("info-circle"),
							" Only the last N nt of each sequence will be compared")
					)
				),
				# Bouton
				column(4,
					div(style = "padding-top:20px;",
						actionButton("btn_run_comparison", "🚀 Run comparison",
							icon = icon("play"), class = "btn-success btn-lg",
							width = "100%", style = "height:65px; font-size:15px;")
					)
				)
			),
			fluidRow(
				column(12,
					uiOutput("lit_removed_allele_info")
				)
			)
		)
	),
	
	# ──  Résultats ───────────────────────────────────────────────
	fluidRow(
		box(
			width = 12,
			title = tagList(
				tags$span(
					style = "background:#7F1D1D; color:white; border-radius:50%; width:26px; height:26px;
                   display:inline-flex; align-items:center; justify-content:center;
                   font-weight:800; font-size:13px; margin-right:8px;",
					"3"),
				"Results"
			),
			status = "info", solidHeader = TRUE,
			
			uiOutput("compare_status"),
			hr(),
			
			tabsetPanel(
				
				# Tab 1 : Heatmap identité
				tabPanel("🔥 Identity Heatmap",
					br(),
					div(
						style = "background:#EFF6FF; padding:10px; border-radius:6px;
                     border-left:4px solid #1565C0; margin-bottom:12px; font-size:12px;",
						icon("info-circle"),
						strong(" Reading guide:"),
						" Dark red = high identity (> 95%) | Blue = low identity (< 50%). ",
						"Sequences are ordered by hierarchical clustering for visual grouping."
					),
					uiOutput("compare_heatmap_ui"),
					br(),
					div(style = "text-align:left;",
						downloadButton("btn_dl_compare_heatmap_png",
							"📥 PNG — Identity heatmap", class = "btn-primary")
					)
				),
				
				# Tab 2 : Distribution identités
				tabPanel("📊 Identity distribution",
					br(),
					div(
						style = "background:#F0FDF4; padding:10px; border-radius:6px;
                     border-left:4px solid #16A34A; margin-bottom:12px; font-size:12px;",
						icon("info-circle"),
						" Per-sequence identity profile: each dot is a sequence's mean identity against all others, ",
						"the line shows its min–max range. Dashed navy line = overall mean identity. ",
						tags$strong("✕ markers"), " flag sequences compared on a short fragment (unreliable identity)."
					),
					fluidRow(
						column(7,
							uiOutput("compare_dist_plot_ui"),
							br(),
							div(style = "text-align:left;",
								downloadButton("btn_dl_compare_dist_png",
									"📥 PNG — Identity profile", class = "btn-primary")
							)
						),
						column(5,
							br(),
							h5(strong("Summary statistics"), style = "color:#374151;"),
							DTOutput("compare_summary_table"),
							br(),
							div(
								style = "background:#FFF9C4; padding:10px; border-radius:6px; font-size:12px;",
								icon("lightbulb"),
								strong(" Interpretation guide:"),
								tags$ul(style = "margin:6px 0 0 0;",
									tags$li("> 95% — Near-identical sequences"),
									tags$li("85–95% — High conservation"),
									tags$li("70–85% — Moderate conservation"),
									tags$li("< 70% — Divergent sequences")
								)
							)
						)
					)
				),
				
				# Tab 3 : Table des séquences
				tabPanel("📋 Sequence table",
					br(),
					div(
						style = "background:#FEF3C7; padding:10px; border-radius:6px;
                     border-left:4px solid #D97706; margin-bottom:12px; font-size:12px;",
						icon("info-circle"),
						" The alignment column (last N nt) is hidden in the view but ",
						strong("included in CSV/Excel export"), "."
					),
					DTOutput("compare_alignment_table"),
					br(),
					div(style = "text-align:left;",
						downloadButton("btn_dl_compare_csv",
							"📥 CSV — Sequences + alignment", class = "btn-warning")
					)
				)
			)
		)
	)
)

