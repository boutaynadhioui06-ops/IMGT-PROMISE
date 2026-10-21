# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  tabItem (datates)  ████████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════
# ════════════════════════════════════════════════════════
# ONGLET IMPORT 5'UTR
# ════════════════════════════════════════════════════════
ui_tab_import_utr <- function()
tabItem(tabName = "import_utr",
	tab_header("file-code", "5′UTR File Import",
		"Which 5\u2032UTR sequences and headers are available for downstream analysis?",
		"#0D2137", "#1565C0"),
	hr(),
	mod_import_ui("import_module")
)

# ════════════════════════════════════════════════════════
# ONGLET IMPORT MOTIFS
# ════════════════════════════════════════════════════════
ui_tab_import_motifs  <- function()
tabItem(tabName = "import_motifs",
	tab_header("puzzle-piece", "Motif File Import",
		"Which regulatory motifs have been detected and annotated across the imported alleles?",
		"#1565C0", "#06A77D"),
	hr(),
	
	# Avertissement si 5'UTR pas encore importé
	uiOutput("warn_no_utr"),
	mod_motifs_ui("motifs_module")
)

# ════════════════════════════════════════════════════════
# ONGLET DATASET FINAL
# ════════════════════════════════════════════════════════
ui_tab_dataset <- function()
tabItem(tabName = "dataset",
	tab_header("table", "Final Dataset",
		"What does the final integrated IG/TR dataset look like once 5\u2032UTR and motif data are joined?",
		"#001f54", "#2E75B6"),
	
	# Value boxes
	fluidRow(
		valueBoxOutput("vbox_total",    width = 3),
		valueBoxOutput("vbox_loci",     width = 3),
		valueBoxOutput("vbox_species",  width = 3),
		valueBoxOutput("vbox_func",     width = 3)
	),
	
	# Bouton construire
	fluidRow(
		box(width = 12, status = "primary", solidHeader = FALSE,
			fluidRow(
				column(6,
					uiOutput("dataset_build_status")
				),
				column(6,
					actionButton("btn_build_dataset", "🔨 Build Dataset",
						icon  = icon("cogs"),
						class = "btn-primary btn-lg",
						width = "100%")
				)
			)
		)
	),
	
	# Tableau interactif — IG et TR séparés
	fluidRow(
		box(width = 12, title = "Interactive table",
			status = "info", solidHeader = TRUE,
			p(icon("lightbulb"), "Filter, sort, and search within the table"),
			div(class = "alert alert-info", style = "font-size:12px; margin-bottom:10px;",
				icon("info-circle"),
				strong(" Note: "),
				"The ", strong("Copy / CSV / Excel"), " buttons above the table only export ",
				"the current view (rows currently visible after your filters/search), ",
				"not the entire dataset. To export the full dataset, use the ",
				strong("⬇️ Export"), " section below (Excel / CSV)."
			),
			tabsetPanel(
				tabPanel("🧬 IG (IGHV / IGKV / IGLV)",
					br(), DTOutput("main_table_ig")
				),
				tabPanel("🧬 TR (TRAV / TRBV / TRDV / TRGV)",
					br(),
					div(class = "alert alert-warning", style = "font-size:12px;",
						icon("exclamation-triangle"), strong(" Note: "),
						"Regulatory motif columns (SEQ_/POS_START_/POS_END_/BIN_/DIST_) are hidden for TR, ",
						"since these motifs are not annotated for TR V-genes."
					),
					DTOutput("main_table_tr")
				)
			)
		)
	),
	
	# ── Export — IG et TR séparés ─────────────────────────────
	fluidRow(
		box(width = 12, title = "⬇️ Export", status = "success",
			solidHeader = FALSE,
			div(class = "alert alert-info", style = "margin-bottom:14px;",
				icon("info-circle"),
				" Export the full final dataset — separately for IG and TR — in Excel or CSV format."
			),
			fluidRow(
				column(6,
					div(style = "text-align:center; padding:18px; background:#F0FBF7;
                           border-radius:10px; border:2px dashed #06A77D;",
						icon("dna", class = "fa-2x", style = "color:#06A77D;"),
						h5(strong("IG dataset")),
						fluidRow(
							column(6, downloadButton("btn_dl_excel_ig", "⬇ Excel",
								class = "btn-success", style = "width:100%;")),
							column(6, downloadButton("btn_dl_csv_ig", "⬇ CSV",
								class = "btn-info", style = "width:100%;"))
						)
					)
				),
				column(6,
					div(style = "text-align:center; padding:18px; background:#EBF5FF;
                           border-radius:10px; border:2px dashed #2E75B6;",
						icon("dna", class = "fa-2x", style = "color:#2E75B6;"),
						h5(strong("TR dataset")),
						fluidRow(
							column(6, downloadButton("btn_dl_excel_tr", "⬇ Excel",
								class = "btn-success", style = "width:100%;")),
							column(6, downloadButton("btn_dl_csv_tr", "⬇ CSV",
								class = "btn-info", style = "width:100%;"))
						)
					)
				)
			)
		)
	)
)
