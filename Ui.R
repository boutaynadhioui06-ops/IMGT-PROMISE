# ═══════════════════════════════════════════════════════════════
# IMGT/PROMISE - Interface Utilisateur
# ═══════════════════════════════════════════════════════════════
ui <- dashboardPage(
	skin = "blue",
	
	# ══ HEADER ══════════════════════════════════════════════════
	dashboardHeader(
		title = tags$div(
			style = "display: flex; align-items: center; padding: 5px;",
			tags$div(
				style = "font-size: 20px; margin-right: 8px;",
				"🦠️🔭 "
			),
			tags$div(
				style = "font-family: 'Arial', sans-serif; font-weight: bold; font-size: 18px;",
				"IMGT-PROMISE"
			)
		),
		titleWidth = 300
	),
	
	# ══ SIDEBAR ═════════════════════════════════════════════════
	dashboardSidebar(
		width = 280,
		
		# Sous-titre
		div(style = "padding:12px 18px; 
		             background:#0D2137; 
		             color:#90CAF9;
                 font-size:15px; 
                 font-weight:600;
                 text-align:center;
                 line-height:1.5;
			           border-bottom:1px solid #E3F2FD;",
			"Statistical analysis of IG and TR", br(),
			"V-gene regulatory regions", br(), br(),
			
		),
		
		
		sidebarMenu( id = "tabs",
			menuItem("🏠 Home", tabName = "home", icon = icon("home")),
			
			menuItem("📊 Dataset Construction", 
				icon = icon("database"),
				menuSubItem("📂 Import 5'UTR", tabName = "import_utr", icon = icon("file-upload")),
				menuSubItem("🧬 Import Motifs", tabName = "import_motifs", icon = icon("dna")),
				menuSubItem("📋 Final Dataset", tabName = "dataset", icon = icon("table"))
			),
			
			menuItem("🧩 Partial Genes",
				icon = icon("cut"),
				menuSubItem("Descriptive Statistics",
					tabName = "partial_desc",
					icon    = icon("table")),
				menuSubItem("Literature Comparison",
					tabName = "partial_compare",
					icon    = icon("align-left"))
			),
			
			menuItem("📈 Sequences & Functionality",
				icon = icon("chart-bar"),
				menuSubItem("Descriptive Table", 
					tabName = "stat_desc_tableau", 
					icon = icon("table")),
				menuSubItem("Sequence Distribution", 
					tabName = "stat_desc_sequences", 
					icon = icon("chart-bar")),
				menuSubItem("Functionality Distribution", 
					tabName = "stat_desc_functionality", 
					icon = icon("dna")),
				menuSubItem("Alluvial Flow",
					tabName = "adv_alluvial",
					icon = icon("stream"))
			),
			
			menuItem("🔬 Motif Conservation",
				icon = icon("microscope"),
				menuSubItem("Presence Frequencies", 
					tabName = "conserv_frequencies", 
					icon = icon("percent")),
				menuSubItem("Sequence Logos", 
					tabName = "conserv_logos", 
					icon = icon("dna")),
				menuSubItem("Conservation / Functionality", 
					tabName = "conserv_functionality", 
					icon = icon("balance-scale")),
				menuSubItem("Logistic Regression", 
					tabName = "conserv_regression", 
					icon = icon("chart-line"))
			),
			
			menuItem("📍 Motif Positioning",
				icon = icon("map-marked-alt"),
				menuSubItem("Genomic Map",
					tabName = "adv_gmap",
					icon = icon("map-marked-alt")),
				menuSubItem("Distance Distributions",
					tabName = "adv_ridgeline",
					icon = icon("mountain"))
			),
			menuItem("🌳 Motif Similarity",
				icon = icon("sitemap"),
				menuSubItem("Hierarchical Clustering", 
					tabName = "conserv_clustering", 
					icon = icon("project-diagram")),
				menuSubItem("Sequence Identity", 
					tabName = "conserv_identity", 
					icon = icon("th")),
				menuSubItem("Motif Phylogeny",
					tabName = "adv_tsne",
					icon = icon("braille"))
			),
			
			menuItem("🧪 Motif Diversity",
				icon = icon("flask"),
				menuSubItem("Motif Diversity",
					tabName = "diversity",
					icon = icon("cubes")),
				menuSubItem("PERMANOVA & PCoA",
					tabName = "adv_permanova",
					icon = icon("vials"))
			),
			
			menuItem("🔗 Motif Co-occurrence",
				icon = icon("share-alt"),
				menuSubItem("Co-occurrence & Network",
					tabName = "adv_cooccurrence",
					icon = icon("project-diagram"))
			),
			menuItem(" 🔎 Motif Discovery", icon = icon("dna"),
				menuSubItem("🧩 Consensus-based",         tabName = "scanner_consensus",  icon = icon("puzzle-piece")),
				menuSubItem("🧬 PWM-based probabilistic", tabName = "scanner_pwm",        icon = icon("chart-line")),
				menuSubItem("📈 Comparison",              tabName = "scanner_comparison", icon = icon("balance-scale")),
				menuSubItem("✅ Final dataset", tabName = "scanner_final", icon = icon("check-double"))
			),
			menuItem("📖 Resources", icon = icon("book"),
				menuSubItem("📘 Documentation",   tabName = "resources_documentation", icon = icon("book-open")),
				menuSubItem("📄 Analysis Report", tabName = "resources_report",        icon = icon("file-alt")),
				menuSubItem("💻 Source Code",     tabName = "resources_source",        icon = icon("code")),
				menuSubItem("📚 References",      tabName = "resources_references",    icon = icon("book")),
				menuSubItem("🕘 Changelog",       tabName = "resources_changelog",     icon = icon("clock-rotate-left"))
			),
			menuItem("ℹ️ About", tabName = "about_page", icon = icon("circle-info"))
			
		
		)
	),
	
	# ══ BODY ════════════════════════════════════════════════════
	dashboardBody(
		tags$head(
			tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
			tags$link(rel = "stylesheet", type = "text/css", href = "resources.css")
		),
		tabItems(
		ui_tab_home(),
		ui_tab_import_utr(),
		ui_tab_import_motifs(),
		ui_tab_dataset(),
		ui_tab_partial_desc(),
		ui_tab_partial_compare(),
		ui_tab_stat_desc_tableau(),
		ui_tab_stat_desc_sequences(),
		ui_tab_stat_desc_functionality(),
		ui_tab_conserv_frequencies(),
		ui_tab_conserv_regression(),
		ui_tab_conserv_logos(),
		ui_tab_conserv_functionality(),
		ui_tab_conserv_clustering(),
		ui_tab_diversity(),
		ui_tab_conserv_identity(),
		ui_tab_adv_cooccurrence(),
		ui_tab_adv_alluvial(),
		ui_tab_adv_ridgeline(),
		ui_tab_adv_permanova(),
		ui_tab_adv_gmap(),
		ui_tab_adv_tsne(),
		ui_tab_scanner_consensus(),
		ui_tab_scanner_pwm(),
		ui_tab_scanner_comparison(),
		ui_tab_scanner_final(),
		ui_tab_resources_documentation(),
		ui_tab_resources_report(),
		ui_tab_resources_source(),
		ui_tab_resources_references(),
		ui_tab_resources_changelog(),
		ui_tab_about_page()
		)	
	) # end dashboardBody
)