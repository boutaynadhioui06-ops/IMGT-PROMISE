# ═════════════════════════════════════════════════════════════════════════════════════════════════
# ██████████████████  tabItem (cooccurence alluvial, distance...)  ████████████████████████████████
# ═════════════════════════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════
# Advanced Analyses 
# ══════════════════════════════════════════════════════════════

# Tab A — Co-occurrence & Network
ui_tab_adv_cooccurrence <- function()
tabItem(tabName = "adv_cooccurrence",
	tab_header("share-alt", "Co-occurrence & Network",
		"Which motifs tend to appear together, and which are mutually exclusive?",
		"#0B2A5B", "#1B7F3B", "#F4C300"),
	uiOutput("warn_no_dataset"),
	hr(),
	imgt_ligm_note(),
	mod_cooccurrence_ui("adv_cooc_mod")
)

# Tab B — Alluvial Flow
ui_tab_adv_alluvial <- function()
tabItem(tabName = "adv_alluvial",
	tab_header("water", "Alluvial Flow",
		"How are sequences and gene functionality distributed across species, loci and subgroups?",
		"#1B7F3B", "#0B2A5B"),
	uiOutput("warn_no_dataset_alluvial"),
	hr(),
	imgt_ligm_note(),
	mod_alluvial_ui("adv_alluvial_mod")
)

# Tab D — Distance Distributions
ui_tab_adv_ridgeline <- function()
tabItem(tabName = "adv_ridgeline",
	tab_header("mountain", "Distance Distributions",
		"How are motifs positioned relative to each other within the 5\u2032UTR?",
		"#92400E", "#1B7F3B"),
	uiOutput("warn_no_dataset_ridgeline"),
	hr(),
	imgt_ligm_note(),
	mod_ridgeline_ui("adv_ridgeline_mod")
)
# ══════════════════════════════════════════════════════════════
# Tab C — Carte Génomique des Motifs
# ══════════════════════════════════════════════════════════════
ui_tab_adv_gmap <- function()
tabItem(tabName = "adv_gmap",
	tab_header("map-marked-alt", "Genomic Map of Motifs",
		"How are motifs positioned relative to each other within the 5\u2032UTR?",
		"#0284C7", "#0B2A5B"),
	uiOutput("warn_no_dataset_gmap"),
	hr(),
	imgt_ligm_note(),
	mod_gmap_ui("adv_gmap_mod")
)
