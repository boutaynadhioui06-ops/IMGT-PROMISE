# ══════════════════════════════════════════════════════════════════════════════
# ██████████████████  tabItem (scanner motifs)  ████████████████████████████████
# ══════════════════════════════════════════════════════════════════════════════
# ─────────────────────────────────────────────────────────────────────────────
# 5UTR SCANNER
# ─────────────────────────────────────────────────────────────────────────────
ui_tab_scanner_consensus <- function()
tabItem(tabName = "scanner_consensus",
	tab_header("dna", "Motif Scanner — Consensus (baseline)",
		"Can regulatory motifs be automatically detected in species lacking curated IMGT motif annotations?",
		"#1B7F3B", "#0B2A5B", "#7B2D8B"),
	mod_scanner_ui_consensus("scanner_mod")
)
ui_tab_scanner_pwm <- function()
tabItem(tabName = "scanner_pwm",
	tab_header("dna", "Motif Scanner — PWM (probabiliste)",
		"Does a PWM-based probabilistic approach improve motif detection over the deterministic consensus baseline?",
		"#1B7F3B", "#0B2A5B", "#7B2D8B"),
	mod_scanner_ui_pwm("scanner_mod")
)
ui_tab_scanner_comparison <-  function()
tabItem(tabName = "scanner_comparison",
	tab_header("dna", "Motif Scanner — Method comparison",
		"Consensus vs PWM, evaluated against confirmed reference annotations.",
		"#1B7F3B", "#0B2A5B", "#7B2D8B"),
	mod_scanner_ui_comparison("scanner_mod")
)

ui_tab_scanner_final <- function()
	tabItem(tabName = "scanner_final",
		tab_header("sliders-h", "Motif Scanner — Final dataset",
			"Which detection method should be used for each motif to build the dataset analysed downstream?",
			"#1B7F3B", "#0B2A5B", "#7B2D8B"),
		mod_scanner_ui_final("scanner_mod")
	)