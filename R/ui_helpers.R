# ═══════════════════════════════════════════════════════════════
# IMGT-PROMISE — Éléments d'interface réutilisables
# ═══════════════════════════════════════════════════════════════

# Dépend de : rien
# Packages  : shiny

# Bandeau jaune Data source note,
imgt_ligm_note <- function() {
	div(
		style = "background-color:#FEF3C7; border-left:4px solid #D97706; padding:10px 14px;
                 border-radius:6px; margin-bottom:15px; font-size:12.5px; color:#374151;",
		icon("database"), strong(" Data source note: "),
		"All results in this section depend on the motif annotations currently available in ",
		strong("IMGT-LIGM-DB"),
		". Motifs not annotated in IMGT-LIGM-DB cannot be detected or analyzed here."
	)
}

# ── Bandeau de titre  + question de recherche  ───
tab_header <- function(icon_name, title, question, color1, color2, color3 = NULL) {
	grad <- if (!is.null(color3)) {
		sprintf("linear-gradient(135deg, %s 0%%, %s 55%%, %s 100%%)", color1, color2, color3)
	} else {
		sprintf("linear-gradient(135deg, %s 0%%, %s 100%%)", color1, color2)
	}
	fluidRow(
		column(12,
			div(
				style = paste0("background:", grad, "; border-radius:12px; padding:16px 24px;
                 margin-bottom:18px; box-shadow:0 4px 14px rgba(0,0,0,0.15);"),
				div(style = "display:flex; align-items:center; gap:16px;",
					div(style = "font-size:40px; color:white; flex-shrink:0;", icon(icon_name)),
					h2(style = "color:white; margin:0; font-size:22px; font-weight:800;", title)
				),
				div(
					style = "background:rgba(255,255,255,0.16); border-left:4px solid #FFD700;
                   border-radius:8px; padding:8px 14px; margin-top:12px; font-size:13.5px; color:#FFFFFF;",
					icon("compass", style = "margin-right:8px;"),
					tags$strong("Research question: "), tags$em(question)
				)
			)
		)
	)
}

