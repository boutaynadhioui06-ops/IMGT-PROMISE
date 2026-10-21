# ══════════════════════════════════════════════════════════════
# PHASE 4 — ADVANCED ANALYSES
# ══════════════════════════════════════════════════════════════

# ── Warnings dataset manquant ─────────────────────────────────
output$warn_no_dataset <- renderUI({
	if (is.null(store$df_final))
		div(class = "alert alert-warning",
			icon("exclamation-triangle"), strong(" Warning : "),
			"Build the dataset first (onglet 📊 Final Dataset ).")
})
output$warn_no_dataset_alluvial <- renderUI({
	if (is.null(store$df_final))
		div(class = "alert alert-warning",
			icon("exclamation-triangle"), strong(" Warning : "),
			"Build the dataset first (Dataset Final tab).")
})
output$warn_no_dataset_ridgeline <- renderUI({
	if (is.null(store$df_final))
		div(class = "alert alert-warning",
			icon("exclamation-triangle"), strong(" Warning : "),
			"Build the dataset first (Dataset Final tab).")
})
output$warn_no_dataset_permanova <- renderUI({
	if (is.null(store$df_final))
		div(class = "alert alert-warning",
			icon("exclamation-triangle"), strong(" Warning : "),
			"Build the dataset first (Dataset Final tab).")
})

# ── warnings pour les onglets C et E ───────────────
output$warn_no_dataset_gmap <- renderUI({
	if (is.null(store$df_final))
		div(class = "alert alert-warning",
			icon("exclamation-triangle"), strong(" Warning : "),
			"Build the dataset first (Dataset Final tab).")
})
output$warn_no_dataset_tsne <- renderUI({
	if (is.null(store$df_final))
		div(class = "alert alert-warning",
			icon("exclamation-triangle"), strong(" Warning : "),
			"Build the dataset first (Dataset Final tab).")
})


mod_cooccurrence_server("adv_cooc_mod",      df_final = reactive(store$df_final))
mod_alluvial_server(    "adv_alluvial_mod",  df_final = reactive(store$df_final))
mod_ridgeline_server(   "adv_ridgeline_mod", df_final = reactive(store$df_final))
mod_permanova_server(   "adv_permanova_mod", df_final = reactive(store$df_final))

mod_gmap_server("adv_gmap_mod", df_final = reactive(store$df_final))
mod_tsne_server("adv_tsne_mod", df_final = reactive(store$df_final))

# ══════════════════════════════════════════════════════════════
# MODULE MACHINE LEARNING
# ══════════════════════════════════════════════════════════════



scanner_res <- mod_scanner_server("scanner_mod", df_final = reactive(store$df_final))

observeEvent(scanner_res$go_analyses(), {
	df <- scanner_res$merged_df()
	req(df)
	store$df_final      <- df
	store$df_no_motif   <- identify_no_motif_sequences(df)
	store$df_with_motif <- identify_with_motif_sequences(df)
	store$df_partiels   <- store$df_no_motif
	showNotification("✅ Analyses now use the new species dataset.", type = "message", duration = 6)
	shinydashboard::updateTabItems(session, "tabs", "stat_desc_tableau")
}, ignoreInit = TRUE, ignoreNULL = TRUE)