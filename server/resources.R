# ═══════════════════════════════════════════════════════════════════════════
#  Serveur de la section Resources
# ═══════════════════════════════════════════════════════════════════════════

.res_dt <- function(df, page_length = 10) {
	DT::datatable(df, rownames = FALSE, class = "cell-border stripe hover compact",
		options = list(pageLength = page_length, scrollX = TRUE, dom = "frtip"))
}

# Documentation
output$doc_data_dictionary <- DT::renderDT(.res_dt(DATA_DICTIONARY, 15))
output$doc_glossary        <- DT::renderDT(.res_dt(GLOSSARY, 25))

output$doc_motif_defs <- DT::renderDT({
	df <- if (exists("CONSENSUS_DEFINITIONS", inherits = TRUE)) {
		defs <- get("CONSENSUS_DEFINITIONS", inherits = TRUE)
		out <- purrr::imap_dfr(defs, function(motifs, loc)
			purrr::imap_dfr(motifs, function(d, m)
				data.frame(Locus = loc, Motif = m, Pattern = paste(d$patterns, collapse = "  |  "),
					`Max mismatches` = as.integer(d$max_mismatch), check.names = FALSE,
					stringsAsFactors = FALSE)))
		rbind(out, data.frame(Locus = "IGHV / IGKV / IGLV", Motif = "TATA_BOX",
			Pattern = "5–8 nt, A/T-rich composition rule",
			`Max mismatches` = 1L, check.names = FALSE, stringsAsFactors = FALSE))
	} else {
		data.frame(Info = "Motif definitions are available once the scanner module is loaded.")
	}
	.res_dt(df, 15)
})

# Citation
output$cite_text <- renderText(promise_citation_text())
output$btn_dl_cite_bibtex <- downloadHandler(
	filename = function() "IMGT-PROMISE.bib",
	content  = function(file) writeLines(promise_bibtex(), file, useBytes = TRUE))
output$btn_dl_cite_ris <- downloadHandler(
	filename = function() "IMGT-PROMISE.ris",
	content  = function(file) writeLines(promise_ris(), file, useBytes = TRUE))

#  About : environnement logiciel
output$about_r_info <- renderUI({
	tags$p(icon("r-project", class = "fab"), " ", R.version.string, " · ", R.version$platform)
})
output$about_pkg_table <- DT::renderDT({
	pk <- unique(c(REQUIRED_PACKAGES))
	ver <- vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed"),
		character(1))
	.res_dt(data.frame(Package = pk, Version = unname(ver), stringsAsFactors = FALSE), 20)
})

# Rapport : disponibilité des sections
report_scanner_cmp <- function() {
	scn <- get0("scanner_res", inherits = TRUE)
	if (is.list(scn) && is.function(scn$comparison_df)) scn$comparison_df() else NULL
}

report_avail <- reactive({
	df <- store$df_final
	has_df <- !is.null(df) && nrow(df) > 0
	list(
		overview    = has_df,
		descriptive = has_df,
		freq        = has_df && any(grepl("^BIN_", names(df))),
		diversity   = has_df && any(grepl("^BIN_", names(df))),
		distances   = has_df && any(grepl("^DIST_", names(df))),
		partial     = !is.null(store$df_no_motif) && nrow(store$df_no_motif) > 0,
		scanner     = !is.null(report_scanner_cmp())
	)
})

.refresh_report_choices <- function(selected = NULL) {
	av <- report_avail()
	vals <- unname(REPORT_SECTIONS)
	labels <- lapply(seq_along(vals), function(i) {
		ok <- isTRUE(av[[vals[i]]])
		HTML(sprintf("%s <span class='res-status %s'>%s</span>", names(REPORT_SECTIONS)[i],
			if (ok) "ok" else "na", if (ok) "ready" else "not available"))
	})
	updateCheckboxGroupInput(session, "report_sections", choiceNames = labels, choiceValues = vals,
		selected = selected)
}

observeEvent(report_avail(), {
	av <- report_avail()
	keep <- intersect(input$report_sections, names(av)[unlist(av)])
	if (is.null(input$report_sections)) keep <- names(av)[unlist(av)]
	.refresh_report_choices(keep)
}, ignoreNULL = FALSE)

observeEvent(input$report_select_all, {
	av <- report_avail(); .refresh_report_choices(names(av)[unlist(av)])
})
observeEvent(input$report_select_none, .refresh_report_choices(character(0)))

observe({
	df <- store$df_final
	if (is.null(df)) return()
	updateSelectizeInput(session, "report_species", choices = sort(unique(na.omit(df$SPECIES))), server = TRUE)
	updateSelectizeInput(session, "report_locus",   choices = sort(unique(na.omit(df$LOCI))),    server = TRUE)
})

.pdf_available <- function() {
	nzchar(Sys.which("pdflatex")) || nzchar(Sys.which("xelatex")) ||
		(requireNamespace("tinytex", quietly = TRUE) && isTRUE(tinytex::is_tinytex()))
}

output$report_format_warning <- renderUI({
	if (identical(input$report_format, "pdf") && !.pdf_available())
		div(class = "alert alert-warning", style = "font-size:12px;", icon("exclamation-triangle"),
			" PDF export needs a LaTeX distribution (TinyTeX) on the server. Use HTML or Word instead.")
})

output$report_summary <- renderUI({
	if (is.null(store$df_final))
		return(div(class = "alert alert-warning", style = "font-size:12px;", icon("exclamation-triangle"),
			" Build the dataset first (Dataset Construction → Final Dataset)."))
	secs <- input$report_sections
	div(class = "res-summary",
		tags$b(length(secs)), " section(s) selected · ",
		tags$b(if (length(input$report_species)) length(input$report_species) else "all"), " species · ",
		tags$b(if (length(input$report_locus)) length(input$report_locus) else "all"), " loci")
})

# Rapport : génération
output$btn_dl_dynamic_report <- downloadHandler(
	filename = function() {
		ext <- switch(input$report_format, html = "html", word = "docx", pdf = "pdf", "html")
		sprintf("IMGT-PROMISE_report_%s.%s", format(Sys.Date(), "%Y%m%d"), ext)
	},
	content = function(file) {
		fail <- function(msg) { showNotification(paste0("❌ ", msg), type = "error", duration = 12); stop(msg, call. = FALSE) }
		if (is.null(store$df_final)) fail("Build the dataset first.")
		if (length(input$report_sections) == 0) fail("Select at least one report section.")
		fmt_key <- input$report_format
		if (fmt_key == "pdf" && !.pdf_available()) fail("PDF export requires LaTeX (TinyTeX). Choose HTML or Word.")
		
		withProgress(message = "Rendering report...", value = 0.1, {
			tmp <- tempfile("promise_report_"); dir.create(tmp)
			on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
			tpl <- "report/report_template.Rmd"
			if (!file.exists(tpl)) fail("Report template not found: report/report_template.Rmd")
			rmd <- file.path(tmp, "report.Rmd"); file.copy(tpl, rmd)
			
			fmt <- switch(fmt_key,
				html = rmarkdown::html_document(toc = TRUE, toc_float = TRUE, number_sections = TRUE,
					theme = "flatly", df_print = "paged", self_contained = TRUE),
				word = rmarkdown::word_document(toc = TRUE),
				pdf  = rmarkdown::pdf_document(toc = TRUE, number_sections = TRUE))
			
			params <- list(
				title = input$report_title, author = input$report_author, notes = input$report_notes,
				data = store$df_final, partial = store$df_no_motif,
				scanner_cmp = report_scanner_cmp(),
				sections = input$report_sections,
				species = input$report_species, locus = input$report_locus,
				include_repro = isTRUE(input$report_include_repro),
				app_name = APP_META$name, app_version = APP_META$version,
				app_release = APP_META$release_date, app_citation = promise_citation_text())
			incProgress(0.3, detail = "Computing statistics and figures")
			out <- tryCatch(
				rmarkdown::render(rmd, output_format = fmt, params = params, quiet = TRUE,
					envir = new.env(parent = globalenv())),
				error = function(e) fail(paste("Report rendering failed:", conditionMessage(e))))
			incProgress(0.9)
			file.copy(out, file, overwrite = TRUE)
		})
	}
)