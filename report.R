# ══ RESOURCES — Analysis Report download ═══════════════════════
output$btn_dl_dynamic_report <- downloadHandler(
	filename = function() {
		ext <- switch(input$report_format, word = "docx", pdf = "pdf", md = "md")
		paste0("IMGT_PROMISE_Report_", Sys.Date(), ".", ext)
	},
	content = function(file) {
		req(store$df_final)
		showNotification("⏳ Generating report…", type = "message", duration = NULL, id = "notif_report")
		
		tryCatch({
			tmp_rmd <- file.path(tempdir(), "report.Rmd")
			file.copy("report_templates/report.Rmd", tmp_rmd, overwrite = TRUE)
			
			report_params <- isolate(list(
				species         = input$tableau_species,
				sections        = input$report_sections,
				dataset_summary = list(
					n    = nrow(store$df_final),
					loci = unique(store$df_final$LOCI)
				),
				tableau    = tableau_results$data,
				freq       = list(data = freq_results$data, stat_tests = freq_results$stat_tests),
				diversity  = list(data = div_results$data,  tests      = div_results$tests),
				clustering = list(hc = clust_results$hc,    assignments = clust_results$assignments),
				identity   = identity_results$data,
				regression = list(data = regr_results$data, cv = regr_valid_results$data),
				cons_func  = list(data = cons_func_results$data, chi = cons_func_results$chi),
				scanner    = NULL
			))
			
			out_format <- switch(input$report_format,
				word = "word_document",
				pdf  = "pdf_document",
				md   = "md_document")
			
			rmarkdown::render(
				input         = tmp_rmd,
				output_format = out_format,
				output_file   = file,
				params        = report_params,
				envir         = new.env(parent = globalenv()),
				quiet         = TRUE
			)
			
			removeNotification("notif_report")
			showNotification("✅ Report generated!", type = "message", duration = 5)
		}, error = function(e) {
			removeNotification("notif_report")
			showNotification(paste0("❌ Report error: ", e$message), type = "error", duration = 10)
		})
	}
)


output$btn_dl_analysis_report <- downloadHandler(
	filename = function() paste0("IMGT_PROMISE_Analysis_Report_", Sys.Date(), ".pdf"),
	content  = function(file) {
		report_path <- "www/analysis_report.pdf"   # je dois placer le pdf ici  
		if (file.exists(report_path)) {
			file.copy(report_path, file)
		} else {
			showNotification("❌ Report file not found (www/analysis_report.pdf).", type = "error", duration = 8)
		}
	}
)