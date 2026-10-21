# ══════════════════════════════════════════════════════════════
# 2. RÉGRESSION LOGISTIQUE 
# ══════════════════════════════════════════════════════════════

regr_results <- reactiveValues(data=NULL, model=NULL, done=FALSE, ref_levels=NULL)
regr_valid_results <- reactiveValues(data = NULL, done = FALSE)
observe({
	req(store$df_final)
	isolate({
		bin_cols <- names(store$df_final)[str_starts(names(store$df_final),"BIN_")]
		motif_names <- str_remove(bin_cols, "^BIN_")
		updateSelectInput(session, "regr_motif", choices=motif_names, selected=motif_names[1])
	})
})

observeEvent(input$btn_launch_regr, {
	req(store$df_final, input$regr_motif, nchar(input$regr_motif) > 0)
	isolate({
		withProgress(message = "Running logistic regression + CV validation...", value = 0.2, {
			tryCatch({
				# 1. Modèle complet 
				results <- logistic_regression_motif(store$df_final, input$regr_motif, "IG")
				regr_results$data    <- results$results
				regr_results$model   <- results$model
				regr_results$summary <- results$summary
				regr_results$done    <- TRUE
				regr_results$ref_levels <- results$ref_levels
				incProgress(0.4, detail = "Cross-validation...")
				
				# 2. Validation croisée 
				valid_res <- validate_logistic_model(
					df        = store$df_final,
					motif_name = input$regr_motif,
					type_data  = "IG",
					n_folds    = 5,
					seed       = 42
				)
				regr_valid_results$data <- valid_res
				regr_valid_results$done <- TRUE
				incProgress(1)
				
			}, error = function(e) {
				showNotification(paste0("❌ Regression: ", e$message), type="error", duration=8)
			})
		})
	})
})

output$regr_status <- renderUI({
	if (!regr_results$done)
		div(class="alert alert-info", icon("info-circle"),
			" Select a motif and click ", strong("🚀 Run model"), ".")
	else
		div(class="alert alert-success", icon("check-circle"), " Model computed for ", strong(input$regr_motif))
})


output$regr_table <- renderDT({
	req(regr_results$data)
	datatable(regr_results$data,
		options = list(pageLength = 25, scrollX = TRUE),
		rownames = FALSE, class = 'cell-border stripe hover') %>%
		formatRound(columns = c('estimate','std.error','OR','CI_lower','CI_upper'), digits = 3) %>%
		formatStyle('Significant', backgroundColor = styleEqual(c('***','**','*','ns'),
			c('#D4EDDA','#D4EDDA','#FFF3CD','#F8D7DA'))) %>%
		formatStyle('Reference', fontStyle = 'italic', color = '#6B7280')
})


output$regr_reference_box <- renderUI({
	req(regr_results$done, regr_results$ref_levels)
	rl <- regr_results$ref_levels
	div(
		style = "background:#EFF6FF; border-left:4px solid #1565C0;
	             border-radius:6px; padding:12px 16px; font-size:13px;",
		icon("info-circle"),
		strong(" Reference levels used in this model (baseline, not shown in the table): "),
		tags$ul(style = "margin-top:6px; margin-bottom:0;",
			tags$li(strong("LOCI: "), rl$LOCI),
			tags$li(strong("SPECIES: "), rl$SPECIES),
			tags$li(strong("FUNCTIONALITY: "), rl$FUNCTIONALITY)
		),
		tags$small(style = "color:#374151;",
			"Every OR / coefficient in the table below is a comparison against these baseline levels.")
	)
})

output$regr_summary <- renderPrint({ req(regr_results$summary); print(regr_results$summary) })

output$btn_download_regr_csv <- downloadHandler(
	filename = function() paste0("Regression_", input$regr_motif, "_", Sys.Date(), ".csv"),
	content  = function(file) { req(regr_results$data); write.csv(regr_results$data, file, row.names=FALSE) }
)
output$regr_interpretation <- renderUI({
	req(regr_valid_results$done)
	text <- format_logistic_interpretation(regr_valid_results$data)
	div(style="background:#F0FDF4;border-left:4px solid #16A34A;
    border-radius:6px;padding:14px;font-size:13px;",
		HTML(gsub("\n","<br>",
			gsub("\\*\\*(.*?)\\*\\*","<strong>\\1</strong>", text))))
})

output$regr_roc_plot <- renderPlotly({
	req(regr_valid_results$done)
	plot_roc_logistic(regr_valid_results$data, interactive = TRUE)
	
})

output$regr_cv_plot <- renderPlotly({
	req(regr_valid_results$done)
	plot_cv_metrics_logistic(regr_valid_results$data, interactive = TRUE)
})


output$regr_cv_wide_table <- renderDT({
	req(regr_valid_results$done)
	df <- format_cv_metrics_wide(regr_valid_results$data$cv_metrics)
	req(df)
	fold_cols <- setdiff(names(df), "Metric")
	datatable(df,
		options = list(dom = "Bfrtip", scrollX = TRUE, buttons = c("copy","csv"), paging = FALSE),
		extensions = "Buttons", rownames = FALSE,
		class = "cell-border stripe hover compact") %>%
		formatRound(columns = fold_cols, digits = 2) %>%
		formatStyle("Metric", fontWeight = "bold", backgroundColor = "#EFF6FF", color = "#1E40AF") %>%
		formatStyle("Mean", fontWeight = "bold", backgroundColor = "#FEF9C3")
})

output$btn_dl_regr_cv_wide_csv <- downloadHandler(
	filename = function() paste0("CV_Metrics_AllFolds_", input$regr_motif, "_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(regr_valid_results$done)
		df <- format_cv_metrics_wide(regr_valid_results$data$cv_metrics)
		write.csv(df, file, row.names = FALSE)
		showNotification("✅ CSV downloaded!", type = "message", duration = 3)
	}
)

output$regr_calibration <- renderPlotly({
	req(regr_valid_results$done)
	plot_calibration_curve(regr_valid_results$data, interactive = TRUE)
})

output$regr_confusion <- renderPlotly({
	req(regr_valid_results$done)
	plot_confusion_logistic(regr_valid_results$data, interactive = TRUE)
})

# ── Downloads validation régression logistique ─────────────────────────
output$btn_dl_regr_roc_png <- downloadHandler(
	filename = function() paste0("ROC_", input$regr_motif, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(regr_valid_results$done)
		p <- plot_roc_logistic(regr_valid_results$data, interactive = FALSE)
		ggplot2::ggsave(file, plot = p, width = 9, height = 7, dpi = 300, bg = "white")
		showNotification("✅ ROC PNG downloaded!", type = "message", duration = 3)
	}
)

output$btn_dl_regr_cv_png <- downloadHandler(
	filename = function() paste0("CV_Metrics_", input$regr_motif, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(regr_valid_results$done)
		p <- plot_cv_metrics_logistic(regr_valid_results$data, interactive = FALSE)
		ggplot2::ggsave(file, plot = p, width = 11, height = 7, dpi = 300, bg = "white")
		showNotification("✅ CV metrics PNG downloaded!", type = "message", duration = 3)
	}
)

output$btn_dl_regr_calib_png <- downloadHandler(
	filename = function() paste0("Calibration_", input$regr_motif, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(regr_valid_results$done)
		p <- plot_calibration_curve(regr_valid_results$data, interactive = FALSE)
		ggplot2::ggsave(file, plot = p, width = 9, height = 7, dpi = 300, bg = "white")
		showNotification("✅ Calibration PNG downloaded!", type = "message", duration = 3)
	}
)

output$btn_dl_regr_conf_png <- downloadHandler(
	filename = function() paste0("Confusion_Matrix_", input$regr_motif, "_", Sys.Date(), ".png"),
	content  = function(file) {
		req(regr_valid_results$done)
		p <- plot_confusion_logistic(regr_valid_results$data, interactive = FALSE)
		ggplot2::ggsave(file, plot = p, width = 7, height = 6, dpi = 300, bg = "white")
		showNotification("✅ Confusion matrix PNG downloaded!", type = "message", duration = 3)
	}
)

output$btn_dl_regr_cv_csv <- downloadHandler(
	filename = function() paste0("CV_Metrics_", input$regr_motif, "_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(regr_valid_results$done)
		write.csv(regr_valid_results$data$cv_metrics, file, row.names = FALSE)
		showNotification("✅ CV metrics CSV downloaded!", type = "message", duration = 3)
	}
)

output$btn_dl_regr_interp_csv <- downloadHandler(
	filename = function() paste0("Logistic_Summary_", input$regr_motif, "_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(regr_valid_results$done)
		m <- regr_valid_results$data$metrics_global
		df_summary <- data.frame(
			Motif       = input$regr_motif,
			N_valid     = m$N_valid,
			N_total     = m$N_total,
			N_folds     = regr_valid_results$data$n_folds,
			Accuracy    = m$Accuracy,
			Precision   = m$Precision,
			Recall      = m$Recall,
			F1          = m$F1,
			Specificity = m$Specificity,
			AUC         = m$AUC
		)
		write.csv(df_summary, file, row.names = FALSE)
		showNotification("✅ Summary CSV downloaded!", type = "message", duration = 3)
	}
)
