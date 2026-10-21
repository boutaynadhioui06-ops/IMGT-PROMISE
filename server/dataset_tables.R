#  Dataset final : tableaux IG / TR et exports
# ══ TABLEAU INTERACTIF IG et TR SÉPARÉS ═══════════════════════
dataset_split <- reactive({
	req(store$df_final)
	split_dataset_ig_tr(store$df_final)
})

output$main_table_ig <- renderDT({
	req(dataset_split()$ig)
	df_display <- rename_loci_display(dataset_split()$ig)
	datatable(
		df_display,
		options = list(
			pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv","excel"),
			fixedHeader = TRUE
		),
		extensions = c("Buttons", "FixedHeader"),
		rownames = FALSE, filter = "top",
		class = "cell-border stripe hover"
	) %>%
		formatStyle("FUNCTIONALITY",
			backgroundColor = styleEqual(c("F","P","ORF"), c("#D4EDDA","#F8D7DA","#FFF3CD"))) %>%
		formatStyle("GROUP", fontWeight = "bold")
})

output$main_table_tr <- renderDT({
	req(dataset_split()$tr)
	df_display <- rename_loci_display(dataset_split()$tr)
	datatable(
		df_display,
		options = list(
			pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv","excel"),
			fixedHeader = TRUE
		),
		extensions = c("Buttons", "FixedHeader"),
		rownames = FALSE, filter = "top",
		class = "cell-border stripe hover"
	) %>%
		formatStyle("FUNCTIONALITY",
			backgroundColor = styleEqual(c("F","P","ORF"), c("#D4EDDA","#F8D7DA","#FFF3CD"))) %>%
		formatStyle("GROUP", fontWeight = "bold")
})

# ══ EXPORT 
output$btn_dl_excel_ig <- downloadHandler(
	filename = function() paste0("IMGT_PROMISE_IG_", Sys.Date(), ".xlsx"),
	content  = function(file) {
		req(dataset_split()$ig)
		df <- dataset_split()$ig
		showNotification("⏳ Generating Excel (IG)...", type="message", duration=NULL, id="notif_xl_ig")
		tryCatch({
			library(openxlsx)
			wb <- createWorkbook()
			header_style <- createStyle(textDecoration="bold", fgFill="#1F4E79", fontColour="#FFFFFF")
			loci_list <- unique(df$LOCI)
			loci_list <- loci_list[!is.na(loci_list)]
			if (length(loci_list) == 0) {
				addWorksheet(wb, "IG_Data")
				writeData(wb, "IG_Data", df, rowNames=FALSE)
				addStyle(wb, "IG_Data", header_style, rows=1, cols=1:ncol(df), gridExpand=TRUE)
			} else {
				for (locus in loci_list) {
					df_l <- df[df$LOCI == locus, ]
					sn   <- substr(locus, 1, 31)
					addWorksheet(wb, sn)
					writeData(wb, sn, df_l, rowNames=FALSE)
					addStyle(wb, sn, header_style, rows=1, cols=1:ncol(df_l), gridExpand=TRUE)
					setColWidths(wb, sn, cols=1:ncol(df_l), widths="auto")
				}
			}
			saveWorkbook(wb, file, overwrite=TRUE)
			removeNotification("notif_xl_ig")
			showNotification("✅ Excel (IG) created!", type="message", duration=4)
		}, error = function(e) {
			removeNotification("notif_xl_ig")
			showNotification(paste("❌", e$message), type="error", duration=8)
		})
	}
)

output$btn_dl_csv_ig <- downloadHandler(
	filename = function() paste0("IMGT_PROMISE_IG_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(dataset_split()$ig)
		write.csv(dataset_split()$ig, file, row.names=FALSE, fileEncoding="UTF-8")
		showNotification("✅ CSV (IG) created!", type="message", duration=4)
	}
)

output$btn_dl_excel_tr <- downloadHandler(
	filename = function() paste0("IMGT_PROMSIE_TR_", Sys.Date(), ".xlsx"),
	content  = function(file) {
		req(dataset_split()$tr)
		df <- dataset_split()$tr
		showNotification("⏳ Generating Excel (TR)...", type="message", duration=NULL, id="notif_xl_tr")
		tryCatch({
			library(openxlsx)
			wb <- createWorkbook()
			header_style <- createStyle(textDecoration="bold", fgFill="#1F4E79", fontColour="#FFFFFF")
			loci_list <- unique(df$LOCI)
			loci_list <- loci_list[!is.na(loci_list)]
			if (length(loci_list) == 0) {
				addWorksheet(wb, "TR_Data")
				writeData(wb, "TR_Data", df, rowNames=FALSE)
				addStyle(wb, "TR_Data", header_style, rows=1, cols=1:ncol(df), gridExpand=TRUE)
			} else {
				for (locus in loci_list) {
					df_l <- df[df$LOCI == locus, ]
					sn   <- substr(locus, 1, 31)
					addWorksheet(wb, sn)
					writeData(wb, sn, df_l, rowNames=FALSE)
					addStyle(wb, sn, header_style, rows=1, cols=1:ncol(df_l), gridExpand=TRUE)
					setColWidths(wb, sn, cols=1:ncol(df_l), widths="auto")
				}
			}
			saveWorkbook(wb, file, overwrite=TRUE)
			removeNotification("notif_xl_tr")
			showNotification("✅ Excel (TR) created!", type="message", duration=4)
		}, error = function(e) {
			removeNotification("notif_xl_tr")
			showNotification(paste("❌", e$message), type="error", duration=8)
		})
	}
)

output$btn_dl_csv_tr <- downloadHandler(
	filename = function() paste0("IMGT_PROMISE_TR_", Sys.Date(), ".csv"),
	content  = function(file) {
		req(dataset_split()$tr)
		write.csv(dataset_split()$tr, file, row.names=FALSE, fileEncoding="UTF-8")
		showNotification("✅ CSV (TR) created!", type="message", duration=4)
	}
)

# Nettoie la fonctionnalité : garde uniquement F / P / ORF
clean_functionality <- function(raw) {
	f <- gsub("[\\[\\]\\(\\)]", "", raw, perl = TRUE)
	f <- trimws(f)
	f <- toupper(f)
	if (!f %in% c("F", "P", "ORF")) f <- NA_character_
	return(f)
}