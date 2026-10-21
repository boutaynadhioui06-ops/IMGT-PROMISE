# ══════════════════════════════════════════════════════════════
# 5. CLUSTERING HIÉRARCHIQUE 
# ══════════════════════════════════════════════════════════════

clust_results <- reactiveValues(hc=NULL, done=FALSE, assignments=NULL)

observe({
	req(store$df_final)
	isolate({
		df_filtered <- store$df_final %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
		updateSelectInput(session, "clust_species", choices=unique(df_filtered$SPECIES), selected=unique(df_filtered$SPECIES)[1])
		seq_cols <- names(store$df_final)[str_starts(names(store$df_final),"SEQ_")]
		updateSelectInput(session, "clust_motif", choices=str_remove(seq_cols,"^SEQ_"), selected=str_remove(seq_cols,"^SEQ_")[1])
		# Update k slider max based on number of subgroups
		n_sg <- length(unique(df_filtered$SUBGROUP))
		k_max <- max(2, min(n_sg - 1, 15))
		updateSliderInput(session, "clust_k", max = k_max, value = min(3, k_max))
	})
})

observeEvent(input$btn_launch_clust, {
	req(store$df_final, input$clust_species)
	isolate({
		withProgress(message = "Running hierarchical clustering…", value = 0.5, {
			tryCatch({
				df      <- store$df_final
				method  <- input$clust_method
				species <- input$clust_species
				if (method == "presence") {
					hc <- hclust_motif_presence(df, species, "IG")
				} else {
					req(input$clust_motif)
					hc <- hclust_sequence_similarity(df, species, input$clust_motif, "IG")
				}
				k_val <- input$clust_k
				cluster_assignments <- tryCatch(
					cutree(hc, k = k_val),
					error = function(e) NULL
				)
				clust_results$hc          <- hc
				clust_results$assignments <- cluster_assignments
				clust_results$done        <- TRUE
				incProgress(1)
			}, error=function(e) {
				showNotification(paste0("❌ Clustering: ", e$message), type="error", duration=8)
			})
		})
	})
})

# Reactive re-couper when k changes 
observeEvent(input$clust_k, {
	req(clust_results$hc, clust_results$done)
	k_val <- input$clust_k
	cluster_assignments <- tryCatch(
		cutree(clust_results$hc, k = k_val),
		error = function(e) NULL
	)
	clust_results$assignments <- cluster_assignments
}, ignoreInit = TRUE)

output$clust_status <- renderUI({
	if (!clust_results$done)
		div(class="alert alert-info", icon("info-circle"),
			" Select a species and click ", strong("🚀 Generate dendrogram"), ".")
	else
		div(class="alert alert-success", icon("check-circle"), " Dendrogram generated for ", strong(input$clust_species))
})

output$clust_k_display <- renderText({ as.character(input$clust_k) })

output$clust_k_summary <- renderUI({
	req(clust_results$assignments)
	tbl <- table(clust_results$assignments)
	div(
		style = "font-size: 12px;",
		lapply(names(tbl), function(k) {
			div(style = paste0("display:inline-block; margin:3px; padding:4px 10px; ",
				"border-radius:4px; background:#E3F2FD; color:#1565C0; font-weight:bold;"),
				paste0("Cluster ", k, ": ", tbl[k], " subgroups"))
		})
	)
})

output$clust_plot <- renderPlot({
	req(clust_results$hc)
	k_val <- input$clust_k
	# couleures différentes pour k clusters 
	k_colors <- grDevices::hcl.colors(k_val, palette = "Dark 3")
	plot(clust_results$hc,
		main  = paste0("Hierarchical Clustering — ", input$clust_species,
			"  (k = ", k_val, " clusters)"),
		xlab  = "Subgroups", ylab = "Height (dissimilarity)", sub = "", cex = 0.8)
	rect.hclust(clust_results$hc, k = k_val, border = k_colors)
})

output$clust_dist_table <- renderDT({
	req(clust_results$hc)
	
	
	dist_mat <- clust_results$hc$dist_matrix
	
	if (is.null(dist_mat) || nrow(dist_mat) == 0) {
		return(datatable(
			data.frame(Message = "Distance matrix not available."),
			options = list(dom = "t"), rownames = FALSE
		))
	}
	
	dist_df           <- as.data.frame(round(dist_mat, 3))
	dist_df$Subgroup  <- rownames(dist_df)
	dist_df           <- dist_df %>% dplyr::select(Subgroup, dplyr::everything())
	
	n_cols <- ncol(dist_df)
	
	datatable(
		dist_df,
		options = list(
			pageLength = 25, scrollX = TRUE,
			dom        = "Bfrtip",
			buttons    = c("copy","csv")
		),
		extensions = "Buttons",
		rownames   = FALSE,
		class      = "cell-border stripe hover compact"
	) %>%
		formatStyle("Subgroup",
			fontWeight       = "bold",
			backgroundColor  = "#EFF6FF",
			color            = "#1E40AF") %>%
		formatRound(columns = 2:n_cols, digits = 3)
})

output$clust_assignments_table <- renderDT({
	req(clust_results$assignments)
	df_assign <- data.frame(
		Subgroup = names(clust_results$assignments),
		Cluster  = as.integer(clust_results$assignments)
	) %>% dplyr::arrange(Cluster, Subgroup)
	datatable(df_assign,
		options = list(pageLength = 50, scrollX = TRUE, dom = "Bfrtip",
			buttons = c("copy","csv")),
		extensions = "Buttons", rownames = FALSE,
		class = "cell-border stripe hover compact") %>%
		formatStyle("Cluster",
			backgroundColor = styleInterval(
				seq(1.5, max(df_assign$Cluster) - 0.5, by = 1),
				colorRampPalette(c("#DBEAFE","#DCFCE7","#FEF9C3","#FCE7F3",
					"#FFE4E6","#E0F2FE","#F3E8FF","#FEF3C7","#ECFDF5"))(
						max(df_assign$Cluster))),
			fontWeight = "bold")
})

output$btn_download_clust_plot <- downloadHandler(
	filename = function() paste0("Clustering_",input$clust_species,"_k",input$clust_k,"_",Sys.Date(),".png"),
	content  = function(file) {
		req(clust_results$hc)
		
		k_val     <- as.integer(input$clust_k)   
		k_colors  <- tryCatch(
			grDevices::hcl.colors(k_val, palette = "Dark 3"),
			error = function(e) rep("#E63946", k_val)
		)
		n_leaves  <- length(clust_results$hc$labels)
		w_px      <- as.integer(min(max(1800, n_leaves * 28 + 400), 8000))
		h_px      <- as.integer(min(max(1000, n_leaves * 8  + 600), 4000))
		cex_label <- max(0.45, min(0.90, 14 / n_leaves))
		
		grDevices::png(filename = file, width = w_px, height = h_px, res = 150)
		par(mar = c(max(6, n_leaves * 0.12 + 4), 5, 4, 2), xpd = TRUE)
		
		#  1 : tracer le dendrogramme 
		plot(clust_results$hc,
			main     = paste0("Hierarchical Clustering — ", input$clust_species,
				"  (k = ", k_val, " clusters)"),
			xlab     = "Subgroups",
			ylab     = "Height (dissimilarity)",
			sub      = "",
			cex      = cex_label,
			cex.main = 1.4,
			cex.lab  = 1.1)
		
		#  2 : tracer les rectangles 
		rect_ok <- tryCatch({
			stats::rect.hclust(clust_results$hc, k = k_val, border = k_colors)
			TRUE
		}, error = function(e) {
			message("rect.hclust failed : ", conditionMessage(e))
			FALSE
		})
		
		grDevices::dev.off()
		
		if (!rect_ok)
			showNotification("⚠️ Dendrogram downloaded without clade boxes (rect.hclust error).", type = "warning", duration = 6)
		else
			showNotification("✅ Dendrogram PNG downloaded!", type = "message", duration = 3)
	}
)

output$btn_download_clust_assignments <- downloadHandler(
	filename = function() paste0("ClusterAssignments_",input$clust_species,"_k",input$clust_k,"_",Sys.Date(),".csv"),
	content  = function(file) {
		req(clust_results$assignments)
		df_assign <- data.frame(
			Subgroup = names(clust_results$assignments),
			Cluster  = as.integer(clust_results$assignments)
		) %>% dplyr::arrange(Cluster, Subgroup)
		write.csv(df_assign, file, row.names = FALSE)
		showNotification("✅ Assignments CSV downloaded!", type="message", duration=3)
	}
)
