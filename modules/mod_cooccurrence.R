# ═══════════════════════════════════════════════════════════════
# mod_cooccurrence
# ═══════════════════════════════════════════════════════════════

# Dépend de : cooccurrence_analysis.R, cooccurrence_network.R, visualization_advanced.R
# Packages  : shiny, plotly, DT, htmlwidgets, webshot2

mod_cooccurrence_ui <- function(id) {
	ns <- shiny::NS(id)
	shiny::tagList(
		
		shiny::div(class = "result-info-bar",
			shiny::icon("info-circle"),
			shiny::strong(" Principle: "),
			"For each pair of motifs, computes the % of sequences carrying both simultaneously, ",
			"the phi (\u03c6) correlation, and tests independence with Fisher's exact test. ",
			"The network visualizes significant associations."
		),
		
		# Note sur Plotly
		shiny::div(
			style = "background:#EFF6FF; border-left:4px solid #1565C0;
                 border-radius:6px; padding:8px 14px; font-size:12px;
                 margin-bottom:10px;"),
		
		shiny::fluidRow(
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("globe"), shiny::strong("Filters")),
					shiny::selectInput(ns("cooc_species"), "Species",
						choices = NULL, selected = NULL),
					shiny::selectInput(ns("cooc_locus"), "Locus",
						choices = c("All" = ""), selected = ""),
					shiny::selectInput(ns("cooc_func"), "Functionality",
						choices = c("All" = "", "F" = "F", "P" = "P", "ORF" = "ORF"),
						selected = "")
				)
			),
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("sliders-h"), shiny::strong("Network options")),
					shiny::sliderInput(ns("cooc_phi_threshold"),
						"Minimum |\u03c6| threshold",
						min = 0, max = 0.8, value = 0.01, step = 0.01),
					shiny::checkboxInput(ns("cooc_sig_only"),
						"Significant pairs only (p < 0.05)", value = FALSE),
					shiny::radioButtons(ns("cooc_heatmap_mode"),
						"Heatmap: display",
						choices = c("\u03c6 correlation" = "phi", "% co-occurrence" = "pct"),
						selected = "phi", inline = TRUE)
				)
			),
			shiny::column(4,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("info-circle"), shiny::strong("Interpretation")),
					shiny::tags$ul(
						shiny::tags$li(shiny::tags$span(style = "color:#16A34A;font-weight:bold;", "\u03c6 > 0"),
							" : the two motifs tend to co-occur"),
						shiny::tags$li(shiny::tags$span(style = "color:#DC2626;font-weight:bold;", "\u03c6 < 0"),
							" : the motifs mutually exclude each other"),
						shiny::tags$li(shiny::tags$span(style = "color:#F59E0B;font-weight:bold;", "OR > 1"),
							" : odds ratio — the presence of one increases the odds of the other"),
						shiny::tags$li("|\u03c6| \u2265 0.2 recommended for biological relevance")
					)
				)
			),
			shiny::column(2,
				shiny::div(style = "padding-top: 30px;",
					shiny::actionButton(ns("btn_run_cooc"),
						"🔬 Run analysis",
						icon  = shiny::icon("play"),
						class = "btn-success btn-lg",
						style = "width:100%; height:70px; font-size:14px;")
				)
			)
		),
		
		shiny::hr(),
		shiny::uiOutput(ns("cooc_status")),
		
		shiny::tabsetPanel(
			
			# Tab 1 : Réseau Plotly 
			shiny::tabPanel(
				shiny::tagList(shiny::icon("share-alt"), " Interactive network"),
				shiny::br(),
				shiny::div(
					style = "background:#F8FAFF; border-radius:6px; padding:8px 12px;
                     font-size:11px; color:#374151; margin-bottom:8px;
                     display:flex; gap:18px; flex-wrap:wrap;",
					shiny::span(
						shiny::span(style = "display:inline-block;width:10px;height:10px;
                                    border-radius:50%;background:#16A34A;margin-right:4px;"),
						"Positive co-occurrence"),
					shiny::span(
						shiny::span(style = "display:inline-block;width:10px;height:10px;
                                    border-radius:50%;background:#DC2626;margin-right:4px;"),
						"Mutual exclusion"),
					shiny::span(
						shiny::span(style = "display:inline-block;width:10px;height:10px;
                                    border-radius:50%;background:#CBD5E1;margin-right:4px;"),
						"Not significant"),
					shiny::em(style = "color:#6B7280;",
						"Node size \u221d frequency | Edge width \u221d |\u03c6|")
				),
				plotly::plotlyOutput(ns("cooc_network"), height = "580px"),
				shiny::br(),
				shiny::downloadButton(ns("dl_cooc_network_png"), "📥 PNG — Network",
					class = "btn-primary btn-block", style = "width:100%; margin-top:8px;"),
				shiny::br()
			),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("th"), " \u03c6 / Co-occurrence heatmap"),
				shiny::br(),
				plotly::plotlyOutput(ns("cooc_heatmap"), height = "500px"),
				shiny::br(),
				shiny::downloadButton(ns("dl_cooc_png"), "📥 PNG — Heatmap",
					class = "btn-primary btn-block", style = "width:100%; margin-top:8px;"),
				shiny::br()
			),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("chart-bar"), " Fisher tests"),
				shiny::br(),
				shiny::div(
					style = "background:#FEF9C3;padding:10px;border-radius:6px;
                   border-left:4px solid #EAB308;margin-bottom:12px;font-size:12px;",
					shiny::icon("exclamation-triangle"),
					shiny::strong(" Bonferroni correction recommended"),
					" for multiple comparisons. Adjusted threshold = 0.05 / N pairs."
				),
				plotly::plotlyOutput(ns("cooc_fisher_plot"), height = "500px"),
				shiny::br(),
				shiny::downloadButton(ns("dl_fisher_plot_png"), "📥 PNG — Fisher plot",
					class = "btn-primary btn-block", style = "width:100%; margin-top:8px;"),
				shiny::br()
			),
			
			shiny::tabPanel(
				shiny::tagList(shiny::icon("table"), " Association table"),
				shiny::br(),
				DT::DTOutput(ns("cooc_fisher_table")),
				shiny::br(),
				shiny::downloadButton(ns("dl_fisher_csv"), "📥 CSV — Fisher tests",
					class = "btn-warning btn-block", style = "width:100%;"),
				shiny::br()
			)
		)
	)
}


mod_cooccurrence_server <- function(id, df_final) {
	shiny::moduleServer(id, function(input, output, session) {
		
		ns         <- session$ns
		cooc_store <- shiny::reactiveValues(
			cooc_result = NULL, fisher_df = NULL, graph = NULL, done = FALSE
		)
		
		shiny::observe({
			shiny::req(df_final())
			df      <- df_final() %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
			sp_list <- sort(unique(df$SPECIES))
			shiny::isolate({
				shiny::updateSelectInput(session, "cooc_species",
					choices  = c("All" = "", setNames(sp_list, sp_list)), selected = "")
				shiny::updateSelectInput(session, "cooc_locus",
					choices  = c("All" = "", "IGHV" = "IGHV", "IGKV" = "IGKV", "IGLV" = "IGLV"),
					selected = "")
			})
		})
		
		# Auto-launch au chargement
		shiny::observeEvent(df_final(), { run_cooc(reset = TRUE) }, ignoreNULL = TRUE)
		shiny::observeEvent(input$btn_run_cooc, { run_cooc() })
		
		run_cooc <- function(reset = FALSE) {
			shiny::req(df_final())
			shiny::withProgress(message = "Computing co-occurrences...", value = 0, {
				tryCatch({
					sp  <- if (reset || is.null(input$cooc_species) || input$cooc_species == "") NULL else input$cooc_species
					loc <- if (reset || is.null(input$cooc_locus)   || input$cooc_locus   == "") NULL else input$cooc_locus
					fn  <- if (reset || is.null(input$cooc_func)    || input$cooc_func    == "") NULL else input$cooc_func
					
					shiny::incProgress(0.1, detail = "Co-occurrence matrix...")
					cooc_res <- calculate_cooccurrence(df_final(), type_data = "IG",
						species = sp, locus = loc, func_filter = fn)
					
					shiny::incProgress(0.4, detail = "Fisher's tests...")
					fisher_res <- fisher_pairwise_tests(cooc_res)
					
					shiny::incProgress(0.6, detail = "Building the graph...")
					phi_thresh <- if (!is.null(input$cooc_phi_threshold)) input$cooc_phi_threshold else 0.01
					sig_only   <- isTRUE(input$cooc_sig_only)
					
					graph_res <- tryCatch(
						build_cooccurrence_graph(fisher_res, cooc_res,
							phi_threshold = phi_thresh, sig_only = sig_only),
						error = function(e) NULL)
					
					shiny::incProgress(1, detail = "Done.")
					cooc_store$cooc_result <- cooc_res
					cooc_store$fisher_df   <- fisher_res
					cooc_store$graph       <- graph_res
					cooc_store$done        <- TRUE
					
				}, error = function(e) {
					shiny::showNotification(paste0("❌ Co-occurrence: ", e$message),
						type = "error", duration = 8)
				})
			})
		}
		
		output$cooc_status <- shiny::renderUI({
			if (!cooc_store$done)
				shiny::div(class = "alert alert-info", shiny::icon("info-circle"),
					" Configure the filters and click ",
					shiny::strong("🔬 Run analysis"), ".")
			else {
				n_sig   <- if (!is.null(cooc_store$fisher_df)) sum(cooc_store$fisher_df$Significant != "ns") else 0
				n_total <- if (!is.null(cooc_store$fisher_df)) nrow(cooc_store$fisher_df) else 0
				shiny::div(class = "alert alert-success", shiny::icon("check-circle"),
					shiny::strong(" Analysis complete — "),
					cooc_store$cooc_result$n_seq, " sequences | ",
					n_sig, " / ", n_total, " significant pairs (p < 0.05)")
			}
		})
		
		# ── Réseau PLOTLY 
		output$cooc_network <- plotly::renderPlotly({
			shiny::req(cooc_store$done, cooc_store$cooc_result)
			
			phi_thresh <- if (!is.null(input$cooc_phi_threshold)) input$cooc_phi_threshold else 0.01
			sig_only   <- isTRUE(input$cooc_sig_only)
			
			tryCatch(
				render_cooc_network_plotly(
					cooc_result   = cooc_store$cooc_result,
					fisher_df     = cooc_store$fisher_df,
					phi_threshold = phi_thresh,
					sig_only      = sig_only
				),
				error = function(e) {
					plotly::plot_ly() %>%
						plotly::layout(title = paste("Network error:", e$message),
							plot_bgcolor = "#FEF2F2")
				}
			)
		})
		
		output$cooc_heatmap <- plotly::renderPlotly({
			shiny::req(cooc_store$done, cooc_store$cooc_result)
			tryCatch(
				plot_cooccurrence_heatmap(cooc_store$cooc_result,
					mode = input$cooc_heatmap_mode %||% "phi", interactive = TRUE),
				error = function(e) plotly::plot_ly() %>%
					plotly::layout(title = paste("Error:", e$message)))
		})
		
		output$cooc_fisher_plot <- plotly::renderPlotly({
			shiny::req(cooc_store$done, cooc_store$fisher_df)
			tryCatch(
				plot_fisher_associations(cooc_store$fisher_df, interactive = TRUE),
				error = function(e) plotly::plot_ly() %>%
					plotly::layout(title = paste("Error:", e$message)))
		})
		
		output$cooc_fisher_table <- DT::renderDT({
			shiny::req(cooc_store$fisher_df)
			DT::datatable(cooc_store$fisher_df,
				options = list(pageLength = 20, scrollX = TRUE, dom = "Bfrtip",
					buttons = c("copy","csv")),
				extensions = "Buttons", rownames = FALSE,
				class = "cell-border stripe hover compact") %>%
				DT::formatRound(
					columns = c("Freq_Motif1","Freq_Motif2","CoOcc_Pct","OddsRatio"),
					digits = 2) %>%
				DT::formatStyle("Significant",
					backgroundColor = DT::styleEqual(
						c("***","**","*","ns"),
						c("#D4EDDA","#D4EDDA","#FFF3CD","#F8D7DA")),
					fontWeight = "bold") %>%
				DT::formatStyle("Association",
					color = DT::styleEqual(
						c("Co-occurrence +","Exclusion mutuelle","Non significatif",
							"Positive co-occurrence","Mutual exclusion","Not significant"),
						c("#16A34A","#DC2626","#6B7280",
							"#16A34A","#DC2626","#6B7280")),
					fontWeight = "bold")
		})
		
		output$dl_fisher_csv <- shiny::downloadHandler(
			filename = function() paste0("Fisher_Cooccurrence_", Sys.Date(), ".csv"),
			content  = function(file) {
				shiny::req(cooc_store$fisher_df)
				write.csv(cooc_store$fisher_df, file, row.names = FALSE)
				shiny::showNotification("✅ CSV downloaded!", type = "message", duration = 3)
			})
		
		output$dl_cooc_png <- shiny::downloadHandler(
			filename = function() paste0("Cooccurrence_Heatmap_", Sys.Date(), ".png"),
			content  = function(file) {
				shiny::req(cooc_store$cooc_result)
				p <- plot_cooccurrence_heatmap(cooc_store$cooc_result,
					mode = input$cooc_heatmap_mode %||% "phi", interactive = FALSE)
				ggplot2::ggsave(file, plot = p, width = 10, height = 8,
					dpi = 300, bg = "white")
				shiny::showNotification("✅ PNG downloaded!", type = "message", duration = 3)
			})
	
		output$dl_cooc_network_png <- shiny::downloadHandler(
			filename = function() paste0("Network_Cooccurrence_", Sys.Date(), ".png"),
			content  = function(file) {
				shiny::req(cooc_store$done, cooc_store$cooc_result)
				shiny::withProgress(message = "Exporting network as PNG…", value = 0.3, {
					tryCatch({
						if (!requireNamespace("webshot2", quietly = TRUE))
							stop("Package 'webshot2' is required. Install with install.packages('webshot2').")
						
						phi_thresh <- if (!is.null(input$cooc_phi_threshold)) input$cooc_phi_threshold else 0.01
						p_net <- render_cooc_network_plotly(
							cooc_result   = cooc_store$cooc_result,
							fisher_df     = cooc_store$fisher_df,
							phi_threshold = phi_thresh,
							sig_only      = isTRUE(input$cooc_sig_only)
						)
						
						tmp_html <- tempfile(fileext = ".html")
						htmlwidgets::saveWidget(plotly::as_widget(p_net), tmp_html, selfcontained = TRUE)
						shiny::incProgress(0.5, detail = "Rendering image…")
						webshot2::webshot(tmp_html, file = file, vwidth = 1200, vheight = 800, delay = 1)
						shiny::incProgress(1)
						shiny::showNotification("✅ Network PNG downloaded!", type = "message", duration = 3)
					}, error = function(e) {
						shiny::showNotification(paste0("❌ Network export: ", e$message),
							type = "error", duration = 8)
					})
				})
			}
		)
		
		output$dl_fisher_plot_png <- shiny::downloadHandler(
			filename = function() paste0("Fisher_Plot_", Sys.Date(), ".png"),
			content  = function(file) {
				shiny::req(cooc_store$done, cooc_store$fisher_df)
				tryCatch({
					p <- plot_fisher_associations(cooc_store$fisher_df, interactive = FALSE)
					ggplot2::ggsave(file, plot = p, width = 11, height = 8, dpi = 300, bg = "white")
					shiny::showNotification("✅ Fisher plot PNG downloaded!", type = "message", duration = 3)
				}, error = function(e) {
					shiny::showNotification(paste0("❌ Fisher plot export: ", e$message),
						type = "error", duration = 8)
				})
			}
		)
	})
}