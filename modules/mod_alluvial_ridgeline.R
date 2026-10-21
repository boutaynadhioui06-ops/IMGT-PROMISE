# ═══════════════════════════════════════════════════════════════
#Alluvial analysis
# ═══════════════════════════════════════════════════════════════

# Dépend de : visualization_advanced.R
# Packages  : shiny, plotly, dplyr, stringr

# ══════════════════════════════════════════════════════════════
# ALLUVIAL FLOW
# ══════════════════════════════════════════════════════════════

mod_alluvial_ui <- function(id) {
	ns <- shiny::NS(id)
	shiny::tagList(
		shiny::div(class = "result-info-bar",
			shiny::icon("info-circle"), shiny::strong(" Principle: "),
			"Visualises the flow of sequences between ",
			shiny::strong("SPECIES \u2192 LOCI \u2192 FUNCTIONALITY"),
			". Each ribbon represents sequences sharing the same path. ",
			"Colour = functionality (Green = F, Red = P, Yellow = ORF)."
		),
		shiny::fluidRow(
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("filter"), shiny::strong("Filters")),
					shiny::selectInput(ns("alluvial_species"), "Species (optional)",
						choices = c("All species" = ""), selected = ""),
					shiny::tags$small(shiny::icon("info-circle"), " Empty = all species")
				)
			),
			shiny::column(3,
				shiny::div(style = "padding-top: 30px;",
					shiny::actionButton(ns("btn_run_alluvial"),
						"🌊 Generate flow",
						icon  = shiny::icon("play"),
						class = "btn-success btn-lg",
						style = "width:100%; height:70px; font-size:14px;")
				)
			)
		),
		shiny::hr(),
		shiny::uiOutput(ns("alluvial_status")),
		plotly::plotlyOutput(ns("alluvial_plot"), height = "540px")
	)
}

mod_alluvial_server <- function(id, df_final) {
	shiny::moduleServer(id, function(input, output, session) {
		ns             <- session$ns
		alluvial_store <- shiny::reactiveValues(plot = NULL, done = FALSE)
		
		shiny::observe({
			shiny::req(df_final())
			df      <- df_final() %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
			sp_list <- sort(unique(df$SPECIES))
			shiny::isolate({
				shiny::updateSelectInput(session, "alluvial_species",
					choices  = c("All species" = "", setNames(sp_list, sp_list)),
					selected = "")
			})
		})
		
		shiny::observeEvent(input$btn_run_alluvial, { run_alluvial() })
		
		run_alluvial <- function() {
			shiny::req(df_final())
			tryCatch({
				sp <- if (is.null(input$alluvial_species) || input$alluvial_species == "")
					NULL else input$alluvial_species
				p <- plot_alluvial_parcats(df_final(), type_data = "IG", species = sp)
				alluvial_store$plot <- p
				alluvial_store$done <- TRUE
			}, error = function(e) {
				shiny::showNotification(paste0("❌ Alluvial: ", e$message),
					type = "error", duration = 8)
			})
		}
		
		output$alluvial_status <- shiny::renderUI({
			if (!alluvial_store$done)
				shiny::div(class = "alert alert-info",
					shiny::icon("info-circle"),
					" Select optional filters and click ",
					shiny::strong("🌊 Generate flow"), " to start.")
			else
				shiny::div(class = "alert alert-success",
					shiny::icon("check-circle"), " Alluvial flow generated.")
		})
		
		output$alluvial_plot <- plotly::renderPlotly({
			shiny::req(alluvial_store$done, alluvial_store$plot)
			alluvial_store$plot
		})
	})
}


# ══════════════════════════════════════════════════════════════
# DISTANCE DISTRIBUTIONS
# ══════════════════════════════════════════════════════════════

mod_ridgeline_ui <- function(id) {
	ns <- shiny::NS(id)
	shiny::tagList(
		shiny::div(class = "result-info-bar",
			shiny::icon("info-circle"), shiny::strong(" Principle: "),
			"Visualises the distribution of distances ",
			shiny::strong("motif \u2192 LPART1"),
			" (in nucleotides) by species, locus or subgroup. ",
			"Helps detect bimodality, inter-species/locus variance and conserved positions."
		),
		shiny::fluidRow(
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("dna"), shiny::strong("Motif")),
					shiny::selectInput(ns("ridge_motif"), NULL,
						choices = c("All motifs" = ""), selected = ""),
					shiny::tags$small(shiny::icon("info-circle"),
						" Empty = all motifs displayed as facets")
				)
			),
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("layer-group"), shiny::strong("Grouping level")),
					shiny::radioButtons(ns("ridge_group"), NULL,
						choices = c(
							"By species"  = "SPECIES",
							"By locus"    = "LOCI",
							"By subgroup" = "SUBGROUP"
						),
						selected = "SPECIES")
				)
			),
			shiny::column(4,
				shiny::conditionalPanel(
					condition = paste0("input['", ns("ridge_group"), "'] == 'LOCI'"),
					shiny::div(class = "param-block",
						shiny::h5(shiny::icon("globe"), shiny::strong("Filters")),
						shiny::selectInput(ns("ridge_species_locus"),
							"Species (optional)",
							choices  = c("All species" = ""),
							selected = ""),
						shiny::tags$small(shiny::icon("info-circle"),
							" Empty = all species combined")
					)
				),
				shiny::conditionalPanel(
					condition = paste0("input['", ns("ridge_group"), "'] == 'SUBGROUP'"),
					shiny::div(class = "param-block",
						shiny::h5(shiny::icon("globe"), shiny::strong("Filters")),
						shiny::selectInput(ns("ridge_species_sg"),
							"Species",
							choices  = NULL,
							selected = NULL),
						shiny::selectInput(ns("ridge_locus_sg"),
							"Locus (optional)",
							choices  = c("All loci" = ""),
							selected = ""),
						shiny::tags$small(shiny::icon("info-circle"),
							" Species is required for subgroup mode")
					)
				),
				shiny::conditionalPanel(
					condition = paste0("input['", ns("ridge_group"), "'] == 'SPECIES'"),
					shiny::div(
						style = "background:#F0FDF4; border-left:4px solid #22C55E;
                       border-radius:6px; padding:12px 14px;",
						shiny::icon("info-circle"),
						shiny::tags$span(style = "color:#166534; font-size:12px;",
							" All IG species are compared simultaneously.")
					)
				)
			),
			shiny::column(2,
				shiny::div(style = "padding-top: 30px;",
					shiny::actionButton(ns("btn_run_ridge"),
						"📊 Run analysis",
						icon  = shiny::icon("play"),
						class = "btn-success btn-lg",
						style = "width:100%; height:70px; font-size:14px;")
				)
			)
		),
		shiny::hr(),
		shiny::uiOutput(ns("ridge_status")),
		shiny::uiOutput(ns("ridge_plot_ui")),
		shiny::br(),
		shiny::fluidRow(
			shiny::column(4,
				shiny::div(
					style = "text-align:center; padding:10px; background:#EFF6FF;
                   border-radius:8px; border:1px solid #93C5FD;",
					shiny::icon("image",
						style = "font-size:20px; color:#1D4ED8; display:block; margin-bottom:6px;"),
					shiny::p(style = "font-size:11px; color:#374151; margin:0 0 6px 0;",
						"Full figure (PNG 300 dpi)"),
					shiny::downloadButton(ns("dl_ridge_png"),
						"📥 PNG — Distributions",
						class = "btn-primary btn-block", style = "width:100%;")
				)
			),
			shiny::column(4,
				shiny::div(
					style = "text-align:center; padding:10px; background:#FFFBEB;
                   border-radius:8px; border:1px solid #FCD34D;",
					shiny::icon("file-csv",
						style = "font-size:20px; color:#D97706; display:block; margin-bottom:6px;"),
					shiny::p(style = "font-size:11px; color:#374151; margin:0 0 6px 0;",
						"Raw data table (CSV)"),
					shiny::downloadButton(ns("dl_ridge_csv"),
						"📥 CSV — Data",
						class = "btn-warning btn-block", style = "width:100%;")
				)
			)
		)
	)
}

mod_ridgeline_server <- function(id, df_final) {
	shiny::moduleServer(id, function(input, output, session) {
		ns          <- session$ns
		ridge_store <- shiny::reactiveValues(plot = NULL, data = NULL, done = FALSE)
		
		# Populate selectors
		shiny::observe({
			shiny::req(df_final())
			df        <- df_final() %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
			sp_list   <- sort(unique(df$SPECIES))
			loci_list <- sort(unique(df$LOCI))
			dist_cols <- names(df_final())[stringr::str_starts(names(df_final()), "DIST_")]
			motif_names <- stringr::str_remove(dist_cols, "^DIST_")
			
			shiny::isolate({
				shiny::updateSelectInput(session, "ridge_motif",
					choices  = c("All motifs" = "", setNames(motif_names, motif_names)),
					selected = "")
				# Locus mode species filter
				shiny::updateSelectInput(session, "ridge_species_locus",
					choices  = c("All species" = "", setNames(sp_list, sp_list)),
					selected = "")
				# Subgroup mode species filter 
				shiny::updateSelectInput(session, "ridge_species_sg",
					choices  = setNames(sp_list, sp_list),
					selected = sp_list[1])
				# Subgroup mode locus filter
				shiny::updateSelectInput(session, "ridge_locus_sg",
					choices  = c("All loci" = "", setNames(loci_list, loci_list)),
					selected = "")
			})
		})
		shiny::observeEvent(input$btn_run_ridge, { run_ridge() })
				run_ridge <- function() {
			shiny::req(df_final())
			tryCatch({
				grp <- input$ridge_group %||% "SPECIES"
				mot <- if (is.null(input$ridge_motif) || input$ridge_motif == "") NULL else input$ridge_motif
				
				if (grp == "LOCI") {
					sp_val  <- if (is.null(input$ridge_species_locus) || input$ridge_species_locus == "")
						NULL else input$ridge_species_locus
					loc_val <- NULL
					
				} else if (grp == "SUBGROUP") {
					shiny::req(input$ridge_species_sg)
					sp_val  <- input$ridge_species_sg
					loc_val <- if (is.null(input$ridge_locus_sg) || input$ridge_locus_sg == "")
						NULL else input$ridge_locus_sg
					
				} else {
					sp_val  <- NULL
					loc_val <- NULL
				}
				
				df_use <- df_final()
				if (!is.null(loc_val) && loc_val != "")
					df_use <- df_use %>% dplyr::filter(LOCI == loc_val)
				
				p <- plot_distance_ridgeline(df_use,
					motif       = mot,
					group_var   = grp,
					species     = sp_val,
					type_data   = "IG",
					interactive = TRUE)
				
				df_filt <- df_use %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
				if (!is.null(sp_val) && sp_val != "")
					df_filt <- df_filt %>% dplyr::filter(SPECIES == sp_val)
				dist_cols <- names(df_filt)[stringr::str_starts(names(df_filt), "DIST_")]
				if (!is.null(mot) && mot != "") {
					dc <- paste0("DIST_", mot)
					if (dc %in% dist_cols) dist_cols <- dc
				}
				ridge_store$data <- df_filt %>%
					dplyr::select(dplyr::all_of(c(grp, dist_cols)))
				ridge_store$plot <- p
				ridge_store$done <- TRUE
				
			}, error = function(e) {
				shiny::showNotification(paste0("❌ Distance distributions: ", e$message),
					type = "error", duration = 8)
			})
		}
		
		#  Outputs
		output$ridge_status <- shiny::renderUI({
			if (!ridge_store$done)
				shiny::div(class = "alert alert-info",
					shiny::icon("info-circle"),
					" Select grouping level and filters, then click ",
					shiny::strong("📊 Run analysis"), " to start.")
			else
				shiny::div(class = "alert alert-success",
					shiny::icon("check-circle"), " Distribution plot generated.")
		})
		
		output$ridge_plot_ui <- shiny::renderUI({
			shiny::req(ridge_store$done)
			grp      <- input$ridge_group %||% "SPECIES"
			n_groups <- if (!is.null(ridge_store$data) && grp %in% names(ridge_store$data))
				length(unique(ridge_store$data[[grp]])) else 4
			h <- max(480, n_groups * 45)
			plotly::plotlyOutput(ns("ridge_plot"), height = paste0(min(h, 850), "px"))
		})
		
		output$ridge_plot <- plotly::renderPlotly({
			shiny::req(ridge_store$done, ridge_store$plot)
			ridge_store$plot
		})
		
		# Downloads
		output$dl_ridge_png <- shiny::downloadHandler(
			filename = function() paste0("Distances_", input$ridge_group, "_", Sys.Date(), ".png"),
			content  = function(file) {
				shiny::req(df_final())
				grp <- input$ridge_group %||% "SPECIES"
				mot <- if (is.null(input$ridge_motif) || input$ridge_motif == "") NULL else input$ridge_motif
				
				if (grp == "LOCI") {
					sp_val  <- if (is.null(input$ridge_species_locus) || input$ridge_species_locus == "")
						NULL else input$ridge_species_locus
					df_use  <- df_final()
				} else if (grp == "SUBGROUP") {
					sp_val  <- input$ridge_species_sg
					loc_val <- if (is.null(input$ridge_locus_sg) || input$ridge_locus_sg == "")
						NULL else input$ridge_locus_sg
					df_use  <- df_final()
					if (!is.null(loc_val)) df_use <- df_use %>% dplyr::filter(LOCI == loc_val)
				} else {
					sp_val <- NULL
					df_use <- df_final()
				}
				
				p <- plot_distance_ridgeline(df_use,
					motif = mot, group_var = grp, species = sp_val,
					type_data = "IG", interactive = FALSE)
				n_groups <- if (!is.null(ridge_store$data) && grp %in% names(ridge_store$data))
					length(unique(ridge_store$data[[grp]])) else 4
				h <- max(6, n_groups * 0.5 + 2)
				ggplot2::ggsave(file, plot = p, width = 12,
					height = min(h, 40), dpi = 300, bg = "white")
				shiny::showNotification("✅ PNG downloaded!", type = "message", duration = 3)
			})
		
		output$dl_ridge_csv <- shiny::downloadHandler(
			filename = function() paste0("Distances_", input$ridge_group, "_", Sys.Date(), ".csv"),
			content  = function(file) {
				shiny::req(ridge_store$data)
				write.csv(ridge_store$data, file, row.names = FALSE)
				shiny::showNotification("✅ CSV downloaded!", type = "message", duration = 3)
			})
	})
}

