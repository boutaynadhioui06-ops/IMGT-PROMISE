# ═══════════════════════════════════════════════════════════════
#  Motif Phylogeny
# ═══════════════════════════════════════════════════════════════

# Dépend de : utils_core.R
# Packages  : shiny, ape, phangorn, ggtree, ggplot2, viridis, DT

# ══════════════════════════════════════════════════════════════
#  MOTIF SEQUENCE EXTRACTION
# ══════════════════════════════════════════════════════════════

.motif_seq_columns <- function(df, motifs = NULL) {
	seq_cols <- names(df)[stringr::str_starts(names(df), "SEQ_")]
	seq_cols <- seq_cols[!stringr::str_detect(seq_cols, "LPART")]
	motif_names <- stringr::str_remove(seq_cols, "^SEQ_")
	if (!is.null(motifs)) {
		keep <- motif_names %in% motifs
		seq_cols    <- seq_cols[keep]
		motif_names <- motif_names[keep]
	}
	has_bin <- paste0("BIN_", motif_names) %in% names(df)
	data.frame(
		motif    = motif_names[has_bin],
		seq_col  = seq_cols[has_bin],
		bin_col  = paste0("BIN_", motif_names[has_bin]),
		stringsAsFactors = FALSE
	)
}
.consensus_motif_seq <- function(seqs) {
	seqs <- seqs[!is.na(seqs) & nchar(seqs) > 0]
	if (length(seqs) == 0) return(NA_character_)
	if (length(seqs) == 1) return(toupper(seqs))
	
	lens     <- nchar(seqs)
	ref_len  <- as.integer(names(sort(table(lens), decreasing = TRUE))[1])
	seqs_up  <- toupper(seqs)
	
	chars_mat <- matrix(NA_character_, nrow = length(seqs_up), ncol = ref_len)
	for (i in seq_along(seqs_up)) {
		s <- seqs_up[i]
		n <- nchar(s)
		use_n <- min(n, ref_len)
		if (use_n > 0) chars_mat[i, seq_len(use_n)] <- strsplit(substr(s, 1, use_n), "")[[1]]
	}
	
	cons <- vapply(seq_len(ref_len), function(j) {
		col <- chars_mat[, j]
		col <- col[!is.na(col)]
		if (length(col) == 0) return("N")
		names(sort(table(col), decreasing = TRUE))[1]
	}, character(1))
	
	paste(cons, collapse = "")
}

.build_entity_sequences <- function(df, group_by, motif_tbl) {
	
	if (group_by == "ALLELE") {
		entities <- unique(df$ALLELE)
		seq_mat  <- matrix(NA_character_, nrow = length(entities),
			ncol = nrow(motif_tbl), dimnames = list(entities, motif_tbl$motif))
		
		df_one <- df %>% dplyr::group_by(ALLELE) %>% dplyr::slice(1) %>% dplyr::ungroup()
		meta <- df_one %>%
			dplyr::select(dplyr::any_of(c("ALLELE","LOCI","SUBGROUP","SPECIES","FUNCTIONALITY"))) %>%
			dplyr::rename(Entity = ALLELE)
		
		for (i in seq_len(nrow(motif_tbl))) {
			m  <- motif_tbl$motif[i]
			sc <- motif_tbl$seq_col[i]
			bc <- motif_tbl$bin_col[i]
			row_match <- match(rownames(seq_mat), df_one$ALLELE)
			present <- !is.na(df_one[[bc]][row_match]) & df_one[[bc]][row_match] == 1
			vals <- df_one[[sc]][row_match]
			vals[!present] <- NA
			seq_mat[, m] <- toupper(as.character(vals))
		}
		
	} else { 
		entities <- unique(df$SUBGROUP)
		seq_mat  <- matrix(NA_character_, nrow = length(entities),
			ncol = nrow(motif_tbl), dimnames = list(entities, motif_tbl$motif))
		
		meta <- df %>%
			dplyr::group_by(SUBGROUP) %>%
			dplyr::summarise(
				LOCI    = dplyr::first(LOCI),
				SPECIES = dplyr::first(SPECIES),
				.groups = "drop") %>%
			dplyr::rename(Entity = SUBGROUP)
		
		for (i in seq_len(nrow(motif_tbl))) {
			m  <- motif_tbl$motif[i]
			sc <- motif_tbl$seq_col[i]
			bc <- motif_tbl$bin_col[i]
			for (sg in entities) {
				sub <- df[df$SUBGROUP == sg, ]
				present <- !is.na(sub[[bc]]) & sub[[bc]] == 1
				vals <- sub[[sc]][present]
				seq_mat[sg, m] <- .consensus_motif_seq(vals)
			}
		}
	}
	
	list(seq_matrix = seq_mat, meta = meta)
}


# ══════════════════════════════════════════════════════════════
#  MOTIF-BY-MOTIF SEQUENCE DISTANCE
# ══════════════════════════════════════════════════════════════

.pdist_pair <- function(s1, s2) {
	if (is.na(s1) && is.na(s2)) return(0)        
	if (is.na(s1) || is.na(s2)) return(1)        
	if (s1 == s2) return(0)
	n1 <- nchar(s1); n2 <- nchar(s2)
	n_min <- min(n1, n2); n_max <- max(n1, n2)
	if (n_min == 0) return(1)
	c1 <- strsplit(substr(s1, 1, n_min), "")[[1]]
	c2 <- strsplit(substr(s2, 1, n_min), "")[[1]]
	n_diff <- sum(c1 != c2) + (n_max - n_min)
	n_diff / n_max
}

.build_phylo_distance <- function(seq_matrix) {
	n <- nrow(seq_matrix)
	d <- matrix(0, n, n, dimnames = list(rownames(seq_matrix), rownames(seq_matrix)))
	for (i in seq_len(n - 1)) {
		for (j in (i + 1):n) {
			per_motif <- vapply(seq_len(ncol(seq_matrix)), function(k)
				.pdist_pair(seq_matrix[i, k], seq_matrix[j, k]), numeric(1))
			dij <- mean(per_motif, na.rm = TRUE)
			d[i, j] <- dij
			d[j, i] <- dij
		}
	}
	as.dist(d)
}


# ══════════════════════════════════════════════════════════════
# MAIN PIPELINE
# ══════════════════════════════════════════════════════════════

compute_motif_phylogeny <- function(df,
	type_data   = "IG",
	species     = NULL,
	locus       = NULL,
	func_filter = NULL,
	group_by    = "SUBGROUP",
	min_motifs  = 2) {
	
	if (type_data == "IG") {
		df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	} else {
		df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	}
	if (!is.null(species) && species != "" && species != "All")
		df <- df %>% dplyr::filter(SPECIES == species)
	if (!is.null(locus) && locus != "" && locus != "All")
		df <- df %>% dplyr::filter(LOCI == locus)
	if (!is.null(func_filter) && func_filter != "" && func_filter != "All")
		df <- df %>% dplyr::filter(FUNCTIONALITY == func_filter)
	if (nrow(df) == 0) stop("No sequences left after filtering.")
	
	motif_tbl <- .motif_seq_columns(df)
	if (nrow(motif_tbl) < min_motifs)
		stop(paste0("Fewer than ", min_motifs, " motifs with SEQ_* columns available."))
	
	built <- .build_entity_sequences(df, group_by, motif_tbl)
	seq_mat <- built$seq_matrix
	meta    <- built$meta
	
	n_motifs_per_entity <- rowSums(!is.na(seq_mat))
	keep <- n_motifs_per_entity >= min_motifs
	if (sum(keep) < 3)
		stop(paste0("Fewer than 3 entities with at least ", min_motifs,
			" sequenced motifs — phylogeny impossible. Try lowering the threshold or changing the grouping level."))
	
	seq_mat <- seq_mat[keep, , drop = FALSE]
	meta    <- meta[match(rownames(seq_mat), meta$Entity), ]
	
	dist_obj <- .build_phylo_distance(seq_mat)
	
	tree <- tryCatch(ape::nj(dist_obj), error = function(e) NULL)
	if (is.null(tree) || any(tree$edge.length < 0, na.rm = TRUE)) {
		hc   <- hclust(dist_obj, method = "average")
		tree <- ape::as.phylo(hc)
	}
	
	tree <- tryCatch({
		if (!ape::is.rooted(tree)) phangorn::midpoint(tree) else tree
	}, error = function(e) tree)

	if (!is.null(tree$edge.length)) {
		bad <- is.na(tree$edge.length) | tree$edge.length < 0
		if (any(bad)) tree$edge.length[bad] <- 0
	}
	
	list(
		tree        = tree,
		dist        = dist_obj,
		seq_matrix  = seq_mat,
		meta        = meta,
		motifs_used = motif_tbl$motif,
		group_by    = group_by,
		n_entities  = nrow(seq_mat)
	)
}


# ══════════════════════════════════════════════════════════════
# VISUALIZATION 
# ══════════════════════════════════════════════════════════════

plot_motif_phylo_tree <- function(phylo_result,
	color_by     = "LOCI",
	layout       = "rectangular",
	interactive  = FALSE,
	branch_style = "cladogram") {
	
	tree <- phylo_result$tree
	meta <- phylo_result$meta
	if (!color_by %in% names(meta)) color_by <- "LOCI"
	
	group_vals <- unique(na.omit(meta[[color_by]]))
	if (color_by == "LOCI") {
		pal <- tryCatch({
			lc <- get("LOCI_COLORS", envir = .GlobalEnv)
			lc[group_vals]
		}, error = function(e) NULL)
		if (is.null(pal) || any(is.na(pal)))
			pal <- setNames(viridis::viridis(length(group_vals)), group_vals)
	} else if (color_by == "FUNCTIONALITY") {
		pal <- c("F" = "#16A34A", "P" = "#DC2626", "ORF" = "#F59E0B")[group_vals]
		pal[is.na(pal)] <- "#94A3B8"
		names(pal) <- group_vals
	} else {
		pal <- tryCatch(get_species_palette(group_vals),
			error = function(e) setNames(viridis::viridis(length(group_vals)), group_vals))
	}
	
	branch_len_arg <- if (identical(branch_style, "real")) "branch.length" else "none"
	
	p <- ggtree::ggtree(tree, layout = layout, branch.length = branch_len_arg,
		color = "#64748B", linewidth = 0.5, open.angle = 0) %<+% meta +
		ggtree::geom_tippoint(ggplot2::aes(color = .data[[color_by]]), size = 2.6) +
		ggtree::geom_tiplab(size = 3, offset = 0.01, color = "#1F2937") +
		ggplot2::scale_color_manual(values = pal, name = color_by, na.value = "#94A3B8") +
		ggplot2::labs(
			title    = paste0("Motif Phylogeny (", phylo_result$group_by, ")"),
			subtitle = paste0(
				"Mean per-motif p-distance | ",
				phylo_result$n_entities, " entities | ",
				length(phylo_result$motifs_used), " motif(s): ",
				paste(phylo_result$motifs_used, collapse = ", "),
				if (identical(branch_len_arg, "none")) " | equal-spacing (topology only)" else "")
		) +
		ggplot2::theme(
			plot.title    = ggplot2::element_text(face = "bold", size = 13, color = "#0B2A5B"),
			plot.subtitle = ggplot2::element_text(size = 9, color = "#4B5563"),
			legend.position = "right",
			plot.margin     = ggplot2::margin(6, 6, 6, 6)
		)
	
	if (layout == "rectangular")
		p <- p + ggtree::theme_tree2()
	
	if (interactive) {
		p <- tryCatch(
			plotly::ggplotly(p, tooltip = c("label")) %>%
				plotly::layout(height = max(450, phylo_result$n_entities * 22 + 100)),
			error = function(e) p   
		)
	}
	return(p)
}


# ══════════════════════════════════════════════════════════════
# UI 
# ══════════════════════════════════════════════════════════════

mod_tsne_ui <- function(id) {
	ns <- shiny::NS(id)
	shiny::tagList(
		
		shiny::div(
			style = "background: linear-gradient(90deg,#EFF6FF,#F0FDF4);
                 border-radius:6px; padding:10px 14px; font-size:12px;
                 color:#374151; margin-bottom:12px;",
			shiny::icon("info-circle"),
			shiny::strong(" Motif Phylogeny — "),
			"Unlike Hierarchical Clustering, which compares one motif at a time, ",
			"this builds a single evolutionary tree from the ",
			shiny::strong("combined sequences of all motifs"), "."
		),
		
		shiny::fluidRow(
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("layer-group"), shiny::strong("Grouping level")),
					shiny::radioButtons(ns("phy_group_by"), NULL,
						choices = c(
							"Subgroup (consensus)" = "SUBGROUP",
							"Allele (individual)"  = "ALLELE"
						),
						selected = "SUBGROUP"),
					shiny::numericInput(ns("phy_min_motifs"),
						"Min. sequenced motifs required",
						value = 2, min = 1, max = 10, step = 1)
				)
			),
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("filter"), shiny::strong("Filters")),
					shiny::selectInput(ns("phy_species"), "Species (optional)",
						choices = c("All" = ""), selected = ""),
					shiny::selectInput(ns("phy_locus"), "Locus (optional)",
						choices = c("All" = ""), selected = ""),
					shiny::selectInput(ns("phy_func"), "Functionality (optional)",
						choices = c("All" = "", "F" = "F", "P" = "P", "ORF" = "ORF"),
						selected = "")
				)
			),
			shiny::column(3,
				shiny::div(class = "param-block",
					shiny::h5(shiny::icon("palette"), shiny::strong("Visualization")),
					shiny::selectInput(ns("phy_color_by"), "Color by",
						choices = c("Locus" = "LOCI", "Species" = "SPECIES"),
						selected = "LOCI"),
					shiny::radioButtons(ns("phy_layout"), "Layout",
						choices = c("Rectangular" = "rectangular", "Circular" = "circular"),
						selected = "rectangular", inline = TRUE),
					shiny::radioButtons(ns("phy_branch_style"), "Branch lengths",
						choices = c(
							"Equal spacing (recommended)" = "cladogram",
							"Real p-distance"             = "real"
						),
						selected = "cladogram")
				)
			),
			shiny::column(3,
				shiny::div(style = "padding-top:30px;",
					shiny::actionButton(ns("btn_run_phylo"),
						"🌳 Build tree",
						icon  = shiny::icon("play"),
						class = "btn-success btn-lg",
						style = "width:100%; height:60px;"),
					shiny::br(), shiny::br(),
					shiny::tags$small(style = "color:#6B7280;",
						shiny::icon("clock"),
						" NJ (ape) on mean per-motif p-distance")
				)
			)
		),
		
		shiny::hr(),
		shiny::uiOutput(ns("phy_status")),
		
		shiny::tabsetPanel(
			shiny::tabPanel(
				shiny::tagList(shiny::icon("sitemap"), " Phylogenetic Tree"),
				shiny::br(),
				shiny::uiOutput(ns("phy_plot_ui")),
				shiny::br()
			),
			shiny::tabPanel(
				shiny::tagList(shiny::icon("table"), " Distance Matrix"),
				shiny::br(),
				DT::DTOutput(ns("phy_dist_table")),
				shiny::br()
			),
			shiny::tabPanel(
				shiny::tagList(shiny::icon("download"), " Exports"),
				shiny::br(),
				shiny::fluidRow(
					shiny::column(4,
						shiny::downloadButton(ns("dl_phy_png"), "📥 PNG — Tree",
							class = "btn-primary btn-block", style = "width:100%;")
					),
					shiny::column(4,
						shiny::downloadButton(ns("dl_phy_newick"), "📥 Newick — Tree",
							class = "btn-info btn-block", style = "width:100%;")
					),
					shiny::column(4,
						shiny::downloadButton(ns("dl_phy_csv"), "📥 CSV — Distances",
							class = "btn-warning btn-block", style = "width:100%;")
					)
				)
			)
		)
	)
}


# ══════════════════════════════════════════════════════════════
# SERVER 
# ══════════════════════════════════════════════════════════════

mod_tsne_server <- function(id, df_final) {
	shiny::moduleServer(id, function(input, output, session) {
		
		ns        <- session$ns
		phy_store <- shiny::reactiveValues(result = NULL, done = FALSE)
		
		shiny::observe({
			shiny::req(df_final())
			df        <- df_final() %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
			sp_list   <- sort(unique(df$SPECIES))
			loci_list <- sort(unique(df$LOCI))
			shiny::isolate({
				shiny::updateSelectInput(session, "phy_species",
					choices = c("All" = "", setNames(sp_list, sp_list)), selected = "")
				shiny::updateSelectInput(session, "phy_locus",
					choices = c("All" = "", setNames(loci_list, loci_list)), selected = "")
			})
		})
		
		shiny::observeEvent(input$btn_run_phylo, {
			shiny::req(df_final())
			shiny::withProgress(message = "Building phylogeny…", value = 0.2, {
				tryCatch({
					sp  <- if (is.null(input$phy_species) || input$phy_species == "") NULL else input$phy_species
					loc <- if (is.null(input$phy_locus)   || input$phy_locus   == "") NULL else input$phy_locus
					fn  <- if (is.null(input$phy_func)    || input$phy_func    == "") NULL else input$phy_func
					
					shiny::incProgress(0.3, detail = "Extracting motif sequences…")
					result <- compute_motif_phylogeny(
						df_final(), type_data = "IG",
						species     = sp, locus = loc, func_filter = fn,
						group_by    = input$phy_group_by,
						min_motifs  = input$phy_min_motifs
					)
					shiny::incProgress(0.9, detail = "NJ tree computed…")
					phy_store$result <- result
					phy_store$done   <- TRUE
				}, error = function(e) {
					shiny::showNotification(paste0("❌ Phylogeny: ", e$message),
						type = "error", duration = 8)
				})
			})
		})
		
		output$phy_status <- shiny::renderUI({
			if (!phy_store$done)
				shiny::div(class = "alert alert-info", shiny::icon("info-circle"),
					" Configure the filters and click ", shiny::strong("🌳 Build tree"), ".")
			else {
				res <- phy_store$result
				shiny::div(class = "alert alert-success", shiny::icon("check-circle"),
					shiny::strong(" Tree built — "),
					res$n_entities, " entities (", tolower(res$group_by), ") | ",
					length(res$motifs_used), " motif(s): ", paste(res$motifs_used, collapse = ", "))
			}
		})
		
		output$phy_plot_ui <- shiny::renderUI({
			shiny::req(phy_store$done, phy_store$result)
			h <- max(450, phy_store$result$n_entities * 22 + 120)
			h <- min(h, 1000)
			
			
			if (identical(input$phy_layout, "circular")) {
				w <- h + 180   
				shiny::div(
					style = "display:flex; justify-content:center;",
					shiny::plotOutput(ns("phy_plot"),
						height = paste0(h, "px"), width = paste0(w, "px"))
				)
			} else {
				shiny::plotOutput(ns("phy_plot"), height = paste0(h, "px"))
			}
		})
		
		output$phy_plot <- shiny::renderPlot({
			shiny::req(phy_store$done, phy_store$result)
			plot_motif_phylo_tree(phy_store$result,
				color_by     = input$phy_color_by,
				layout       = input$phy_layout,
				interactive  = FALSE,
				branch_style = input$phy_branch_style %||% "cladogram")
		})
		
		output$phy_dist_table <- DT::renderDT({
			shiny::req(phy_store$result)
			m <- as.matrix(phy_store$result$dist)
			df_tbl <- as.data.frame(round(m, 3)) %>% tibble::rownames_to_column("Entity")
			DT::datatable(df_tbl,
				options = list(pageLength = 20, scrollX = TRUE, dom = "Bfrtip",
					buttons = c("copy","csv","excel")),
				extensions = "Buttons", rownames = FALSE,
				class = "cell-border stripe hover compact"
			) %>%
				DT::formatStyle("Entity", fontWeight = "bold", color = "#0B2A5B")
		})
		
		output$dl_phy_png <- shiny::downloadHandler(
			filename = function() paste0("Motif_Phylogeny_",
				input$phy_group_by, "_", Sys.Date(), ".png"),
			content = function(file) {
				shiny::req(phy_store$result)
				p <- plot_motif_phylo_tree(phy_store$result,
					color_by     = input$phy_color_by,
					layout       = input$phy_layout,
					interactive  = FALSE,
					branch_style = input$phy_branch_style %||% "cladogram")
				h <- max(6, phy_store$result$n_entities * 0.25 + 2)
				h <- min(h, 50)
				
				w <- if (identical(input$phy_layout, "circular")) h + 3 else 12
				ggplot2::ggsave(file, plot = p, width = w, height = h,
					dpi = 300, bg = "white")
				shiny::showNotification("✅ PNG downloaded!", type = "message", duration = 3)
			}
		)
		
		output$dl_phy_newick <- shiny::downloadHandler(
			filename = function() paste0("Motif_Phylogeny_",
				input$phy_group_by, "_", Sys.Date(), ".nwk"),
			content = function(file) {
				shiny::req(phy_store$result)
				ape::write.tree(phy_store$result$tree, file)
				shiny::showNotification("✅ Newick downloaded!", type = "message", duration = 3)
			}
		)
		
		output$dl_phy_csv <- shiny::downloadHandler(
			filename = function() paste0("Motif_Phylogeny_Distances_",
				input$phy_group_by, "_", Sys.Date(), ".csv"),
			content = function(file) {
				shiny::req(phy_store$result)
				m <- as.matrix(phy_store$result$dist)
				write.csv(as.data.frame(m) %>% tibble::rownames_to_column("Entity"),
					file, row.names = FALSE)
				shiny::showNotification("✅ CSV downloaded!", type = "message", duration = 3)
			}
		)
	})
}