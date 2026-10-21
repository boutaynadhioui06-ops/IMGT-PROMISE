# ═══════════════════════════════════════════════════════════════
# Analyse Statistique Descriptive
# ═══════════════════════════════════════════════════════════════

# Dépend de : utils_core.R
# Packages  : tidyverse, DT, scales, plotly

library(tidyverse)
library(DT)
library(scales)

# ══════════════════════════════════════════════════════════════
# TABLEAU DESCRIPTIF COMPLET 
# ══════════════════════════════════════════════════════════════

generate_descriptive_table_species <- function(df, type_data = "IG", species) {
	
	# ── Filtre type + espèce
	if (type_data == "IG") {
		df <- df %>% dplyr::filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	} else {
		df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	}
	df <- df %>% dplyr::filter(SPECIES == species)
	
	if (nrow(df) == 0) stop(paste("No sequence for species:", species))
	
	# Pour TR : pas de motifs → on n'inclut ni Taux ni Dist ni BIN colonnes
	bin_cols  <- if (type_data == "IG") names(df)[stringr::str_starts(names(df), "BIN_")]  else character(0)
	dist_cols <- if (type_data == "IG") names(df)[stringr::str_starts(names(df), "DIST_")] else character(0)
	
	id_base_cols <- intersect(c("SPECIES", "LOCI", "SUBGROUP", "GENE", "ALLELE"), names(df))
	
	# Bloc Identification + Effectifs + Taille 5'UTR
	result <- df %>%
		dplyr::group_by(SPECIES, LOCI, SUBGROUP) %>%
		dplyr::summarise(
			Nb_alleles      = dplyr::n_distinct(ALLELE),
			Fonctionnels    = sum(FUNCTIONALITY == "F",   na.rm = TRUE),
			Pseudogenes     = sum(FUNCTIONALITY == "P",   na.rm = TRUE),
			ORF             = sum(FUNCTIONALITY == "ORF", na.rm = TRUE),
			Taille_5UTR_moy = as.integer(round(mean(LENGTH_NT,   na.rm = TRUE), 0)),
			Taille_5UTR_med = as.integer(round(median(LENGTH_NT, na.rm = TRUE), 0)),
			.groups = "drop"
		)
	
	# ── Taux de présence des motifs (%) 
	if (type_data == "IG" && length(bin_cols) > 0) {
		taux_df <- df %>%
			dplyr::group_by(SPECIES, LOCI, SUBGROUP) %>%
			dplyr::summarise(
				dplyr::across(
					dplyr::all_of(bin_cols),
					~ ifelse(sum(!is.na(.x)) == 0, NA_real_,
						round(mean(.x == 1, na.rm = TRUE) * 100, 1)),
					.names = "Taux_{.col}"
				),
				.groups = "drop"
			) %>%
			dplyr::rename_with(
				~ stringr::str_replace(.x, "Taux_BIN_", "Taux_"),
				dplyr::starts_with("Taux_BIN_")
			)
		
		result <- dplyr::left_join(result, taux_df, by = c("SPECIES", "LOCI", "SUBGROUP"))
	}
	
	# ── Distances moyennes au LPART1 
	if (type_data == "IG" && length(dist_cols) > 0) {
		dist_df <- df %>%
			dplyr::group_by(SPECIES, LOCI, SUBGROUP) %>%
			dplyr::summarise(
				dplyr::across(
					dplyr::all_of(dist_cols),
					~ ifelse(all(is.na(.x)), NA_integer_,
						as.integer(round(mean(.x, na.rm = TRUE), 0))),
					.names = "Mean_{.col}"
				),
				.groups = "drop"
			) %>%
			dplyr::rename_with(
				~ stringr::str_replace(.x, "Mean_DIST_", "Dist_"),
				dplyr::starts_with("Mean_DIST_")
			)
		
		result <- dplyr::left_join(result, dist_df, by = c("SPECIES", "LOCI", "SUBGROUP"))
	}
	
	# ── Ordre des colonnes
	id_cols      <- c("SPECIES", "LOCI", "SUBGROUP")
	eff_cols     <- c("Nb_alleles", "Fonctionnels", "Pseudogenes", "ORF")
	taux_present <- names(result)[stringr::str_starts(names(result), "Taux_")]
	dist_present <- names(result)[stringr::str_starts(names(result), "Dist_")]
	taille_cols  <- c("Taille_5UTR_moy", "Taille_5UTR_med")
	
	col_order <- c(id_cols, eff_cols, taux_present, dist_present, taille_cols)
	col_order <- col_order[col_order %in% names(result)]
	
	result <- result %>%
		dplyr::select(dplyr::all_of(col_order)) %>%
		dplyr::arrange(LOCI, SUBGROUP)
	
	return(result)
}

get_tableau_col_groups <- function(df) {
	list(
		id      = intersect(c("SPECIES", "LOCI", "SUBGROUP"), names(df)),
		eff     = intersect(c("Nb_alleles", "Fonctionnels", "Pseudogenes", "ORF"), names(df)),
		taux    = names(df)[stringr::str_starts(names(df), "Taux_")],
		dist    = names(df)[stringr::str_starts(names(df), "Dist_")],
		taille  = intersect(c("Taille_5UTR_moy", "Taille_5UTR_med"), names(df))
	)
}
# ══════════════════════════════════════════════════════════════
# SEQUENCE DISTRIBUTION
# ══════════════════════════════════════════════════════════════

plot_seq_by_species <- function(df, type_data = "IG", interactive = TRUE) {
	if (type_data == "IG") df <- df %>% filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	else df <- df %>% filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	
	data_plot <- df %>% group_by(SPECIES) %>%
		summarise(Count = n(), .groups = "drop") %>% arrange(desc(Count))
	sp_colors <- get_species_palette(data_plot$SPECIES)
	n_sp <- nrow(data_plot)
	
	p <- ggplot(data_plot, aes(x = reorder(SPECIES, Count), y = Count, fill = SPECIES)) +
		geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
		geom_text(aes(label = scales::comma(Count)),
			hjust = -0.15, size = 4.5, fontface = "bold", color = "grey20") +
		coord_flip() +
		scale_fill_manual(values = sp_colors) +
		scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
		labs(title = paste("Total sequences per species —", type_data),
			x = NULL, y = "Number of sequences") +
		theme_imgt(base_size = 13) +
		theme(
			axis.text.y = element_text(size = 12, face = "bold"),
			plot.margin = margin(10, 55, 10, 10)
		)
	
	h <- max(400, n_sp * 55)
	if (interactive) p <- ggplotly(p, tooltip = c("x","y"), height = min(h, 700)) %>%
		plotly::layout(margin = list(r = 110))
	return(p)
}


plot_seq_by_locus <- function(df, type_data = "IG", interactive = TRUE) {
	if (type_data == "IG") df <- df %>% filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	else df <- df %>% filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	
	seq_loci <- df %>% group_by(LOCI, SPECIES) %>% summarise(Count = n(), .groups = "drop")
	n_sp  <- length(unique(seq_loci$SPECIES))
	n_col  <- if (n_sp <= 3) 1 else 2
	n_rows <- ceiling(n_sp / n_col)
	
	p <- ggplot(seq_loci, aes(x = LOCI, y = Count, fill = LOCI)) +
		geom_bar(stat = "identity", width = 0.50) +
		geom_text(aes(label = scales::comma(Count)),
			vjust = -0.5, size = 4.5, fontface = "bold", color = "grey20") +
		scale_fill_manual(values = LOCI_COLORS, guide = "none") +
		scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
		labs(title = paste("Sequence distribution by locus —", type_data),
			subtitle = if (n_sp == 1) paste("Species:", unique(seq_loci$SPECIES))
			else paste(n_sp, "species — each panel = one species"),
			x = "Locus", y = "Number of sequences") +
		facet_wrap(~ SPECIES, scales = "free_y", ncol = n_col) +
		theme_imgt(base_size = 12) +
		theme(
			strip.text    = element_text(size = 12, face = "bold"),
			axis.text.x   = element_text(size = 11, face = "bold"),
			axis.text.y   = element_text(size = 10),
			panel.spacing = unit(1.8, "lines"),
			plot.margin   = margin(10, 20, 10, 10)
		)
		h <- if (n_sp == 1) 420 else max(420, n_rows * 340)
	if (interactive) p <- ggplotly(p, tooltip = c("x","y"), height = h)
	return(p)
}




# ══════════════════════════════════════════════════════════════
# FUNCTIONALITY DISTRIBUTION
# ══════════════════════════════════════════════════════════════

plot_func_by_species <- function(df, type_data = "IG", interactive = TRUE) {
	if (type_data == "IG") df <- df %>% filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	else df <- df %>% filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	
	data_plot <- df %>% group_by(SPECIES, FUNCTIONALITY) %>%
		summarise(n = n(), .groups = "drop") %>% group_by(SPECIES) %>%
		mutate(Total = sum(n), Proportion = n / Total) %>% ungroup() %>%
		filter(!is.na(FUNCTIONALITY))
	n_sp <- length(unique(data_plot$SPECIES))
	sp_order <- data_plot %>% group_by(SPECIES) %>%
		summarise(tot = sum(n), .groups = "drop") %>%
		arrange(tot) %>% pull(SPECIES)
	data_plot$SPECIES <- factor(data_plot$SPECIES, levels = sp_order)
	
	p <- ggplot(data_plot, aes(x = SPECIES, y = Proportion, fill = FUNCTIONALITY)) +
		geom_bar(stat = "identity", position = "stack", width = 0.65) +
		geom_text(aes(
			label = ifelse(n >= 3, n, ""),
			color = FUNCTIONALITY),
			position = position_stack(vjust = 0.5),
			size = 4.5, fontface = "bold", show.legend = FALSE) +
		coord_flip() +
		scale_color_manual(values = FUNC_LABEL_COLORS, guide = "none") +
		scale_fill_manual(values = FUNC_COLORS, breaks = c("F","P","ORF"),
			labels = c("Functional (F)", "Pseudogene (P)", "ORF")) +
		scale_y_continuous(labels = percent_format(accuracy = 1),
			expand = expansion(mult = c(0, 0.04))) +
		labs(title = paste("Functionality distribution by species —", type_data),
			x = NULL, y = "Proportion", fill = "Gene type") +
		theme_imgt(base_size = 13) +
		theme(axis.text.y = element_text(size = 12, face = "bold"))
	
	h <- max(450, n_sp * 55)
	if (interactive) p <- ggplotly(p, tooltip = c("x","y","fill"), height = min(h, 700)) %>%
		plotly::layout(
			legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.12),
			margin = list(b = 60, l = 130, r = 40, t = 60)
		)
	return(p)
}


plot_func_by_locus <- function(df, type_data = "IG", interactive = TRUE) {
	if (type_data == "IG") df <- df %>% filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	else df <- df %>% filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	
	n_sp <- length(unique(df$SPECIES))
	n_col  <- if (n_sp <= 3) 1 else 2
	n_rows <- ceiling(n_sp / n_col)
	
	data_plot <- df %>% group_by(LOCI, SPECIES, FUNCTIONALITY) %>%
		summarise(n = n(), .groups = "drop") %>% group_by(LOCI, SPECIES) %>%
		mutate(Total = sum(n), Proportion = n / Total) %>% ungroup() %>%
		filter(!is.na(FUNCTIONALITY))
	
	p <- ggplot(data_plot, aes(x = LOCI, y = Proportion, fill = FUNCTIONALITY)) +
		geom_bar(stat = "identity", position = "stack", width = 0.50) +
		geom_text(aes(
			label = ifelse(n >= 3, n, ""),
			color = FUNCTIONALITY),
			position = position_stack(vjust = 0.5),
			size = 4.5, fontface = "bold", show.legend = FALSE) +
		scale_color_manual(values = FUNC_LABEL_COLORS, guide = "none") +
		scale_y_continuous(labels = percent_format(accuracy = 1),
			expand = expansion(mult = c(0, 0.04))) +
		scale_fill_manual(values = FUNC_COLORS, breaks = c("F","P","ORF"),
			labels = c("Functional (F)","Pseudogene (P)","ORF")) +
		labs(title = paste("Functionality distribution by locus —", type_data),
			subtitle = if (n_sp == 1) paste("Species:", unique(df$SPECIES))
			else paste(n_sp, "species — each panel = one species"),
			x = "Locus", y = "Proportion", fill = "Gene type") +
		facet_wrap(~ SPECIES, ncol = n_col, scales = "free_x") +
		theme_imgt(base_size = 12) +
		theme(
			strip.text    = element_text(size = 12, face = "bold"),
			axis.text.x   = element_text(size = 11, face = "bold"),
			axis.text.y   = element_text(size = 10),
			panel.spacing = unit(1.8, "lines"),
			plot.margin   = margin(10, 20, 10, 10)
		)
	
	h <- if (n_sp == 1) 420 else max(420, n_rows * 340)
	if (interactive) p <- ggplotly(p, tooltip = c("x","y","fill"), height = h + 60) %>%
		plotly::layout(
			legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.15),
			margin = list(b = 90)
		)
	return(p)
}


