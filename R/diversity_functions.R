# ═══════════════════════════════════════════════════════════════
# Analyse de Diversity des Motifs
# ═══════════════════════════════════════════════════════════════

# Dépend de : utils_core.R
# Packages  : tidyverse, vegan, plotly, ggrepel

library(tidyverse)
library(vegan)

# ══════════════════════════════════════════════════════════════
#  CALCUL DES INDICES DE DIVERSITÉ
# ══════════════════════════════════════════════════════════════

calculate_diversity_indices <- function(df, type_data = "IG", grouping = "species",
	species_filter = NULL, locus_filter = NULL) {
	
	if (type_data == "IG") {
		df <- df %>% filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	} else {
		df <- df %>% filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	}
	
	# Filtre espèce 
	if (!is.null(species_filter) && length(species_filter) == 1 &&
			species_filter != "" && species_filter != "All") {
		df <- df %>% filter(SPECIES == species_filter)
	}
	
	# Filtre locus
	if (!is.null(locus_filter) && locus_filter != "" && locus_filter != "All") {
		df <- df %>% filter(LOCI == locus_filter)
	}
	
	if (nrow(df) == 0) {
		warning("calculate_diversity_indices: no data after filtering.")
		return(NULL)
	}
	
	bin_cols <- names(df)[str_starts(names(df), "BIN_")]
	if (length(bin_cols) == 0) stop("No BIN_* column found")
	
	if (grouping == "species") {
		group_vars <- "SPECIES"
	} else if (grouping == "locus") {
		group_vars <- c("LOCI", "SPECIES")
	} else {
		group_vars <- c("SUBGROUP", "SPECIES")
	}
	
	result <- df %>%
		group_by(across(all_of(group_vars))) %>%
		summarise(
			N_sequences = n(),
			Richness    = sum(across(all_of(bin_cols), ~ sum(.x == 1, na.rm = TRUE)) > 0),
			.groups = "drop"
		)
	
	diversity_list <- df %>%
		group_by(across(all_of(group_vars))) %>%
		group_split()
	
	diversity_values <- lapply(diversity_list, function(subdf) {
		counts <- sapply(bin_cols, function(col) sum(subdf[[col]] == 1, na.rm = TRUE))
		counts <- counts[counts > 0]
		if (length(counts) == 0) {
			return(data.frame(Shannon = 0, Simpson = 0, Evenness = 0))
		}
		shannon  <- vegan::diversity(counts, "shannon")
		simpson  <- vegan::diversity(counts, "simpson")
		richness <- length(counts)
		evenness <- if (richness > 1) shannon / log(richness) else 0
		data.frame(Shannon = round(shannon, 4), Simpson = round(simpson, 4),
			Evenness = round(evenness, 4))
	})
	
	diversity_df <- bind_rows(diversity_values)
	result <- cbind(result, diversity_df)
	return(result)
}


# ══════════════════════════════════════════════════════════════
# TESTS STATISTIQUES 
# ══════════════════════════════════════════════════════════════

# Tests Kruskal-Wallis  
diversity_statistical_tests <- function(df, type_data = "IG", grouping = "species",
	species_filter = NULL, locus_filter = NULL) {
	
	if (type_data == "IG") {
		df <- df %>% filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	} else {
		df <- df %>% filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	}
	
	if (!is.null(species_filter) && species_filter != "" && species_filter != "All")
		df <- df %>% filter(SPECIES == species_filter)
	if (!is.null(locus_filter) && locus_filter != "" && locus_filter != "All")
		df <- df %>% filter(LOCI == locus_filter)
	
	if (nrow(df) == 0) return(NULL)
	
	bin_cols <- names(df)[str_starts(names(df), "BIN_")]
	if (length(bin_cols) == 0) return(NULL)
	
	group_col <- switch(grouping,
		species  = "SPECIES",
		locus    = "LOCI",
		subgroup = "SUBGROUP"
	)
	
	if (!group_col %in% names(df)) return(NULL)
	n_groups <- length(unique(na.omit(df[[group_col]])))
	
	if (n_groups < 2) {
		return(data.frame(
			Variable    = "—",
			Test        = "—",
			N_sequences = nrow(df),
			Statistic   = NA_real_,
			df_test     = NA_real_,
			p.value     = NA_real_,
			Significant = paste0("Only 1 group (", group_col, " = 1 level)")
		))
	}
	
	# Calcul de la richesse 
	df_seq <- df %>%
		rowwise() %>%
		mutate(
			N_motifs_present = sum(c_across(all_of(bin_cols)) == 1, na.rm = TRUE),
			Pct_motifs       = N_motifs_present / length(bin_cols) * 100
		) %>%
		ungroup()
	
	results <- list()
	
	# Test global 
	for (var in c("N_motifs_present", "Pct_motifs")) {
		label <- if (var == "N_motifs_present") "N motifs per sequence" else "% motifs per sequence"
		test <- tryCatch(
			kruskal.test(as.formula(paste(var, "~", group_col)), data = df_seq),
			error = function(e) NULL
		)
		if (!is.null(test)) {
			results[[var]] <- data.frame(
				Variable    = label,
				Test        = "Kruskal-Wallis",
				N_sequences = nrow(df_seq),
				Statistic   = round(test$statistic, 3),
				df_test     = test$parameter,
				p.value     = round(test$p.value, 4),
				Significant = ifelse(test$p.value < 0.001, "***",
					ifelse(test$p.value < 0.01,  "**",
						ifelse(test$p.value < 0.05,  "*", "ns")))
			)
		}
	}
	
	# Tests par motif individuel
	for (bin_col in bin_cols) {
		motif_name <- str_remove(bin_col, "^BIN_")
		test <- tryCatch(
			kruskal.test(as.formula(paste(bin_col, "~", group_col)), data = df_seq),
			error = function(e) NULL
		)
		if (!is.null(test)) {
			results[[motif_name]] <- data.frame(
				Variable    = paste0("Motif — ", motif_name),
				Test        = "Kruskal-Wallis",
				N_sequences = nrow(df_seq),
				Statistic   = round(test$statistic, 3),
				df_test     = test$parameter,
				p.value     = round(test$p.value, 4),
				Significant = ifelse(test$p.value < 0.001, "***",
					ifelse(test$p.value < 0.01,  "**",
						ifelse(test$p.value < 0.05,  "*", "ns")))
			)
		}
	}
	
	if (length(results) > 0) bind_rows(results) else NULL
}


# ══════════════════════════════════════════════════════════════
#  VISUALISATIONS
# ══════════════════════════════════════════════════════════════

plot_diversity_barplot <- function(data, grouping = "species", interactive = TRUE,
	index_filter = "all") {
	
	if (is.null(data) || nrow(data) == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "No data available.", size = 6, color = "grey50") +
				ggplot2::theme_void())
	}
	
	# Indices disponibles et sélectionnés
	all_indices <- c("Shannon", "Simpson", "Richness", "Evenness")
	indices_use <- if (!is.null(index_filter) && index_filter != "all" &&
			index_filter %in% all_indices) index_filter else all_indices
	single_idx  <- length(indices_use) == 1
	
	if (grouping == "species") {
		data_long <- data %>%
			select(SPECIES, all_of(all_indices)) %>%
			pivot_longer(cols = all_of(all_indices), names_to = "Index", values_to = "Value") %>%
			filter(Index %in% indices_use)
		sp_colors <- get_species_palette(unique(data_long$SPECIES))
		
		if (single_idx) {
			p <- ggplot(data_long, aes(x = SPECIES, y = Value, fill = SPECIES)) +
				geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
				geom_text(aes(label = round(Value, 3)),
					vjust = -0.4, size = 3.5, color = "grey20", fontface = "bold") +
				scale_fill_manual(values = sp_colors) +
				scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
				labs(title = paste0(indices_use, " by species"), x = NULL, y = indices_use) +
				theme_imgt(base_size = 12) +
				theme(axis.text.x = element_text(angle = 35, hjust = 1))
		} else {
			p <- ggplot(data_long, aes(x = SPECIES, y = Value, fill = SPECIES)) +
				geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
				geom_text(aes(label = round(Value, 2)), vjust = -0.4, size = 3, color = "grey20") +
				facet_wrap(~ Index, scales = "free_y") +
				scale_fill_manual(values = sp_colors) +
				scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
				labs(title = "Diversity indices by species", x = NULL, y = "Value") +
				theme_imgt(base_size = 12) +
				theme(axis.text.x = element_text(angle = 35, hjust = 1))
		}
		
	} else if (grouping == "locus") {
		data_long <- data %>%
			select(LOCI, SPECIES, all_of(all_indices)) %>%
			pivot_longer(cols = all_of(all_indices), names_to = "Index", values_to = "Value") %>%
			filter(Index %in% indices_use) %>%
			mutate(Label = paste0(SPECIES, "\n[", LOCI, "]"))
		loci_unique <- unique(data_long$LOCI)
		loci_pal    <- LOCI_COLORS[loci_unique]
		loci_pal[is.na(loci_pal)] <- "#999999"
		
		if (single_idx) {
			p <- ggplot(data_long, aes(x = Label, y = Value, fill = LOCI)) +
				geom_bar(stat = "identity", width = 0.7) +
				geom_text(aes(label = round(Value, 3)),
					vjust = -0.4, size = 3.2, color = "grey20", fontface = "bold") +
				scale_fill_manual(values = loci_pal, name = "Locus") +
				scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
				labs(title = paste0(indices_use, " by locus"), x = NULL, y = indices_use) +
				theme_imgt(base_size = 12) +
				theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9))
		} else {
			p <- ggplot(data_long, aes(x = Label, y = Value, fill = LOCI)) +
				geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
				geom_text(aes(label = round(Value, 2)),
					position = position_dodge(width = 0.8),
					vjust = -0.4, size = 2.8, color = "grey20") +
				facet_wrap(~ Index, scales = "free_y") +
				scale_fill_manual(values = loci_pal, name = "Locus") +
				scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
				labs(title = "Diversity indices by locus", x = NULL, y = "Value") +
				theme_imgt(base_size = 12) +
				theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9))
		}
		
	} else {
		n_subgrp  <- length(unique(data$SUBGROUP))
		txt_size  <- max(5, min(8, 100 / n_subgrp))
		sp_colors <- get_species_palette(unique(data$SPECIES))
		
		subgrp_order <- data %>% arrange(Shannon) %>% pull(SUBGROUP)
		data_ordered <- data %>%
			mutate(SUBGROUP = factor(SUBGROUP, levels = subgrp_order))
		
		data_long <- data_ordered %>%
			select(SUBGROUP, SPECIES, all_of(all_indices)) %>%
			pivot_longer(cols = all_of(all_indices), names_to = "Index", values_to = "Value") %>%
			filter(Index %in% indices_use)
		
		if (single_idx) {
			max_val <- max(data_long$Value, na.rm = TRUE)
			p <- ggplot(data_long, aes(x = SUBGROUP, y = Value, fill = SPECIES)) +
				geom_bar(stat = "identity", width = 0.80) +
				geom_text(
					aes(label = ifelse(Value > 0, round(Value, 3), "")),
					hjust = -0.15, size = 2.5, color = "grey20", fontface = "bold") +
				coord_flip() +
				scale_fill_manual(values = sp_colors) +
				scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
				labs(title    = paste0(indices_use, " by subgroup"),
					subtitle = "Sorted by index value (ascending)",
					x = NULL, y = indices_use, fill = "Species") +
				theme_imgt(base_size = 9) +
				theme(axis.text.y = element_text(size = txt_size), legend.position = "right")
		} else {
			p <- ggplot(data_long, aes(x = SUBGROUP, y = Value, fill = SPECIES)) +
				geom_bar(stat = "identity", width = 0.80) +
				geom_text(
					aes(label = ifelse(Value > 0, round(Value, 2), "")),
					hjust = -0.1, size = 1.8, color = "grey20") +
				coord_flip() +
				facet_wrap(~ Index, scales = "free_x", nrow = 1) +
				scale_fill_manual(values = sp_colors) +
				scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
				labs(title    = "Diversity indices by subgroup",
					subtitle = "Tip: select a single index above for detailed values",
					x = NULL, y = "Value", fill = "Species") +
				theme_imgt(base_size = 8) +
				theme(
					axis.text.y     = element_text(size = txt_size),
					legend.position = "bottom",
					strip.text      = element_text(size = 9, face = "bold"),
					panel.spacing   = unit(0.8, "lines")
				)
		}
	}
	
	if (interactive) {
		if (grouping == "subgroup") {
			h <- if (single_idx)
				min(max(400, n_subgrp * 22), 1000)
			else
				min(max(500, n_subgrp * 18), 1200)
			p <- ggplotly(p, tooltip = c("x", "y", "fill"), height = h)
		} else {
			h <- if (single_idx) 480 else 560
			p <- ggplotly(p, tooltip = c("x", "y", "fill"), height = h)
		}
	}
	return(p)
}

#  Radar chart
plot_diversity_radar <- function(data, grouping = "species") {
	
	if (is.null(data) || nrow(data) == 0)
		return(plotly::plot_ly() %>% plotly::layout(title = "No data available"))
	
	data_norm <- data %>%
		mutate(
			Shannon_norm  = if (max(Shannon,  na.rm=TRUE) > 0) Shannon  / max(Shannon,  na.rm=TRUE) else 0,
			Simpson_norm  = if (max(Simpson,  na.rm=TRUE) > 0) Simpson  / max(Simpson,  na.rm=TRUE) else 0,
			Richness_norm = if (max(Richness, na.rm=TRUE) > 0) Richness / max(Richness, na.rm=TRUE) else 0,
			Evenness_norm = if (max(Evenness, na.rm=TRUE) > 0) Evenness / max(Evenness, na.rm=TRUE) else 0
		)
	
	fig <- plot_ly(type = 'scatterpolar', fill = 'toself', mode = 'lines+markers')
	
	n_show <- if (grouping == "subgroup") min(nrow(data_norm), 12) else nrow(data_norm)
	
	for (i in seq_len(n_show)) {
		name_label <- switch(grouping,
			species  = data_norm$SPECIES[i],
			locus    = paste(data_norm$SPECIES[i], data_norm$LOCI[i], sep = " - "),
			subgroup = paste(data_norm$SPECIES[i], data_norm$SUBGROUP[i], sep = " - ")
		)
		fig <- fig %>% add_trace(
			r     = c(data_norm$Shannon_norm[i], data_norm$Simpson_norm[i],
				data_norm$Richness_norm[i], data_norm$Evenness_norm[i]),
			theta = c('Shannon', 'Simpson', 'Richness', 'Evenness'),
			name  = name_label
		)
	}
	
	fig %>% layout(
		polar = list(radialaxis = list(visible = TRUE, range = c(0, 1))),
		title = list(text = "Normalized diversity profile (radar)",
			font = list(size = 14))
	)
}


