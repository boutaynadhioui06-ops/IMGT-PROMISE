# ═══════════════════════════════════════════════════════════════
# Analyse de Conservation des Motifs
# ═══════════════════════════════════════════════════════════════

# Dépend de : utils_core.R, subgroup_display.R
# Packages  : tidyverse, scales, broom, plotly

library(tidyverse)
library(ggplot2)
library(scales)
library(broom)


# ══════════════════════════════════════════════════════════════
# FRÉQUENCES DE PRÉSENCE DES MOTIFS
# ══════════════════════════════════════════════════════════════

calculate_motif_frequencies <- function(df, grouping = "species", type_data= "IG") {
	
	# Conservation des motifs : IG
	df <- df %>% filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	
	bin_cols <- names(df)[str_starts(names(df), "BIN_")]
	
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
			N_total = n(),
			across(all_of(bin_cols),
				~ sum(.x == 1, na.rm = TRUE),
				.names = "N_{.col}"),
			.groups = "drop"
		)
	
	for (col in bin_cols) {
		n_col   <- paste0("N_", col)
		pct_col <- paste0("Pct_", str_remove(col, "^BIN_"))
		result[[pct_col]] <- round(result[[n_col]] / result$N_total * 100, 1)
	}
	
	return(result)
}


plot_motif_frequencies <- function(data, grouping = "species", interactive = TRUE, type_data = "IG") {
	
	pct_cols <- names(data)[str_starts(names(data), "Pct_")]
	
	if (grouping == "species") {
		
		data_long <- data %>%
			select(SPECIES, all_of(pct_cols)) %>%
			pivot_longer(cols = all_of(pct_cols),
				names_to = "Motif", values_to = "Percentage") %>%
			mutate(Motif = str_remove(Motif, "^Pct_"))
		
		motifs_present <- unique(data_long$Motif)
		motif_pal      <- get_motif_palette(motifs_present)
		n_motifs       <- length(motifs_present)
		
		p <- ggplot(data_long,
			aes(x = SPECIES, y = Percentage, fill = Motif)) +
			geom_bar(stat = "identity",
				position = position_dodge(width = 0.85), width = 0.75) +
			{ if (n_motifs <= 5)
				geom_text(aes(label = ifelse(Percentage > 0, paste0(Percentage, "%"), "")),
					position = position_dodge(width = 0.85),
					vjust = -0.4, size = 2.8, fontface = "bold", color = "grey20")
				else
					geom_blank()
			} +
			scale_fill_manual(values = motif_pal, name = "Motif") +
			scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
			labs(title    = "Motif presence frequencies by species",
				x = NULL, y = "Presence (%)") +
			theme_imgt() +
			theme(axis.text.x = element_text(angle = 35, hjust = 1))
		
		n_sp     <- length(unique(data_long$SPECIES))
		n_subgrp <- 0
		
	} else if (grouping == "locus") {
		
		data_long <- data %>%
			select(LOCI, SPECIES, all_of(pct_cols)) %>%
			pivot_longer(cols = all_of(pct_cols),
				names_to = "Motif", values_to = "Percentage") %>%
			mutate(Motif = str_remove(Motif, "^Pct_"))
		
		if (nrow(data_long) == 0) {
			return(ggplot2::ggplot() +
					ggplot2::annotate("text", x = 0.5, y = 0.5,
						label = "No data available for locus view.", size = 5, color = "grey50") +
					ggplot2::theme_void())
		}
		
		motifs_present <- unique(data_long$Motif)
		motif_pal      <- get_motif_palette(motifs_present)
		n_sp           <- length(unique(data_long$SPECIES))
		n_loci         <- length(unique(data_long$LOCI))
		n_subgrp       <- 0
		
		
		data_long <- data_long %>%
			mutate(
				SPECIES = as.character(SPECIES),
				LOCI    = as.character(LOCI),
				Panel   = factor(paste0(SPECIES, "\n[", LOCI, "]"))
			)
		
		n_panels     <- length(unique(data_long$Panel))
		n_cols_facet <- min(n_loci, 3)
		
		p <- ggplot(data_long,
			aes(x = Motif, y = Percentage, fill = Motif)) +
			geom_bar(stat = "identity", width = 0.72) +
			geom_text(aes(label = ifelse(Percentage > 3, paste0(round(Percentage, 0), "%"), "")),
				vjust = -0.4, size = 2.5, fontface = "bold", color = "grey20") +
			scale_fill_manual(values = motif_pal, guide = "none") +
			scale_y_continuous(expand = expansion(mult = c(0, 0.22)), limits = c(0, 108)) +
			facet_wrap(~ Panel, ncol = n_cols_facet, scales = "free_x") +
			labs(title    = "Motif presence frequencies by locus",
				subtitle = paste0(n_sp, " species | ", n_loci, " loci | Each panel = species × locus"),
				x = "Motif", y = "Presence (%)") +
			theme_imgt(base_size = 11) +
			theme(
				axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
				strip.text       = element_text(face = "bold", color = "white", size = 8,
					lineheight = 1.3, margin = margin(t = 6, b = 6)),
				strip.background = element_rect(fill = "#1565C0", color = NA),
				panel.spacing    = unit(1.2, "lines")
			)
		
	} else {  
		
		data_long <- data %>%
			select(SUBGROUP, SPECIES, all_of(pct_cols)) %>%
			pivot_longer(cols = all_of(pct_cols),
				names_to = "Motif", values_to = "Percentage") %>%
			mutate(Motif = str_remove(Motif, "^Pct_"))
		
		motifs_present <- unique(data_long$Motif)
		motif_pal      <- get_motif_palette(motifs_present)
		n_subgrp       <- length(unique(data_long$SUBGROUP))
		n_motifs       <- length(motifs_present)
		sp_label       <- if (length(unique(data_long$SPECIES)) == 1) unique(data_long$SPECIES) else ""
		
		subgrp_order <- data_long %>%
			group_by(SUBGROUP) %>%
			summarise(total = sum(Percentage, na.rm = TRUE), .groups = "drop") %>%
			arrange(total) %>%
			pull(SUBGROUP)
		data_long$SUBGROUP <- factor(data_long$SUBGROUP, levels = subgrp_order)
		
		txt_size <- max(6, min(9, 120 / n_subgrp))
		
		p <- ggplot(data_long,
			aes(x = SUBGROUP, y = Percentage, fill = Motif)) +
			geom_bar(stat = "identity", position = "stack", width = 0.80) +
			geom_text(aes(label = ifelse(Percentage > 8, paste0(round(Percentage, 0), "%"), "")),
				position = position_stack(vjust = 0.5),
				size = 2.5, fontface = "bold", color = "white") +
			coord_flip() +
			scale_fill_manual(values = motif_pal, name = "Motif") +
			scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
			labs(title    = "Motif presence frequencies by subgroup",
				subtitle = if (sp_label != "")
					paste("Species:", sp_label, "| Stacked = cumulative % per subgroup")
				else "Each color = one motif",
				x = NULL, y = "Cumulative presence (%)") +
			theme_imgt(base_size = 10) +
			theme(
				axis.text.y     = element_text(size = txt_size),
				legend.position = "right"
			)
	}
	
	if (interactive) {
		h <- if (grouping == "subgroup") max(400, n_subgrp * 22)
		else if (grouping == "locus") min(900, max(400, n_panels * 220))
		else 520
		p <- ggplotly(p, tooltip = c("x", "y", "fill"), height = h)
		
		if (grouping %in% c("species", "subgroup")) {
			p <- p %>% plotly::layout(
				legend = list(orientation = "v", x = 1.02, y = 1, xanchor = "left"),
				margin = list(r = 140)
			)
		}
	}
	
	return(p)
}


# ══════════════════════════════════════════════════════════════
# PAIRWISE STATISTICAL TESTS BETWEEN LOCI (post-hoc)
# ══════════════════════════════════════════════════════════════

calculate_locus_pairwise_tests <- function(df) {
	df <- df %>% dplyr::filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	if (nrow(df) == 0) return(NULL)
	
	bin_cols <- names(df)[stringr::str_starts(names(df), "BIN_")]
	if (length(bin_cols) == 0) return(NULL)
	
	loci_vals <- sort(unique(na.omit(df$LOCI)))
	if (length(loci_vals) < 2) return(NULL)
	
	pairs   <- utils::combn(loci_vals, 2, simplify = FALSE)
	n_tests <- length(bin_cols) * length(pairs)
	results <- list()
	
	for (bin_col in bin_cols) {
		motif_name <- stringr::str_remove(bin_col, "^BIN_")
		
		for (pair in pairs) {
			l1 <- pair[1]; l2 <- pair[2]
			df_pair <- df %>% dplyr::filter(LOCI %in% c(l1, l2))
			
			ct <- tryCatch(table(df_pair$LOCI, df_pair[[bin_col]]),
				error = function(e) NULL)
			if (is.null(ct) || nrow(ct) < 2 || ncol(ct) < 2) next
			
			expected   <- tryCatch(chisq.test(ct)$expected, error = function(e) NULL)
			use_fisher <- !is.null(expected) && any(expected < 5, na.rm = TRUE)
			
			if (use_fisher) {
				test_res  <- tryCatch(fisher.test(ct), error = function(e) NULL)
				test_name <- "Fisher exact"
				stat_val  <- NA_real_
				pval      <- if (!is.null(test_res)) test_res$p.value else NA_real_
			} else {
				test_res  <- tryCatch(chisq.test(ct, correct = FALSE), error = function(e) NULL)
				test_name <- "Chi-square"
				stat_val  <- if (!is.null(test_res)) round(test_res$statistic, 3) else NA_real_
				pval      <- if (!is.null(test_res)) test_res$p.value else NA_real_
			}
			
			if (is.null(test_res) || is.na(pval)) next
			
			pval_adj <- min(pval * n_tests, 1)
			sig <- dplyr::case_when(
				pval_adj < 0.001 ~ "***",
				pval_adj < 0.01  ~ "**",
				pval_adj < 0.05  ~ "*",
				TRUE             ~ "ns"
			)
			
			rates <- df_pair %>%
				dplyr::group_by(LOCI) %>%
				dplyr::summarise(
					pct = round(mean(.data[[bin_col]] == 1, na.rm = TRUE) * 100, 1),
					n   = dplyr::n(), .groups = "drop"
				) %>%
				dplyr::mutate(label = paste0(LOCI, ": ", pct, "% (n=", n, ")")) %>%
				dplyr::pull(label)
			
			results[[paste0(motif_name, "_", l1, "_vs_", l2)]] <- data.frame(
				Motif       = motif_name,
				Locus_1     = l1,
				Locus_2     = l2,
				Test        = test_name,
				Statistic   = stat_val,
				p.value     = round(pval, 5),
				p.adjusted  = round(pval_adj, 5),
				Correction  = paste0("Bonferroni x", n_tests),
				Significant = sig,
				Rates       = paste(rates, collapse = " | "),
				stringsAsFactors = FALSE
			)
		}
	}
	
	if (length(results) == 0) return(NULL)
	dplyr::bind_rows(results) %>% dplyr::arrange(Motif, p.adjusted)
}


# ══════════════════════════════════════════════════════════════
# STATISTICAL TESTS FOR MOTIF PRESENCE FREQUENCIES
# ══════════════════════════════════════════════════════════════

calculate_frequency_stat_tests <- function(df, grouping = "species",
	group_filter = NULL) {
	df <- df %>% dplyr::filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	
	if (!is.null(group_filter) && group_filter != "" && group_filter != "All" &&
			group_filter != "All species") {
		if (grouping %in% c("locus", "subgroup")) {
			df <- df %>% dplyr::filter(SPECIES == group_filter)
		}
	}
	
	if (nrow(df) == 0) return(NULL)
	
	bin_cols <- names(df)[stringr::str_starts(names(df), "BIN_")]
	if (length(bin_cols) == 0) return(NULL)
	
	group_col <- switch(grouping,
		species  = "SPECIES",
		locus    = "LOCI",
		subgroup = "SUBGROUP"
	)
	if (!group_col %in% names(df)) return(NULL)
	
	group_vals <- unique(na.omit(df[[group_col]]))
	if (length(group_vals) < 2) {
		return(data.frame(
			Motif = "—", Group_variable = group_col,
			N_groups = 1, N_total = nrow(df),
			Test = "—", Statistic = NA_real_, df_test = NA_real_,
			p.value = NA_real_, p.adjusted = NA_real_,
			Significant = paste0("Only 1 level in ", group_col),
			stringsAsFactors = FALSE
		))
	}
	
	n_tests <- length(bin_cols)
	results <- list()
	
	for (bin_col in bin_cols) {
		motif_name <- stringr::str_remove(bin_col, "^BIN_")
		
		
		ct <- tryCatch(
			table(df[[group_col]], df[[bin_col]]),
			error = function(e) NULL
		)
		if (is.null(ct) || nrow(ct) < 2 || ncol(ct) < 2) next
		
		
		expected <- tryCatch(chisq.test(ct)$expected, error = function(e) NULL)
		use_fisher <- !is.null(expected) && any(expected < 5, na.rm = TRUE)
		
		if (use_fisher && nrow(ct) == 2) {
			test_res <- tryCatch(
				fisher.test(ct, simulate.p.value = TRUE, B = 2000),
				error = function(e) NULL
			)
			test_name <- "Fisher (simulated)"
			stat_val  <- NA_real_
			df_val    <- NA_real_
			pval      <- if (!is.null(test_res)) test_res$p.value else NA_real_
		} else {
			test_res <- tryCatch(
				chisq.test(ct, correct = FALSE),
				error = function(e) NULL
			)
			test_name <- "Chi-square"
			stat_val  <- if (!is.null(test_res)) round(test_res$statistic, 3) else NA_real_
			df_val    <- if (!is.null(test_res)) test_res$parameter else NA_real_
			pval      <- if (!is.null(test_res)) test_res$p.value else NA_real_
		}
		
		if (is.null(test_res) || is.na(pval)) next
		
		
		pval_adj <- min(pval * n_tests, 1)
		
		sig <- dplyr::case_when(
			pval_adj < 0.001 ~ "***",
			pval_adj < 0.01  ~ "**",
			pval_adj < 0.05  ~ "*",
			TRUE             ~ "ns"
		)
		
		
		group_rates <- df %>%
			dplyr::group_by(.data[[group_col]]) %>%
			dplyr::summarise(
				pct = round(mean(.data[[bin_col]] == 1, na.rm = TRUE) * 100, 1),
				n   = dplyr::n(),
				.groups = "drop"
			) %>%
			dplyr::mutate(label = paste0(.data[[group_col]], ": ", pct, "% (n=", n, ")")) %>%
			dplyr::pull(label)
		
		results[[motif_name]] <- data.frame(
			Motif          = motif_name,
			Group_variable = group_col,
			N_groups       = length(group_vals),
			N_total        = nrow(df),
			Test           = test_name,
			Statistic      = stat_val,
			df_test        = df_val,
			p.value        = round(pval, 5),
			p.adjusted     = round(pval_adj, 5),
			Correction     = paste0("Bonferroni ×", n_tests),
			Significant    = sig,
			Group_rates    = paste(group_rates, collapse = " | "),
			stringsAsFactors = FALSE
		)
	}
	
	if (length(results) == 0) return(NULL)
	dplyr::bind_rows(results) %>%
		dplyr::arrange(p.adjusted)
}

# ══════════════════════════════════════════════════════════════
# RÉGRESSION LOGISTIQUE
# ══════════════════════════════════════════════════════════════

logistic_regression_motif <- function(df, motif_name, type_data = "IG") {
	
	df <- df %>% dplyr::filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	
	bin_col <- paste0("BIN_", motif_name)
	if (!bin_col %in% names(df))
		stop(paste("Column", bin_col, "not found in dataset"))
	
	df_model <- df %>%
		dplyr::select(dplyr::all_of(bin_col), LOCI, SPECIES, FUNCTIONALITY) %>%
		dplyr::mutate(
			Presence      = as.factor(.data[[bin_col]]),
			LOCI          = as.factor(LOCI),
			SPECIES       = as.factor(SPECIES),
			FUNCTIONALITY = as.factor(FUNCTIONALITY)
		) %>%
		na.omit()
	df_model <- droplevels(df_model)
	facs <- c("LOCI", "SPECIES", "FUNCTIONALITY")
	facs <- facs[sapply(facs, function(f) nlevels(df_model[[f]]) >= 2)]
	if (length(facs) == 0) stop("Aucun facteur avec au moins 2 niveaux : régression impossible.")
	fml <- as.formula(paste("Presence ~", paste(facs, collapse = " + ")))
	
	if (nlevels(df_model$Presence) < 2) stop("Motif is constant (always present or always absent): regression impossible.")
	lf <- .logit_terms(df_model)
	
	model <- glm(lf$formula, data = df_model, family = binomial(link = "logit"))
	# Niveaux de référence utilisés 
	ref_loci  <- levels(df_model$LOCI)[1]
	ref_sp    <- levels(df_model$SPECIES)[1]
	ref_func  <- levels(df_model$FUNCTIONALITY)[1]
	
	results_tidy <- broom::tidy(model) %>%
		dplyr::mutate(
			OR        = exp(estimate),
			CI_lower  = exp(estimate - 1.96 * std.error),
			CI_upper  = exp(estimate + 1.96 * std.error),
			Reference = dplyr::case_when(
				stringr::str_starts(term, "LOCI")          ~ paste0("vs ", ref_loci),
				stringr::str_starts(term, "SPECIES")       ~ paste0("vs ", ref_sp),
				stringr::str_starts(term, "FUNCTIONALITY") ~ paste0("vs ", ref_func),
				term == "(Intercept)"                       ~ "Baseline (all references)",
				TRUE ~ ""
			),
			p_value_num = p.value,
			p.value   = format.pval(p.value, digits = 3),
			Significant = dplyr::case_when(
				p_value_num < 0.001 ~ "***",
				p_value_num < 0.01  ~ "**",
				p_value_num < 0.05  ~ "*",
				TRUE                 ~ "ns"
			)
		) %>%
		dplyr::select(term, Reference, estimate, std.error, OR,
			CI_lower, CI_upper, statistic, p.value, Significant)
	
	return(list(
		model      = model,
		results    = results_tidy,
		summary    = summary(model),
		ref_levels = list(LOCI = ref_loci, SPECIES = ref_sp, FUNCTIONALITY = ref_func)
	))
}


# ══════════════════════════════════════════════════════════════
# LOGOS DE SÉQUENCES
# ══════════════════════════════════════════════════════════════

prepare_sequences_for_logo <- function(df, motif_name, type_data= "IG") {
	
	
	df <- df %>% filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	
	seq_col <- paste0("SEQ_", motif_name)
	if (!seq_col %in% names(df))
		stop(paste("Column", seq_col, "not found in dataset"))
	
	df %>%
		select(SPECIES, LOCI, all_of(seq_col)) %>%
		filter(!is.na(.data[[seq_col]]), .data[[seq_col]] != "") %>%
		rename(Sequence = all_of(seq_col))
}


# ══════════════════════════════════════════════════════════════
# CONSERVATION SELON FONCTIONNALITÉ
# ══════════════════════════════════════════════════════════════

conservation_by_functionality <- function(df, grouping = "species", type_data = "IG",
	motif_filter = NULL) {
	
	
	df <- df %>% filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	
	bin_cols <- names(df)[str_starts(names(df), "BIN_")]
	
	
	if (!is.null(motif_filter) && motif_filter != "") {
		bin_col_single <- paste0("BIN_", motif_filter)
		if (bin_col_single %in% bin_cols) {
			bin_cols <- bin_col_single
		}
	}
	
	if (grouping == "species") {
		group_vars <- c("SPECIES", "FUNCTIONALITY")
	} else if (grouping == "locus") {
		group_vars <- c("LOCI", "SPECIES", "FUNCTIONALITY")
	} else {
		group_vars <- c("SUBGROUP", "SPECIES", "FUNCTIONALITY")
	}
	
	
	if (nrow(df) == 0) {
		warning("conservation_by_functionality: no data left after filtering.")
		return(list(data = data.frame(), chi_tests = NULL))
	}
	
	result <- df %>%
		group_by(across(all_of(group_vars))) %>%
		summarise(
			N_total = n(),
			across(all_of(bin_cols),
				~ sum(.x == 1, na.rm = TRUE),
				.names = "N_{.col}"),
			.groups = "drop"
		)
	
	for (col in bin_cols) {
		n_col   <- paste0("N_", col)
		pct_col <- paste0("Pct_", str_remove(col, "^BIN_"))
		result[[pct_col]] <- round(result[[n_col]] / result$N_total * 100, 1)
	}
	
	
	chi_results <- list()
	func_levels <- unique(na.omit(df$FUNCTIONALITY))
	if (length(func_levels) >= 2) {
		for (bin_col in bin_cols) {
			motif_name <- str_remove(bin_col, "^BIN_")
			contingency <- tryCatch(
				table(df$FUNCTIONALITY, df[[bin_col]]),
				error = function(e) NULL
			)
			if (is.null(contingency) || nrow(contingency) < 2 || ncol(contingency) < 2) next
			test <- tryCatch(chisq.test(contingency), error = function(e) NULL)
			if (!is.null(test)) {
				chi_results[[motif_name]] <- data.frame(
					Motif       = motif_name,
					Test        = "Chi2",
					Chi2        = test$statistic,
					df          = test$parameter,
					p.value     = test$p.value,
					Significant = ifelse(test$p.value < 0.05, "***", "ns")
				)
			}
		}
	}
	
	chi_table <- if (length(chi_results) > 0) bind_rows(chi_results) else NULL
	
	# TESTS STATISTIQUES
	
	stat_tests <- tryCatch({
		compute_conservation_stat_tests(df, bin_cols, grouping)
	}, error = function(e) {
		warning(paste0("compute_conservation_stat_tests: ", e$message))
		NULL
	})
	
	return(list(data = result, chi_tests = chi_table, stat_tests = stat_tests))
}

compute_conservation_stat_tests <- function(df, bin_cols, grouping) {
		group_col <- switch(grouping,
		species  = "FUNCTIONALITY",
		locus    = "FUNCTIONALITY",
		subgroup = "FUNCTIONALITY"
	)
	
	if (!group_col %in% names(df)) return(NULL)
	
	df_clean <- df %>%
		dplyr::filter(!is.na(.data[[group_col]]))
	
	group_vals <- unique(na.omit(df_clean[[group_col]]))
	n_groups   <- length(group_vals)
	if (n_groups < 2) return(NULL)
	
	n_total  <- nrow(df_clean)
	n_tests  <- length(bin_cols)
	results  <- list()
	
	for (bin_col in bin_cols) {
		motif_name <- str_remove(bin_col, "^BIN_")
		vals       <- df_clean[[bin_col]]
		if (all(is.na(vals)) || var(vals, na.rm = TRUE) == 0) next
		
		# Test de Kruskal-Wallis global
		kw <- tryCatch(
			kruskal.test(as.formula(paste(bin_col, "~", group_col)), data = df_clean),
			error = function(e) NULL
		)
		if (is.null(kw)) next
		
		pv_kw <- kw$p.value
		# Correction Bonferroni 
		pv_adj <- min(pv_kw * n_tests, 1)
		
		sig_label <- dplyr::case_when(
			pv_adj < 0.001 ~ "***",
			pv_adj < 0.01  ~ "**",
			pv_adj < 0.05  ~ "*",
			TRUE           ~ "ns"
		)
		
		results[[paste0("KW_", motif_name)]] <- data.frame(
			Motif         = motif_name,
			Test          = "Kruskal-Wallis",
			Comparison    = paste(group_vals, collapse = " vs "),
			N             = n_total,
			Statistic     = round(kw$statistic, 3),
			df_test       = kw$parameter,
			p.value       = round(pv_kw,  4),
			p.adjusted    = round(pv_adj, 4),
			Correction    = paste0("Bonferroni (×", n_tests, ")"),
			Effect_size   = NA_real_,
			Significant   = sig_label,
			stringsAsFactors = FALSE
		)
		
		# Comparisons pairées Mann-Whitney si KW significatif 
		if (pv_adj < 0.05 && n_groups >= 2) {
			pairs <- combn(group_vals, 2, simplify = FALSE)
			n_pairs <- length(pairs)
			
			for (pair in pairs) {
				g1 <- pair[1]; g2 <- pair[2]
				v1 <- df_clean[[bin_col]][df_clean[[group_col]] == g1]
				v2 <- df_clean[[bin_col]][df_clean[[group_col]] == g2]
				
				v1 <- v1[!is.na(v1)]; v2 <- v2[!is.na(v2)]
				if (length(v1) < 3 || length(v2) < 3) next
				if (var(v1) == 0 && var(v2) == 0) next
				
				mw <- tryCatch(
					wilcox.test(v1, v2, exact = FALSE),
					error = function(e) NULL
				)
				if (is.null(mw)) next
				
				pv_mw <- mw$p.value
				pv_mw_adj <- min(pv_mw * n_tests * n_pairs, 1)
				n_pair <- length(v1) + length(v2)
				z_val  <- tryCatch({
					qnorm(pv_mw / 2) * sign(mean(v1, na.rm=TRUE) - mean(v2, na.rm=TRUE))
				}, error = function(e) NA_real_)
				r_eff  <- if (!is.na(z_val)) round(abs(z_val) / sqrt(n_pair), 3) else NA_real_
				
				sig_mw <- dplyr::case_when(
					pv_mw_adj < 0.001 ~ "***",
					pv_mw_adj < 0.01  ~ "**",
					pv_mw_adj < 0.05  ~ "*",
					TRUE              ~ "ns"
				)
				
				results[[paste0("MW_", motif_name, "_", g1, "_vs_", g2)]] <- data.frame(
					Motif         = motif_name,
					Test          = "Mann-Whitney",
					Comparison    = paste0(g1, " vs ", g2),
					N             = n_pair,
					Statistic     = round(mw$statistic, 1),
					df_test       = NA_real_,
					p.value       = round(pv_mw,     4),
					p.adjusted    = round(pv_mw_adj, 4),
					Correction    = paste0("Bonferroni (×", n_tests * n_pairs, ")"),
					Effect_size   = r_eff,
					Significant   = sig_mw,
					stringsAsFactors = FALSE
				)
			}
		}
	}
	
	if (length(results) == 0) return(NULL)
	dplyr::bind_rows(results) %>%
		dplyr::arrange(Motif, Test, p.adjusted)
}


plot_conservation_functionality <- function(data, grouping = "species", type_data = "IG",
	interactive = TRUE, motif_filter = NULL, force_mode = "auto") {
	
	pct_cols <- names(data)[str_starts(names(data), "Pct_")]
	
	
	single_motif_mode <- (length(pct_cols) == 1) ||
		(!is.null(motif_filter) && motif_filter != "" && length(pct_cols) >= 1)
	if (single_motif_mode && is.null(motif_filter)) {
		motif_filter <- str_remove(pct_cols[1], "^Pct_")
	}
	
	# ═══════════════════════════════════════════════════════════
	# MODE MOTIF UNIQUE FUNCTIONALITY
	# ═══════════════════════════════════════════════════════════
	if (single_motif_mode) {
		
		pct_col    <- pct_cols[1]
		motif_name <- str_remove(pct_col, "^Pct_")
		func_pal   <- c("F" = "#28A745", "P" = "#DC3545", "ORF" = "#FFC107")
		
		if (grouping == "species") {
			data_plot <- data %>%
				select(SPECIES, FUNCTIONALITY, all_of(pct_col)) %>%
				rename(Percentage = all_of(pct_col)) %>%
				filter(!is.na(FUNCTIONALITY))
			
			n_sp <- length(unique(data_plot$SPECIES))
			
			p <- ggplot(data_plot,
				aes(x = SPECIES, y = Percentage, fill = FUNCTIONALITY)) +
				geom_bar(stat = "identity",
					position = position_dodge(width = 0.85), width = 0.75) +
				geom_text(aes(label = ifelse(Percentage > 0, paste0(Percentage, "%"), "")),
					position = position_dodge(width = 0.85),
					vjust = -0.4, size = 3, fontface = "bold", color = "grey20") +
				scale_fill_manual(values = func_pal,
					breaks = c("F", "P", "ORF"),
					labels = c("Functional(F)", "Pseudogene (P)", "ORF"),
					name   = "Functionality") +
				scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
				labs(title    = paste0("Conservation of motif ", motif_name, " by functionality"),
					subtitle = "Each bar = species × functionality combination",
					x = NULL, y = "Presence (%)") +
				theme_imgt(base_size = 11) +
				theme(axis.text.x = element_text(angle = 35, hjust = 1))
			
		} else if (grouping == "locus") {
			data_plot <- data %>%
				select(LOCI, SPECIES, FUNCTIONALITY, all_of(pct_col)) %>%
				rename(Percentage = all_of(pct_col)) %>%
				filter(!is.na(FUNCTIONALITY))
			
			n_sp   <- length(unique(data_plot$SPECIES))
			n_loci <- length(unique(data_plot$LOCI))
			n_col  <- if (n_sp == 1) 1 else 2
			n_rows <- ceiling(n_sp / n_col)
			
			p <- ggplot(data_plot,
				aes(x = LOCI, y = Percentage, fill = FUNCTIONALITY)) +
				geom_bar(stat = "identity",
					position = position_dodge(width = 0.85), width = 0.65) +
				geom_text(aes(label = ifelse(Percentage > 0, paste0(Percentage, "%"), "")),
					position = position_dodge(width = 0.85),
					vjust = -0.5, size = 4, fontface = "bold", color = "grey20") +
				scale_fill_manual(values = func_pal,
					breaks = c("F", "P", "ORF"),
					labels = c("Functional (F)", "Pseudogene (P)", "ORF"),
					name   = "Functionality") +
				scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
				facet_wrap(~ SPECIES, scales = "free_x", ncol = n_col) +
				labs(title    = paste0("Conservation of motif ", motif_name, " by locus — Functionality"),
					subtitle = if (n_sp == 1) paste("Species:", unique(data_plot$SPECIES))
					else "Each panel = one species | Bars colored by functionality",
					x = "Locus", y = "Presence (%)") +
				theme_imgt(base_size = 13) +
				theme(
					strip.text    = element_text(size = 13, face = "bold"),
					axis.text.x   = element_text(size = 12, face = "bold"),
					panel.spacing = unit(1.5, "lines")
				)
			
		} else {  
			data_plot <- data %>%
				select(SUBGROUP, SPECIES, FUNCTIONALITY, all_of(pct_col)) %>%
				rename(Percentage = all_of(pct_col)) %>%
				filter(!is.na(FUNCTIONALITY))
			
			n_subgrp <- length(unique(data_plot$SUBGROUP))
			sp_label <- if (length(unique(data_plot$SPECIES)) == 1) unique(data_plot$SPECIES) else ""
			
			use_heatmap <- (force_mode == "heatmap") ||
				(force_mode == "auto" && sg_use_heatmap(n_subgrp))
			
			if (use_heatmap) {
				sg_order <- data_plot %>%
					group_by(SUBGROUP) %>%
					summarise(mean_pct = mean(Percentage, na.rm = TRUE), .groups = "drop") %>%
					arrange(dplyr::desc(mean_pct)) %>% pull(SUBGROUP)
				data_plot$SUBGROUP <- factor(data_plot$SUBGROUP, levels = rev(sg_order))
				
				txt_size        <- sg_text_size(n_subgrp)
				tile_label_size <- sg_tile_label_size(n_subgrp)
				
				p <- ggplot(data_plot,
					aes(x = FUNCTIONALITY, y = SUBGROUP, fill = Percentage)) +
					geom_tile(color = "white", linewidth = 0.4) +
					geom_text(aes(
						label = ifelse(Percentage > 0, paste0(Percentage, "%"), ""),
						color = ifelse(Percentage <= 55, "white", "grey20")),
						size = tile_label_size, fontface = "bold") +
					scale_color_identity() +
					scale_fill_viridis_c(
						option = "viridis",
						limits = c(0, 100),
						name   = "Presence (%)",
						na.value = "#F9FAFB"
					) +
					scale_x_discrete(limits = c("F", "P", "ORF")) +
					labs(
						title    = paste0("Conservation of motif ", motif_name, " — Heatmap by subgroup"),
						subtitle = if (sp_label != "") paste("Species:", sp_label) else paste(n_subgrp, "subgroups"),
						x = "Functionality", y = NULL
					) +
					theme_imgt(base_size = 10) +
					theme(
						axis.text.y       = element_text(size = txt_size),
						axis.text.x       = element_text(size = 11, face = "bold"),
						panel.grid        = element_blank(),
						legend.position   = "right",
						legend.key.height = unit(1.2, "cm")
					)
				if (interactive) {
					h  <- sg_plot_height_px(n_subgrp, min_h = 420)
					ml <- sg_margin_left_px(unique(data_plot$SUBGROUP))
					p <- ggplotly(p, tooltip = c("x","y","fill"), height = h) %>%
						plotly::layout(margin = list(l = ml))
				}
				return(p)
			}
			
			subgrp_order <- data_plot %>%
				group_by(SUBGROUP) %>%
				summarise(mean_pct = mean(Percentage, na.rm = TRUE), .groups = "drop") %>%
				arrange(mean_pct) %>% pull(SUBGROUP)
			data_plot$SUBGROUP <- factor(data_plot$SUBGROUP, levels = subgrp_order)
			
			txt_size <- sg_text_size(n_subgrp)
			
			p <- ggplot(data_plot,
				aes(x = SUBGROUP, y = Percentage, fill = FUNCTIONALITY)) +
				geom_bar(stat = "identity", position = "stack", width = 0.80) +
				geom_text(aes(
					label = ifelse(Percentage > 8, paste0(round(Percentage, 0), "%"), ""),
					color = FUNCTIONALITY),
					position = position_stack(vjust = 0.5),
					size = 2.5, fontface = "bold", show.legend = FALSE) +
				scale_color_manual(values = FUNC_LABEL_COLORS, guide = "none") +
				coord_flip() +
				scale_fill_manual(values = func_pal,
					breaks = c("F", "P", "ORF"),
					labels = c("Functional (F)", "Pseudogene (P)", "ORF"),
					name   = "Functionality") +
				scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
				labs(title    = paste0("Conservation of motif ", motif_name, " by subgroup"),
					subtitle = if (sp_label != "")
						paste("Species:", sp_label, "| Stacked = cumulative % per subgroup")
					else "Bars colored by functionality",
					x = NULL, y = "Cumulative presence (%)") +
				theme_imgt(base_size = 9) +
				theme(
					axis.text.y     = element_text(size = txt_size),
					legend.position = "right"
				)
		}
		
		if (interactive) {
			if (grouping == "subgroup") {
				h  <- sg_plot_height_px(n_subgrp, min_h = 400)
				ml <- sg_margin_left_px(unique(data_plot$SUBGROUP))
				p  <- ggplotly(p, tooltip = c("x", "y", "fill"), height = h) %>%
					plotly::layout(margin = list(l = ml))
			} else if (grouping == "locus") {
				n_rows <- ceiling(n_sp / (if (n_sp == 1) 1 else 2))
				h <- if (n_sp == 1) 500 else max(500, n_rows * 500)
				p <- ggplotly(p, tooltip = c("x", "y", "fill"), height = min(h, 900))
			} else {
				h <- max(450, n_sp * 120)
				p <- ggplotly(p, tooltip = c("x", "y", "fill"), height = min(h, 650))
			}
		}
		return(p)
	}
	
	# ═══════════════════════════════════════════════════════════
	# MODE GLOBAL : taux de conservation  de tous motifs
	# ═══════════════════════════════════════════════════════════
	
	func_pal_global <- c("F" = "#28A745", "P" = "#DC3545", "ORF" = "#FFC107")
	
	if (grouping == "species") {
		
		# Taux moyen de conservation
		data_summary <- data %>%
			select(SPECIES, FUNCTIONALITY, all_of(pct_cols)) %>%
			filter(!is.na(FUNCTIONALITY)) %>%
			pivot_longer(cols = all_of(pct_cols), names_to = "Motif", values_to = "Pct") %>%
			group_by(SPECIES, FUNCTIONALITY) %>%
			summarise(MeanConservation = round(mean(Pct, na.rm = TRUE), 1), .groups = "drop")
		
		n_sp     <- length(unique(data_summary$SPECIES))
		n_subgrp <- 0
		
		p <- ggplot(data_summary,
			aes(x = SPECIES, y = MeanConservation, fill = FUNCTIONALITY)) +
			geom_bar(stat = "identity",
				position = position_dodge(width = 0.85), width = 0.75) +
			geom_text(aes(label = ifelse(MeanConservation > 0, paste0(MeanConservation, "%"), "")),
				position = position_dodge(width = 0.85),
				vjust = -0.4, size = 3, fontface = "bold", color = "grey20") +
			scale_fill_manual(values = func_pal_global,
				breaks = c("F", "P", "ORF"),
				labels = c("Functional (F)", "Pseudogene (P)", "ORF"),
				name   = "Functionality") +
			scale_y_continuous(expand = expansion(mult = c(0, 0.20)), limits = c(0, 105)) +
			labs(title    = "Global motif conservation by species and functionality",
				subtitle  = "Mean presence rate across all motifs — higher = better conserved",
				x = NULL, y = "Mean conservation rate (%)") +
			theme_imgt(base_size = 11) +
			theme(axis.text.x = element_text(angle = 35, hjust = 1))
		
	} else if (grouping == "locus") {
		
		data_summary <- data %>%
			select(LOCI, SPECIES, FUNCTIONALITY, all_of(pct_cols)) %>%
			filter(!is.na(FUNCTIONALITY)) %>%
			pivot_longer(cols = all_of(pct_cols), names_to = "Motif", values_to = "Pct") %>%
			group_by(LOCI, SPECIES, FUNCTIONALITY) %>%
			summarise(MeanConservation = round(mean(Pct, na.rm = TRUE), 1), .groups = "drop")
		
		n_sp     <- length(unique(data_summary$SPECIES))
		n_subgrp <- 0
		n_col    <- if (n_sp == 1) 1 else 2
		n_rows   <- ceiling(n_sp / n_col)
		
		p <- ggplot(data_summary,
			aes(x = LOCI, y = MeanConservation, fill = FUNCTIONALITY)) +
			geom_bar(stat = "identity",
				position = position_dodge(width = 0.85), width = 0.65) +
			geom_text(aes(label = ifelse(MeanConservation > 0, paste0(MeanConservation, "%"), "")),
				position = position_dodge(width = 0.85),
				vjust = -0.5, size = 4, fontface = "bold", color = "grey20") +
			scale_fill_manual(values = func_pal_global,
				breaks = c("F", "P", "ORF"),
				labels = c("Functional (F)", "Pseudogene (P)", "ORF"),
				name   = "Functionality") +
			scale_y_continuous(expand = expansion(mult = c(0, 0.25)), limits = c(0, 108)) +
			facet_wrap(~ SPECIES, scales = "free_x", ncol = n_col) +
			labs(title    = "Global motif conservation by locus and functionality",
				subtitle  = if (n_sp == 1) paste("Species:", unique(data_summary$SPECIES))
				else "Mean presence rate across all motifs per locus × species",
				x = "Locus", y = "Mean conservation rate (%)") +
			theme_imgt(base_size = 13) +
			theme(
				strip.text    = element_text(size = 13, face = "bold"),
				axis.text.x   = element_text(size = 12, face = "bold"),
				axis.text.y   = element_text(size = 11),
				panel.spacing = unit(1.5, "lines")
			)
		
	} else {  
		
		data_long_raw <- data %>%
			select(SUBGROUP, SPECIES, FUNCTIONALITY, all_of(pct_cols)) %>%
			filter(!is.na(FUNCTIONALITY)) %>%
			pivot_longer(cols = all_of(pct_cols), names_to = "Motif", values_to = "Pct")
		
		if (nrow(data_long_raw) == 0 || length(unique(data_long_raw$SUBGROUP)) == 0) {
			return(ggplot2::ggplot() +
					ggplot2::annotate("text", x = 0.5, y = 0.5,
						label = "No data available for this selection.", size = 5, color = "grey50") +
					ggplot2::theme_void())
		}
		
		# Taux moyen tous motifs  par sous-groupe × fonctionnalité
		data_summary <- data_long_raw %>%
			group_by(SUBGROUP, SPECIES, FUNCTIONALITY) %>%
			summarise(MeanConservation = round(mean(Pct, na.rm = TRUE), 1), .groups = "drop")
		
		n_subgrp <- length(unique(data_summary$SUBGROUP))
		n_sp     <- length(unique(data_summary$SPECIES))
		sp_label <- if (n_sp == 1) unique(data_summary$SPECIES) else ""
		
		# Dispatch heatmap / barplot
		use_heatmap <- (force_mode == "heatmap") ||
			(force_mode == "auto" && sg_use_heatmap(n_subgrp))
		
		if (use_heatmap) {
			# Heatmap globale : SUBGROUP × FUNCTIONALITY
			sg_order <- data_summary %>%
				group_by(SUBGROUP) %>%
				summarise(m = mean(MeanConservation, na.rm = TRUE), .groups = "drop") %>%
				arrange(dplyr::desc(m)) %>% pull(SUBGROUP)
			data_summary$SUBGROUP <- factor(data_summary$SUBGROUP, levels = rev(sg_order))
			
			txt_size        <- sg_text_size(n_subgrp)
			tile_label_size <- sg_tile_label_size(n_subgrp)
			
			p <- ggplot(data_summary,
				aes(x = FUNCTIONALITY, y = SUBGROUP, fill = MeanConservation)) +
				geom_tile(color = "white", linewidth = 0.4) +
				geom_text(aes(
					label = ifelse(!is.na(MeanConservation) & MeanConservation > 0,
						paste0(MeanConservation, "%"), ""),
					color = ifelse(is.na(MeanConservation) | MeanConservation <= 55,
						"white", "grey20")),
					size = tile_label_size, fontface = "bold") +
				scale_color_identity() +
				scale_fill_viridis_c(
					option = "viridis",
					limits = c(0, 100),
					name   = "Mean conservation (%)",
					na.value = "#F9FAFB"
				) +
				scale_x_discrete(limits = c("F", "P", "ORF")) +
				labs(
					title    = "Global motif conservation — Heatmap by subgroup",
					subtitle = if (sp_label != "")
						paste("Species:", sp_label, "| Mean presence rate across all motifs")
					else "Mean presence rate across all motifs",
					x = "Functionality", y = NULL
				) +
				theme_imgt(base_size = 10) +
				theme(
					axis.text.y       = element_text(size = txt_size),
					axis.text.x       = element_text(size = 11, face = "bold"),
					panel.grid        = element_blank(),
					legend.position   = "right",
					legend.key.height = unit(1.2, "cm")
				)
			if (interactive) {
				h  <- sg_plot_height_px(n_subgrp, min_h = 420)
				ml <- sg_margin_left_px(unique(data_summary$SUBGROUP))
				p <- ggplotly(p, tooltip = c("x","y","fill"), height = h) %>%
					plotly::layout(margin = list(l = ml))
			}
			return(p)
		}
		
		# Barplot STACKED global 
		subgrp_order <- data_summary %>%
			group_by(SUBGROUP) %>%
			summarise(tot = sum(MeanConservation, na.rm = TRUE), .groups = "drop") %>%
			arrange(tot) %>% pull(SUBGROUP)
		data_summary$SUBGROUP <- factor(data_summary$SUBGROUP, levels = subgrp_order)
		
		txt_size <- sg_text_size(n_subgrp)
		
		p <- ggplot(data_summary,
			aes(x = SUBGROUP, y = MeanConservation, fill = FUNCTIONALITY)) +
			geom_bar(stat = "identity", position = "stack", width = 0.80) +
			geom_text(aes(
				label = ifelse(MeanConservation > 5,
					paste0(round(MeanConservation, 0), "%"), ""),
				color = FUNCTIONALITY),
				position = position_stack(vjust = 0.5),
				size = 2.5, fontface = "bold", show.legend = FALSE) +
			scale_color_manual(values = FUNC_LABEL_COLORS, guide = "none") +
			coord_flip() +
			scale_fill_manual(values = func_pal_global,
				breaks = c("F", "P", "ORF"),
				labels = c("Functional (F)", "Pseudogene (P)", "ORF"),
				name   = "Functionality") +
			scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
			labs(title    = "Global motif conservation by subgroup",
				subtitle  = if (sp_label != "")
					paste("Species:", sp_label, "| Mean rate across all motifs | Fill = functionality")
				else "Mean conservation rate across all motifs | Fill = functionality",
				x = NULL, y = "Mean conservation rate (%)") +
			theme_imgt(base_size = 9) +
			theme(
				axis.text.y     = element_text(size = txt_size),
				legend.position = "right"
			)
	}
	
	if (interactive) {
		if (grouping == "subgroup") {
			h  <- sg_plot_height_px(n_subgrp, min_h = 400)
			ml <- sg_margin_left_px(unique(data_summary$SUBGROUP))
			p  <- ggplotly(p, tooltip = c("x", "y", "fill"), height = h) %>%
				plotly::layout(margin = list(l = ml))
		} else if (grouping == "locus") {
			n_rows_g <- ceiling(n_sp / n_col)
			h <- if (n_sp == 1) 500 else max(500, n_rows_g * 500)
			p <- ggplotly(p, tooltip = c("x", "y", "fill"), height = min(h, 900))
		} else {
			h <- max(450, n_sp * 120)
			p <- ggplotly(p, tooltip = c("x", "y", "fill"), height = min(h, 650))
		}
	}
	
	return(p)
}

# ══════════════════════════════════════════════════════════════
# CLUSTERING HIÉRARCHIQUE
# ══════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════
# Fonction commune 
# ══════════════════════════════════════════════════════════════

get_consensus_sequences <- function(df, seq_col, group_col) {
	
	dat <- df %>%
		select(all_of(c(group_col, seq_col))) %>%
		rename(Group = all_of(group_col), Sequence = all_of(seq_col)) %>%
		filter(!is.na(Sequence), Sequence != "", Sequence != "absent")
	
	groups <- unique(dat$Group)
	
	consensus_list <- lapply(groups, function(g) {
		
		seqs <- dat$Sequence[dat$Group == g]
		n_seqs <- length(seqs)
		
		# Longueur de référence = longueur la plus fréquente 
		lens <- nchar(seqs)
		ref_len <- as.integer(names(sort(table(lens), decreasing = TRUE))[1])
		seqs <- seqs[lens == ref_len]  
		
		if (length(seqs) == 0) {
			return(data.frame(Group = g, Sequence = NA_character_, n_seqs = n_seqs))
		}
		
		# Matrice caractère x position
		char_mat <- do.call(rbind, strsplit(seqs, ""))
		
		# Consensus position par position 
		consensus_seq <- apply(char_mat, 2, function(col) {
			tab <- sort(table(col), decreasing = TRUE)
			if (length(tab) > 1 && tab[1] == tab[2]) "N" else names(tab)[1]
		})
		
		data.frame(Group = g,
			Sequence = paste(consensus_seq, collapse = ""),
			n_seqs = n_seqs)
	})
	
	out <- bind_rows(consensus_list)
	names(out)[1] <- group_col
	return(out)
}

hclust_motif_presence <- function(df, species, type_data= "IG") {
		df <- df %>% filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	
	df <- df %>% filter(SPECIES == species)
	bin_cols <- names(df)[str_starts(names(df), "BIN_")]
	
	mat <- df %>%
		group_by(SUBGROUP) %>%
		summarise(across(all_of(bin_cols), ~ mean(.x, na.rm = TRUE)),
			.groups = "drop") %>%
		column_to_rownames("SUBGROUP") %>%
		as.matrix()
	
	dist_obj       <- dist(mat, method = "euclidean")
	hc             <- hclust(dist_obj, method = "ward.D2")
	hc$dist_matrix <- as.matrix(dist_obj)
	
	return(hc)
}


hclust_sequence_similarity <- function(df, species, motif_name, type_data = "IG") {
	
	df <- df %>% filter(LOCI %in% c("IGHV", "IGKV", "IGLV"), SPECIES == species)
	seq_col <- paste0("SEQ_", motif_name)
	if (!seq_col %in% names(df)) stop(paste("Column", seq_col, "not found"))
	
	seqs_df <- get_consensus_sequences(df, seq_col, "SUBGROUP") %>%
		filter(!is.na(Sequence))
	
	subgroups <- seqs_df$SUBGROUP
	sequences <- seqs_df$Sequence
	n <- length(sequences)
	dist_mat <- matrix(0, n, n, dimnames = list(subgroups, subgroups))
	
	for (i in seq_len(n - 1)) {
		for (j in (i + 1):n) {
			s1 <- strsplit(sequences[i], "")[[1]]
			s2 <- strsplit(sequences[j], "")[[1]]
			l  <- min(length(s1), length(s2))
			d  <- sum(s1[1:l] != s2[1:l]) + abs(length(s1) - length(s2))
			dist_mat[i, j] <- dist_mat[j, i] <- d
		}
	}
	
	dist_obj <- as.dist(dist_mat)
	hc <- hclust(dist_obj, method = "ward.D2")
	hc$dist_matrix <- dist_mat
	return(hc)
}


# ══════════════════════════════════════════════════════════════
#  HEATMAP D'IDENTITÉ DES SÉQUENCES ENTRE MOTIFS / ESPÈCES
# ══════════════════════════════════════════════════════════════

calculate_sequence_identity <- function(df, motif_name, type_data = "IG") {
	
	df <- df %>% filter(LOCI %in% c("IGHV", "IGKV", "IGLV"))
	
	seq_col <- paste0("SEQ_", motif_name)
	if (!seq_col %in% names(df))
		stop(paste("Column", seq_col, "not found"))
	
	consensus <- get_consensus_sequences(df, seq_col, "SPECIES") %>%
		filter(!is.na(Sequence))
	
	if (nrow(consensus) < 2)
		stop("Not enough species with valid sequences (minimum 2)")
	
	species_list <- consensus$SPECIES
	seqs         <- consensus$Sequence
	n            <- length(species_list)
	
	mat <- matrix(NA_real_, n, n,
		dimnames = list(species_list, species_list))
	
	pct_identity <- function(s1, s2) {
		s1 <- strsplit(toupper(s1), "")[[1]]
		s2 <- strsplit(toupper(s2), "")[[1]]
		l  <- min(length(s1), length(s2))
		if (l == 0) return(0)
		round(sum(s1[1:l] == s2[1:l]) / max(length(s1), length(s2)) * 100, 1)
	}
	
	for (i in seq_len(n)) {
		for (j in seq_len(n)) {
			mat[i, j] <- pct_identity(seqs[i], seqs[j])
		}
	}
	
	return(as.data.frame(mat) %>%
			rownames_to_column("Species_1") %>%
			pivot_longer(-Species_1,
				names_to  = "Species_2",
				values_to = "Identity"))
}



plot_sequence_identity_heatmap <- function(df, motif_name, type_data= "IG",
	interactive = TRUE) {
	
	identity_data <- calculate_sequence_identity(df, motif_name)
	
	# Seuil 75% : couleur foncée au-dessus
	p <- ggplot(identity_data,
		aes(x = Species_2, y = Species_1, fill = Identity)) +
		geom_tile(color = "white", linewidth = 0.6) +
		geom_text(aes(
			label = paste0(Identity, "%"),
			color = ifelse(Identity <= 60, "white", "grey20")
		),
			size = 3.5, fontface = "bold") +
		scale_color_identity() +
		scale_fill_viridis_c(
			option = "cividis",
			limits = c(0, 100),
			name   = "Identity (%)",
			breaks = c(0, 25, 50, 75, 100),
			labels = c("0%", "25%", "50%", "75%", "100%")
		) +
		# Ligne de seuil 75% — annotation visuelle
		annotate("text",
			x     = Inf, y = Inf,
			label = "■ > 75% = high identity",
			hjust = 1.05, vjust = 1.5,
			size  = 3, color = "#264E52", fontface = "italic") +
		labs(
			title    = paste0("Sequence identity heatmap — Motif: ", motif_name),
			subtitle = "Pairwise identity (%) between consensus sequences per species",
			x = NULL, y = NULL
		) +
		theme_imgt(base_size = 11) +
		theme(
			axis.text.x      = element_text(angle = 40, hjust = 1, size = 10),
			axis.text.y      = element_text(size = 10),
			panel.grid       = element_blank(),
			legend.position  = "right",
			legend.key.height = unit(1.2, "cm")
		)
	
	if (interactive) {
		p <- ggplotly(p, tooltip = c("x", "y", "fill"), height = 550)
	}
	
	return(p)
}


# ══════════════════════════════════════════════════════════════
# HEATMAP D'IDENTITÉ ENTRE LOCUS 
# ══════════════════════════════════════════════════════════════

calculate_locus_identity <- function(df, motif_name, species) {
	
	df <- df %>% filter(LOCI %in% c("IGHV", "IGKV", "IGLV"),
		SPECIES == species)
	
	seq_col <- paste0("SEQ_", motif_name)
	if (!seq_col %in% names(df))
		stop(paste("Column", seq_col, "not found"))
	
	consensus <- get_consensus_sequences(df, seq_col, "LOCI") %>%
		filter(!is.na(Sequence))
	
	if (nrow(consensus) < 2)
		stop("Not enough loci with valid sequences (minimum 2)")
	
	loci_list <- consensus$LOCI
	seqs      <- consensus$Sequence
	n         <- length(loci_list)
	
	mat <- matrix(NA_real_, n, n, dimnames = list(loci_list, loci_list))
	
	pct_identity <- function(s1, s2) {
		s1 <- strsplit(toupper(s1), "")[[1]]
		s2 <- strsplit(toupper(s2), "")[[1]]
		l  <- min(length(s1), length(s2))
		if (l == 0) return(0)
		round(sum(s1[1:l] == s2[1:l]) / max(length(s1), length(s2)) * 100, 1)
	}
	
	for (i in seq_len(n))
		for (j in seq_len(n))
			mat[i, j] <- pct_identity(seqs[i], seqs[j])
	
	as.data.frame(mat) %>%
		rownames_to_column("Locus_1") %>%
		pivot_longer(-Locus_1, names_to = "Locus_2", values_to = "Identity")
}

# Heatmap d'identité entre locus d'une espèce
plot_locus_identity_heatmap <- function(df, motif_name, species, interactive = TRUE) {
	
	identity_data <- calculate_locus_identity(df, motif_name, species)
	
	p <- ggplot(identity_data,
		aes(x = Locus_2, y = Locus_1, fill = Identity)) +
		geom_tile(color = "white", linewidth = 0.8) +
		geom_text(aes(
			label = paste0(Identity, "%"),
			color = ifelse(Identity <= 60, "white", "grey20")
		), size = 4.5, fontface = "bold") +
		scale_color_identity() +
		scale_fill_viridis_c(
			option = "cividis",
			limits = c(0, 100),
			name   = "Identity (%)",
			breaks = c(0, 25, 50, 75, 100),
			labels = c("0%", "25%", "50%", "75%", "100%")
		) +
		annotate("text", x = Inf, y = Inf,
			label = "■ > 75% = high identity",
			hjust = 1.05, vjust = 1.5,
			size = 3, color = "#264E52", fontface = "italic") +
		labs(
			title    = paste0("Sequence identity between loci — Motif: ", motif_name),
			subtitle = paste0("Species: ", species, " | Consensus sequences per locus"),
			x = NULL, y = NULL
		) +
		theme_imgt(base_size = 12) +
		theme(
			axis.text.x      = element_text(size = 12),
			axis.text.y      = element_text(size = 12),
			panel.grid       = element_blank(),
			legend.position  = "right",
			legend.key.height = unit(1.2, "cm")
		)
	
	if (interactive)
		p <- ggplotly(p, tooltip = c("x", "y", "fill"), height = 420)
	
	return(p)
}


# ══════════════════════════════════════════════════════════════
# HEATMAP D'IDENTITÉ ENTRE SOUS-GROUPES 
# ══════════════════════════════════════════════════════════════

calculate_subgroup_identity <- function(df, motif_name, species) {
	
	df <- df %>% filter(LOCI %in% c("IGHV", "IGKV", "IGLV"),
		SPECIES == species)
	
	seq_col <- paste0("SEQ_", motif_name)
	if (!seq_col %in% names(df))
		stop(paste("Column", seq_col, "not found"))
	
	# Séquence consensus par sous-groupe (la plus fréquente)
	consensus <- df %>%
		select(SUBGROUP, all_of(seq_col)) %>%
		filter(!is.na(.data[[seq_col]]),
			.data[[seq_col]] != "",
			.data[[seq_col]] != "absent") %>%
		rename(Sequence = all_of(seq_col)) %>%
		group_by(SUBGROUP, Sequence) %>%
		summarise(n = n(), .groups = "drop") %>%
		group_by(SUBGROUP) %>%
		slice_max(n, n = 1, with_ties = FALSE) %>%
		ungroup() %>%
		select(SUBGROUP, Sequence)
	
	if (nrow(consensus) < 2)
		stop("Not enough subgroups with valid sequences (minimum 2)")
	
	sg_list <- consensus$SUBGROUP
	seqs    <- consensus$Sequence
	n       <- length(sg_list)
	
	mat <- matrix(NA_real_, n, n, dimnames = list(sg_list, sg_list))
	
	pct_identity <- function(s1, s2) {
		s1 <- strsplit(toupper(s1), "")[[1]]
		s2 <- strsplit(toupper(s2), "")[[1]]
		l  <- min(length(s1), length(s2))
		if (l == 0) return(0)
		round(sum(s1[1:l] == s2[1:l]) / max(length(s1), length(s2)) * 100, 1)
	}
	
	for (i in seq_len(n))
		for (j in seq_len(n))
			mat[i, j] <- pct_identity(seqs[i], seqs[j])
	
	as.data.frame(mat) %>%
		rownames_to_column("Subgroup_1") %>%
		pivot_longer(-Subgroup_1, names_to = "Subgroup_2", values_to = "Identity")
}

# Heatmap d'identité entre sous-groupes d'une espèce
plot_subgroup_identity_heatmap <- function(df, motif_name, species, interactive = TRUE) {
	
	identity_data <- calculate_subgroup_identity(df, motif_name, species)
	n_sg <- length(unique(identity_data$Subgroup_1))
	h    <- max(450, n_sg * 28)
	
	p <- ggplot(identity_data,
		aes(x = Subgroup_2, y = Subgroup_1, fill = Identity)) +
		geom_tile(color = "white", linewidth = 0.5) +
		geom_text(aes(
			label = paste0(Identity, "%"),
			color = ifelse(Identity <= 60, "white", "grey20")
		), size = 2.8, fontface = "bold") +
		scale_color_identity() +
		scale_fill_viridis_c(
			option = "cividis",
			limits = c(0, 100),
			name   = "Identity (%)",
			breaks = c(0, 25, 50, 75, 100),
			labels = c("0%", "25%", "50%", "75%", "100%")
		) +
		annotate("text", x = Inf, y = Inf,
			label = "■ > 75% = high identity",
			hjust = 1.05, vjust = 1.5,
			size = 3, color = "#264E52", fontface = "italic") +
		labs(
			title    = paste0("Sequence identity between subgroups — Motif: ", motif_name),
			subtitle = paste0("Species: ", species, " | Consensus sequences per subgroup"),
			x = NULL, y = NULL
		) +
		theme_imgt(base_size = 10) +
		theme(
			axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
			axis.text.y      = element_text(size = 8),
			panel.grid       = element_blank(),
			legend.position  = "right",
			legend.key.height = unit(1.2, "cm")
		)
	
	if (interactive)
		p <- ggplotly(p, tooltip = c("x", "y", "fill"), height = h)
	
	return(p)
}
