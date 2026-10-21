# ═══════════════════════════════════════════════════════════════
#  Validation des Modèles de Régression Logistique
# ═══════════════════════════════════════════════════════════════

# Dépend de : utils_core.R
# Packages  : tidyverse, broom, pROC, plotly, scales

library(tidyverse)
library(broom)

# ══════════════════════════════════════════════════════════════
#  VALIDATION CROISÉE k-fold STRATIFIÉE
# ══════════════════════════════════════════════════════════════


validate_logistic_model <- function(df,
	motif_name,
	type_data = "IG",
	n_folds   = 5,
	seed      = 42) { 

	
	if (!requireNamespace("pROC", quietly = TRUE))
		stop("Package 'pROC' is required. Install it with: install.packages('pROC')")
	
	#  Filtrage et préparation
	if (type_data == "IG") {
		df <- df %>% dplyr::filter(LOCI %in% c("IGHV","IGKV","IGLV"))
	} else {
		df <- df %>% dplyr::filter(LOCI %in% c("TRAV","TRBV","TRDV","TRGV"))
	}
	
	bin_col <- paste0("BIN_", motif_name)
	if (!bin_col %in% names(df))
		stop(paste0("Column '", bin_col, "' missing from the dataset."))
	
	df_model <- df %>%
		dplyr::select(dplyr::all_of(bin_col), LOCI, SPECIES, FUNCTIONALITY) %>%
		dplyr::mutate(
			Presence      = as.factor(.data[[bin_col]]),
			LOCI          = as.factor(LOCI),
			SPECIES       = as.factor(SPECIES),
			FUNCTIONALITY = as.factor(FUNCTIONALITY)
		) %>%
		dplyr::filter(!is.na(Presence), !is.na(LOCI),
			!is.na(SPECIES), !is.na(FUNCTIONALITY)) %>%
		dplyr::select(-dplyr::all_of(bin_col))
	
	if (nrow(df_model) < 20)
		stop("Fewer than 20 valid observations — validation impossible.")
	
	n_levels <- length(unique(df_model$Presence))
	if (n_levels < 2)
		stop(paste0("Motif '", motif_name, "' has only one class — regression impossible."))
	
	df_model <- droplevels(df_model)
	lf_full  <- .logit_terms(df_model)
	
	n <- nrow(df_model)
	set.seed(seed)
	
	# Folds STRATIFIÉS par classe
	idx_1 <- which(df_model$Presence == "1")
	idx_0 <- which(df_model$Presence == "0")
	
	folds <- integer(n)
	if (length(idx_1) >= n_folds)
		folds[idx_1] <- sample(rep(1:n_folds, length.out = length(idx_1)))
	else
		folds[idx_1] <- sample(1:n_folds, length(idx_1), replace = TRUE)
	
	if (length(idx_0) >= n_folds)
		folds[idx_0] <- sample(rep(1:n_folds, length.out = length(idx_0)))
	else
		folds[idx_0] <- sample(1:n_folds, length(idx_0), replace = TRUE)
	
	# Boucle CV
	cv_probs  <- numeric(n)
	cv_preds  <- factor(rep(NA, n), levels = c("0","1"))
	fold_metrics <- list()
	
	for (k in seq_len(n_folds)) {
		idx_train <- which(folds != k)
		idx_test  <- which(folds == k)
		
		df_train <- df_model[idx_train, ]
		df_train <- droplevels(df_train)
		df_test  <- df_model[idx_test,  ]
		
		# Vérifier présence des deux classes dans train
		if (length(unique(df_train$Presence)) < 2) {
			message(paste0("Fold ", k, " : only one class — fold skipped."))
			next
		}
		
		# Gérer les niveaux absent du fold de test dans les facteurs
		for (col_fac in c("LOCI","SPECIES","FUNCTIONALITY")) {
			df_test[[col_fac]] <- factor(df_test[[col_fac]],
				levels = levels(df_train[[col_fac]]))
		}
		
		fit_k <- tryCatch(
			suppressWarnings(
				glm(
					.logit_terms(df_train)$formula,
					data = df_train,
					family = binomial(link = "logit")
				)
			),
			error = function(e) NULL
		)
		
		if (is.null(fit_k)) next
		
		prob_k <- tryCatch(
			predict(fit_k, newdata = df_test, type = "response"),
			error = function(e) rep(NA_real_, nrow(df_test))
		)
		
		cv_probs[idx_test] <- prob_k
		cv_preds[idx_test] <- factor(
			ifelse(prob_k >= 0.5, "1", "0"), levels = c("0","1")
		)
		
		y_test <- df_test$Presence
		tp <- sum(cv_preds[idx_test] == "1" & y_test == "1", na.rm = TRUE)
		tn <- sum(cv_preds[idx_test] == "0" & y_test == "0", na.rm = TRUE)
		fp <- sum(cv_preds[idx_test] == "1" & y_test == "0", na.rm = TRUE)
		fn <- sum(cv_preds[idx_test] == "0" & y_test == "1", na.rm = TRUE)
		
		fold_metrics[[k]] <- data.frame(
			Fold      = k,
			N         = length(idx_test),
			Accuracy  = round((tp + tn) / length(idx_test) * 100, 2),
			Precision = round(ifelse(tp + fp > 0, tp / (tp + fp) * 100, NA_real_), 2),
			Recall    = round(ifelse(tp + fn > 0, tp / (tp + fn) * 100, NA_real_), 2),
			F1        = round(ifelse((tp + fp > 0) && (tp + fn > 0),
				2 * tp / (2 * tp + fp + fn) * 100, NA_real_), 2),
			stringsAsFactors = FALSE
		)
	}
	
	cv_metrics <- dplyr::bind_rows(fold_metrics)
	
	# Métriques globales
	valid_idx <- !is.na(cv_preds)
	y_true    <- df_model$Presence[valid_idx]
	y_pred    <- cv_preds[valid_idx]
	y_prob    <- cv_probs[valid_idx]
	
	tp_g <- sum(y_pred == "1" & y_true == "1", na.rm = TRUE)
	tn_g <- sum(y_pred == "0" & y_true == "0", na.rm = TRUE)
	fp_g <- sum(y_pred == "1" & y_true == "0", na.rm = TRUE)
	fn_g <- sum(y_pred == "0" & y_true == "1", na.rm = TRUE)
	
	metrics_global <- list(
		Accuracy  = round((tp_g + tn_g) / sum(valid_idx) * 100, 2),
		Precision = round(ifelse(tp_g + fp_g > 0, tp_g / (tp_g + fp_g) * 100, NA_real_), 2),
		Recall    = round(ifelse(tp_g + fn_g > 0, tp_g / (tp_g + fn_g) * 100, NA_real_), 2),
		F1        = round(ifelse((tp_g + fp_g > 0) && (tp_g + fn_g > 0),
			2 * tp_g / (2 * tp_g + fp_g + fn_g) * 100, NA_real_), 2),
		Specificity = round(ifelse(tn_g + fp_g > 0, tn_g / (tn_g + fp_g) * 100, NA_real_), 2),
		N_valid   = sum(valid_idx),
		N_total   = n
	)
	
	# ROC
	roc_data <- tryCatch({
		roc_obj <- pROC::roc(
			response  = as.integer(y_true == "1"),
			predictor = y_prob,
			quiet     = TRUE
		)
		auc_val <- as.numeric(pROC::auc(roc_obj))
		data.frame(
			FPR = 1 - roc_obj$specificities,
			TPR = roc_obj$sensitivities,
			AUC = auc_val
		)
	}, error = function(e) NULL)
	
	if (!is.null(roc_data))
		metrics_global$AUC <- round(roc_data$AUC[1] * 100, 2)
	else
		metrics_global$AUC <- NA_real_
	
	# Courbe de calibration
	calibration_data <- tryCatch({
		df_cal <- data.frame(prob = y_prob, obs = as.integer(y_true == "1")) %>%
			dplyr::filter(!is.na(prob)) %>%
			dplyr::mutate(bin = cut(prob, breaks = seq(0, 1, 0.1),
				include.lowest = TRUE)) %>%
			dplyr::group_by(bin) %>%
			dplyr::summarise(
				Mean_Predicted = mean(prob, na.rm = TRUE),
				Mean_Observed  = mean(obs,  na.rm = TRUE),
				N              = dplyr::n(),
				.groups = "drop"
			) %>%
			dplyr::filter(!is.na(Mean_Predicted))
		df_cal
	}, error = function(e) NULL)
	
	# Matrice de confusion globale
	conf_mat <- table(Predicted = y_pred, Actual = y_true)
	
	# Modèle final sur toutes les données
	fit_full <- tryCatch(
		suppressWarnings(
			glm(
				lf_full$formula,
				data = df_model,
				family = binomial(link = "logit")
			)
		),
		error = function(e) NULL
	)
	
	coef_df <- NULL
	if (!is.null(fit_full)) {
		coef_df <- broom::tidy(fit_full) %>%
			dplyr::mutate(
				OR        = exp(estimate),
				CI_lower  = exp(estimate - 1.96 * std.error),
				CI_upper  = exp(estimate + 1.96 * std.error),
				Significant = dplyr::case_when(
					is.na(p.value) ~ NA_character_,
					p.value < 0.001 ~ "***",
					p.value < 0.01  ~ "**",
					p.value < 0.05  ~ "*",
					TRUE ~ "ns"
				)
			) %>%
			dplyr::select(term, estimate, std.error, OR,
				CI_lower, CI_upper, statistic, p.value, Significant)
	}
	
	return(list(
		cv_metrics      = cv_metrics,
		metrics_global  = metrics_global,
		roc_data        = roc_data,
		calibration     = calibration_data,
		conf_mat        = conf_mat,
		coef_df         = coef_df,
		fit_full        = fit_full,
		motif           = motif_name,
		n_folds         = n_folds,
		n               = n
	))
}


# ══════════════════════════════════════════════════════════════
# COURBE ROC
# ══════════════════════════════════════════════════════════════

plot_roc_logistic <- function(validation_result, interactive = TRUE) {
	
	if (is.null(validation_result$roc_data)) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "ROC curve unavailable.", size = 5, color = "grey50") +
				ggplot2::theme_void())
	}
	
	roc_df <- validation_result$roc_data
	auc    <- round(roc_df$AUC[1] * 100, 1)
	
	auc_label <- dplyr::case_when(
		auc >= 90 ~ "Excellent",
		auc >= 80 ~ "Very good",
		auc >= 70 ~ "Good",
		auc >= 60 ~ "Moderate",
		TRUE      ~ "Poor"
	)
	
	p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = FPR, y = TPR)) +
		ggplot2::geom_ribbon(ggplot2::aes(ymin = FPR, ymax = TPR),
			fill = "#1565C0", alpha = 0.15) +
		ggplot2::geom_line(color = "#1565C0", linewidth = 1.8) +
		ggplot2::geom_abline(slope = 1, intercept = 0,
			color = "#DC2626", linetype = "dashed", linewidth = 1.0) +
		ggplot2::annotate("text", x = 0.62, y = 0.12,
			label = paste0("AUC = ", auc, "% (", auc_label, ")"),
			size = 5.2, fontface = "bold", color = "#1565C0") +
		ggplot2::annotate("text", x = 0.78, y = 0.06,
			label = "Random classifier",
			size = 3.2, color = "#DC2626", fontface = "italic") +
		ggplot2::scale_x_continuous(limits = c(0, 1),
			labels = scales::percent_format(accuracy = 1),
			expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
		ggplot2::scale_y_continuous(limits = c(0, 1),
			labels = scales::percent_format(accuracy = 1),
			expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
		ggplot2::coord_fixed(ratio = 1) +  
		ggplot2::labs(
			title    = paste0("ROC Curve — Logistic Regression | Motif: '",
				validation_result$motif, "'"),
			subtitle = paste0(validation_result$n_folds,
				"-fold stratified CV | N = ", validation_result$n, " sequences"),
			x = "False Positive Rate (1 - Specificity)",
			y = "True Positive Rate (Sensitivity)"
		) +
		theme_imgt(base_size = 12)
	
	if (interactive) {
		p <- plotly::ggplotly(p, height = 620, width = 620) %>%
			plotly::layout(
				xaxis  = list(scaleanchor = "y", scaleratio = 1, constrain = "domain"),
				yaxis  = list(constrain = "domain"),
				margin = list(t = 60)
			)
	}
	
	return(p)
}


# ══════════════════════════════════════════════════════════════
# COURBE DE CALIBRATION
# ══════════════════════════════════════════════════════════════

plot_calibration_curve <- function(validation_result, interactive = TRUE) {
	
	cal <- validation_result$calibration
	
	if (is.null(cal) || nrow(cal) == 0) {
		return(ggplot2::ggplot() +
				ggplot2::annotate("text", x = 0.5, y = 0.5,
					label = "Calibration data unavailable.", size = 5, color = "grey50") +
				ggplot2::theme_void())
	}
	
	p <- ggplot2::ggplot(cal,
		ggplot2::aes(x = Mean_Predicted, y = Mean_Observed,
			size = N, color = N,
			text = paste0("Predicted: ", round(Mean_Predicted, 3),
				"\nObserved: ", round(Mean_Observed, 3),
				"\nN: ", N))) +
		ggplot2::geom_abline(slope = 1, intercept = 0,
			color = "#16A34A", linetype = "dashed", linewidth = 0.9) +
		ggplot2::geom_segment(ggplot2::aes(
			x = Mean_Predicted, xend = Mean_Predicted,
			y = Mean_Predicted, yend = Mean_Observed
		), color = "grey70", linewidth = 0.6) +
		ggplot2::geom_point(alpha = 0.85) +
		ggplot2::scale_color_viridis_c(option = "viridis", name = "N seq.") +
		ggplot2::scale_size_continuous(range = c(3, 12), guide = "none") +
		ggplot2::scale_x_continuous(limits = c(0, 1),
			labels = scales::percent_format(accuracy = 1)) +
		ggplot2::scale_y_continuous(limits = c(0, 1),
			labels = scales::percent_format(accuracy = 1)) +
		ggplot2::labs(
			title    = paste0("Calibration Curve — '", validation_result$motif, "'"),
			subtitle = "Perfect calibration = points on the green dashed diagonal | Point size ∝ N",
			x = "Mean Predicted Probability",
			y = "Mean Observed Frequency"
		) +
		theme_imgt(base_size = 12)
	
	if (interactive)
		p <- plotly::ggplotly(p, tooltip = "text", height = 440)
	
	return(p)
}


# ══════════════════════════════════════════════════════════════
#  CV PAR FOLD
# ══════════════════════════════════════════════════════════════

plot_cv_metrics_logistic <- function(validation_result, interactive = TRUE) {
	
	metric_cols <- c("Accuracy", "Precision", "Recall", "F1")
	
	summary_df <- validation_result$cv_metrics %>%
		tidyr::pivot_longer(cols = dplyr::all_of(metric_cols),
			names_to = "Metric", values_to = "Value") %>%
		dplyr::filter(!is.na(Value)) %>%
		dplyr::group_by(Metric) %>%
		dplyr::summarise(
			Mean = mean(Value, na.rm = TRUE),
			SD   = sd(Value, na.rm = TRUE),
			N    = dplyr::n(),
			.groups = "drop"
		) %>%
		dplyr::mutate(Metric = factor(Metric, levels = metric_cols))
	
	metric_colors <- c(
		"Accuracy"  = "#1565C0",
		"Precision" = "#166534",
		"Recall"    = "#854D0E",
		"F1"        = "#7C3AED"
	)
	
	p <- ggplot2::ggplot(summary_df,
		ggplot2::aes(x = Metric, y = Mean, fill = Metric,
			text = paste0(Metric, ": ", round(Mean, 1), "% \u00b1 ", round(SD, 1), "%"))) +
		ggplot2::geom_col(width = 0.55) +
		ggplot2::geom_errorbar(
			ggplot2::aes(ymin = pmax(0, Mean - SD), ymax = pmin(100, Mean + SD)),
			width = 0.18, linewidth = 0.8, color = "grey30"
		) +
		ggplot2::geom_text(
			ggplot2::aes(label = paste0(round(Mean, 1), "%")),
			vjust = -1.6, size = 4, fontface = "bold"
		) +
		ggplot2::scale_fill_manual(values = metric_colors, guide = "none") +
		ggplot2::scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 20),
			labels = function(x) paste0(x, "%")) +
		ggplot2::labs(
			title    = paste0("Cross-validation metrics (mean \u00b1 SD) \u2014 Motif: '",
				validation_result$motif, "'"),
			subtitle = paste0(validation_result$n_folds,
				"-fold stratified CV | Error bars = \u00b1 1 SD across folds"),
			x = NULL, y = "Score (%)"
		) +
		theme_imgt(base_size = 12)
	
	if (interactive)
		p <- plotly::ggplotly(p, tooltip = "text", height = 440)
	
	return(p)
}

# ══════════════════════════════════════════════════════════════
#  TABLEAU CV 
# ══════════════════════════════════════════════════════════════

format_cv_metrics_wide <- function(cv_metrics) {
	
	if (is.null(cv_metrics) || nrow(cv_metrics) == 0) return(NULL)
	
	metric_cols <- intersect(c("N", "Accuracy", "Precision", "Recall", "F1"),
		names(cv_metrics))
	if (length(metric_cols) == 0) return(NULL)
	
	df_long <- cv_metrics %>%
		dplyr::select(Fold, dplyr::all_of(metric_cols)) %>%
		tidyr::pivot_longer(cols = dplyr::all_of(metric_cols),
			names_to = "Metric", values_to = "Value")
	
	df_wide <- df_long %>%
		tidyr::pivot_wider(id_cols = Metric, names_from = Fold,
			names_prefix = "Fold_", values_from = Value)
	
	mean_col <- df_long %>%
		dplyr::group_by(Metric) %>%
		dplyr::summarise(Mean = round(mean(Value, na.rm = TRUE), 2), .groups = "drop")
	
	df_wide <- df_wide %>%
		dplyr::left_join(mean_col, by = "Metric") %>%
		dplyr::mutate(Metric = factor(Metric, levels = metric_cols)) %>%
		dplyr::arrange(Metric) %>%
		dplyr::mutate(Metric = as.character(Metric))
	
	df_wide
}


# ══════════════════════════════════════════════════════════════
# INTERPRÉTATION AUTOMATIQUE
# ══════════════════════════════════════════════════════════════

format_logistic_interpretation <- function(validation_result) {
	
	m      <- validation_result$metrics_global
	motif  <- validation_result$motif
	n_folds <- validation_result$n_folds
	
	auc_label <- dplyr::case_when(
		is.na(m$AUC)  ~ "not computed",
		m$AUC >= 90   ~ "excellent (≥ 90%)",
		m$AUC >= 80   ~ "very good (80–90%)",
		m$AUC >= 70   ~ "good (70–80%)",
		m$AUC >= 60   ~ "moderate (60–70%)",
		TRUE           ~ "poor (< 60%)"
	)
	
	# Top prédicteur (coefficient significatif)
	top_predictor <- ""
	if (!is.null(validation_result$coef_df)) {
		sig_coefs <- validation_result$coef_df %>%
			dplyr::filter(Significant %in% c("***","**","*")) %>%
			dplyr::arrange(dplyr::desc(abs(estimate)))
		if (nrow(sig_coefs) > 0) {
			top_predictor <- paste0(
				"\n\nThe most influential predictor is **",
				sig_coefs$term[1], "** (OR = ",
				round(sig_coefs$OR[1], 2), ", ",
				sig_coefs$Significant[1],
				"). An OR > 1 indicates that this factor ",
				"increases the probability of the motif being present."
			)
		}
	}
	
	# Déséquilibre
	balance_note <- ""
	if (!is.null(validation_result$conf_mat)) {
		tab <- validation_result$conf_mat
		n1  <- if ("1" %in% colnames(tab)) sum(tab[, "1"]) else NA
		n0  <- if ("0" %in% colnames(tab)) sum(tab[, "0"]) else NA
		if (!is.na(n1) && !is.na(n0) && min(n1, n0) > 0) {
			ratio <- max(n1, n0) / min(n1, n0)
			if (ratio > 3)
				balance_note <- paste0(
					"\n\n⚠️ **Class imbalance detected** (ratio ",
					round(ratio, 1), ":1). Accuracy may be misleading — ",
					"prioritise AUC and Recall to assess real performance."
				)
		}
	}
	
	paste0(
		"The logistic regression model predicts the presence of motif **", motif,
		"** with an AUC of **", m$AUC, "%** (", auc_label, ") ",
		"and an overall Accuracy of **", m$Accuracy, "%** ",
		"in **", n_folds, "-fold stratified cross-validation** ",
		"on ", m$N_valid, " valid sequences.\n\n",
		"**Recall = ", m$Recall, "%**: proportion of truly present sequences ",
		"correctly identified by the model. ",
		"**Specificity = ", m$Specificity, "%**: true negative rate. ",
		"**F1 = ", m$F1, "%**: harmonic mean of Precision and Recall.",
		top_predictor,
		balance_note,
		"\n\n",
		"**Calibration**: a well-calibrated model has predicted probabilities ",
		"close to observed frequencies (points near the diagonal). ",
		"Deviations indicate systematic over- or under-confidence."
	)
}


# ══════════════════════════════════════════════════════════════
# CONFUSION MATRIX PLOT
# ══════════════════════════════════════════════════════════════


plot_confusion_logistic <- function(validation_result, interactive = TRUE) {
	
	conf  <- as.data.frame(validation_result$conf_mat)
	motif <- validation_result$motif
	
	conf_pct <- conf %>%
		dplyr::group_by(Actual) %>%
		dplyr::mutate(Pct = round(Freq / sum(Freq) * 100, 1)) %>%
		dplyr::ungroup() %>%
		dplyr::mutate(
			
			CellType = dplyr::case_when(
				Predicted == "1" & Actual == "1" ~ "TP",  # Vrai Positif
				Predicted == "0" & Actual == "0" ~ "TN",  # Vrai Négatif
				Predicted == "1" & Actual == "0" ~ "FP",  # Faux Positif
				Predicted == "0" & Actual == "1" ~ "FN"   # Faux Négatif
			),
			Label = paste0(Freq, " (", CellType, ")\n", Pct, "%"),
			Fill  = dplyr::case_when(
				Predicted == Actual ~ Pct,
				TRUE                ~ -Pct
			)
		)
	
	p <- ggplot2::ggplot(conf_pct,
		ggplot2::aes(x = Actual, y = Predicted, fill = Fill,
			text = paste0("Actual: ", Actual,
				"\nPredicted: ", Predicted,
				"\n", CellType, ": ", Freq, " (", Pct, "%)"))) +
		ggplot2::geom_tile(color = "white", linewidth = 2.0) +
		ggplot2::geom_text(ggplot2::aes(label = Label),
			size = 5.5, fontface = "bold", color = "#111827", lineheight = 0.9) +
		ggplot2::scale_fill_gradient2(
			low = "#F5B7B1", mid = "#F9FAFB", high = "#ABEBC6",
			midpoint = 0, guide = "none"
		) +
		ggplot2::scale_x_discrete(labels = c("0" = "Absent", "1" = "Present")) +
		ggplot2::scale_y_discrete(labels = c("0" = "Absent", "1" = "Present")) +
		ggplot2::labs(
			title    = paste0("Confusion Matrix — Logistic Regression | Motif: '", motif, "'"),
			subtitle = "VP = True Positive | VN = True Negative | FP = False Positive | FN = False Negative",
			x = "Actual class",
			y = "Predicted class"
		) +
		theme_imgt(base_size = 13) +
		ggplot2::theme(
			panel.grid = ggplot2::element_blank(),
			plot.subtitle = ggplot2::element_text(size = 10, color = "#374151")
		)
	
	if (interactive)
		p <- plotly::ggplotly(p, tooltip = "text", height = 380)
	
	return(p)
}
