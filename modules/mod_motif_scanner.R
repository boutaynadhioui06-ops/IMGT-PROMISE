# ═══════════════════════════════════════════════════════════════════════════
#  Module : Détection des motifs régulateurs 5'UTR
# ═══════════════════════════════════════════════════════════════════════════

# Dépend de : dataset_helpers.R (clean_subgroup, optionnel)
# Packages  : shiny, shinydashboard, DT, dplyr, tidyr, purrr, stringr, plotly, openxlsx

# ══════════════════════════════════════════════════════════════════════════
# 1. PÉRIMÈTRE
# ══════════════════════════════════════════════════════════════════════════

SCANNER_LOCI <- c("IGHV", "IGKV", "IGLV")


PWM_MOTIFS_BY_LOCUS <- list(
	IGHV = c("OCTAMER", "HEPTANUCLEOTIDE"),
	IGKV = c("DECAMER", "PENTADECAMER"),
	IGLV = c("DECAMER", "PENTADECAMER")
)

RULE_MOTIFS_BY_LOCUS <- list(
	IGHV = c("TATA_BOX", "E_BOX"),
	IGKV = c("TATA_BOX", "E_BOX"),
	IGLV = c("TATA_BOX", "E_BOX")
)

ALL_MOTIFS_BY_LOCUS <- purrr::imap(PWM_MOTIFS_BY_LOCUS, function(v, loc) c(v, RULE_MOTIFS_BY_LOCUS[[loc]]))

# Mismatch maximal toléré par la méthode PWM 
HARD_MISMATCH_MAX <- c(
	OCTAMER         = 1L,
	DECAMER         = 1L,
	HEPTANUCLEOTIDE = 2L,
	PENTADECAMER    = 2L
)

# Correspondance motif du scanner 
APP_MOTIF_COLUMN_MAP <- c(
	TATA_BOX        = "TATA_BOX",
	E_BOX           = "E_BOX",
	OCTAMER         = "OCTA_DECA",
	DECAMER         = "OCTA_DECA",
	HEPTANUCLEOTIDE = "HEPTA_PENTA",
	PENTADECAMER    = "HEPTA_PENTA"
)

SCANNER_SCANNER_MOTIF_COLORS <- c(
	LPART_1         = "#7B2D8B",
	TATA_BOX        = "#E63946",
	E_BOX           = "#457B9D",
	OCTAMER         = "#2A9D8F",
	HEPTANUCLEOTIDE = "#F4A261",
	DECAMER         = "#1B7F3B",
	PENTADECAMER    = "#F9A825"
)

CONSENSUS_DEFINITIONS <- list(
	IGHV = list(
		OCTAMER          = list(patterns = "ATGCAAAT", max_mismatch = 1L),
		HEPTANUCLEOTIDE  = list(patterns = "CTCATGA",  max_mismatch = 1L),
		E_BOX            = list(patterns = "CANNTG",   max_mismatch = 0L)
	),
	IGKV = list(
		DECAMER          = list(patterns = "TNATTTGCAT", max_mismatch = 1L),
		PENTADECAMER     = list(patterns = c("TGCAGNNNNNNCCAG", "CCTGTGNNNNNCCAG"), max_mismatch = 3L),
		E_BOX            = list(patterns = "CANNTG", max_mismatch = 0L)
	),
	IGLV = list(
		DECAMER          = list(patterns = "TNATTTGCAT", max_mismatch = 1L),
		PENTADECAMER     = list(patterns = c("TGCAGNNNNNNCCAG", "CCTGTGNNNNNCCAG"), max_mismatch = 3L),
		E_BOX            = list(patterns = "CANNTG", max_mismatch = 0L)
	)
)


# ══════════════════════════════════════════════════════════════════════════
# 3. MOTEUR DE MATCHING IUPAC 
# ══════════════════════════════════════════════════════════════════════════

.scan_pattern_all_windows <- function(seq_chars, valid_acgt, len_s, pattern, max_mismatch = 0L) {
	len_p <- nchar(pattern)
	if (len_s < len_p) return(NULL)
	n_win <- len_s - len_p + 1L
	pattern_chars <- strsplit(pattern, "")[[1]]
	pattern_is_n  <- pattern_chars == "N"
	
	idx <- outer(seq_len(len_p), seq_len(n_win), FUN = "+") - 1L
	window_mat <- matrix(seq_chars[idx], nrow = len_p, ncol = n_win)
	
	mismatch_mat <- window_mat != pattern_chars              
	if (any(pattern_is_n)) mismatch_mat[pattern_is_n, ] <- FALSE
	mismatches <- .colSums(mismatch_mat, m = len_p, n = n_win)
	
	valid_mat  <- matrix(valid_acgt[idx], nrow = len_p, ncol = n_win)
	all_valid  <- .colSums(valid_mat, m = len_p, n = n_win) == len_p
	
	keep <- which(mismatches <= max_mismatch & all_valid)
	if (length(keep) == 0L) return(NULL)
	
	data.frame(
		pos_start    = keep,
		pos_end      = keep + len_p - 1L,
		matched_seq  = vapply(keep, function(i) paste(window_mat[, i], collapse = ""), character(1)),
		mismatch     = mismatches[keep],
		pattern_used = pattern,
		stringsAsFactors = FALSE
	)
}


scan_consensus_pattern <- function(seq_upstream, pattern, max_mismatch = 0L, db = NULL) {
	seq_upstream <- toupper(seq_upstream)
	pattern      <- toupper(pattern)
	len_s <- nchar(seq_upstream)
	if (len_s < nchar(pattern)) return(NULL)
	seq_chars  <- strsplit(seq_upstream, "")[[1]]
	valid_acgt <- seq_chars %in% c("A", "C", "G", "T")
	res <- .scan_pattern_all_windows(seq_chars, valid_acgt, len_s, pattern, max_mismatch)
	if (is.null(res)) return(NULL)
	res <- filter_hits_by_distance(res, len_s, db)
	if (is.null(res)) return(NULL)
	res[order(res$mismatch, -res$pos_start), ][1, ]
}

scan_consensus_motif <- function(seq_upstream, patterns, max_mismatch = 0L, db = NULL) {
	seq_upstream <- toupper(seq_upstream)
	len_s <- nchar(seq_upstream)
	seq_chars  <- strsplit(seq_upstream, "")[[1]]
	valid_acgt <- seq_chars %in% c("A", "C", "G", "T")
	all_hits <- purrr::map_dfr(toupper(patterns), function(p)
		.scan_pattern_all_windows(seq_chars, valid_acgt, len_s, p, max_mismatch))
	if (is.null(all_hits) || nrow(all_hits) == 0) return(NULL)
	all_hits <- filter_hits_by_distance(all_hits, len_s, db)
	if (is.null(all_hits)) return(NULL)
	all_hits[order(all_hits$mismatch, -all_hits$pos_start), ][1, ]
}


# ══════════════════════════════════════════════════════════════════════════
# 4. EXCLUSION DES ALLÈLES  À NUMÉROTATION ROMAINE
# ══════════════════════════════════════════════════════════════════════════

is_roman_clan_allele <- function(allele_or_gene) {
	x <- toupper(trimws(as.character(allele_or_gene)))
	std_form   <- grepl("^(IGHV|IGKV|IGLV)[IVXLCM]+([-*]|$)", x)
	paren_form <- grepl("^(IGHV|IGKV|IGLV)\\([IVXLCM]+\\)", x)
	std_form | paren_form
}

is_roman_clan_subgroup <- function(subgroup) {
	x <- toupper(trimws(as.character(subgroup)))
	suffix <- sub("^(IGHV|IGKV|IGLV)", "", x)
	suffix <- gsub("[()]", "", suffix)
	nchar(suffix) > 0 & grepl("^[IVXLCM]+$", suffix)
}

exclude_roman_clan_alleles <- function(df, id_col = "ALLELE_ID") {
	if (is.null(df) || nrow(df) == 0 || !id_col %in% names(df)) return(df)
	excluded <- is_roman_clan_allele(df[[id_col]])
	if (any(excluded)) attr(df, "n_excluded_roman") <- sum(excluded)
	df[!excluded, , drop = FALSE]
}


# ══════════════════════════════════════════════════════════════════════════
# 5. JEU DE RÉFÉRENCE EMBARQUÉ
# ══════════════════════════════════════════════════════════════════════════

load_bundled_reference <- function(data_path = "data/motif_reference_calibration_UPDATED.csv") {
	if (!file.exists(data_path)) {
		return(list(data = NULL, version = NA_character_, has_upstream_seq = FALSE,
			error = sprintf("Reference file not found: %s", data_path)))
	}
	df <- tryCatch(utils::read.csv(data_path, stringsAsFactors = FALSE), error = function(e) NULL)
	if (is.null(df) || nrow(df) == 0) {
		return(list(data = NULL, version = NA_character_, has_upstream_seq = FALSE,
			error = "Reference file is empty or unreadable."))
	}
	df <- df[df$LOCI %in% SCANNER_LOCI, , drop = FALSE]
	if ("ALLELE" %in% names(df)) df <- df[!is_roman_clan_allele(df$ALLELE), , drop = FALSE]
	version <- format(file.info(data_path)$mtime, "%Y-%m-%d")
	list(data = df, version = version, has_upstream_seq = "SEQ_UPSTREAM" %in% names(df), error = NULL)
}


# ══════════════════════════════════════════════════════════════════════════
# 6. MODÈLE DE FOND — fréquences A/C/G/T locales (méthode PWM uniquement)
# ══════════════════════════════════════════════════════════════════════════

compute_background_freqs <- function(seqs) {
	seqs <- toupper(seqs); seqs <- seqs[nchar(seqs) > 0]
	if (length(seqs) == 0) return(c(A = 0.25, C = 0.25, G = 0.25, T = 0.25))
	all_chars <- unlist(strsplit(paste(seqs, collapse = ""), ""))
	tab   <- table(factor(all_chars, levels = c("A", "C", "G", "T")))
	freqs <- as.numeric(tab) / sum(tab)
	names(freqs) <- c("A", "C", "G", "T")
	freqs[freqs < 0.01] <- 0.01
	freqs / sum(freqs)
}


# ══════════════════════════════════════════════════════════════════════════
# 7. PWM — construction, scoring, consensus empirique, score maximum
# ══════════════════════════════════════════════════════════════════════════

build_pwm <- function(seqs, bg_freqs, pseudocount = 0.5, min_n = 15L) {
	seqs <- toupper(trimws(seqs))
	seqs <- seqs[grepl("^[ACGT]+$", seqs)]
	if (length(seqs) == 0) return(NULL)
	
	len_tab  <- table(nchar(seqs))
	len_mode <- as.integer(names(len_tab)[which.max(len_tab)])
	seqs     <- seqs[nchar(seqs) == len_mode]
	if (length(seqs) < min_n) return(NULL)
	
	bases <- c("A", "C", "G", "T")
	mat_chars <- do.call(rbind, strsplit(seqs, ""))
	counts <- matrix(0, nrow = 4, ncol = len_mode, dimnames = list(bases, NULL))
	for (j in seq_len(len_mode)) {
		tb <- table(factor(mat_chars[, j], levels = bases))
		counts[, j] <- as.numeric(tb)
	}
	n       <- length(seqs)
	freqs   <- (counts + pseudocount) / (n + 4 * pseudocount)
	logodds <- log2(sweep(freqs, 1, bg_freqs[bases], "/"))
	
	list(freqs = freqs, logodds = logodds, length = len_mode, n_seqs = n, bg = bg_freqs[bases])
}

empirical_consensus_from_pwm <- function(pwm) {
	if (is.null(pwm)) return(NULL)
	bases <- rownames(pwm$freqs)
	paste(bases[apply(pwm$freqs, 2, which.max)], collapse = "")
}





max_pwm_score_over_sequence <- function(seq, pwm) {
	seq <- toupper(seq)
	len_c <- pwm$length; len_s <- nchar(seq)
	if (len_s < len_c) return(list(max_score = NA_real_, best_pos_start = NA_integer_))
	n_windows <- len_s - len_c + 1L
	
	seq_chars <- strsplit(seq, "")[[1]]
	bases     <- rownames(pwm$logodds)   # "A","C","G","T", dans cet ordre
	
	idx <- outer(seq_len(len_c), seq_len(n_windows), FUN = "+") - 1L   
	char_mat <- matrix(seq_chars[idx], nrow = len_c, ncol = n_windows)
	base_idx <- matrix(match(char_mat, bases), nrow = len_c, ncol = n_windows)  # NA si hors ACGT
	
	valid_cols <- .colSums(is.na(base_idx), m = len_c, n = n_windows) == 0L
	scores <- rep(NA_real_, n_windows)
	if (any(valid_cols)) {
		pos_mat  <- matrix(rep(seq_len(len_c), n_windows), nrow = len_c, ncol = n_windows)
		# pwm$logodds a pour dimensions (base, position) = (4, len_c) ; on veut
		# logodds[base_idx[j,i], j] pour chaque (j = position dans le motif, i = fenêtre)
		# -> indice linéaire (matrice column-major, nrow = 4 = length(bases))
		lin_idx  <- base_idx + (pos_mat - 1L) * length(bases)
		vals <- matrix(NA_real_, nrow = len_c, ncol = n_windows)
		# NB : `lin_idx[, valid_cols]` est une matrice à 2 dimensions ;
		# (ligne, colonne) et non comme des indices linéairesd'où
		# as.vector() avant l'extraction, puis reconstruction de la forme
		vals[, valid_cols] <- matrix(pwm$logodds[as.vector(lin_idx[, valid_cols, drop = FALSE])],
			nrow = len_c, ncol = sum(valid_cols))
		scores[valid_cols] <- .colSums(vals[, valid_cols, drop = FALSE], m = len_c, n = sum(valid_cols))
	}
	
	best <- which.max(scores)
	if (length(best) == 0 || is.na(scores[best]))
		return(list(max_score = NA_real_, best_pos_start = NA_integer_))
	list(max_score = scores[best], best_pos_start = best)
}

hamming_distance <- function(s1, s2) {
	if (nchar(s1) != nchar(s2)) return(Inf)
	sum(strsplit(toupper(s1), "")[[1]] != strsplit(toupper(s2), "")[[1]])
}


# ══════════════════════════════════════════════════════════════════════════
# 8. DISTRIBUTION NULLE DU SCORE MAXIMUM (méthode PWM) — corrige les tests
#    multiples  
# ══════════════════════════════════════════════════════════════════════════

null_max_score_distribution <- function(pwm, seq_length, bg_freqs, n_sim = 5000L, seed = 1L) {
	if (seq_length < pwm$length) return(numeric(0))
	old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
	set.seed(seed)
	bases <- c("A", "C", "G", "T")
	probs <- bg_freqs[bases]
	sims <- vapply(seq_len(n_sim), function(i) {
		s <- paste(sample(bases, seq_length, replace = TRUE, prob = probs), collapse = "")
		max_pwm_score_over_sequence(s, pwm)$max_score
	}, numeric(1))
	if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
	sims[!is.na(sims)]
}

make_null_distribution_cache <- function(pwm, bg_freqs, n_sim = 5000L, bucket_width = 10L, seed = 1L) {
	cache <- new.env(parent = emptyenv())
	function(seq_length) {
		bucket <- as.character(round(seq_length / bucket_width) * bucket_width)
		if (!exists(bucket, envir = cache)) {
			bucket_len <- max(pwm$length, as.integer(bucket))
			assign(bucket, null_max_score_distribution(pwm, bucket_len, bg_freqs, n_sim, seed), envir = cache)
		}
		get(bucket, envir = cache)
	}
}

empirical_pvalue <- function(observed_score, null_dist) {
	if (is.na(observed_score) || length(null_dist) == 0) return(NA_real_)
	(sum(null_dist >= observed_score) + 1) / (length(null_dist) + 1)
}


# ══════════════════════════════════════════════════════════════════════════
# 9. RÈGLES DE COMPOSITION PARTAGÉES (E_BOX, TATA_BOX) — identiques dans les
#    deux méthodes, car ce ne sont pas des motifs modélisés statistiquement.
# ══════════════════════════════════════════════════════════════════════════

# E-BOX : CA-N-N-TG, longueur 6, 0 mismatch Implémenté  via le moteur
# IUPAC générique (pattern "CANNTG") réutilisé tel quel comme règle "PWM"/production.
find_ebox_rule <- function(seq_upstream, db = NULL) {
	scan_consensus_pattern(seq_upstream, "CANNTG", max_mismatch = 0L, db = db)
}

# TATA-BOX : région de 5 à 8 nt majoritairement A/T, avec jusqu'à
# `max_mismatch` base(s) non-A/T tolérée(s) dans la fenêtre (ex.
# "TACTATAA" reste valide malgré le C — cas réel observé dans le jeu de
# référence). Filtre par plausibilité de distance  On
# retient la fenêtre valide la plus longue, puis la moins mismatchée puis la plus proche du L-PART1

find_tata_box_rule <- function(seq_upstream, min_w = 5L, max_w = 8L, max_mismatch = 1L, db = NULL) {
	seq_upstream <- toupper(seq_upstream)
	len_s <- nchar(seq_upstream)
	if (len_s < min_w) return(NULL)
	
	seq_chars  <- strsplit(seq_upstream, "")[[1]]
	is_at      <- seq_chars %in% c("A", "T")
	is_valid   <- seq_chars %in% c("A", "C", "G", "T")
	cs_at      <- c(0L, cumsum(is_at))
	cs_valid   <- c(0L, cumsum(as.integer(is_valid)))
	
	hits <- list(); k <- 0L
	for (w in max_w:min_w) {
		if (len_s < w) next
		starts <- seq_len(len_s - w + 1L)
		ends   <- starts + w - 1L
		at_count    <- cs_at[ends + 1L]    - cs_at[starts]
		valid_count <- cs_valid[ends + 1L] - cs_valid[starts]
		mm   <- w - at_count
		keep <- which(mm <= max_mismatch & valid_count == w)
		if (length(keep) > 0L) {
			k <- k + 1L
			hits[[k]] <- data.frame(
				pos_start   = starts[keep], pos_end = ends[keep],
				matched_seq = vapply(keep, function(j) paste(seq_chars[starts[j]:ends[j]], collapse = ""), character(1)),
				width = w, mismatch = mm[keep], stringsAsFactors = FALSE)
		}
	}
	if (k == 0L) return(NULL)
	res <- dplyr::bind_rows(hits)
	res <- filter_hits_by_distance(res, len_s, db)
	if (is.null(res)) return(NULL)
	res[order(-res$width, res$mismatch, -res$pos_start), ][1, ]
}

# Filtre dur de mismatch méthode PWM uniquement S'applique EN PLUS du critère de q-value
passes_hard_mismatch <- function(matched_seq, empirical_consensus, motif_name) {
	max_mm <- HARD_MISMATCH_MAX[[motif_name]]
	if (is.null(max_mm) || is.null(empirical_consensus)) return(TRUE)
	if (nchar(matched_seq) != nchar(empirical_consensus)) return(FALSE)
	hamming_distance(matched_seq, empirical_consensus) <= max_mm
}


# ══════════════════════════════════════════════════════════════════════════
# 9bis. PLAUSIBILITÉ POSITIONNELLE — bornes de distance au L-PART1
#
# Un candidat qui satisfait la composition  ou le score mais se trouve à une
# distance du L-PART1 totalement hors de ce qui est observé dans le jeu de
# référence n'est très probablement PAS le bon motif (élément régulateur
# homologue trouvé par hasard ailleurs dans la région amont). On dérive donc,
# par (locus, motif), un intervalle de distance plausible 
# à partir des instances confirmées du jeu de référence avec une marge de
# sécurité Ce filtre s'applique EN plus des critères existants (mismatch/
# q-value) — un candidat doit satisfaire tous les critères pertinents
# ══════════════════════════════════════════════════════════════════════════

DISTANCE_MARGIN_NT <- 15L  # marge de sécurité ajoutée de part et d'autre du min/max observé

.motif_dist_col <- function(ref_data, motif_name) {
	app_col  <- APP_MOTIF_COLUMN_MAP[[motif_name]]
	dist_col <- paste0("DIST_", motif_name)
	if (!dist_col %in% names(ref_data)) dist_col <- paste0("DIST_", app_col)
	dist_col
}

# Calcule les bornes de distance plausibles par (locus, motif) à partir d'un
# sous-ensembl du jeu de référence

compute_distance_bounds <- function(ref_data, margin_nt = DISTANCE_MARGIN_NT) {
	bounds <- list()
	for (locus in SCANNER_LOCI) {
		bounds[[locus]] <- list()
		motifs <- ALL_MOTIFS_BY_LOCUS[[locus]]
		if (is.null(motifs)) next
		for (motif_name in motifs) {
			dist_col <- .motif_dist_col(ref_data, motif_name)
			if (!dist_col %in% names(ref_data)) {
				bounds[[locus]][[motif_name]] <- list(available = FALSE, dist_min = NA_real_, dist_max = NA_real_, n = 0L,
					warning = sprintf("Distance %s missing from the reference dataset for %s — positional plausibility filter disabled for this motif.", dist_col, locus))
				next
			}
			d <- suppressWarnings(as.numeric(ref_data[[dist_col]][ref_data$LOCI == locus]))
			d <- d[!is.na(d)]
			if (length(d) < 3L) {
				bounds[[locus]][[motif_name]] <- list(available = FALSE, dist_min = NA_real_, dist_max = NA_real_, n = length(d),
					warning = sprintf("Distance %s (%s): only %d confirmed value(s), not enough to derive reliable bounds — filter disabled.", dist_col, locus, length(d)))
				next
			}
			bounds[[locus]][[motif_name]] <- list(available = TRUE,
				dist_min = max(0, min(d) - margin_nt), dist_max = max(d) + margin_nt, n = length(d), warning = NULL)
		}
	}
	bounds
}

collect_distance_bound_warnings <- function(dist_bounds) {
	unlist(lapply(dist_bounds, function(loc) lapply(loc, function(x) x$warning)), use.names = FALSE)
}

#  Filtre un data.frame de candidats par plausibilité de distance au L-PART1, calculée directement à partir de la
# longueur de la région amont  distance = len_s - pos_end + 1

filter_hits_by_distance <- function(hits_df, len_s, db) {
	if (is.null(hits_df) || nrow(hits_df) == 0) return(NULL)
	if (is.null(db) || !isTRUE(db$available)) return(hits_df)   # pas de bornes -> pas de filtrage
	dist <- len_s - hits_df$pos_end + 1L
	kept <- hits_df[dist >= db$dist_min & dist <= db$dist_max, , drop = FALSE]
	if (nrow(kept) == 0) return(NULL)
	kept
}

#  Fenêtre de recherche restreinte à la zone de distance plausible, pour la
#  méthode PWM : la région amont est réduite à la sous-séquence contenant
#  TOUTES les fenêtres candidates de longueur  dont la distance
# au L-PART1 tombe dans [dist_min, dist_max]. Le score maximum ET la
#' distribution nulle  sont alors calculés sur cette zone réduite
#' — pas sur toute la région amont — ce qui rend le test statistique
#' cohérent avec l'espace de recherche biologiquement plausible 
restrict_to_distance_zone <- function(seq_upstream, motif_length, db) {
	len_s <- nchar(seq_upstream)
	if (is.null(db) || !isTRUE(db$available)) return(list(zone_seq = seq_upstream, offset = 0L))
	pe_min <- max(motif_length, len_s - db$dist_max + 1L)
	pe_max <- min(len_s, len_s - db$dist_min + 1L)
	if (pe_min > pe_max) return(NULL)   # aucune fenêtre plausible pour cet allèle
	ps_min <- max(1L, pe_min - motif_length + 1L)
	list(zone_seq = substr(seq_upstream, ps_min, pe_max), offset = ps_min - 1L)
}


# ══════════════════════════════════════════════════════════════════════════
# 10. ANCRAGE DU L-PART1 
# ══════════════════════════════════════════════════════════════════════════

locate_lpart1_exact <- function(seq_utr5, seq_lpart1) {
	seq_utr5   <- toupper(seq_utr5)
	seq_lpart1 <- toupper(seq_lpart1)
	if (nchar(seq_lpart1) == 0 || nchar(seq_utr5) < nchar(seq_lpart1))
		return(list(found = FALSE, pos_start = NA_integer_, pos_end = NA_integer_))
	all_pos <- gregexpr(seq_lpart1, seq_utr5, fixed = TRUE)[[1]]
	if (all_pos[1] == -1L) return(list(found = FALSE, pos_start = NA_integer_, pos_end = NA_integer_))
	best_start <- max(all_pos)
	list(found = TRUE, pos_start = best_start, pos_end = best_start + nchar(seq_lpart1) - 1L)
}
locate_lpart1_tolerant <- function(seq_utr5, seq_lpart1, max_mm_frac = 0.10) {
	ex <- locate_lpart1_exact(seq_utr5, seq_lpart1)
	if (ex$found) return(c(ex, list(mismatch = 0L)))
	
	nf <- list(found = FALSE, pos_start = NA_integer_, pos_end = NA_integer_, mismatch = NA_integer_)
	s <- strsplit(toupper(seq_utr5), "")[[1]]
	p <- strsplit(toupper(seq_lpart1), "")[[1]]
	lp <- length(p); ls <- length(s)
	if (lp == 0L || ls < lp) return(nf)
	
	n_win <- ls - lp + 1L
	idx <- outer(seq_len(lp), seq_len(n_win), "+") - 1L
	mm  <- .colSums(matrix(s[idx], lp, n_win) != p, lp, n_win)
	best <- max(which(mm == min(mm)))                 
	if (mm[best] > floor(lp * max_mm_frac)) return(nf)
	list(found = TRUE, pos_start = best, pos_end = best + lp - 1L, mismatch = as.integer(mm[best]))
}

locate_lpart1_by_coords <- function(row, max_ed_frac = 0.15) {
	nf <- list(found = FALSE)
	if (anyNA(c(row$FLANK5, row$GENOMIC_START, row$GENOMIC_START_LPART1, row$GENOMIC_END_LPART1))) return(nf)
	ps <- row$FLANK5 + (row$GENOMIC_START_LPART1 - row$GENOMIC_START) + 1
	pe <- ps + (row$GENOMIC_END_LPART1 - row$GENOMIC_START_LPART1)
	pe <- min(pe, nchar(row$SEQ_FULL))
	if (ps < 1 || ps > pe) return(nf)
	ed <- as.integer(utils::adist(toupper(substr(row$SEQ_FULL, ps, pe)), toupper(row$SEQ_LPART1)))
	if (ed > floor(nchar(row$SEQ_LPART1) * max_ed_frac)) return(nf)
	list(found = TRUE, pos_start = ps, pos_end = pe, mismatch = ed)
}

locate_lpart1_seed <- function(seq_utr5, seq_lpart1, seed_len = 20L, max_ed_frac = 0.15) {
	nf <- list(found = FALSE)
	s <- toupper(seq_utr5); p <- toupper(seq_lpart1)
	if (nchar(p) < seed_len) return(nf)
	pos <- gregexpr(substr(p, 1L, seed_len), s, fixed = TRUE)[[1]]
	if (pos[1] == -1L) return(nf)
	ps <- max(pos)
	pe <- min(ps + nchar(p) - 1L, nchar(s))
	ed <- as.integer(utils::adist(substr(s, ps, pe), p))
	if (ed > floor(nchar(p) * max_ed_frac)) return(nf)
	list(found = TRUE, pos_start = ps, pos_end = pe, mismatch = ed)
}
# ══════════════════════════════════════════════════════════════════════════
# 11. CALIBRATION PWM PAR LOCUS ET PAR SOUS-GROUPE 
# ══════════════════════════════════════════════════════════════════════════

compute_ref_subgroup <- function(ref_data) {
	if ("SUBGROUP" %in% names(ref_data)) return(as.character(ref_data$SUBGROUP))
	if ("GENE" %in% names(ref_data)) {
		if (exists("clean_subgroup", mode = "function")) return(clean_subgroup(ref_data$GENE))
		return(sub("^([A-Za-z]+[0-9]+).*$", "\\1", ref_data$GENE))
	}
	rep(NA_character_, nrow(ref_data))
}

build_motif_subgroup_calibration <- function(df_scope, seq_col, bg_freqs, min_n_pwm = 15L, min_n_subgroup = 8L) {
	get_clean_seqs <- function(mask) {
		if (!seq_col %in% names(df_scope)) return(character(0))
		s <- df_scope[[seq_col]][mask]
		toupper(s[!is.na(s) & s != "" & s != "absent"])
	}
	
	locus_seqs <- get_clean_seqs(rep(TRUE, nrow(df_scope)))
	locus_pwm  <- build_pwm(locus_seqs, bg_freqs, min_n = min_n_pwm)
	
	subgroup_cal <- list()
	sgs <- unique(df_scope$SUBGROUP[!is.na(df_scope$SUBGROUP) & df_scope$SUBGROUP != ""])
	for (sg in sgs) {
		sg_seqs <- get_clean_seqs(df_scope$SUBGROUP == sg)
		sg_pwm  <- build_pwm(sg_seqs, bg_freqs, min_n = min_n_subgroup)
		if (!is.null(sg_pwm))
			subgroup_cal[[sg]] <- list(pwm = sg_pwm, empirical_consensus = empirical_consensus_from_pwm(sg_pwm),
				n_training = sg_pwm$n_seqs)
	}
	
	list(
		locus_pwm = locus_pwm,
		locus_consensus  = if (!is.null(locus_pwm)) empirical_consensus_from_pwm(locus_pwm) else NULL,
		locus_n_training = if (!is.null(locus_pwm)) locus_pwm$n_seqs else 0L,
		subgroup_cal = subgroup_cal,
		n_total_seqs = length(locus_seqs)
	)
}

build_locus_calibration <- function(ref_data, locus, bg_freqs, min_n_pwm = 15L, min_n_subgroup = 8L) {
	motif_names <- PWM_MOTIFS_BY_LOCUS[[locus]]
	ref_loc <- ref_data[ref_data$LOCI == locus, , drop = FALSE]
	ref_loc$SUBGROUP <- compute_ref_subgroup(ref_loc)
	
	out <- list()
	for (motif_name in motif_names) {
		app_col <- APP_MOTIF_COLUMN_MAP[[motif_name]]
		seq_col <- paste0("SEQ_", motif_name)
		if (!seq_col %in% names(ref_loc)) seq_col <- paste0("SEQ_", app_col)
		
		cal <- build_motif_subgroup_calibration(ref_loc, seq_col, bg_freqs, min_n_pwm, min_n_subgroup)
		n_sg_ok <- length(cal$subgroup_cal)
		
		if (is.null(cal$locus_pwm) && n_sg_ok == 0) {
			cal$warning <- sprintf(
				"Motif '%s' (%s): %d confirmed instance(s) in total, minimum required %d — no PWM could be built (neither by subgroup nor at locus level).",
				motif_name, locus, cal$n_total_seqs, min_n_pwm)
		} else if (is.null(cal$locus_pwm) && n_sg_ok > 0) {
			cal$warning <- sprintf(
				"Motif '%s' (%s): locus-level pooled PWM could not be built (no fallback available for uncalibrated subgroups); %d subgroup(s) calibrated individually.",
				motif_name, locus, n_sg_ok)
		} else {
			cal$warning <- NULL
		}
		out[[motif_name]] <- cal
	}
	out
}


get_effective_calibration <- function(cal_motif, subgroup) {
	if (is.null(cal_motif)) return(NULL)
	if (!is.null(subgroup) && !is.na(subgroup) && subgroup %in% names(cal_motif$subgroup_cal)) {
		sg <- cal_motif$subgroup_cal[[subgroup]]
		return(list(pwm = sg$pwm, empirical_consensus = sg$empirical_consensus, source = "subgroup", source_id = subgroup))
	}
	if (!is.null(cal_motif$locus_pwm))
		return(list(pwm = cal_motif$locus_pwm, empirical_consensus = cal_motif$locus_consensus, source = "locus", source_id = "locus"))
	NULL
}

collect_calibration_warnings <- function(calib_by_locus) {
	unlist(lapply(calib_by_locus, function(loc_cal) lapply(loc_cal, function(x) x$warning)), use.names = FALSE)
}


# ══════════════════════════════════════════════════════════════════════════
# 12. FASTA — parseur et jointure 5'UTR <-> L-PART1
# ══════════════════════════════════════════════════════════════════════════

.split_fasta_blocks <- function(fasta_text) {
	if (is.null(fasta_text) || nchar(trimws(fasta_text)) == 0) return(NULL)
	lines <- trimws(strsplit(fasta_text, "\n")[[1]])
	lines <- lines[lines != ""]
	is_header <- startsWith(lines, ">")
	records <- list(); cur_hdr <- NULL; cur_seq <- character(0)
	flush <- function() {
		if (is.null(cur_hdr)) return(NULL)
		list(header = cur_hdr, seq = toupper(paste(cur_seq, collapse = "")))
	}
	for (i in seq_along(lines)) {
		if (is_header[i]) {
			rec <- flush(); if (!is.null(rec)) records[[length(records) + 1L]] <- rec
			cur_hdr <- lines[i]; cur_seq <- character(0)
		} else {
			clean <- gsub("[^ACGTNacgtn]", "", lines[i])
			if (nchar(clean) > 0) cur_seq <- c(cur_seq, clean)
		}
	}
	rec <- flush(); if (!is.null(rec)) records[[length(records) + 1L]] <- rec
	if (length(records) == 0) return(NULL)
	records
}

.extract_genomic_coords <- function(header) {
	m <- regmatches(header, regexpr("[0-9]{3,}\\.\\.[0-9]{3,}", header))
	if (length(m) == 0 || nchar(m) == 0) return(list(start = NA_real_, end = NA_real_))
	parts <- as.numeric(strsplit(m, "\\.\\.")[[1]])
	list(start = parts[1], end = parts[2])
}

parse_fasta <- function(fasta_text) {
	records <- .split_fasta_blocks(fasta_text)
	if (is.null(records)) return(NULL)
	
	df <- purrr::map_dfr(records, function(r) {
		hdr   <- sub("^>", "", r$header)
		parts <- trimws(strsplit(hdr, "\\|")[[1]])
		allele_idx <- which(grepl("^(IGHV|IGKV|IGLV)", parts))[1]
		if (is.na(allele_idx)) allele_idx <- min(2, length(parts))
		allele_id     <- parts[allele_idx]
		species       <- if (length(parts) >= allele_idx + 1) parts[allele_idx + 1] else NA_character_
		species       <- sub("_.*$", "", species)
		functionality <- if (length(parts) >= allele_idx + 2) parts[allele_idx + 2] else NA_character_
		functionality <- trimws(gsub("[\\(\\)\\[\\]]", "", functionality %||% ""))
		loci <- stringr::str_extract(allele_id, "^(IGHV|IGKV|IGLV)")
		gene <- sub("\\*.*$", "", allele_id)
		gcoord <- .extract_genomic_coords(hdr)
		fl <- regmatches(hdr, regexec("([0-9]+)nt\\+[0-9]+\\.\\.[0-9]+", hdr))[[1]]
		flank5 <- if (length(fl) == 2) as.numeric(fl[2]) else NA_real_
		data.frame(
			ALLELE_ID = allele_id, SPECIES = species %||% "?", GENE = gene,
			LOCI = loci %||% NA_character_, FUNCTIONALITY = functionality,
			SEQ_FULL = r$seq, LEN_FULL = nchar(r$seq),
			GENOMIC_START = gcoord$start, GENOMIC_END = gcoord$end,
			GENOMIC_START = gcoord$start, GENOMIC_END = gcoord$end,
			FLANK5 = flank5,
			stringsAsFactors = FALSE
		)
	})
	df <- df[!is.na(df$LOCI) & df$LOCI %in% SCANNER_LOCI, , drop = FALSE]
	exclude_roman_clan_alleles(df, "ALLELE_ID")
}

join_utr5_lpart1 <- function(df_utr5, df_lpart1) {
	norm <- function(x) toupper(trimws(x))
	df_utr5$JOIN_KEY   <- paste(norm(df_utr5$ALLELE_ID),   norm(df_utr5$SPECIES))
	df_lpart1$JOIN_KEY <- paste(norm(df_lpart1$ALLELE_ID), norm(df_lpart1$SPECIES))
	df_lpart1_sel <- df_lpart1 %>%
		dplyr::select(JOIN_KEY, SEQ_LPART1 = SEQ_FULL, LEN_LPART1 = LEN_FULL,
			GENOMIC_START_LPART1 = GENOMIC_START, GENOMIC_END_LPART1 = GENOMIC_END)
	matched <- dplyr::inner_join(df_utr5, df_lpart1_sel, by = "JOIN_KEY")
	list(
		matched = matched,
		unmatched_utr5   = df_utr5[!df_utr5$JOIN_KEY %in% df_lpart1$JOIN_KEY, , drop = FALSE],
		unmatched_lpart1 = df_lpart1[!df_lpart1$JOIN_KEY %in% df_utr5$JOIN_KEY, , drop = FALSE]
	)
}


# ══════════════════════════════════════════════════════════════════════════
# 13. ANCRAGE COMMUN — jointure + ancrage L-PART1 + extraction de la région
# ══════════════════════════════════════════════════════════════════════════

MAX_UPSTREAM_SCAN_NT <- 400L  

anchor_alleles <- function(df_utr5, df_lpart1) {
	join_res <- join_utr5_lpart1(df_utr5, df_lpart1)
	matched  <- join_res$matched
	
	diag_rows <- list()
	add_diag <- function(id, sp, gene, loci, status) {
		diag_rows[[length(diag_rows) + 1L]] <<- data.frame(
			ALLELE_ID = id, SPECIES = sp, GENE = gene, LOCI = loci,
			STATUS = status, stringsAsFactors = FALSE)
	}
	
	um <- join_res$unmatched_utr5
	if (nrow(um) > 0) for (i in seq_len(nrow(um)))
		add_diag(um$ALLELE_ID[i], um$SPECIES[i], um$GENE[i], um$LOCI[i], "1. Not in L-PART1 file (ALLELE+SPECIES)")
	ul <- join_res$unmatched_lpart1
	if (nrow(ul) > 0) for (i in seq_len(nrow(ul)))
		add_diag(ul$ALLELE_ID[i], ul$SPECIES[i], ul$GENE[i], ul$LOCI[i], "1b. Only in L-PART1 file")
	
	if (nrow(matched) > 0) {
		is_rom <- is_roman_clan_allele(matched$ALLELE_ID) | is_roman_clan_allele(matched$GENE)
		for (i in which(is_rom))
			add_diag(matched$ALLELE_ID[i], matched$SPECIES[i], matched$GENE[i], matched$LOCI[i], "2. Roman clan (excluded)")
		matched <- matched[!is_rom, , drop = FALSE]
	}
	if (nrow(matched) == 0)
		stop("No ALLELE+SPECIES match between the two FASTA files (or all alleles belong to Roman-numeral clans, which are excluded from scope).")
	
	anchored_env <- new.env(parent = emptyenv())
	anchored_env$rows <- list()
	
	anchored_rows_meta <- purrr::map_dfr(seq_len(nrow(matched)), function(i) {
		row <- matched[i, ]
		subgroup_val <- if (exists("clean_subgroup", mode = "function")) clean_subgroup(row$GENE)
		else sub("^([A-Za-z]+[0-9]+).*$", "\\1", row$GENE)
		if (is_roman_clan_subgroup(subgroup_val)) {
			add_diag(row$ALLELE_ID, row$SPECIES, row$GENE, row$LOCI, "2. Roman clan (excluded)"); return(NULL) }
		anchor <- locate_lpart1_by_coords(row)
		if (!anchor$found) anchor <- locate_lpart1_seed(row$SEQ_FULL, row$SEQ_LPART1)
		if (!anchor$found) anchor <- locate_lpart1_tolerant(row$SEQ_FULL, row$SEQ_LPART1, 0.15)
		if (!anchor$found) {
			add_diag(row$ALLELE_ID, row$SPECIES, row$GENE, row$LOCI, "3. L-PART1 not found in 5'UTR (exact match)"); return(NULL) }
		seq_upstream_full <- if (anchor$pos_start > 1L) substr(row$SEQ_FULL, 1L, anchor$pos_start - 1L) else ""
		if (nchar(seq_upstream_full) == 0) {
			add_diag(row$ALLELE_ID, row$SPECIES, row$GENE, row$LOCI, "4. Empty upstream region"); return(NULL) }
		
		add_diag(row$ALLELE_ID, row$SPECIES, row$GENE, row$LOCI,
			if (anchor$mismatch > 0L) "OK (anchored, with mismatches)" else "OK (anchored)")
		
		len_full <- nchar(seq_upstream_full)
		if (len_full > MAX_UPSTREAM_SCAN_NT) {
			trunc_offset <- len_full - MAX_UPSTREAM_SCAN_NT
			seq_upstream <- substr(seq_upstream_full, trunc_offset + 1L, len_full)
		} else {
			trunc_offset <- 0L
			seq_upstream <- seq_upstream_full
		}
		
		seq_lpart1_full <- substr(row$SEQ_FULL, anchor$pos_start, anchor$pos_end)
		seq_downstream  <- if (anchor$pos_end < nchar(row$SEQ_FULL))
			substr(row$SEQ_FULL, anchor$pos_end + 1L, nchar(row$SEQ_FULL)) else ""
		
		# ── tableau d'AFFICHAGE 
		anchored_env$rows[[length(anchored_env$rows) + 1L]] <- data.frame(
			ALLELE_ID = row$ALLELE_ID, SPECIES = row$SPECIES, GENE = row$GENE,
			LOCI = row$LOCI, SUBGROUP = subgroup_val,
			SEQ_UPSTREAM = seq_upstream, SEQ_LPART1 = seq_lpart1_full, SEQ_DOWNSTREAM = seq_downstream,
			UPSTREAM_TRUNCATED = len_full > MAX_UPSTREAM_SCAN_NT,
			TRUNC_OFFSET = trunc_offset,
			stringsAsFactors = FALSE
		)
		
		# ── tableau de SCAN 
		data.frame(ALLELE_ID = row$ALLELE_ID, SPECIES = row$SPECIES, GENE = row$GENE,
			LOCI = row$LOCI, SUBGROUP = subgroup_val, FUNCTIONALITY = row$FUNCTIONALITY,
			SEQ_UPSTREAM = seq_upstream,
			LPART1_POS_START = anchor$pos_start, LPART1_POS_END = anchor$pos_end,
			UPSTREAM_TRUNCATED = len_full > MAX_UPSTREAM_SCAN_NT,
			TRUNC_OFFSET = trunc_offset,
			stringsAsFactors = FALSE)
	})
	
	anchored_data <- if (length(anchored_env$rows) > 0) dplyr::bind_rows(anchored_env$rows) else NULL
	list(join_report = join_res, anchored_meta = anchored_rows_meta, anchored_data = anchored_data,
		diagnostic = dplyr::bind_rows(diag_rows))
}


# ══════════════════════════════════════════════════════════════════════════
# 13b SOUS-ONGLET 1 — PIPELINE DE PRODUCTION CONSENSUS 
# ══════════════════════════════════════════════════════════════════════════

scan_one_allele_consensus <- function(seq_upstream, locus, dist_bounds_locus = NULL) {
	if (nchar(seq_upstream) == 0) return(NULL)
	defs <- CONSENSUS_DEFINITIONS[[locus]]
	if (is.null(defs)) return(NULL)
	
	rows <- purrr::imap_dfr(defs, function(def, motif_name) {
		db  <- dist_bounds_locus[[motif_name]]
		hit <- scan_consensus_motif(seq_upstream, def$patterns, def$max_mismatch, db = db)
		if (is.null(hit)) return(NULL)
		data.frame(MOTIF = motif_name, pos_start = hit$pos_start, pos_end = hit$pos_end,
			matched_seq = hit$matched_seq, MISMATCH = hit$mismatch, DETECTED = TRUE,
			DETECTION_MODE = "CONSENSUS", stringsAsFactors = FALSE)
	})
	
	tata_hit <- find_tata_box_rule(seq_upstream, db = dist_bounds_locus[["TATA_BOX"]])
	if (!is.null(tata_hit)) {
		rows <- dplyr::bind_rows(rows, data.frame(MOTIF = "TATA_BOX",
			pos_start = tata_hit$pos_start, pos_end = tata_hit$pos_end,
			matched_seq = tata_hit$matched_seq, MISMATCH = tata_hit$mismatch, DETECTED = TRUE,
			DETECTION_MODE = "CONSENSUS", stringsAsFactors = FALSE))
	}
	rows
}

# méthode consensus Aucune calibration requise en entrée  applique le filtre de plausibilité positionnelle
# un candidat trop loin du L-PART1 par rapport à ce qui est observé dans le jeu de référence n'est pas retenu

run_consensus_scan <- function(df_utr5, df_lpart1, dist_bounds = NULL) {
	anch <- anchor_alleles(df_utr5, df_lpart1)
	meta <- anch$anchored_meta
	if (is.null(meta) || nrow(meta) == 0)
		return(list(join_report = anch$join_report, hits_detected = NULL, anchored_data = anch$anchored_data))
	
	all_hits <- purrr::map_dfr(seq_len(nrow(meta)), function(i) {
		row  <- meta[i, ]
		hits <- scan_one_allele_consensus(row$SEQ_UPSTREAM, row$LOCI, dist_bounds_locus = dist_bounds[[row$LOCI]])
		if (is.null(hits) || nrow(hits) == 0) return(NULL)
		# fonctions de scan travaillent sur la région amont tronquée
		hits$pos_start <- hits$pos_start + row$TRUNC_OFFSET
		hits$pos_end   <- hits$pos_end   + row$TRUNC_OFFSET
		hits$ALLELE_ID  <- row$ALLELE_ID
		hits$SPECIES    <- row$SPECIES
		hits$GENE       <- row$GENE
		hits$LOCI       <- row$LOCI
		hits$SUBGROUP   <- row$SUBGROUP
		hits$FUNCTIONALITY    <- row$FUNCTIONALITY
		hits$LPART1_POS_START <- row$LPART1_POS_START
		hits$LPART1_POS_END   <- row$LPART1_POS_END
		# POS_START_LPART1 - POS_END_motif, sans -1
		hits$DIST_TO_LPART1 <- row$LPART1_POS_START - hits$pos_end
		hits
	})
	
	list(join_report = anch$join_report, hits_detected = all_hits, anchored_data = anch$anchored_data)
}


# ══════════════════════════════════════════════════════════════════════════
# 14. SOUS-ONGLET 2 — PWM (probabiliste)
# ══════════════════════════════════════════════════════════════════════════

scan_one_allele_pwm <- function(seq_upstream, locus, subgroup, calib_locus, null_caches_env, dist_bounds_locus = NULL, n_sim = 5000L, bg_freqs = NULL) {
	if (nchar(seq_upstream) == 0) return(NULL)
	
	rule_rows <- purrr::map_dfr(RULE_MOTIFS_BY_LOCUS[[locus]], function(motif_name) {
		db  <- dist_bounds_locus[[motif_name]]
		hit <- if (motif_name == "E_BOX") find_ebox_rule(seq_upstream, db = db) else find_tata_box_rule(seq_upstream, db = db)
		if (is.null(hit)) return(NULL)
		data.frame(MOTIF = motif_name, pos_start = hit$pos_start, pos_end = hit$pos_end,
			matched_seq = hit$matched_seq, PWM_SCORE = NA_real_, P_VALUE = NA_real_,
			Q_VALUE = NA_real_, MISMATCH = as.integer(hit$mismatch %||% NA_integer_), DETECTED = TRUE,
			DETECTION_MODE = "RULE", HARD_FILTER_PASS = NA, stringsAsFactors = FALSE)
	})
	
	pwm_rows <- purrr::map_dfr(PWM_MOTIFS_BY_LOCUS[[locus]], function(motif_name) {
		cal_motif <- calib_locus[[motif_name]]
		eff <- get_effective_calibration(cal_motif, subgroup)
		if (is.null(eff)) return(NULL)
		
		db <- dist_bounds_locus[[motif_name]]
		zi <- restrict_to_distance_zone(seq_upstream, eff$pwm$length, db)
		if (is.null(zi)) return(NULL)   # aucune fenêtre plausible pour cet allèle
		
		r <- max_pwm_score_over_sequence(zi$zone_seq, eff$pwm)
		if (is.na(r$max_score)) return(NULL)
		
		cache_key <- paste(motif_name, eff$source_id, sep = "::")
		if (!exists(cache_key, envir = null_caches_env)) {
			assign(cache_key, make_null_distribution_cache(eff$pwm, bg_freqs, n_sim = n_sim), envir = null_caches_env)
		}
		null_cache_fn <- get(cache_key, envir = null_caches_env)
		pv <- empirical_pvalue(r$max_score, null_cache_fn(nchar(zi$zone_seq)))
		
		matched <- substr(zi$zone_seq, r$best_pos_start, r$best_pos_start + eff$pwm$length - 1L)
		mm <- if (!is.null(eff$empirical_consensus) && nchar(matched) == nchar(eff$empirical_consensus))
			hamming_distance(matched, eff$empirical_consensus) else NA_integer_
		hard_ok <- passes_hard_mismatch(matched, eff$empirical_consensus, motif_name)
		global_pos_start <- r$best_pos_start + zi$offset
		data.frame(MOTIF = motif_name, pos_start = global_pos_start,
			pos_end = global_pos_start + eff$pwm$length - 1L, matched_seq = matched,
			PWM_SCORE = round(r$max_score, 3), P_VALUE = pv, Q_VALUE = NA_real_, MISMATCH = mm,
			DETECTED = NA, DETECTION_MODE = "PWM", HARD_FILTER_PASS = hard_ok,
			CALIBRATION_SOURCE = eff$source, stringsAsFactors = FALSE)
	})
	
	dplyr::bind_rows(rule_rows, pwm_rows)
}

run_motif_scan_pwm <- function(df_utr5, df_lpart1, calibration_by_locus, q_threshold = 0.05, n_sim = 5000L, dist_bounds = NULL) {
	anch <- anchor_alleles(df_utr5, df_lpart1)
	meta <- anch$anchored_meta
	if (is.null(meta) || nrow(meta) == 0)
		return(list(join_report = anch$join_report, hits_all = NULL, hits_detected = NULL, anchored_data = anch$anchored_data))
	
	
	bg_by_locus <- list()
	null_caches_env_by_locus <- list()
	for (loc in SCANNER_LOCI) {
		rows_loc <- meta[meta$LOCI == loc, , drop = FALSE]
		if (nrow(rows_loc) == 0) next
		bg_by_locus[[loc]] <- compute_background_freqs(rows_loc$SEQ_UPSTREAM)
		null_caches_env_by_locus[[loc]] <- new.env(parent = emptyenv())
	}
	
	all_hits <- purrr::map_dfr(seq_len(nrow(meta)), function(i) {
		row  <- meta[i, ]
		hits <- scan_one_allele_pwm(row$SEQ_UPSTREAM, row$LOCI, row$SUBGROUP,
			calibration_by_locus[[row$LOCI]], null_caches_env_by_locus[[row$LOCI]],
			dist_bounds_locus = dist_bounds[[row$LOCI]], n_sim = n_sim, bg_freqs = bg_by_locus[[row$LOCI]])
		if (is.null(hits) || nrow(hits) == 0) return(NULL)
		hits$pos_start <- hits$pos_start + row$TRUNC_OFFSET
		hits$pos_end   <- hits$pos_end   + row$TRUNC_OFFSET
		hits$ALLELE_ID  <- row$ALLELE_ID
		hits$SPECIES    <- row$SPECIES
		hits$GENE       <- row$GENE
		hits$LOCI       <- row$LOCI
		hits$SUBGROUP   <- row$SUBGROUP
		hits$FUNCTIONALITY    <- row$FUNCTIONALITY
		hits$LPART1_POS_START <- row$LPART1_POS_START
		hits$LPART1_POS_END   <- row$LPART1_POS_END
		hits$DIST_TO_LPART1   <- row$LPART1_POS_START - hits$pos_end
		hits
	})
	
	if (is.null(all_hits) || nrow(all_hits) == 0)
		return(list(join_report = anch$join_report, hits_all = NULL, hits_detected = NULL, anchored_data = anch$anchored_data))
	
	pwm_rows  <- all_hits[all_hits$DETECTION_MODE == "PWM", , drop = FALSE]
	rule_rows <- all_hits[all_hits$DETECTION_MODE == "RULE", , drop = FALSE]
	
	if (nrow(pwm_rows) > 0) {
		pwm_rows <- pwm_rows %>%
			dplyr::mutate(match_len = pos_end - pos_start + 1L) %>%
			dplyr::group_by(LOCI, MOTIF, match_len) %>%
			dplyr::mutate(Q_VALUE = stats::p.adjust(P_VALUE, method = "BH")) %>%
			dplyr::ungroup() %>%
			dplyr::mutate(DETECTED = !is.na(Q_VALUE) & Q_VALUE <= q_threshold & (HARD_FILTER_PASS %in% TRUE)) %>%
			dplyr::select(-match_len)
	}
	if (nrow(rule_rows) > 0) rule_rows$Q_VALUE <- NA_real_
	
	all_hits <- dplyr::bind_rows(pwm_rows, rule_rows)
	hits_detected <- all_hits %>% dplyr::filter(is.na(DETECTED) | DETECTED)
	
	list(join_report = anch$join_report, hits_all = all_hits, hits_detected = hits_detected,
		anchored_data = anch$anchored_data)
}
truncate_upstream_tail <- function(seqs, max_nt = MAX_UPSTREAM_SCAN_NT) {
	n <- nchar(seqs)
	ifelse(n > max_nt, substr(seqs, n - max_nt + 1L, n), seqs)
}


# ══════════════════════════════════════════════════════════════════════════
# 15. SOUS-ONGLET 3 — COMPARAISON DES MÉTHODES / ÉVALUATION DE PERFORMANCE
#
# Les deux méthodes sont évaluées contre les instances confirmées du jeu de
# référence, PAS l'une contre l'autre
# Définition des VP/FP/VN/FN, identique aux deux méthodes pour rester
# comparable :
#   - VP (vrai positif) : instance confirmée du jeu de référence -> la
#     méthode détecte un hit dans sa région amont réelle 
#   - FN (faux négatif) : instance confirmée -> la méthode NE détecte PAS
#     de hit.
#   - contrôle négatif : même région amont, nucléotides permutés (mélange
#     aléatoire, préserve la composition locale) -> VN si aucun hit, FP si
#     un hit est trouvé 
# Méthode PWM : validation croisée k=5 stricte — à chaque fold la PWM est
# reconstruite UNIQUEMENT sur les k-1 folds d'entraînement (le fond A/C/G/T
# est également recalculé par fold sur les séquences d'entraînement, pour
# éviter toute fuite d'information), puis testée sur le fold laissé de côté
#
# Méthode consensus : AUCUN paramètre n'est appris sur les données (le
# pattern est fixe, issu de la littérature) — il n'y a donc structurellement
# aucun risque de fuite, et pas besoin de k-fold : la performance est
# évaluée directement sur toutes les instances disponibles du jeu de
# référence, en une seule passe. C'est un choix délibéré pour rester simple
# et transparent, pas un raccourci méthodologique : un k-fold sur une
# méthode sans paramètre entraîné donnerait k fois la même réponse
# ══════════════════════════════════════════════════════════════════════════

.motif_seq_col <- function(ref_data, motif_name) {
	app_col <- APP_MOTIF_COLUMN_MAP[[motif_name]]
	seq_col <- paste0("SEQ_", motif_name)
	if (!seq_col %in% names(ref_data)) seq_col <- paste0("SEQ_", app_col)
	seq_col
}

.compute_confusion_metrics <- function(tp, fn, tn, fp) {
	sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
	spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
	prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
	f1   <- if (!is.na(prec) && !is.na(sens) && (prec + sens) > 0) 2 * prec * sens / (prec + sens) else NA_real_
	data.frame(sensitivity = sens, specificity = spec, precision = prec, f1_score = f1,
		tp = tp, fn = fn, tn = tn, fp = fp)
}

# Détection consensus sur une séquence amont isolée 
.consensus_call <- function(seq_upstream, locus, motif_name, db = NULL) {
	if (motif_name == "TATA_BOX") return(!is.null(find_tata_box_rule(seq_upstream, db = db)))
	if (motif_name == "E_BOX")    return(!is.null(find_ebox_rule(seq_upstream, db = db)))
	def <- CONSENSUS_DEFINITIONS[[locus]][[motif_name]]
	if (is.null(def)) return(FALSE)
	!is.null(scan_consensus_motif(seq_upstream, def$patterns, def$max_mismatch, db = db))
}

# Évaluation de la méthode CONSENSUS 

run_consensus_evaluation <- function(ref_data, locus, motif_name, dist_bounds = NULL, seed = 42L) {
	seq_col <- .motif_seq_col(ref_data, motif_name)
	if (!"SEQ_UPSTREAM" %in% names(ref_data))
		return(list(done = FALSE, reason = "SEQ_UPSTREAM column missing from the reference dataset (see calibration status)."))
	if (!seq_col %in% names(ref_data))
		return(list(done = FALSE, reason = sprintf("Column %s missing from the reference dataset.", seq_col)))
	
	df_loc <- ref_data[ref_data$LOCI == locus, , drop = FALSE]
	df_loc <- df_loc[!is.na(df_loc[[seq_col]]) & df_loc[[seq_col]] != "" &
			!is.na(df_loc$SEQ_UPSTREAM) & df_loc$SEQ_UPSTREAM != "", , drop = FALSE]
	if (nrow(df_loc) == 0)
		return(list(done = FALSE, reason = sprintf("No usable confirmed instance for (%s, %s) — check SEQ_UPSTREAM and %s in the reference dataset.", locus, motif_name, seq_col)))
	df_loc$SEQ_UPSTREAM <- truncate_upstream_tail(df_loc$SEQ_UPSTREAM)
	
	db <- dist_bounds[[locus]][[motif_name]]
	
	old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) .Random.seed else NULL
	set.seed(seed)
	
	pos_calls <- vapply(seq_len(nrow(df_loc)), function(i)
		.consensus_call(toupper(df_loc$SEQ_UPSTREAM[i]), locus, motif_name, db = db), logical(1))
	neg_calls <- vapply(seq_len(nrow(df_loc)), function(i) {
		s <- toupper(df_loc$SEQ_UPSTREAM[i])
		shuffled <- paste(sample(strsplit(s, "")[[1]]), collapse = "")
		.consensus_call(shuffled, locus, motif_name, db = db)
	}, logical(1))
	
	if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
	
	tp <- sum(pos_calls); fn <- sum(!pos_calls); tn <- sum(!neg_calls); fp <- sum(neg_calls)
	metrics <- .compute_confusion_metrics(tp, fn, tn, fp)
	list(done = TRUE, method = "Consensus", locus = locus, motif = motif_name,
		n_test = nrow(df_loc), metrics = metrics)
}

# Évaluation de la méthode PWM 

run_pwm_evaluation <- function(ref_data, locus, motif_name, k = 5L, q_threshold = 0.05, n_sim = 2000L, seed = 42L) {
	seq_col <- .motif_seq_col(ref_data, motif_name)
	if (!"SEQ_UPSTREAM" %in% names(ref_data))
		return(list(done = FALSE, reason = "SEQ_UPSTREAM column missing from the reference dataset (see calibration status)."))
	if (!seq_col %in% names(ref_data))
		return(list(done = FALSE, reason = sprintf("Column %s missing from the reference dataset.", seq_col)))
	
	df_loc <- ref_data[ref_data$LOCI == locus, , drop = FALSE]
	df_loc <- df_loc[!is.na(df_loc[[seq_col]]) & df_loc[[seq_col]] != "" &
			!is.na(df_loc$SEQ_UPSTREAM) & df_loc$SEQ_UPSTREAM != "", , drop = FALSE]
	df_loc$SEQ_UPSTREAM <- truncate_upstream_tail(df_loc$SEQ_UPSTREAM)
	df_loc$SUBGROUP <- compute_ref_subgroup(df_loc)
	n_total <- nrow(df_loc)
	if (n_total < k)
		return(list(done = FALSE, reason = sprintf("n=%d usable instances for (%s, %s), insufficient for k=%d folds — check SEQ_UPSTREAM and %s.", n_total, locus, motif_name, k, seq_col)))
	
	set.seed(seed)
	fold_id <- sample(rep_len(seq_len(k), n_total))
	
	tp_tot <- fn_tot <- tn_tot <- fp_tot <- 0L
	n_folds_ok <- 0L
	diag_last <- NULL   # diagnostic du dernier fold tenté, pour un message d'échec actionnable
	for (fold in seq_len(k)) {
		train <- df_loc[fold_id != fold, , drop = FALSE]
		test  <- df_loc[fold_id == fold, , drop = FALSE]
		
		bg_freqs <- compute_background_freqs(train$SEQ_UPSTREAM)
		
		cal <- build_motif_subgroup_calibration(train, seq_col, bg_freqs, min_n_pwm = 15L, min_n_subgroup = 8L)
		if (is.null(cal$locus_pwm) && length(cal$subgroup_cal) == 0) {
			seqs_avail <- train[[seq_col]][!is.na(train[[seq_col]]) & train[[seq_col]] != "" & train[[seq_col]] != "absent"]
			seqs_acgt  <- toupper(seqs_avail)[grepl("^[ACGT]+$", toupper(seqs_avail))]
			mode_len   <- if (length(seqs_acgt) > 0) as.integer(names(sort(table(nchar(seqs_acgt)), decreasing = TRUE))[1]) else NA_integer_
			n_mode     <- if (length(seqs_acgt) > 0) sum(nchar(seqs_acgt) == mode_len) else 0L
			diag_last <- sprintf("train n=%d ; %d with %s filled in, %d valid ACGT, %d at the dominant length (%s nt) — need 15 min.",
				nrow(train), length(seqs_avail), seq_col, length(seqs_acgt), n_mode, mode_len %||% "NA")
			next
		}
		n_folds_ok <- n_folds_ok + 1L
		
		fold_dist_bounds <- list(); fold_dist_bounds[[locus]] <- compute_distance_bounds(train)[[locus]]
		db <- fold_dist_bounds[[locus]][[motif_name]]
		
		null_caches_env <- new.env(parent = emptyenv())
		score_one <- function(seq_up, sg) {
			eff <- get_effective_calibration(cal, sg)
			if (is.null(eff)) return(FALSE)
			zi <- restrict_to_distance_zone(seq_up, eff$pwm$length, db)
			if (is.null(zi)) return(FALSE)
			r <- max_pwm_score_over_sequence(zi$zone_seq, eff$pwm)
			if (is.na(r$max_score)) return(FALSE)
			cache_key <- paste(motif_name, eff$source_id, sep = "::")
			if (!exists(cache_key, envir = null_caches_env))
				assign(cache_key, make_null_distribution_cache(eff$pwm, bg_freqs, n_sim = n_sim, seed = seed + fold), envir = null_caches_env)
			pv <- empirical_pvalue(r$max_score, get(cache_key, envir = null_caches_env)(nchar(zi$zone_seq)))
			matched <- substr(zi$zone_seq, r$best_pos_start, r$best_pos_start + eff$pwm$length - 1L)
			!is.na(pv) && pv <= q_threshold && passes_hard_mismatch(matched, eff$empirical_consensus, motif_name)
		}
		
		pos_calls <- vapply(seq_len(nrow(test)), function(i)
			score_one(toupper(test$SEQ_UPSTREAM[i]), test$SUBGROUP[i]), logical(1))
		neg_calls <- vapply(seq_len(nrow(test)), function(i) {
			s <- toupper(test$SEQ_UPSTREAM[i])
			shuffled <- paste(sample(strsplit(s, "")[[1]]), collapse = "")
			score_one(shuffled, test$SUBGROUP[i])
		}, logical(1))
		
		tp_tot <- tp_tot + sum(pos_calls);  fn_tot <- fn_tot + sum(!pos_calls)
		tn_tot <- tn_tot + sum(!neg_calls); fp_tot <- fp_tot + sum(neg_calls)
	}
	
	if (n_folds_ok == 0) return(list(done = FALSE, reason = sprintf(
		"No fold could build a valid PWM (neither by subgroup nor at locus level). Diagnostic (last fold attempted): %s",
		diag_last %||% "n/a")))
	
	metrics <- .compute_confusion_metrics(tp_tot, fn_tot, tn_tot, fp_tot)
	list(done = TRUE, method = "PWM", locus = locus, motif = motif_name, k = k,
		n_folds_ok = n_folds_ok, n_test = n_total, metrics = metrics)
}

#  E_BOX / TATA_BOX : les deux méthodes utilisent  la même règle de composition ces motifs ne sont modélisés
# statistiquement dans aucune des deux méthodes. 
run_method_comparison <- function(ref_data, k = 5L, q_threshold = 0.05, n_sim = 2000L, seed = 42L) {
	dist_bounds <- compute_distance_bounds(ref_data)
	rows <- list(); reasons <- character(0)
	
	for (locus in SCANNER_LOCI) {
		motifs <- ALL_MOTIFS_BY_LOCUS[[locus]]
		for (motif_name in motifs) {
			cons <- run_consensus_evaluation(ref_data, locus, motif_name, dist_bounds = dist_bounds, seed = seed)
			if (isTRUE(cons$done)) {
				rows[[length(rows) + 1L]] <- cbind(
					data.frame(Method = "Consensus", Locus = locus, Motif = motif_name, N = cons$n_test),
					cons$metrics)
			} else {
				reasons <- c(reasons, sprintf("[Consensus] %s / %s : %s", locus, motif_name, cons$reason))
			}
			
			if (motif_name %in% PWM_MOTIFS_BY_LOCUS[[locus]]) {
				pwm_res <- run_pwm_evaluation(ref_data, locus, motif_name, k = k, q_threshold = q_threshold, n_sim = n_sim, seed = seed)
				if (isTRUE(pwm_res$done)) {
					rows[[length(rows) + 1L]] <- cbind(
						data.frame(Method = "PWM", Locus = locus, Motif = motif_name, N = pwm_res$n_test),
						pwm_res$metrics)
				} else {
					reasons <- c(reasons, sprintf("[PWM] %s / %s : %s", locus, motif_name, pwm_res$reason))
				}
			} else if (isTRUE(cons$done)) {
				# E_BOX / TATA_BOX : même règle, même résultat 
				rows[[length(rows) + 1L]] <- cbind(
					data.frame(Method = "PWM (shared rule)", Locus = locus, Motif = motif_name, N = cons$n_test),
					cons$metrics)
			}
		}
	}
	
	list(data = if (length(rows) > 0) dplyr::bind_rows(rows) else NULL, reasons = reasons)
}


# ══════════════════════════════════════════════════════════════════════════
# 16. CONSTRUCTION DU DATASET FINAL (par méthode)
# ══════════════════════════════════════════════════════════════════════════

build_base_info <- function(matched_df) {
	matched_df %>%
		dplyr::transmute(ALLELE = ALLELE_ID, SPECIES = SPECIES, GENE = GENE,
			LOCI = LOCI, FUNCTIONALITY = FUNCTIONALITY, LENGTH_NT = LEN_FULL) %>%
		dplyr::distinct(ALLELE, SPECIES, .keep_all = TRUE)
}

build_final_dataset <- function(base_info, hits_detected, reference_version) {
	df <- base_info
	if (is.null(hits_detected) || nrow(hits_detected) == 0) {
		df$REFERENCE_VERSION <- reference_version
		return(df)
	}
	
	lpart1_tbl <- hits_detected %>%
		dplyr::distinct(ALLELE_ID, SPECIES, GENE, LPART1_POS_START, LPART1_POS_END) %>%
		dplyr::transmute(ALLELE = ALLELE_ID, SPECIES = SPECIES, GENE = GENE,
			SEQ_LPART_1 = "present", POS_START_LPART_1 = LPART1_POS_START,
			POS_END_LPART_1 = LPART1_POS_END)
	df <- merge(df, lpart1_tbl, by = c("ALLELE", "SPECIES", "GENE"), all.x = TRUE)
	
	has_pval <- "P_VALUE" %in% names(hits_detected)
	has_qval <- "Q_VALUE" %in% names(hits_detected)
	
	for (m in unique(hits_detected$MOTIF)) {
		sub_df <- hits_detected %>% dplyr::filter(MOTIF == m)
		sub <- sub_df %>%
			dplyr::transmute(ALLELE = ALLELE_ID, SPECIES = SPECIES, GENE = GENE,
				!!paste0("SEQ_", m)       := matched_seq,
				!!paste0("POS_START_", m) := pos_start,
				!!paste0("POS_END_", m)   := pos_end,
				!!paste0("MISMATCH_", m)  := MISMATCH,
				!!paste0("DIST_", m)      := DIST_TO_LPART1)
		if (has_pval) sub[[paste0("PVAL_", m)]] <- sub_df$P_VALUE
		if (has_qval) sub[[paste0("QVAL_", m)]] <- sub_df$Q_VALUE
		df <- merge(df, sub, by = c("ALLELE", "SPECIES", "GENE"), all.x = TRUE)
	}
	
	for (m in unique(hits_detected$MOTIF)) {
		seq_col <- paste0("SEQ_", m); bin_col <- paste0("BIN_", m)
		if (seq_col %in% names(df)) {
			df[[bin_col]] <- ifelse(is.na(df[[seq_col]]) | df[[seq_col]] == "", 0L, 1L)
			df[[seq_col]] <- ifelse(is.na(df[[seq_col]]) | df[[seq_col]] == "", "absent", df[[seq_col]])
		}
	}
	df$REFERENCE_VERSION <- reference_version
	df
}

# ══════════════════════════════════════════════════════════════════════════
# 16b. FUSION DES DEUX MÉTHODES -> jeu de données pour les analyses
# ══════════════════════════════════════════════════════════════════════════

MOTIF_COL_PREFIXES <- c("SEQ", "POS_START", "POS_END", "MISMATCH", "DIST", "BIN", "PVAL", "QVAL")

.clean_func_vec <- function(x) {
	f <- toupper(trimws(gsub("[\\[\\]\\(\\)]", "", as.character(x), perl = TRUE)))
	ifelse(f %in% c("F", "P", "ORF"), f, NA_character_)
}

build_merged_dataset <- function(cons_df, pwm_df, choice, anchored_data, overrides = list()) {
	keys <- c("ALLELE", "SPECIES", "GENE")
	if (is.null(anchored_data) || nrow(anchored_data) == 0) stop("No anchored allele available.")
	base_df <- if (!is.null(cons_df)) cons_df else pwm_df
	
	lp <- anchored_data %>%
		dplyr::transmute(ALLELE = ALLELE_ID, SPECIES = SPECIES, GENE = GENE,
			SEQ_LPART_1 = SEQ_LPART1,
			POS_START_LPART_1 = TRUNC_OFFSET + nchar(SEQ_UPSTREAM) + 1L,
			POS_END_LPART_1   = POS_START_LPART_1 + nchar(SEQ_LPART1) - 1L,
			BIN_LPART_1 = 1L) %>%
		dplyr::distinct(ALLELE, SPECIES, GENE, .keep_all = TRUE)
	info <- base_df[, intersect(c(keys, "LOCI", "FUNCTIONALITY", "LENGTH_NT"), names(base_df)), drop = FALSE]
	out  <- merge(info, lp, by = keys)
	out$SUBGROUP <- if (exists("clean_subgroup", mode = "function")) clean_subgroup(out$GENE)
	else sub("^([A-Za-z]+[0-9]+).*$", "\\1", out$GENE)
	out$FUNCTIONALITY <- .clean_func_vec(out$FUNCTIONALITY)
	
	# source retenue pour (motif, sous-groupe) : override sinon choix par défaut du motif
	sg_src <- function(m, sg) {
		o <- overrides[[m]]
		if (!is.null(o) && sg %in% names(o)) o[[sg]] else choice[[m]]
	}
	active <- names(choice)[vapply(names(choice), function(m)
		choice[[m]] != "Exclude" || any(overrides[[m]] != "Exclude"), logical(1))]
	app_of <- unname(APP_MOTIF_COLUMN_MAP[active])
	empty  <- character(0)
	
	for (app_col in unique(app_of)) {
		subs <- list(); analysed <- rep(FALSE, nrow(out))
		for (m in active[app_of == app_col]) {
			loci_ok <- names(ALL_MOTIFS_BY_LOCUS)[vapply(ALL_MOTIFS_BY_LOCUS, function(v) m %in% v, logical(1))]
			src_row <- vapply(seq_len(nrow(out)), function(i)
				if (out$LOCI[i] %in% loci_ok) sg_src(m, out$SUBGROUP[i]) else "Exclude", character(1))
			analysed <- analysed | src_row != "Exclude"
			for (src in c("Consensus", "PWM")) {
				ids <- out[src_row == src, keys, drop = FALSE]
				if (nrow(ids) == 0) next
				df   <- if (src == "PWM") pwm_df else cons_df
				cols <- if (is.null(df)) character(0) else intersect(paste0(MOTIF_COL_PREFIXES, "_", m), names(df))
				if (length(cols) == 0) { empty <- c(empty, paste0(m, " (", src, ")")); next }
				s <- dplyr::semi_join(df[, c(keys, cols), drop = FALSE], ids, by = keys)
				names(s) <- c(keys, sub(paste0("_", m, "$"), paste0("_", app_col), cols))
				s[[paste0("METHOD_", app_col)]] <- src
				subs[[paste(m, src)]] <- s
			}
		}
		if (length(subs) > 0) out <- dplyr::left_join(out, dplyr::bind_rows(subs), by = keys)  # conserve l'ordre
		
		seq_c <- paste0("SEQ_", app_col); bin_c <- paste0("BIN_", app_col)
		if (!seq_c %in% names(out)) out[[seq_c]] <- NA_character_
		out[[seq_c]] <- ifelse(is.na(out[[seq_c]]) | out[[seq_c]] == "", "absent", out[[seq_c]])
		out[[bin_c]] <- ifelse(out[[seq_c]] == "absent", 0L, 1L)
		out[[seq_c]][!analysed] <- NA_character_   # non analysé (exclu) ≠ absent
		out[[bin_c]][!analysed] <- NA_integer_
	}
	
	req_cols <- c("ALLELE", "SPECIES", "GENE", "LOCI", "SUBGROUP", "FUNCTIONALITY", "LENGTH_NT")
	list(data = out, missing_cols = setdiff(req_cols, names(out)), empty_motifs = empty,
		n_not_anchored = nrow(info) - nrow(out), n_func_na = sum(is.na(out$FUNCTIONALITY)))
}
# ══════════════════════════════════════════════════════════════════════════
# 17. VUE ANCRÉE SURLIGNÉE  
# ══════════════════════════════════════════════════════════════════════════

.scanner_seg_span <- function(txt, color) {
	if (nchar(txt) == 0) return("")
	if (is.na(color)) sprintf('<span>%s</span>', htmltools::htmlEscape(txt))
	else sprintf('<span style="background:%s;color:#ffffff;">%s</span>', color, htmltools::htmlEscape(txt))
}

render_anchored_alignment <- function(anchored_data, hits_detected, subgroup_sel, max_downstream = 40L) {
	if (is.null(anchored_data) || nrow(anchored_data) == 0)
		return(shiny::div(class = "alert alert-warning", shiny::icon("exclamation-triangle"),
			" No anchored sequence available."))
	
	df_anchor <- anchored_data[anchored_data$SUBGROUP == subgroup_sel, , drop = FALSE]
	if (nrow(df_anchor) == 0)
		return(shiny::div(class = "alert alert-warning", " No sequence for this subgroup."))
	
	df_anchor <- df_anchor[order(df_anchor$ALLELE_ID), , drop = FALSE]
	df_anchor$SEQ_DOWNSTREAM <- substr(df_anchor$SEQ_DOWNSTREAM, 1L, max_downstream)
	max_up <- max(nchar(df_anchor$SEQ_UPSTREAM))
	
	df_hits_sg <- NULL
	if (!is.null(hits_detected) && nrow(hits_detected) > 0)
		df_hits_sg <- hits_detected[which(hits_detected$SUBGROUP == subgroup_sel), , drop = FALSE]
	
	rows_html <- lapply(seq_len(nrow(df_anchor)), function(i) {
		r   <- df_anchor[i, ]
		pad <- max_up - nchar(r$SEQ_UPSTREAM)
		full_seq <- paste0(strrep("-", pad), r$SEQ_UPSTREAM, r$SEQ_LPART1, r$SEQ_DOWNSTREAM)
		n_chars  <- nchar(full_seq)
		colors   <- rep(NA_character_, n_chars)
		
		lp_start <- pad + nchar(r$SEQ_UPSTREAM) + 1L
		lp_end   <- lp_start + nchar(r$SEQ_LPART1) - 1L
		if (lp_end >= lp_start) colors[lp_start:lp_end] <- unname(SCANNER_SCANNER_MOTIF_COLORS[["LPART_1"]])
		
		if (!is.null(df_hits_sg) && nrow(df_hits_sg) > 0) {
			hits_allele <- df_hits_sg[which(df_hits_sg$ALLELE_ID == r$ALLELE_ID &
					df_hits_sg$SPECIES == r$SPECIES &
					!is.na(df_hits_sg$pos_start) & !is.na(df_hits_sg$pos_end)), , drop = FALSE]
			for (j in seq_len(nrow(hits_allele))) {
				hs  <- pad + hits_allele$pos_start[j] - r$TRUNC_OFFSET
				he  <- pad + hits_allele$pos_end[j]   - r$TRUNC_OFFSET
				col <- SCANNER_SCANNER_MOTIF_COLORS[hits_allele$MOTIF[j]]
				if (isTRUE(hs >= 1L && he <= n_chars && he >= hs) && !is.na(col))
					colors[hs:he] <- unname(col)
			}
		}
		
		chars <- strsplit(full_seq, "")[[1]]
		spans <- character(0)
		cur_color <- colors[1]; cur_txt <- chars[1]
		if (length(chars) > 1) {
			for (k in 2:length(chars)) {
				same <- (is.na(colors[k]) && is.na(cur_color)) ||
					(!is.na(colors[k]) && !is.na(cur_color) && colors[k] == cur_color)
				if (same) cur_txt <- paste0(cur_txt, chars[k])
				else { spans <- c(spans, .scanner_seg_span(cur_txt, cur_color)); cur_color <- colors[k]; cur_txt <- chars[k] }
			}
		}
		spans <- c(spans, .scanner_seg_span(cur_txt, cur_color))
		
		label <- sprintf("%s | %s", r$ALLELE_ID, r$SPECIES)
		htmltools::tags$div(
			style = "white-space:nowrap; font-family:'Courier New',monospace; font-size:12px; line-height:19px; letter-spacing:0.5px;",
			htmltools::tags$span(
				style = "display:inline-block; width:280px; font-weight:bold; color:#374151; overflow:hidden; text-overflow:ellipsis; vertical-align:top;",
				htmltools::htmlEscape(label)),
			htmltools::HTML(paste(spans, collapse = ""))
		)
	})
	
	legend <- htmltools::tags$div(
		style = "margin-bottom:10px; padding:8px 10px; background:#F9FAFB; border-radius:6px;",
		htmltools::tags$strong("Legend: ", style = "margin-right:8px; color:#374151;"),
		lapply(names(SCANNER_SCANNER_MOTIF_COLORS), function(m) {
			htmltools::tags$span(
				style = sprintf(
					"background:%s;color:#ffffff;padding:2px 8px;border-radius:4px;margin-right:6px;font-size:11px;font-weight:bold;display:inline-block;margin-bottom:4px;",
					unname(SCANNER_SCANNER_MOTIF_COLORS[[m]])), m)
		})
	)
	
	htmltools::tagList(
		legend,
		htmltools::tags$div(
			style = "overflow-x:auto; background:#FAFAFA; padding:12px; border-radius:8px; border:1px solid #E5E7EB; max-height:650px; overflow-y:auto;",
			rows_html
		)
	)
}


# ══════════════════════════════════════════════════════════════════════════
# 18. UI
# ══════════════════════════════════════════════════════════════════════════

.html2canvas_dependency <- function() {
	htmltools::singleton(htmltools::tags$head(
		htmltools::tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js")
	))
}

.png_export_button <- function(capture_id, filename, label = "📷 Export anchored view (PNG)") {
	onclick_js <- sprintf(
		"var el = document.getElementById('%s'); if (!el) { return false; }
		 var btn = event.target; var oldTxt = btn.innerHTML; btn.innerHTML = '⏳ Generating PNG...'; btn.disabled = true;
		 html2canvas(el, {backgroundColor: '#ffffff', scale: 2, useCORS: true}).then(function(canvas) {
		   var link = document.createElement('a');
		   link.download = '%s_' + new Date().toISOString().slice(0,10) + '.png';
		   link.href = canvas.toDataURL('image/png');
		   link.click();
		   btn.innerHTML = oldTxt; btn.disabled = false;
		 }).catch(function(err) {
		   btn.innerHTML = oldTxt; btn.disabled = false;
		   alert('PNG export failed: ' + err);
		 });
		 return false;",
		capture_id, filename)
	shiny::tags$button(type = "button", class = "btn btn-info", style = "margin-top:6px;",
		onclick = onclick_js, label)
}

.mod_method_results_ui <- function(ns, prefix, show_candidates_tab = TRUE) {
	capture_id <- ns(paste0(prefix, "_align_capture"))
	shiny::tabsetPanel(
		shiny::tabPanel(shiny::tagList(shiny::icon("dna"), " Anchored view"),
			.html2canvas_dependency(),
			shiny::br(),
			shiny::p(shiny::em("Not a multiple sequence alignment — sequences are simply anchored on the start of L-PART1."), style = "font-size:11px;color:#888;"),
			shiny::fluidRow(
				shiny::column(8, shiny::selectInput(ns(paste0(prefix, "_align_subgroup")), "Subgroup", choices = NULL)),
				shiny::column(4, style = "text-align:right; padding-top:24px;",
					.png_export_button(capture_id, paste0("anchored_view_", prefix)))
			),
			shiny::div(id = capture_id, style = "background:#ffffff; padding:6px;",
				shiny::uiOutput(ns(paste0(prefix, "_alignment_view")))
			),
			shiny::br()
		),
		shiny::tabPanel(shiny::tagList(shiny::icon("star"), " Detected motifs"),
			shiny::br(),
			DT::DTOutput(ns(paste0(prefix, "_detected_table"))), shiny::br(),
			shiny::downloadButton(ns(paste0(prefix, "_dl_detected_csv")), "📥 CSV — Detected motifs", class = "btn-warning")
		),
		if (show_candidates_tab) shiny::tabPanel(shiny::tagList(shiny::icon("list"), " All candidates"),
			shiny::br(),
			DT::DTOutput(ns(paste0(prefix, "_all_table"))), shiny::br()
		),
		shiny::tabPanel(shiny::tagList(shiny::icon("table"), " Final dataset"),
			shiny::br(),
			DT::DTOutput(ns(paste0(prefix, "_final_table"))), shiny::br(),
			shiny::downloadButton(ns(paste0(prefix, "_dl_final_csv")),  "📥 CSV — Final dataset",  class = "btn-warning"),
			shiny::downloadButton(ns(paste0(prefix, "_dl_final_xlsx")), "📊 Excel — Final dataset", class = "btn-primary")
		)
	)
}


mod_scanner_ui_consensus <- function(id) {
	ns <- shiny::NS(id)
	
	shiny::tagList(
		shiny::uiOutput(ns("reference_status")),
		
		shiny::fluidRow(
			shinydashboard::box(
				width = 12, title = "① Import sequences to annotate", status = "primary", solidHeader = TRUE,
				shiny::fluidRow(
					shiny::column(6,
						shiny::div(style = "background:#FFF9C4; padding:14px; border-radius:8px; height:100%;",
							shiny::h5(shiny::icon("dna"), shiny::strong("Full 5'UTR (contains L-PART1)")),
							shiny::p(style = "font-size:11px;color:#555;",
								"Expected header: ", shiny::code(">ACC|ALLELE|Species|Func|...")),
							shiny::fileInput(ns("f_utr5"), "📂 FASTA file", accept = c(".fasta", ".fa", ".txt", ".fna")),
							shiny::textAreaInput(ns("t_utr5"), "— or paste here —", rows = 5, width = "100%"),
							shiny::uiOutput(ns("status_utr5"))
						)
					),
					shiny::column(6,
						shiny::div(style = "background:#E8F5E9; padding:14px; border-radius:8px; height:100%;",
							shiny::h5(shiny::icon("bullseye"), shiny::strong("L-PART1 only")),
							shiny::fileInput(ns("f_lpart1"), "📂 FASTA file", accept = c(".fasta", ".fa", ".txt", ".fna")),
							shiny::textAreaInput(ns("t_lpart1"), "— or paste here —", rows = 5, width = "100%"),
							shiny::uiOutput(ns("status_lpart1"))
						)
					)
				)
			)
		),
		
		shiny::fluidRow(
			shinydashboard::box(
				width = 12, title = "② Detection parameters", status = "primary", solidHeader = TRUE,
				shiny::fluidRow(
					shiny::column(4, shiny::sliderInput(ns("q_threshold"),
						"q-value threshold (BH, PWM method)", min = 0.01, max = 0.10, value = 0.05, step = 0.01)),
					shiny::column(4, shiny::numericInput(ns("n_sim"),
						"Monte Carlo simulations (null distribution, PWM)", value = 5000, min = 1000, max = 20000, step = 1000)),
					shiny::column(4, shiny::checkboxInput(ns("run_comparison"),
						"Also compute method comparison (k=5 for PWM) — slower", value = FALSE))
				),
				shiny::br(),
				shiny::actionButton(ns("btn_scan"), "🔍 Run scan (Consensus + PWM)", icon = shiny::icon("search"),
					class = "btn-success btn-lg", style = "width:100%; height:55px; font-size:14px;")
			)
		),
		
		shiny::uiOutput(ns("scan_status")),
		
		shiny::fluidRow(
			shinydashboard::box(
				width = 12, title = "🔎 Allele diagnostic (where alleles are lost)",
				status = "primary", solidHeader = TRUE,
				collapsible = TRUE, collapsed = TRUE,
				DT::DTOutput(ns("diag_summary")),
				shiny::h5("Alleles not anchored"),
				DT::DTOutput(ns("diag_lost"))
			)
		),
		shiny::hr(),
		
		shiny::h4(shiny::icon("ruler"), " Results — Consensus"),
		.mod_method_results_ui(ns, "cons", show_candidates_tab = FALSE)
	)
}

mod_scanner_ui_pwm <- function(id) {
	ns <- shiny::NS(id)
	shiny::tagList(
		shiny::uiOutput(ns("scan_status")),
		shiny::h4(shiny::icon("chart-line"), " Results — PWM"),
		.mod_method_results_ui(ns, "pwm", show_candidates_tab = TRUE)
	)
}


mod_scanner_ui_comparison <- function(id) {
	ns <- shiny::NS(id)
	shiny::tagList(
		shiny::h4(shiny::icon("scale-balanced"), " Method comparison"),
		shiny::uiOutput(ns("comparison_status")),
		
		shinydashboard::box(
			width = 12, title = "📋 Comparison table (full data, all loci & motifs)",
			status = "primary", solidHeader = TRUE,
			DT::DTOutput(ns("comparison_table")), shiny::br(),
			shiny::downloadButton(ns("dl_comparison_csv"), "📥 CSV — Full comparison data", class = "btn-warning")
		),
		
		shinydashboard::box(
			width = 12, title = "📊 Chart", status = "primary", solidHeader = TRUE,
			shiny::fluidRow(
				shiny::column(4, shiny::selectInput(ns("cmp_view"), "View",
					choices = c("Global (all motifs & loci)" = "global",
						"By motif" = "motif", "By locus" = "locus"))),
				shiny::column(4, shiny::uiOutput(ns("cmp_filter_secondary"))),
				shiny::column(4, style = "padding-top:24px; display:flex; gap:8px;",
					shiny::tags$button(type = "button", class = "btn btn-info", style = "width:100%;",
						onclick = sprintf(
							"var gd = document.getElementById('%s'); if (!gd) return false;
							 Plotly.downloadImage(gd, {format:'png', filename:'method_comparison_' + new Date().toISOString().slice(0,10), height:600, width:1000, scale:2});
							 return false;", ns("comparison_plot")),
						"📷 PNG")
				)
			),
			shiny::br(),
			plotly::plotlyOutput(ns("comparison_plot")), shiny::br(),
			shiny::downloadButton(ns("dl_comparison_view_csv"), "📥 CSV — Data shown in this chart (current view)", class = "btn-warning")
		)
	)
}


mod_scanner_ui_final <- function(id) {
	ns <- shiny::NS(id)
	shiny::tagList(
		shiny::h4(shiny::icon("sliders"), " Final dataset — choose the method per motif"),
		shiny::p(style = "font-size:12px;color:#666;",
			"Run the scan first (Consensus tab). For each motif, pick the method that detected it best. ",
			"TATA_BOX and E_BOX use the same rule in both methods."),
		shiny::actionButton(ns("btn_compute_f1"), "📐 Compute F1 now (reference set)", class = "btn-info"),
		shiny::br(), shiny::br(),
		shiny::uiOutput(ns("f1_status")),
		shiny::uiOutput(ns("choice_table")),
		shiny::br(),
		shiny::actionButton(ns("btn_run_analyses"), "🚀 Build dataset & run analyses",
			class = "btn-success btn-lg", style = "width:100%; height:55px;"),
		shiny::br(), shiny::br(),
		shiny::uiOutput(ns("schema_status")),
		DT::DTOutput(ns("merged_table")), shiny::br(),
		shiny::downloadButton(ns("dl_merged_csv"),  "📥 CSV — Dataset for analyses", class = "btn-warning"),
		shiny::downloadButton(ns("dl_merged_xlsx"), "📊 Excel — Dataset for analyses", class = "btn-primary")
	)
}

# ══════════════════════════════════════════════════════════════════════════
# 19. SERVER
# ══════════════════════════════════════════════════════════════════════════

mod_scanner_server <- function(id, df_final = shiny::reactive(NULL)) {
	shiny::moduleServer(id, function(input, output, session) {
		
		store <- shiny::reactiveValues( diag = NULL,
			done = FALSE, ref_info = NULL, calib_warnings = NULL,
			cons_join_report = NULL, cons_hits_detected = NULL, cons_base_info = NULL, cons_final_df = NULL, cons_anchored_data = NULL,
			pwm_join_report  = NULL, pwm_hits_all = NULL, pwm_hits_detected = NULL, pwm_base_info = NULL, pwm_final_df = NULL, pwm_anchored_data = NULL,
			comparison_df = NULL, comparison_reasons = NULL
		)
		
		store$ref_info <- load_bundled_reference()
		
		output$reference_status <- shiny::renderUI({
			ri <- store$ref_info
			if (!is.null(ri$error))
				return(shiny::div(class = "alert alert-danger", style = "font-size:12px;",
					shiny::icon("exclamation-triangle"), " ", ri$error))
			upstream_note <- if (!ri$has_upstream_seq)
				" (method comparison unavailable without SEQ_UPSTREAM — annotation works normally without it)" else ""
			shiny::div(class = "alert alert-success", style = "font-size:12px;",
				shiny::icon("database"),
				sprintf(" PWM calibration: curated IMGT reference dataset, version %s (%d instances, %d loci).%s",
					ri$version, nrow(ri$data), length(unique(ri$data$LOCI)), upstream_note))
		})
		
		shiny::observeEvent(input$f_utr5, {
			shiny::req(input$f_utr5)
			txt <- paste(readLines(input$f_utr5$datapath, warn = FALSE), collapse = "\n")
			shiny::updateTextAreaInput(session, "t_utr5", value = txt)
		})
		shiny::observeEvent(input$f_lpart1, {
			shiny::req(input$f_lpart1)
			txt <- paste(readLines(input$f_lpart1$datapath, warn = FALSE), collapse = "\n")
			shiny::updateTextAreaInput(session, "t_lpart1", value = txt)
		})
		
		mk_status <- function(txt, label) {
			if (is.null(txt) || nchar(trimws(txt)) == 0) return(NULL)
			n <- length(grep("^>", strsplit(txt, "\n")[[1]]))
			shiny::div(class = "alert alert-success", style = "margin:4px 0 0 0; padding:4px 8px; font-size:11px;",
				shiny::icon("check-circle"), paste0(" ", n, " ", label, " loaded"))
		}
		output$status_utr5   <- shiny::renderUI(mk_status(input$t_utr5,   "5'UTR sequence(s)"))
		output$status_lpart1 <- shiny::renderUI(mk_status(input$t_lpart1, "L-PART1 sequence(s)"))
		
		shiny::observeEvent(input$btn_scan, {
			txt_utr5 <- input$t_utr5; txt_lp <- input$t_lpart1
			if (is.null(txt_utr5) || nchar(trimws(txt_utr5)) == 0 || is.null(txt_lp) || nchar(trimws(txt_lp)) == 0) {
				shiny::showNotification("❌ Both FASTA files (①) are required.", type = "error", duration = 6)
				return()
			}
			ri <- store$ref_info
			if (is.null(ri$data)) {
				shiny::showNotification("❌ Reference dataset unavailable — cannot calibrate the PWM method.", type = "error", duration = 8)
				return()
			}
			
			store$done <- FALSE
			
			shiny::withProgress(message = "Scan in progress...", value = 0, {
				tryCatch({
					shiny::incProgress(0.05, detail = "Parsing FASTA files...")
					df_utr5   <- parse_fasta(txt_utr5)
					df_lpart1 <- parse_fasta(txt_lp)
					if (is.null(df_utr5)   || nrow(df_utr5)   == 0) stop("Invalid or empty 5'UTR file, or out of scope (IG only).")
					if (is.null(df_lpart1) || nrow(df_lpart1) == 0) stop("Invalid or empty L-PART1 file, or out of scope (IG only).")
					store$diag <- tryCatch(anchor_alleles(df_utr5, df_lpart1)$diagnostic, error = function(e) NULL)
				
					shiny::incProgress(0.1, detail = "Computing plausible distance bounds...")
					dist_bounds <- compute_distance_bounds(ri$data)
					dist_warnings <- collect_distance_bound_warnings(dist_bounds)
					
					shiny::incProgress(0.25, detail = "Method ① Consensus (baseline)...")
					res_cons <- run_consensus_scan(df_utr5, df_lpart1, dist_bounds = dist_bounds)
					base_info_cons <- build_base_info(res_cons$join_report$matched)
					final_df_cons  <- build_final_dataset(base_info_cons, res_cons$hits_detected, store$ref_info$version)
					store$cons_join_report   <- res_cons$join_report
					store$cons_hits_detected <- res_cons$hits_detected
					store$cons_base_info     <- base_info_cons
					store$cons_final_df      <- final_df_cons
					store$cons_anchored_data <- res_cons$anchored_data
					
					shiny::incProgress(0.4, detail = "PWM calibration (by subgroup, locus fallback)...")
					calib_by_locus <- list()
					for (loc in SCANNER_LOCI) {
						rows_loc <- df_utr5[df_utr5$LOCI == loc, , drop = FALSE]
						bg <- compute_background_freqs(
							if (nrow(rows_loc) > 0) rows_loc$SEQ_FULL else ri$data$SEQ_UPSTREAM %||% character(0))
						calib_by_locus[[loc]] <- build_locus_calibration(ri$data, loc, bg)
					}
					all_warnings <- collect_calibration_warnings(calib_by_locus)
					store$calib_warnings <- unique(c(
						all_warnings[!vapply(all_warnings, is.null, logical(1))],
						dist_warnings))
					
					shiny::incProgress(0.6, detail = "Method ② PWM — scan + p-value + BH...")
					res_pwm <- run_motif_scan_pwm(df_utr5, df_lpart1, calib_by_locus,
						q_threshold = input$q_threshold, n_sim = input$n_sim, dist_bounds = dist_bounds)
					base_info_pwm <- build_base_info(res_pwm$join_report$matched)
					final_df_pwm  <- build_final_dataset(base_info_pwm, res_pwm$hits_detected, store$ref_info$version)
					store$pwm_join_report   <- res_pwm$join_report
					store$pwm_hits_all      <- res_pwm$hits_all
					store$pwm_hits_detected <- res_pwm$hits_detected
					store$pwm_base_info     <- base_info_pwm
					store$pwm_final_df      <- final_df_pwm
					store$pwm_anchored_data <- res_pwm$anchored_data
					
					store$done <- TRUE
					
					if (isTRUE(input$run_comparison)) {
						shiny::incProgress(0.85, detail = "③ Method comparison (k=5 CV for PWM)...")
						cmp <- run_method_comparison(ri$data, k = 5L, q_threshold = input$q_threshold, n_sim = 2000L)
						store$comparison_df      <- cmp$data
						store$comparison_reasons <- cmp$reasons
					} else {
						store$comparison_df <- NULL
						store$comparison_reasons <- NULL
					}
					
					shiny::incProgress(1, detail = "Done!")
					
					if (length(store$calib_warnings) > 0) {
						shiny::showNotification(
							shiny::HTML(paste0("<strong>⚠ Incomplete calibration:</strong><br>", paste(store$calib_warnings, collapse = "<br>"))),
							type = "warning", duration = 20)
					}
					
					sg_cons <- if (!is.null(res_cons$anchored_data) && nrow(res_cons$anchored_data) > 0) sort(unique(res_cons$anchored_data$SUBGROUP)) else character(0)
					sg_pwm  <- if (!is.null(res_pwm$anchored_data)  && nrow(res_pwm$anchored_data)  > 0) sort(unique(res_pwm$anchored_data$SUBGROUP))  else character(0)
					shiny::updateSelectInput(session, "cons_align_subgroup", choices = sg_cons)
					shiny::updateSelectInput(session, "pwm_align_subgroup",  choices = sg_pwm)
					
				}, error = function(e) {
					shiny::showNotification(paste0("❌ Error: ", e$message), type = "error", duration = 12)
				})
			})
		})
		
		output$scan_status <- shiny::renderUI({
			if (!store$done) {
				shiny::div(class = "alert alert-info", shiny::icon("info-circle"),
					" Import both FASTA files then click ", shiny::strong("🔍 Run scan"), ".")
			} else {
				n_cons <- if (!is.null(store$cons_hits_detected)) length(unique(store$cons_hits_detected$ALLELE_ID)) else 0L
				n_pwm  <- if (!is.null(store$pwm_hits_detected))  length(unique(store$pwm_hits_detected$ALLELE_ID))  else 0L
				shiny::div(class = "alert alert-success", shiny::icon("check-circle"),
					shiny::strong(sprintf(" ✅ Scan complete — Consensus: %d allele(s) with ≥1 motif; PWM: %d allele(s) with ≥1 motif.", n_cons, n_pwm)))
			}
		})
		
	
		.register_method_outputs <- function(prefix, hits_detected_r, hits_all_r, anchored_data_r, final_df_r, dl_prefix) {
			
			output[[paste0(prefix, "_alignment_view")]] <- shiny::renderUI({
				shiny::req(store$done, input[[paste0(prefix, "_align_subgroup")]])
				shiny::validate(shiny::need(nchar(input[[paste0(prefix, "_align_subgroup")]]) > 0, "Select a subgroup."))
				render_anchored_alignment(anchored_data_r(), hits_detected_r(), input[[paste0(prefix, "_align_subgroup")]])
			})
			
			output[[paste0(prefix, "_detected_table")]] <- DT::renderDT({
				shiny::req(store$done)
				df <- hits_detected_r(); shiny::req(df)
				DT::datatable(df, options = list(pageLength = 20, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")),
					extensions = "Buttons", rownames = FALSE, class = "cell-border stripe hover compact")
			})
			
			if (!is.null(hits_all_r)) {
				output[[paste0(prefix, "_all_table")]] <- DT::renderDT({
					shiny::req(store$done)
					df <- hits_all_r(); shiny::req(df)
					DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")),
						extensions = "Buttons", rownames = FALSE, class = "cell-border stripe hover compact")
				})
			}
			
			output[[paste0(prefix, "_final_table")]] <- DT::renderDT({
				shiny::req(store$done)
				df <- final_df_r(); shiny::req(df)
				DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv", "excel")),
					extensions = "Buttons", rownames = FALSE, class = "cell-border stripe hover compact")
			})
			
			output[[paste0(prefix, "_dl_detected_csv")]] <- shiny::downloadHandler(
				filename = function() paste0("DetectedMotifs_", dl_prefix, "_", Sys.Date(), ".csv"),
				content  = function(file) { df <- hits_detected_r(); shiny::req(df); write.csv(df, file, row.names = FALSE) }
			)
			output[[paste0(prefix, "_dl_final_csv")]] <- shiny::downloadHandler(
				filename = function() paste0("FinalDataset_", dl_prefix, "_", Sys.Date(), ".csv"),
				content  = function(file) { df <- final_df_r(); shiny::req(df); write.csv(df, file, row.names = FALSE) }
			)
			output[[paste0(prefix, "_dl_final_xlsx")]] <- shiny::downloadHandler(
				filename = function() paste0("FinalDataset_", dl_prefix, "_", Sys.Date(), ".xlsx"),
				content  = function(file) { df <- final_df_r(); shiny::req(df); openxlsx::write.xlsx(df, file, asTable = FALSE) }
			)
		}
		
		.register_method_outputs("cons",
			hits_detected_r = shiny::reactive(store$cons_hits_detected), hits_all_r = NULL,
			anchored_data_r = shiny::reactive(store$cons_anchored_data), final_df_r = shiny::reactive(store$cons_final_df),
			dl_prefix = "Consensus")
		
		.register_method_outputs("pwm",
			hits_detected_r = shiny::reactive(store$pwm_hits_detected), hits_all_r = shiny::reactive(store$pwm_hits_all),
			anchored_data_r = shiny::reactive(store$pwm_anchored_data), final_df_r = shiny::reactive(store$pwm_final_df),
			dl_prefix = "PWM")
		

		output$comparison_status <- shiny::renderUI({
			if (is.null(store$comparison_df) && is.null(store$comparison_reasons))
				return(shiny::div(class = "alert alert-info", shiny::icon("info-circle"),
					" Check \u00ab Also compute method comparison \u00bb before launching the scan."))
			
			reasons_block <- if (length(store$comparison_reasons) > 0)
				shiny::div(class = "alert alert-warning", style = "font-size:11px; margin-top:6px;",
					shiny::icon("exclamation-triangle"),
					shiny::strong(sprintf(" %d (locus, motif) combination(s) not computed: ", length(store$comparison_reasons))),
					shiny::tags$ul(style = "margin-bottom:0;",
						lapply(store$comparison_reasons, function(r) shiny::tags$li(r)))
				) else NULL
			
			if (is.null(store$comparison_df)) {
				return(shiny::tagList(
					shiny::div(class = "alert alert-danger", shiny::icon("times-circle"),
						" No comparison could be computed — see details below."),
					reasons_block
				))
			}
			
			n_rows <- nrow(store$comparison_df)
			shiny::tagList(
				shiny::div(class = "alert alert-success", shiny::icon("check-circle"),
					sprintf(" Comparison computed: %d row(s) (Method \u00d7 Locus \u00d7 Motif).", n_rows)),
				reasons_block
			)
		})
		
		output$comparison_table <- DT::renderDT({
			shiny::req(store$comparison_df)
			df <- store$comparison_df
			df$sensitivity <- round(df$sensitivity, 3); df$specificity <- round(df$specificity, 3)
			df$precision   <- round(df$precision, 3);   df$f1_score    <- round(df$f1_score, 3)
			DT::datatable(df, options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")),
				extensions = "Buttons", rownames = FALSE, class = "cell-border stripe hover compact")
		})
		
		output$cmp_filter_secondary <- shiny::renderUI({
			shiny::req(store$comparison_df)
			df <- store$comparison_df
			if (identical(input$cmp_view, "motif")) {
				shiny::selectInput(session$ns("cmp_motif_sel"), "Motif", choices = sort(unique(df$Motif)))
			} else if (identical(input$cmp_view, "locus")) {
				shiny::selectInput(session$ns("cmp_locus_sel"), "Locus", choices = sort(unique(df$Locus)))
			} else {
				shiny::div(style = "padding-top:24px; font-size:11px; color:#888;",
					shiny::icon("info-circle"), " All (locus, motif) combinations are aggregated by method.")
			}
		})
		
	
		.aggregate_comparison_global <- function(df) {
			agg <- df %>%
				dplyr::group_by(Method) %>%
				dplyr::summarise(N = sum(N), tp = sum(tp), fn = sum(fn), tn = sum(tn), fp = sum(fp), .groups = "drop")
			purrr::pmap_dfr(agg, function(Method, N, tp, fn, tn, fp) {
				cbind(data.frame(Method = Method, Locus = "All loci", Motif = "All motifs", N = N,
					stringsAsFactors = FALSE), .compute_confusion_metrics(tp, fn, tn, fp))
			})
		}
		
		
		comparison_view_df <- shiny::reactive({
			shiny::req(store$comparison_df)
			df <- store$comparison_df
			view <- input$cmp_view %||% "global"
			if (view == "global") {
				.aggregate_comparison_global(df)
			} else if (view == "motif") {
				shiny::req(input$cmp_motif_sel)
				df[df$Motif == input$cmp_motif_sel, , drop = FALSE]
			} else if (view == "locus") {
				shiny::req(input$cmp_locus_sel)
				df[df$Locus == input$cmp_locus_sel, , drop = FALSE]
			} else df
		})
		
		output$comparison_plot <- plotly::renderPlotly({
			df <- comparison_view_df()
			shiny::req(nrow(df) > 0)
			view <- input$cmp_view %||% "global"
			
			
			df$label <- switch(view, global = df$Method, motif = df$Locus, locus = df$Motif, paste(df$Locus, df$Motif, sep = " / "))
			plot_title <- switch(view,
				global = "Global method comparison — all motifs and loci combined",
				motif  = sprintf("Method comparison — motif %s (by locus)", unique(df$Motif)[1] %||% ""),
				locus  = sprintf("Method comparison — locus %s (by motif)", unique(df$Locus)[1] %||% ""),
				"Method comparison vs confirmed annotations")
			
			dfl <- df %>%
				dplyr::select(label, Method, sensitivity, specificity, precision, f1_score) %>%
				tidyr::pivot_longer(cols = c(sensitivity, specificity, precision, f1_score), names_to = "Metric", values_to = "Value")
			
			p <- ggplot2::ggplot(dfl, ggplot2::aes(x = label, y = Value, fill = Method)) +
				ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.7) +
				ggplot2::facet_wrap(~Metric, ncol = 2) +
				ggplot2::labs(title = plot_title, x = NULL, y = NULL) +
				ggplot2::theme_minimal(base_size = 10) +
				ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))
			
			plotly::ggplotly(p, height = 550) %>% plotly::config(displaylogo = FALSE)
		})
		
		output$dl_comparison_csv <- shiny::downloadHandler(
			filename = function() paste0("MethodComparison_Full_", Sys.Date(), ".csv"),
			content  = function(file) { shiny::req(store$comparison_df); write.csv(store$comparison_df, file, row.names = FALSE) }
		)
		
		
		output$dl_comparison_view_csv <- shiny::downloadHandler(
			filename = function() {
				view <- input$cmp_view %||% "global"
				paste0("MethodComparison_", tools::toTitleCase(view), "View_", Sys.Date(), ".csv")
			},
			content = function(file) { df <- comparison_view_df(); shiny::req(df); write.csv(df, file, row.names = FALSE) }
		)
		.diag_sg <- function(d) ifelse(is.na(d$GENE), NA_character_,
			if (exists("clean_subgroup", mode = "function")) clean_subgroup(d$GENE)
			else sub("^([A-Za-z]+[0-9]+).*$", "\\1", d$GENE))
		
		output$diag_summary <- DT::renderDT({
			shiny::req(store$diag)
			d <- store$diag; d$SUBGROUP <- .diag_sg(d)
			tab <- d %>% dplyr::count(SPECIES, LOCI, SUBGROUP, STATUS) %>%
				tidyr::pivot_wider(names_from = STATUS, values_from = n, values_fill = 0)
			DT::datatable(tab, rownames = FALSE, filter = "top",
				options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")),
				extensions = "Buttons", class = "cell-border stripe hover compact")
		})
		
		output$diag_lost <- DT::renderDT({
			shiny::req(store$diag)
			d <- store$diag; d$SUBGROUP <- .diag_sg(d)
			d <- d[!grepl("^OK", d$STATUS), , drop = FALSE]
			DT::datatable(d, rownames = FALSE, filter = "top",
				options = list(pageLength = 15, scrollX = TRUE, dom = "Bfrtip", buttons = c("copy", "csv")),
				extensions = "Buttons", class = "cell-border stripe hover compact")
		})
		
		
		# ── Jeu de données final
		# serveur du module
		shiny::observeEvent(input$btn_compute_f1, {
			shiny::req(store$ref_info$data)
			shiny::withProgress(message = "Computing F1 on the reference set…", value = 0.3, {
				cmp <- tryCatch(
					run_method_comparison(store$ref_info$data, k = 5L, q_threshold = 0.05, n_sim = 2000L),
					error = function(e) list(data = NULL, reasons = paste("Error:", e$message)))
				store$comparison_df      <- cmp$data
				store$comparison_reasons <- cmp$reasons
			})
			if (is.null(store$comparison_df)) {
				shiny::showNotification(
					shiny::HTML(paste0("<strong>❌ F1 could not be computed:</strong><br>",
						paste(utils::head(store$comparison_reasons, 4), collapse = "<br>"))),
					type = "error", duration = 20)
			} else {
				shiny::showNotification("✅ F1 computed.", type = "message")
			}
		})
		
		detected_motifs <- shiny::reactive({
			shiny::req(store$done)
			sort(unique(c(store$cons_hits_detected$MOTIF, store$pwm_hits_detected$MOTIF)))
		})
		
		.f1_by_method <- function(motif, method) {
			cmp <- store$comparison_df
			if (is.null(cmp)) return(NA_real_)
			d <- cmp[cmp$Motif == motif & grepl(paste0("^", method), cmp$Method), , drop = FALSE]
			if (nrow(d) == 0) return(NA_real_)
			tp <- sum(d$tp); fp <- sum(d$fp); fn <- sum(d$fn)
			if ((2 * tp + fp + fn) == 0) NA_real_ else 2 * tp / (2 * tp + fp + fn)
		}
		
		.sgs_for_motif <- function(m) {
			loci_ok <- names(ALL_MOTIFS_BY_LOCUS)[vapply(ALL_MOTIFS_BY_LOCUS, function(v) m %in% v, logical(1))]
			ad <- store$cons_anchored_data
			sort(unique(ad$SUBGROUP[ad$LOCI %in% loci_ok]))
		}
		
		output$f1_status <- shiny::renderUI({
			if (!is.null(store$comparison_df) || length(store$comparison_reasons) == 0) return(NULL)
			shiny::div(class = "alert alert-danger", style = "font-size:11px;",
				shiny::strong("F1 unavailable: "),
				shiny::tags$ul(lapply(utils::head(store$comparison_reasons, 6), shiny::tags$li)))
		})
		
		output$choice_table <- shiny::renderUI({
			ms <- detected_motifs(); shiny::req(length(ms) > 0)
			fmt <- function(x) if (is.na(x)) "n/a" else sprintf("%.2f", x)
			shiny::tagList(lapply(ms, function(m) {
				n_c <- length(unique(store$cons_hits_detected$ALLELE_ID[store$cons_hits_detected$MOTIF == m]))
				n_p <- length(unique(store$pwm_hits_detected$ALLELE_ID[store$pwm_hits_detected$MOTIF == m]))
				f_c <- .f1_by_method(m, "Consensus"); f_p <- .f1_by_method(m, "PWM")
				default <- if (!is.na(f_p) && !is.na(f_c) && f_p > f_c) "PWM" else "Consensus"
				cid <- session$ns(paste0("custom_", m))
				shiny::div(style = "border-bottom:1px solid #eee; padding:6px 0;",
					shiny::fluidRow(
						shiny::column(3, shiny::strong(m)),
						shiny::column(3, sprintf("Alleles detected: Consensus %d / PWM %d", n_c, n_p)),
						shiny::column(2, sprintf("F1: %s / %s", fmt(f_c), fmt(f_p))),
						shiny::column(4, shiny::radioButtons(session$ns(paste0("src_", m)), NULL,
							choices = c("Consensus", "PWM", "Exclude"), selected = default, inline = TRUE))),
					shiny::checkboxInput(cid, "Customize per subgroup", FALSE),
					shiny::conditionalPanel(sprintf("input['%s']", cid),
						shiny::uiOutput(session$ns(paste0("sg_ui_", m)))))
			}))
		})
		
		shiny::observe({
			ms <- detected_motifs()
			lapply(ms, function(m) local({
				mm <- m
				output[[paste0("sg_ui_", mm)]] <- shiny::renderUI({
					def <- input[[paste0("src_", mm)]] %||% "Consensus"
					shiny::tagList(lapply(.sgs_for_motif(mm), function(sg)
						shiny::fluidRow(
							shiny::column(3, shiny::tags$em(sg)),
							shiny::column(9, shiny::radioButtons(session$ns(paste0("sgsrc_", mm, "__", sg)), NULL,
								choices = c("Consensus", "PWM", "Exclude"), selected = def, inline = TRUE)))))
				})
			}))
		})
		
		shiny::observeEvent(input$btn_run_analyses, {
			ms <- detected_motifs()
			choice <- vapply(ms, function(m) input[[paste0("src_", m)]] %||% "Exclude", character(1))
			names(choice) <- ms
			overrides <- list()
			for (m in ms) if (isTRUE(input[[paste0("custom_", m)]])) {
				sgs <- .sgs_for_motif(m)
				v <- vapply(sgs, function(sg) input[[paste0("sgsrc_", m, "__", sg)]] %||% choice[[m]], character(1))
				names(v) <- sgs; overrides[[m]] <- v
			}
			any_on <- any(choice != "Exclude") || any(unlist(overrides) != "Exclude")
			if (!any_on) { shiny::showNotification("❌ Select at least one motif.", type = "error"); return() }
			tryCatch({
				res <- build_merged_dataset(store$cons_final_df, store$pwm_final_df, choice,
					store$cons_anchored_data, overrides)
				store$merged_df <- res$data; store$merged_check <- res; store$merged_stamp <- Sys.time()
				shiny::showNotification(sprintf("✅ Dataset ready: %d alleles, %d columns.", nrow(res$data), ncol(res$data)), type = "message")
			}, error = function(e) shiny::showNotification(paste0("❌ ", e$message), type = "error", duration = 10))
		})
		
		output$schema_status <- shiny::renderUI({
			r <- store$merged_check; shiny::req(r)
			msgs <- list()
			if (length(r$missing_cols) > 0)
				msgs <- c(msgs, list(shiny::div(class = "alert alert-danger",
					paste("Missing required columns:", paste(r$missing_cols, collapse = ", ")))))
			if (length(r$empty_motifs) > 0)
				msgs <- c(msgs, list(shiny::div(class = "alert alert-warning",
					paste("No detection available for:", paste(r$empty_motifs, collapse = ", "), "— motif left empty."))))
			if (r$n_not_anchored > 0)
				msgs <- c(msgs, list(shiny::div(class = "alert alert-info",
					sprintf("%d allele(s) excluded because L-PART1 could not be anchored.", r$n_not_anchored))))
			if (r$n_func_na > 0)
				msgs <- c(msgs, list(shiny::div(class = "alert alert-warning",
					sprintf("%d allele(s) with unrecognised functionality (NA).", r$n_func_na))))
			if (length(msgs) == 0)
				msgs <- list(shiny::div(class = "alert alert-success", shiny::icon("check-circle"), " Dataset consistent with the application schema."))
			shiny::tagList(msgs)
		})
		
		output$merged_table <- DT::renderDT({
			shiny::req(store$merged_df)
			DT::datatable(store$merged_df, rownames = FALSE, class = "cell-border stripe hover compact",
				options = list(pageLength = 20, scrollX = TRUE))
		})
		output$dl_merged_csv <- shiny::downloadHandler(
			filename = function() paste0("DatasetForAnalyses_", Sys.Date(), ".csv"),
			content  = function(file) { shiny::req(store$merged_df); write.csv(store$merged_df, file, row.names = FALSE) })
		output$dl_merged_xlsx <- shiny::downloadHandler(
			filename = function() paste0("DatasetForAnalyses_", Sys.Date(), ".xlsx"),
			content  = function(file) { shiny::req(store$merged_df); openxlsx::write.xlsx(store$merged_df, file) })
		
		return(list(
			base_info_consensus     = shiny::reactive(store$cons_base_info),
			final_df_consensus      = shiny::reactive(store$cons_final_df),
			hits_detected_consensus = shiny::reactive(store$cons_hits_detected),
			base_info_pwm           = shiny::reactive(store$pwm_base_info),
			final_df_pwm            = shiny::reactive(store$pwm_final_df),
			hits_detected_pwm       = shiny::reactive(store$pwm_hits_detected),
			comparison_df           = shiny::reactive(store$comparison_df),
			merged_df               = shiny::reactive(store$merged_df),
			go_analyses             = shiny::reactive(store$merged_stamp)
		))
	})
}
