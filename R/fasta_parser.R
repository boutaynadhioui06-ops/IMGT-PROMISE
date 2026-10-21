# ═══════════════════════════════════════════════════════════════
#  Parser FASTA
# ═══════════════════════════════════════════════════════════════

# Dépend de : dataset_helpers.R (clean_subgroup, clean_functionality)
# Packages  : base R


# Parser 5'UTR
parse_5utr_fasta <- function(fasta_file) {
	
	lines      <- readLines(fasta_file, warn = FALSE)
	header_idx <- which(startsWith(lines, ">"))
	
	if (length(header_idx) == 0)
		stop("No FASTA header found (line starting with'>')")
	
	headers <- lines[header_idx]
	n       <- length(headers)
	
	df <- data.frame(
		ACCESSION     = character(n),
		LOCI          = character(n),
		SUBGROUP      = character(n),
		GENE          = character(n),
		ALLELE        = character(n),
		SPECIES       = character(n),
		FUNCTIONALITY = character(n),
		POS_START     = integer(n),
		POS_END       = integer(n),
		LENGTH_NT     = integer(n),
		stringsAsFactors = FALSE
	)
	
	for (i in seq_len(n)) {
		
		h      <- sub("^>", "", headers[i])
		fields <- strsplit(h, "\\|")[[1]]
		
		#  ACCESSION ──
		df$ACCESSION[i] <- trimws(fields[1])
		
		#  GENE / ALLELE
		ga <- trimws(fields[2])
		df$LOCI[i]     <- gsub("^([A-Za-z]+).*", "\\1", ga)
		df$SUBGROUP[i] <- gsub("[-*].*$", "", ga)
		df$GENE[i]     <- strsplit(ga, "\\*")[[1]][1]
		df$ALLELE[i]   <- ga
		
		#  SPECIES  
		if (length(fields) >= 3) {
			sp <- trimws(fields[3])
			df$SPECIES[i] <- gsub("_.*$", "", sp)
		}
		
		# FUNCTIONALITY
		if (length(fields) >= 4) {
			df$FUNCTIONALITY[i] <- clean_functionality(fields[4])
		}
		
		#  POSITION 
		if (length(fields) >= 6) {
			pos   <- trimws(fields[6])
			parts <- strsplit(pos, "\\.\\.")[[1]]
			df$POS_START[i] <- suppressWarnings(as.integer(gsub("[^0-9]", "", parts[1])))
			df$POS_END[i]   <- suppressWarnings(as.integer(gsub("[^0-9]", "", parts[2])))
		}
		
		# Champ 7 : LENGTH  
		if (length(fields) >= 7) {
			df$LENGTH_NT[i] <- suppressWarnings(
				as.integer(regmatches(fields[7], regexpr("^[0-9]+", fields[7])))
			)
		}
	}
	
	n_before <- nrow(df)
	df <- df[!grepl("\\([IVX]+\\)", df$GENE), ]
	n_roman  <- n_before - nrow(df)
	
	
	n_before_v <- nrow(df)
	df <- df[grepl("V$", df$LOCI, ignore.case = FALSE), ]
	n_v_removed <- n_before_v - nrow(df)
	
	
	n_before_sg <- nrow(df)
	
	
	exclude_pattern <- paste0(
		"[0-9]+(OR|S[0-9]|H[0-9]|L[0-9]|-D)(-|$)|",
		"[0-9]+/?OR[0-9]*$|",
		"[0-9]+S[0-9]+$|[0-9]+H[0-9]+$|[0-9]+L[0-9]+$"	)
	
	df_sg_before <- nrow(df)
	df <- df[
		is.na(df$SUBGROUP) |
			!grepl(exclude_pattern, df$SUBGROUP, perl = TRUE, ignore.case = TRUE),
	]
	n_sg_removed <- df_sg_before - nrow(df)
	n_after      <- nrow(df)
	
	n_after      <- nrow(df)
	
	df$SUBGROUP <- clean_subgroup(df$SUBGROUP)
	
	cat("   ✓ Genes with Roman numerals removed         :", n_roman, "\n")
	cat("   ✓ Non-V genes removed (J/D/C/E/...)         :", n_v_removed, "\n")
	cat("   ✓ Non-standard subgroups removed (OR/S/H/L) :", n_sg_removed, "\n")
	cat("   ✓ V genes retained                          :", n_after, "\n")
	cat("   ℹ Short 5'UTR (< 100 nt) filtered later,\n")
	cat("     after join with motifs (see Dataset Build step).\n\n")
	
	return(df)
}


# Parser Motif
parse_motif_fasta <- function(fasta_file, motif_name) {
	
	lines      <- readLines(fasta_file, warn = FALSE)
	header_idx <- which(startsWith(lines, ">"))
	
	if (length(header_idx) == 0)
		stop(paste("No header detected in the motif file", motif_name))
	
	n           <- length(header_idx)
	seq_end_idx <- c(header_idx[-1] - 1, length(lines))
	
	seq_col   <- paste0("SEQ_",       motif_name)
	pos_start <- paste0("POS_START_", motif_name)
	pos_end   <- paste0("POS_END_",   motif_name)
	
	df <- data.frame(
		ALLELE        = character(n),
		GENE          = character(n),
		SUBGROUP      = character(n),
		LOCI          = character(n),
		SPECIES       = character(n),
		FUNCTIONALITY = character(n),
		stringsAsFactors = FALSE
	)
	df[[seq_col]]   <- character(n)
	df[[pos_start]] <- integer(n)
	df[[pos_end]]   <- integer(n)
	
	for (i in seq_len(n)) {
		
		h      <- sub("^>", "", lines[header_idx[i]])
		fields <- strsplit(h, "\\|")[[1]]
		
		#  GENE / ALLELE 
		ga <- trimws(fields[2])
		df$LOCI[i]     <- gsub("^([A-Za-z]+).*", "\\1", ga)
		df$SUBGROUP[i] <- clean_subgroup(gsub("[-*].*$", "", ga))
		df$GENE[i]     <- strsplit(ga, "\\*")[[1]][1]
		df$ALLELE[i]   <- ga
		
		#  SPECIES 
		if (length(fields) >= 3) {
			df$SPECIES[i] <- gsub("_.*$", "", trimws(fields[3]))
		}
		
		# FUNCTIONALITY 
		if (length(fields) >= 4) {
			df$FUNCTIONALITY[i] <- clean_functionality(fields[4])
		}
		
		# SÉQUENCE 
		seq_lines        <- lines[(header_idx[i] + 1):seq_end_idx[i]]
		df[[seq_col]][i] <- paste(trimws(seq_lines), collapse = "")
		
		# POSITION 
		if (length(fields) >= 6) {
			parts              <- strsplit(trimws(fields[6]), "\\.\\.")[[1]]
			df[[pos_start]][i] <- suppressWarnings(as.integer(gsub("[^0-9]", "", parts[1])))
			df[[pos_end]][i]   <- suppressWarnings(as.integer(gsub("[^0-9]", "", parts[2])))
		}
	}
	
	# Conserver uniquement les gènes de type V 
	df <- df[grepl("V$", df$LOCI, ignore.case = FALSE), ]
	
	return(df)
}