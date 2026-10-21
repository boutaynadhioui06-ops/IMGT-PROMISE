# inventaire_dependances.R — à lancer depuis la racine du projet
fichiers <- list.files(c("R", "modules"), pattern = "\\.R$", full.names = TRUE)

texte <- lapply(fichiers, function(f) {
	l <- readLines(f, warn = FALSE, encoding = "UTF-8")
	l[!grepl("^\\s*#", l)]
})
names(texte) <- fichiers

# 1. objets définis au premier niveau de chaque fichier
defs <- list()
for (f in fichiers) {
	exprs <- tryCatch(parse(f, encoding = "UTF-8", keep.source = FALSE), error = function(e) NULL)
	if (is.null(exprs)) next
	for (e in exprs) {
		if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") && is.name(e[[2]])) {
			defs[[length(defs) + 1]] <- data.frame(objet = as.character(e[[2]]), fichier = f,
				stringsAsFactors = FALSE)
		}
	}
}
defs <- do.call(rbind, defs)
defs <- defs[defs$objet != "%||%", ]

# 2. pour chaque objet, quels AUTRES fichiers l'utilisent
esc <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
deps <- list()
for (i in seq_len(nrow(defs))) {
	motif <- paste0("(?<![A-Za-z0-9._])", esc(defs$objet[i]), "(?![A-Za-z0-9._])")
	for (g in setdiff(fichiers, defs$fichier[i])) {
		n <- sum(lengths(regmatches(texte[[g]], gregexpr(motif, texte[[g]], perl = TRUE))))
		if (n > 0)
			deps[[length(deps) + 1]] <- data.frame(fichier = g, depend_de = defs$fichier[i],
				utilise = defs$objet[i], stringsAsFactors = FALSE)
	}
}
deps <- do.call(rbind, deps)

# 3. résumé : un fichier -> les fichiers dont il dépend, et ce qu'il y prend
resume <- aggregate(utilise ~ fichier + depend_de, data = deps,
	FUN = function(x) paste(sort(unique(x)), collapse = ", "))
resume <- resume[order(resume$fichier, resume$depend_de), ]

cat("\n=== DÉPENDANCES ENTRE FICHIERS ===\n")
for (f in unique(resume$fichier)) {
	cat("\n", f, "\n", sep = "")
	sous <- resume[resume$fichier == f, ]
	for (k in seq_len(nrow(sous)))
		cat("   <- ", sous$depend_de[k], "  [", sous$utilise[k], "]\n", sep = "")
}

cat("\n=== FICHIERS SANS AUCUNE DÉPENDANCE (autonomes) ===\n")
print(setdiff(fichiers, resume$fichier))

write.csv(resume, "dependances.csv", row.names = FALSE, fileEncoding = "UTF-8")