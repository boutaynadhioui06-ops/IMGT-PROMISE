# detecter_fonctions_mortes.R — à lancer depuis la racine du projet
fichiers <- c(list.files(c("R", "modules"), pattern = "\\.R$", full.names = TRUE),
	"Global.R", "Server.R", "Ui.R")
fichiers <- fichiers[file.exists(fichiers)]

# texte sans les lignes commentées
texte <- lapply(fichiers, function(f) {
	l <- readLines(f, warn = FALSE, encoding = "UTF-8")
	l[!grepl("^\\s*#", l)]
})
names(texte) <- fichiers
tout <- unlist(texte)

# noms des fonctions définies au premier niveau
defs <- list()
for (f in setdiff(fichiers, c("Global.R", "Server.R", "Ui.R"))) {
	exprs <- tryCatch(parse(f, encoding = "UTF-8", keep.source = FALSE), error = function(e) NULL)
	if (is.null(exprs)) next
	for (e in exprs) {
		if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=") &&
				is.call(e[[3]]) && identical(as.character(e[[3]][[1]]), "function")) {
			nom <- as.character(e[[2]])
			if (length(nom) == 1)
				defs[[length(defs) + 1]] <- data.frame(fonction = nom, fichier = f, stringsAsFactors = FALSE)
		}
	}
}
defs <- do.call(rbind, defs)
defs <- defs[defs$fonction != "%||%", ]

# nombre d'occurrences de chaque nom dans tout le code (hors commentaires)
esc <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
defs$occurrences <- vapply(defs$fonction, function(nom) {
	motif <- paste0("(?<![A-Za-z0-9._])", esc(nom), "(?![A-Za-z0-9._])")
	sum(lengths(regmatches(tout, gregexpr(motif, tout, perl = TRUE))))
}, numeric(1))

# 1 occurrence = seulement sa propre définition = jamais appelée
morts <- defs[defs$occurrences <= 1, ]
morts <- morts[order(morts$fichier, morts$fonction), ]

cat("\n=== FONCTIONS JAMAIS APPELÉES (candidates à la suppression) ===\n")
print(morts, row.names = FALSE)
write.csv(morts, "fonctions_mortes.csv", row.names = FALSE, fileEncoding = "UTF-8")