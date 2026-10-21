# detecter_francais.R — à lancer depuis la racine du projet
fichiers <- c("Global.R", "Ui.R", "Server.R",
	list.files(c("R", "modules", "server", "ui"), pattern = "\\.R$", full.names = TRUE))

marqueurs <- paste0("[éèêëàâùûôîïç]|",
	"\\b(aucun|aucune|toutes?|tous|données|séquences?|espèces?|erreur|taille|",
	"comparaison|résultats?|fonctionnel|légende|couleurs?|lancer|télécharg\\w*)\\b")

res <- list()
for (f in fichiers) {
	exprs <- tryCatch(parse(f, encoding = "UTF-8", keep.source = TRUE), error = function(e) NULL)
	if (is.null(exprs)) next
	pd <- utils::getParseData(exprs)
	chaines <- pd[pd$token == "STR_CONST", ]
	if (nrow(chaines) == 0) next
	garde <- grepl(marqueurs, chaines$text, ignore.case = TRUE, perl = TRUE)
	if (any(garde))
		res[[length(res) + 1]] <- data.frame(fichier = f, ligne = chaines$line1[garde],
			texte = substr(chaines$text[garde], 1, 90),
			stringsAsFactors = FALSE)
}
res <- do.call(rbind, res)
res <- res[order(res$fichier, res$ligne), ]
cat("\n=== TEXTES POSSIBLEMENT EN FRANÇAIS ===\n")
print(res, row.names = FALSE)
write.csv(res, "textes_francais.csv", row.names = FALSE, fileEncoding = "UTF-8")