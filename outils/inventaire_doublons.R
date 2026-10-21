# inventaire_doublons.R

# 1. Récupérer la liste des fichiers sourcés, dans l'ordre de Global.R
lignes <- readLines("Global.R", warn = FALSE, encoding = "UTF-8")
lignes <- lignes[!grepl("^\\s*#", lignes)]                 # ignorer les lignes commentées
m <- regmatches(lignes, regexpr('source\\("[^"]+"', lignes))
fichiers <- sub('"$', "", sub('^source\\("', "", m))

# 2. Lister tous les objets définis au premier niveau de chaque fichier
objets <- list()
for (i in seq_along(fichiers)) {
  f <- fichiers[i]
  if (!file.exists(f)) { message("INTROUVABLE : ", f); next }
  exprs <- tryCatch(
    parse(f, encoding = "UTF-8", keep.source = FALSE),
    error = function(e) { message("ERREUR de lecture : ", f, " -> ", conditionMessage(e)); NULL }
  )
  if (is.null(exprs)) next
  for (e in exprs) {
    if (is.call(e) && as.character(e[[1]])[1] %in% c("<-", "=", "<<-")) {
      lhs <- e[[2]]
      nom <- if (is.name(lhs)) as.character(lhs) else if (is.character(lhs)) lhs else NA
      if (!is.na(nom))
        objets[[length(objets) + 1]] <- data.frame(objet = nom, fichier = f, ordre = i,
          stringsAsFactors = FALSE)
    }
  }
}
inv <- do.call(rbind, objets)

# 3. Garder uniquement les objets définis plusieurs fois
dup <- inv[inv$objet %in% inv$objet[duplicated(inv$objet)], ]
dup <- dup[order(dup$objet, dup$ordre), ]
dup$statut <- ifelse(dup$ordre == ave(dup$ordre, dup$objet, FUN = max),
  "ACTIVE", "écrasée")

# 4. Fichiers présents dans R/ et modules/ mais jamais sourcés
tous <- list.files(c("R", "modules"), pattern = "\\.R$", full.names = TRUE)
non_sources <- setdiff(tous, fichiers)

# 5. Résultats
cat("\n=== OBJETS DÉFINIS PLUSIEURS FOIS ===\n")
print(dup, row.names = FALSE)
cat("\n=== FICHIERS JAMAIS SOURCÉS ===\n")
print(non_sources)

write.csv(dup, "inventaire_doublons.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("\nRésultat enregistré dans inventaire_doublons.csv\n")