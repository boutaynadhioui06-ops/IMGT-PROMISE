# outils/detecter_images_inutilisees
fichiers <- list.files("www", recursive = TRUE)
fichiers <- fichiers[!dir.exists(file.path("www", fichiers))]

sources <- c("Global.R", "Ui.R", "Server.R",
	list.files(c("R", "modules", "server", "ui", "report_templates"),
		pattern = "\\.(R|Rmd|css|js)$", full.names = TRUE),
	list.files("www", pattern = "\\.(css|js)$", full.names = TRUE))
texte <- tolower(paste(unlist(lapply(sources, readLines, warn = FALSE, encoding = "UTF-8")),
	collapse = "\n"))

utilise <- vapply(fichiers, function(f) grepl(tolower(basename(f)), texte, fixed = TRUE),
	logical(1))
res <- data.frame(fichier   = fichiers,
	taille_ko = round(file.size(file.path("www", fichiers)) / 1024, 1),
	utilise   = utilise)
inutiles <- res[!res$utilise, ]

cat("\n=== FICHIERS DE www/ UTILISÉS :", sum(res$utilise), "===\n")
cat("\n=== FICHIERS DE www/ JAMAIS RÉFÉRENCÉS ===\n")
print(inutiles, row.names = FALSE)
cat("\nPoids total des fichiers inutilisés :", sum(inutiles$taille_ko), "Ko\n")

if (nrow(inutiles) > 0) {
	cat("\n=== COMMANDES À COLLER DANS LE TERMINAL ===\n")
	cat("mkdir archive\\www_inutilises\n")
	cat(sprintf('move "www\\%s" archive\\www_inutilises\\',
		gsub("/", "\\\\", inutiles$fichier)), sep = "\n")
}