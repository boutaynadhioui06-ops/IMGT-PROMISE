# detecter_outputs_orphelins.R 
srv <- readLines("Server.R", warn = FALSE, encoding = "UTF-8")
ui  <- readLines("Ui.R",     warn = FALSE, encoding = "UTF-8")
srv <- srv[!grepl("^\\s*#", srv)]
ui  <- ui[!grepl("^\\s*#", ui)]

# 1. sorties définies dans Server.R : output$xxx <- ...
m <- regmatches(srv, gregexpr("output\\$[A-Za-z0-9_.]+\\s*<-", srv))
out_srv <- unique(sub("\\s*<-$", "", sub("^output\\$", "", unlist(m))))

# 2. ids d'écran : xxxOutput("id"), downloadButton("id"), downloadLink("id")
ids_depuis <- function(lignes) {
	m <- regmatches(lignes, gregexpr("[A-Za-z]*(Output|downloadButton|downloadLink)\\(\\s*\"[A-Za-z0-9_.]+\"", lignes))
	unique(sub("^.*\\(\\s*\"", "", sub("\"$", "", unlist(m))))
}
ids_ui  <- ids_depuis(ui)
ids_srv <- ids_depuis(srv)   # écrans créés par renderUI() dans Server.R

# 3. références via conditionalPanel : condition = "output.xxx"
m2 <- regmatches(ui, gregexpr("output\\.[A-Za-z0-9_]+", ui))
cond <- unique(sub("^output\\.", "", unlist(m2)))

refs <- unique(c(ids_ui, ids_srv, cond))

orphelines <- sort(setdiff(out_srv, refs))      # code serveur sans écran
manquantes <- sort(setdiff(ids_ui, out_srv))    # écran sans code serveur

cat("\n=== A. SORTIES DU SERVEUR SANS ECRAN (calculent pour rien) ===\n")
print(orphelines)
cat("\n=== B. ECRANS DE Ui.R SANS CODE SERVEUR (vides ou boutons inactifs) ===\n")
print(manquantes)

# write.csv(rbind(
# 	data.frame(type = "A_serveur_sans_ecran", id = orphelines),
# 	data.frame(type = "B_ecran_sans_serveur", id = manquantes)),
# 	"outputs_orphelins.csv", row.names = FALSE, fileEncoding = "UTF-8")