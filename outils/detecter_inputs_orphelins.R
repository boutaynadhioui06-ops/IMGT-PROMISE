# detecter_inputs_orphelins.R — à lancer depuis la racine du projet
lire <- function(f) {
	l <- readLines(f, warn = FALSE, encoding = "UTF-8")
	paste(l[!grepl("^\\s*#", l)], collapse = "\n")
}
esc <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
compte <- function(id, txt) {
	motif <- paste0("(?<![A-Za-z0-9._])", esc(id), "(?![A-Za-z0-9._])")
	sum(lengths(regmatches(txt, gregexpr(motif, txt, perl = TRUE))))
}

fn_input <- "(?:selectInput|selectizeInput|radioButtons|checkboxInput|checkboxGroupInput|numericInput|sliderInput|textInput|textAreaInput|fileInput|actionButton|dateInput|pickerInput)"
motif_def <- paste0(fn_input, "\\(\\s*(?:(?:session\\$)?ns\\()?\"\\K[A-Za-z0-9_.]+(?=\")")

ids_definis <- function(txt) regmatches(txt, gregexpr(motif_def, txt, perl = TRUE))[[1]]
ids_lus     <- function(txt) unique(regmatches(txt, gregexpr("input\\$\\K[A-Za-z0-9_.]+", txt, perl = TRUE))[[1]])

resultats <- list()
analyser <- function(portee, f_ui, f_srv) {
	txt_ui  <- paste(vapply(f_ui,  lire, ""), collapse = "\n")
	txt_srv <- paste(vapply(f_srv, lire, ""), collapse = "\n")
	tout    <- if (identical(f_ui, f_srv)) txt_ui else paste(txt_ui, txt_srv, sep = "\n")
	defs    <- ids_definis(txt_ui)
	lus     <- ids_lus(txt_srv)
	
	doubles   <- unique(defs[duplicated(defs)])
	inutiles  <- unique(defs[vapply(unique(defs), compte, numeric(1), txt = tout) <= 1])
	manquants <- setdiff(lus, defs)
	
	ajouter <- function(type, ids) if (length(ids) > 0)
		resultats[[length(resultats) + 1]] <<- data.frame(portee = portee, type = type, id = sort(ids),
			stringsAsFactors = FALSE)
	ajouter("A_defini_mais_jamais_utilise", inutiles)
	ajouter("B_lu_mais_jamais_defini",      manquants)
	ajouter("C_defini_en_double",           doubles)
}

# Application principale : Ui.R + ui/  <->  Server.R + server/
analyser("APP",
	c("Ui.R",     list.files("ui",     pattern = "\\.R$", full.names = TRUE)),
	c("Server.R", list.files("server", pattern = "\\.R$", full.names = TRUE)))

# Modules : interface et serveur dans le même fichier
for (f in list.files("modules", pattern = "\\.R$", full.names = TRUE))
	analyser(f, f, f)

res <- if (length(resultats) > 0) do.call(rbind, resultats) else
	data.frame(portee = character(), type = character(), id = character())

cat("\n=== INPUTS À EXAMINER ===\n")
print(res, row.names = FALSE)
write.csv(res, "inputs_orphelins.csv", row.names = FALSE, fileEncoding = "UTF-8")