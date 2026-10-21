# ═══════════════════════════════════════════════════════════════════════════
# IMGT/PROMISE — RESOURCES 
# ═══════════════════════════════════════════════════════════════════════════

.res_note <- function(...) div(class = "res-note", icon("info-circle"), ...)

.faq_item <- function(question, ...) {
	tags$details(class = "faq-item", tags$summary(question), div(class = "faq-body", ...))
}

.doi_link <- function(doi) tags$a(href = paste0("https://doi.org/", doi), target = "_blank",
	rel = "noopener", paste0("doi:", doi))

# ═══════════════════════════════════════════════════════════════
# DOCUMENTATION
# ═══════════════════════════════════════════════════════════════
ui_tab_resources_documentation <- function()
	tabItem(tabName = "resources_documentation",
		tab_header("book", "Documentation",
			"How do you use IMGT/PROMISE, from data import to interpretation of results?",
			"#0B2A5B", "#1565C0"),
		
		tabsetPanel(
			
			# ── Overview ─────────────────────────────────────────────
			tabPanel("\U0001F9ED Overview", br(),
				box(width = 12, status = "primary", solidHeader = FALSE,
					h4(icon("bullseye"), strong(" Objective")),
					p("IMGT/PROMISE provides an integrated environment for exploring and analyzing regulatory ",
						"features of 5′UTR regions in IG and TR V-genes available in IMGT/LIGM-DB. It combines ",
						"statistical analyses, visualization, comparative approaches and automatic motif discovery ",
						"for species lacking curated motif annotations."),
					hr(),
					h4(icon("route"), strong(" Workflow")),
					div(class = "res-workflow",
						div(class = "res-phase", span("1"), strong("Dataset"),
							p("Import 5′UTR and motif files, build the final dataset.")),
						div(class = "res-arrow", icon("chevron-right")),
						div(class = "res-phase", span("2"), strong("Descriptive & conservation"),
							p("Counts, functionality, presence frequencies, logos, regression.")),
						div(class = "res-arrow", icon("chevron-right")),
						div(class = "res-phase", span("3"), strong("Advanced analyses"),
							p("Positioning, similarity, diversity, co-occurrence, PERMANOVA.")),
						div(class = "res-arrow", icon("chevron-right")),
						div(class = "res-phase", span("4"), strong("Motif discovery"),
							p("Consensus and PWM scanners, compared against confirmed annotations."))
					),
					hr(),
					h4(icon("users"), strong(" Target audience")),
					tags$ul(
						tags$li("Immunogenetics researchers working on IG/TR V-gene regulatory regions."),
						tags$li("Bioinformaticians curating or extending IMGT/LIGM-DB motif annotations."),
						tags$li("Reviewers and students who need a transparent, reproducible workflow.")
					),
					.res_note(strong(" Scope: "),
						"motif-based analyses are restricted to immunoglobulins (IGHV, IGKV, IGLV); TR genes currently ",
						"lack sufficient curated 5′UTR motif annotations.")
				)
			),
			
			# ── Getting started ──────────────────────────────────────
			tabPanel("\U0001F680 Getting Started", br(),
				box(width = 12, status = "success", solidHeader = FALSE,
					h4(icon("list-ol"), strong(" Six steps to your first analysis")),
					tags$ol(class = "res-steps",
						tags$li(strong("Import the 5′UTR file"), " (\U0001F4CA Dataset Construction → Import 5′UTR)."),
						tags$li(strong("Import the motif files"), " (Dataset Construction → Import Motifs)."),
						tags$li(strong("Build the dataset"), " (Dataset Construction → Final Dataset → Build Dataset)."),
						tags$li(strong("Pick a question"), " in the left menu; each page states its research question at the top."),
						tags$li(strong("Set parameters and click Launch"), "; most pages offer cascading Species → Locus → Subgroup filters."),
						tags$li(strong("Export"), " figures (PNG/PDF) and tables (CSV/Excel) from each page, or generate a full ",
							"report from Resources → Analysis Report.")
					),
					.res_note(strong(" Tip: "),
						"analyses only use what you have already launched; the Analysis Report shows which sections are ready.")
				)
			),
			
			# ── Input data ───────────────────────────────────────────
			tabPanel("\U0001F4C2 Input Data", br(),
				box(width = 12, status = "warning", solidHeader = FALSE,
					h4(icon("file-code"), strong(" FASTA files")),
					p("The 5′UTR file (and, for the scanner, the L-PART1 file) follow the IMGT header convention:"),
					tags$code(">Accession|Allele|Species|Functionality|...|Genomic_coordinates|..."),
					p(style = "margin-top:8px;", "Only the header is used for the dataset; sequences are read by the scanner only."),
					hr(),
					h4(icon("table"), strong(" Data dictionary")),
					p("Columns of the final dataset (",
						tags$code("<MOTIF>"), " is replaced by each imported motif name)."),
					DT::DTOutput("doc_data_dictionary"),
					br(),
					.res_note(strong(" Excluded automatically: "),
						"alleles of Roman-numeral clan subgroups (e.g. IGHVII, IGHVIII) are out of scope.")
				)
			),
			
			# ── Methods ──────────────────────────────────────────────
			tabPanel("\U0001F52C Methods", br(),
				box(width = 12, status = "info", solidHeader = FALSE,
					h4(icon("puzzle-piece"), strong(" Consensus-based detection (baseline)")),
					p("Deterministic reproduction of the IMGT manual annotation procedure: fixed IUPAC patterns from ",
						"the literature with a fixed mismatch tolerance. No parameter is learned from the data, so the ",
						"method serves as an external reference. The search is limited to the 600 nt nearest to L-PART1 ",
						"(~500 bp prescribed by the annotation procedure plus a safety margin)."),
					DT::DTOutput("doc_motif_defs"),
					hr(),
					h4(icon("chart-line"), strong(" PWM-based probabilistic detection")),
					tags$ol(
						tags$li("A PWM is calibrated per subgroup from confirmed instances (fallback to a locus-level PWM)."),
						tags$li("The search is restricted to the plausible distance window to L-PART1, derived from the reference set (with a 15 nt margin)."),
						tags$li("The maximum PWM score is compared with a Monte Carlo null distribution to obtain a p-value."),
						tags$li("p-values are adjusted by Benjamini–Hochberg (q-value); a hard mismatch filter is applied in addition.")
					),
					hr(),
					h4(icon("scale-balanced"), strong(" Method comparison")),
					p("Both methods are evaluated against confirmed reference instances (true positives / false negatives) ",
						"and against nucleotide-shuffled controls (true negatives / false positives). PWM uses stratified ",
						"k-fold cross-validation (k = 5) with calibration rebuilt on training folds only; the consensus ",
						"method has no learned parameter and is evaluated in a single pass."),
					hr(),
					h4(icon("chart-bar"), strong(" Statistical analyses")),
					tags$ul(
						tags$li(strong("Conservation:"), " Chi-square or Fisher's exact test (Bonferroni-corrected), logistic regression with cross-validation and ROC."),
						tags$li(strong("Diversity:"), " Shannon, Simpson, richness and evenness; Kruskal–Wallis tests on individual sequences."),
						tags$li(strong("Similarity:"), " pairwise identity of consensus sequences, hierarchical clustering."),
						tags$li(strong("Multivariate:"), " PERMANOVA and PCoA on motif profiles.")
					)
				)
			),
			
			# ── Interpretation ───────────────────────────────────────
			tabPanel("\U0001F3AF Interpretation", br(),
				box(width = 12, status = "danger", solidHeader = FALSE,
					h4(icon("question-circle"), strong(" How to read the results")),
					tags$ul(class = "res-steps",
						tags$li(strong("p-value:"), " p < 0.05 ★, p < 0.01 ★★, p < 0.001 ★★★, ns = not significant."),
						tags$li(strong("q-value:"), " p-value corrected for multiple testing; compare it, not the raw p-value, with the threshold."),
						tags$li(strong("Odds ratio:"), " relative to the reference level shown in the regression table, not to the average."),
						tags$li(strong("PWM score:"), " only comparable within the same motif and calibration."),
						tags$li(strong("Sensitivity / specificity:"), " a high sensitivity with low precision signals many false detections."),
						tags$li(strong("PERMANOVA:"), " a significant result means group centroids differ, not which groups differ (use pairwise tests)."),
						tags$li(strong("Shannon / Simpson:"), " compare relative values between groups rather than absolute thresholds."),
						tags$li(strong("Small groups:"), " frequencies computed on few sequences are unstable; check the N shown in tables.")
					)
				)
			),
			
			# ── FAQ ──────────────────────────────────────────────────
			tabPanel("❓ FAQ & Troubleshooting", br(),
				box(width = 12, status = "primary", solidHeader = FALSE,
					.faq_item("“No ALLELE+SPECIES match between the two FASTA files”",
						"The scanner joins the 5′UTR and L-PART1 files on allele + species. Check that both files use ",
						"the same allele names and species spelling, and that alleles are not all Roman-numeral clans (excluded)."),
					.faq_item("Some alleles are missing after the scan",
						"Alleles are dropped when L-PART1 cannot be found with 0 mismatch in the 5′UTR, when no upstream ",
						"sequence exists, or when they belong to a Roman-numeral clan."),
					.faq_item("“Method comparison unavailable” / SEQ_UPSTREAM missing",
						"The reference dataset must contain the SEQ_UPSTREAM column and the SEQ_<MOTIF> columns. ",
						"Annotation itself works without it."),
					.faq_item("A motif is not detected for a species or subgroup",
						"The PWM method needs at least 15 confirmed instances at locus level (8 per subgroup). See the ",
						"calibration warnings displayed after the scan."),
					.faq_item("An analysis page is empty",
						"Build the dataset first (Dataset Construction → Final Dataset). TR genes are not covered by motif analyses."),
					.faq_item("The report is missing a section",
						"Sections appear as unavailable until the corresponding data exist (dataset built, scanner comparison computed).")
				)
			),
			
			# ── Glossary ─────────────────────────────────────────────
			tabPanel("\U0001F4D6 Glossary", br(),
				box(width = 12, status = "info", solidHeader = FALSE, DT::DTOutput("doc_glossary"))
			)
		)
	)

# ═══════════════════════════════════════════════════════════════
# ANALYSIS REPORT (dynamic)
# ═══════════════════════════════════════════════════════════════
ui_tab_resources_report <- function()
	tabItem(tabName = "resources_report",
		tab_header("file-alt", "Analysis Report",
			"Can the analyses of this session be exported as a reproducible, dated report?",
			"#1565C0", "#06A77D"),
		
		fluidRow(
			column(7,
				box(width = 12, title = "① Report content", status = "primary", solidHeader = TRUE,
					fluidRow(
						column(8, textInput("report_title", "Title", value = "IMGT/PROMISE analysis report")),
						column(4, textInput("report_author", "Author", value = APP_META$authors[[1]]$name))
					),
					fluidRow(
						column(6, selectizeInput("report_species", "Species filter", choices = NULL, multiple = TRUE,
							options = list(placeholder = "All species (leave empty)"))),
						column(6, selectizeInput("report_locus", "Locus filter", choices = NULL, multiple = TRUE,
							options = list(placeholder = "All loci (leave empty)")))
					),
					h5(strong("Sections")),
					div(style = "margin-bottom:6px; font-size:12px;",
						actionLink("report_select_all", "Select all available"), " · ",
						actionLink("report_select_none", "Clear")),
					checkboxGroupInput("report_sections", label = NULL, choices = REPORT_SECTIONS),
					checkboxInput("report_include_repro", "Include reproducibility information (versions, parameters, session)", TRUE),
					textAreaInput("report_notes", "Analyst notes (optional, added to the report)", rows = 3, width = "100%")
				)
			),
			column(5,
				box(width = 12, title = "② Export", status = "success", solidHeader = TRUE,
					radioButtons("report_format", "Format",
						choices = c("HTML (interactive, recommended)" = "html",
							"Word (.docx)" = "word", "PDF" = "pdf"), selected = "html"),
					uiOutput("report_format_warning"),
					uiOutput("report_summary"),
					br(),
					downloadButton("btn_dl_dynamic_report", "\U0001F4E5 Generate & download report",
						class = "btn-success btn-lg", style = "width:100%; height:55px;")
				),
				.res_note(strong(" Note: "),
					"the report is computed from the current dataset and the filters above. It does not depend on which ",
					"tabs you opened, but sections that need the scanner comparison require it to have been computed.")
			)
		)
	)

# ═══════════════════════════════════════════════════════════════
# SOURCE CODE
# ═══════════════════════════════════════════════════════════════
ui_tab_resources_source <- function() {
	repo <- APP_META$repo_url
	repo_ui <- if (!is.null(repo)) {
		tags$a(href = repo, target = "_blank", rel = "noopener", class = "btn btn-dark btn-lg",
			icon("github", class = "fab"), " View source code on GitHub")
	} else {
		div(class = "res-pending", icon("clock"),
			" The source repository will be made public together with the associated publication.")
	}
	clone_url <- if (!is.null(repo)) repo else "<repository-url>"
	
	tabItem(tabName = "resources_source",
		tab_header("code", "Source Code",
			"Where can the IMGT/PROMISE source code be reviewed, installed and reused?",
			"#06A77D", "#0B2A5B"),
		fluidRow(
			box(width = 12, status = "success", solidHeader = FALSE,
				div(style = "text-align:center; padding:22px;",
					icon("github", class = "fab", style = "font-size:54px; color:#24292E;"),
					h3(strong("Open, transparent and reproducible")),
					p(style = "max-width:700px; margin:0 auto 18px auto; color:#374151;",
						"IMGT/PROMISE is developed to be inspected, installed locally and cited."),
					repo_ui,
					div(style = "margin-top:16px;",
						tags$span(class = "res-badge", paste("version", APP_META$version)),
						tags$span(class = "res-badge", paste("license:", meta_or(APP_META$license, "pending"))),
						if (!is.null(APP_META$doi))
							tags$span(class = "res-badge", paste("DOI", APP_META$doi)))
				)
			)
		),
		fluidRow(
			box(width = 6, title = "⚙️ Local installation", status = "primary", solidHeader = TRUE,
				tags$pre(class = "res-code", paste0(
					"git clone ", clone_url, "\n",
					"cd IMGT-PROMISE\n\n",
					"# restore the exact package versions (renv)\n",
					"Rscript -e 'install.packages(\"renv\"); renv::restore()'\n\n",
					"# launch the application\n",
					"Rscript -e 'shiny::runApp()'")),
				.res_note(" PDF reports additionally require a LaTeX distribution (e.g. TinyTeX).")
			),
			box(width = 6, title = "\U0001F4E6 R packages", status = "info", solidHeader = TRUE,
				p("Main dependencies:"),
				div(lapply(REQUIRED_PACKAGES, function(p) tags$span(class = "res-badge res-badge-pkg", p)))
			)
		),
		fluidRow(
			box(width = 6, title = "\U0001F41B Issues & contributions", status = "warning", solidHeader = TRUE,
				p("Bug reports and suggestions are welcome through the repository issue tracker once public, ",
					"or by e-mail:"),
				tags$a(href = paste0("mailto:", APP_META$contact), APP_META$contact)),
			box(width = 6, title = "\U0001F4C4 License", status = "warning", solidHeader = TRUE,
				p(strong(meta_or(APP_META$license, "To be defined before public release.")))
			)
		)
	)
}

# ═══════════════════════════════════════════════════════════════
# REFERENCES
# ═══════════════════════════════════════════════════════════════
ui_tab_resources_references <- function()
	tabItem(tabName = "resources_references",
		tab_header("book-open", "References",
			"Which scientific sources and tools underpin the methods implemented in IMGT/PROMISE?",
			"#0B2A5B", "#7B2D8B"),
		
		fluidRow(
			box(width = 12, title = "\U0001F4CC How to cite IMGT/PROMISE", status = "danger", solidHeader = TRUE,
				div(class = "res-cite", textOutput("cite_text")),
				br(),
				tags$button(type = "button", class = "btn btn-default btn-sm",
					onclick = "var t=document.querySelector('.res-cite').innerText; navigator.clipboard.writeText(t);",
					icon("copy"), " Copy"),
				downloadButton("btn_dl_cite_bibtex", "BibTeX", class = "btn-primary btn-sm"),
				downloadButton("btn_dl_cite_ris", "RIS", class = "btn-primary btn-sm"),
				if (is.null(APP_META$article))
					div(class = "res-pending", style = "margin-top:10px;",
						icon("clock"), " The associated article is in preparation; please cite the article once published.")
			)
		),
		fluidRow(
			box(width = 12, title = "\U0001F9EC Motif detection", status = "primary", solidHeader = TRUE,
				tags$ul(class = "res-refs",
					tags$li("IMGT® scientific chart and IMGT promoter annotation procedure: reference definitions of ",
						"octamer, heptanucleotide, decamer, pentadecamer, TATA-box and E-box motifs."),
					tags$li("Bailey, T.L. & Gribskov, M. (1998). Combining evidence using p-values: application to sequence ",
						"homology searches. ", tags$em("Bioinformatics"), " 14(1), 48–54. ", .doi_link("10.1093/bioinformatics/14.1.48")),
					tags$li("Grant, C.E., Bailey, T.L. & Noble, W.S. (2011). FIMO: scanning for occurrences of a given motif. ",
						tags$em("Bioinformatics"), " 27(7), 1017–1018. ", .doi_link("10.1093/bioinformatics/btr064"))
				)
			)
		),
		fluidRow(
			box(width = 12, title = "\U0001F5C2️ IMGT", status = "info", solidHeader = TRUE,
				tags$ul(class = "res-refs",
					tags$li("Lefranc, M.-P. et al. (2009). IMGT®, the international ImMunoGeneTics information system®. ",
						tags$em("Nucleic Acids Research"), " 37, D1006–D1012. ", .doi_link("10.1093/nar/gkn838")),
					tags$li("IMGT/LIGM-DB: the IMGT comprehensive database of immunoglobulin and T cell receptor nucleotide sequences. ",
						tags$a(href = APP_META$website, target = "_blank", rel = "noopener", "imgt.org"))
				)
			)
		),
		fluidRow(
			box(width = 12, title = "\U0001F4CA Statistical methods", status = "success", solidHeader = TRUE,
				tags$ul(class = "res-refs",
					tags$li("Anderson, M.J. (2001). A new method for non-parametric multivariate analysis of variance. ",
						tags$em("Austral Ecology"), " 26(1), 32–46. ", .doi_link("10.1111/j.1442-9993.2001.01070.pp.x")),
					tags$li("Shannon, C.E. (1948). A mathematical theory of communication. ",
						tags$em("Bell System Technical Journal"), " 27(3), 379–423. ", .doi_link("10.1002/j.1538-7305.1948.tb01338.x")),
					tags$li("Simpson, E.H. (1949). Measurement of diversity. ", tags$em("Nature"), " 163, 688. ", .doi_link("10.1038/163688a0")),
					tags$li("Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate. ",
						tags$em("Journal of the Royal Statistical Society B"), " 57(1), 289–300. ", .doi_link("10.1111/j.2517-6161.1995.tb02031.x"))
				)
			)
		),
		fluidRow(
			box(width = 12, title = "\U0001F4BB Software", status = "warning", solidHeader = TRUE,
				p("Exact package versions of the running session are listed in ", strong("About → Software environment"), "."),
				div(lapply(c("R", REQUIRED_PACKAGES), function(p) tags$span(class = "res-badge res-badge-pkg", p)))
			)
		)
	)

# ═══════════════════════════════════════════════════════════════
# CHANGELOG
# ═══════════════════════════════════════════════════════════════
ui_tab_resources_changelog <- function()
	tabItem(tabName = "resources_changelog",
		tab_header("clock-rotate-left", "Changelog",
			"What changed between versions of IMGT/PROMISE?",
			"#1565C0", "#0B2A5B"),
		lapply(CHANGELOG, function(rel) {
			box(width = 12, status = "primary", solidHeader = FALSE,
				title = tagList(strong(paste0("v", rel$version)), " — ", rel$date,
					tags$span(class = "res-badge", rel$status)),
				tags$ul(lapply(rel$items, tags$li)))
		})
	)

# ═══════════════════════════════════════════════════════════════
# ABOUT
# ═══════════════════════════════════════════════════════════════
ui_tab_about_page <- function() {
	a <- APP_META$authors[[1]]
	row <- function(label, value) tags$tr(tags$td(class = "res-k", label), tags$td(value))
	
	tabItem(tabName = "about_page",
		tab_header("circle-info", "About",
			"Who develops and maintains IMGT/PROMISE, and how can they be contacted?",
			"#7B2D8B", "#0B2A5B"),
		fluidRow(
			box(width = 12, status = "primary", solidHeader = FALSE,
				div(style = "padding:6px 16px;",
					h3(strong(APP_META$name), style = "color:#0B2A5B; margin-bottom:2px;"),
					h5(em(APP_META$expansion), style = "color:#6B7280; margin-top:0;"),
					p(APP_META$tagline),
					hr(),
					fluidRow(
						column(6, tags$table(class = "table res-table",
							row("Version", paste0(APP_META$version, " (", APP_META$status, ")")),
							row("Release date", APP_META$release_date),
							row("Developer", a$name),
							row("Team / Institute", a$affiliation))),
						column(6, tags$table(class = "table res-table",
							row("Contact", tags$a(href = paste0("mailto:", APP_META$contact), APP_META$contact)),
							row("Website", tags$a(href = APP_META$website, target = "_blank", rel = "noopener", APP_META$website)),
							row("License", meta_or(APP_META$license, "To be defined at public release")),
							row("Source code", if (!is.null(APP_META$repo_url))
								tags$a(href = APP_META$repo_url, target = "_blank", rel = "noopener", APP_META$repo_url)
								else "Available at publication"),
							row("Related article", if (!is.null(APP_META$article))
								.doi_link(APP_META$article$doi) else "In preparation")))
					)
				)
			)
		),
		fluidRow(
			box(width = 6, title = "\U0001F9EA Intended use", status = "warning", solidHeader = TRUE,
				p("IMGT/PROMISE is provided for ", strong("research use only"), ". Results depend on the quality ",
					"and completeness of the imported annotations and must be interpreted by a domain expert.")),
			box(width = 6, title = "\U0001F91D Acknowledgements", status = "info", solidHeader = TRUE,
				p("Developed within the IMGT® team (IGH, CNRS – Université de Montpellier), using data from ",
					"IMGT/LIGM-DB and curated motif annotations from IMGT biocurators."))
		),
		fluidRow(
			box(width = 12, title = "\U0001F4BB Software environment", status = "success", solidHeader = TRUE,
				collapsible = TRUE, collapsed = TRUE,
				uiOutput("about_r_info"),
				DT::DTOutput("about_pkg_table"))
		)
	)
}