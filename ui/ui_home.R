# ═══════════════════════════════════════════════════════════════════════════
# ██████████████████  tabItem (home)  ███████████████████████████████████████
# ═══════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════
# ONGLET ACCUEIL 
# ═══════════════════════════════════════════════════════════════
ui_tab_home <- function()
tabItem(tabName = "home",
	
	# ══ BANNIÈRE COMPACTE ═══════════════════════════════════════
	fluidRow(
		column(12,
			div(class = "imgt-banner",
				
				# ── GAUCHE : panneau animé  ────────────────────
				tags$div(class = "banner-side",
					
					# Hélice ADN animée
					tags$div(id = "dna-helix-left", style = "position:relative;width:130px;height:130px;"),
					
					# Symboles  flottants
					tags$div(class = "gene-letters",
						tags$span(style = "animation-delay:0s;",   "IGH"),
						tags$span(style = "animation-delay:0.4s;", "IGK"),
						tags$span(style = "animation-delay:0.8s;", "IGL"),
						tags$span(style = "animation-delay:1.2s;", "TRB"),
						tags$span(style = "animation-delay:1.2s;", "TRA")
					),
					
					# Carte : années d'expertise IMGT
					tags$div(class = "side-stat-card",
						tags$div(style = "display:flex; align-items:baseline; gap:6px; justify-content:center;",
							tags$div(id = "banner_cnt1", style = "font-size:28px;font-weight:700;color:#fff;line-height:1;text-shadow:1px 2px 4px rgba(0,0,0,0.3);", "0"),
							tags$div(style = "font-size:13px;color:#fff;", "years")
						),
						tags$div(style = "font-size:9.5px;color:rgba(255,255,255,0.85);margin-top:2px;letter-spacing:0.5px;text-align:center;",
							"OF IMGT® EXPERTISE")
					),
					
					tags$div(class = "side-caption", "Genomic Motif Analysis")
				),
				
				# ── CENTRE : logo et titre
				tags$div(class = "imgt-banner-center",
					tags$img(src = "logo.jfif", height = "200px",
						style = "filter: drop-shadow(0 3px 10px rgba(0,0,0,0.45)); margin-bottom: 4px;"),
					h1("IMGT®— PROMISE"),
					p(class = "powered","Powered by IMGT® | Statistical Analysis of Immunoglobulin & T-cell Receptor Genes")
				),
				
				# ── DROITE : panneau animé STATISTIQUE ──────────────────
				tags$div(class = "banner-side",
					
					# Nuage de points réseau animé  
					tags$div(
						HTML('
									<svg class="stats-svg" viewBox="0 0 130 130" width="130" height="130">
										<line x1="8" y1="118" x2="124" y2="118" class="axis-line"/>
										<line x1="8" y1="118" x2="8" y2="10" class="axis-line"/>
										<line x1="30" y1="92" x2="62" y2="55" class="net-line" style="animation-delay:0s;"/>
										<line x1="62" y1="55" x2="96" y2="76" class="net-line" style="animation-delay:0.4s;"/>
										<line x1="62" y1="55" x2="42" y2="30" class="net-line" style="animation-delay:0.8s;"/>
										<line x1="96" y1="76" x2="112" y2="42" class="net-line" style="animation-delay:1.2s;"/>
										<circle class="stat-dot gold" cx="30" cy="92" r="4"  style="animation-delay:0s;"/>
										<circle class="stat-dot"      cx="62" cy="55" r="5"  style="animation-delay:0.3s;"/>
										<circle class="stat-dot gold" cx="96" cy="76" r="4"  style="animation-delay:0.6s;"/>
										<circle class="stat-dot"      cx="42" cy="30" r="3.5" style="animation-delay:0.9s;"/>
										<circle class="stat-dot gold" cx="112" cy="42" r="4.5" style="animation-delay:1.2s;"/>
										<circle class="stat-dot"      cx="76" cy="102" r="3"  style="animation-delay:1.5s;"/>
									</svg>
									')
					),
					
					tags$div(class = "gene-letters",
						tags$span(style = "animation-delay:0s;",   "A"),
						tags$span(style = "animation-delay:0.4s;", "T"),
						tags$span(style = "animation-delay:0.8s;", "G"),
						tags$span(style = "animation-delay:1.2s;", "C")
					),
					
					# Carte : fréquence nucléotidique 
					tags$div(class = "side-stat-card",
						tags$div(style = "font-size:9.5px;color:rgba(255,255,255,0.85);letter-spacing:0.5px;margin-bottom:6px;text-align:center;",
							"NUCLEOTIDE FREQUENCY"),
						tags$div(style = "display:flex;align-items:flex-end;gap:5px;height:30px;",
							tags$div(class = "banner-bb", `data-base` = "A", style = "flex:1;background:#fff;border-radius:2px 2px 0 0;height:10%;opacity:0.95;"),
							tags$div(class = "banner-bb", `data-base` = "T", style = "flex:1;background:#001f54;border-radius:2px 2px 0 0;height:10%;"),
							tags$div(class = "banner-bb", `data-base` = "G", style = "flex:1;background:#fff;border-radius:2px 2px 0 0;height:10%;opacity:0.65;"),
							tags$div(class = "banner-bb", `data-base` = "C", style = "flex:1;background:#001f54;border-radius:2px 2px 0 0;height:10%;opacity:0.65;")
						),
						tags$div(style = "display:flex;gap:5px;margin-top:3px;",
							tags$div(style = "flex:1;text-align:center;font-size:9px;color:rgba(255,255,255,0.75);", "A"),
							tags$div(style = "flex:1;text-align:center;font-size:9px;color:rgba(255,255,255,0.75);", "T"),
							tags$div(style = "flex:1;text-align:center;font-size:9px;color:rgba(255,255,255,0.75);", "G"),
							tags$div(style = "flex:1;text-align:center;font-size:9px;color:rgba(255,255,255,0.75);", "C")
						)
					),
					
					tags$div(class = "side-caption", "Statistical Analysis Engine")
				)
			),
			
			
			tags$script(HTML("
$(document).ready(function() {

  var helix = document.getElementById('dna-helix-left');
  var N = 13, duration = 3.4;
  for (var i = 0; i < N; i++) {
    var rung = document.createElement('div');
    rung.className = 'dna-rung';
    rung.style.top = (i * 9.5) + 'px';
    var delay = -(i * (duration / N));

    var styleTag = document.createElement('style');
    styleTag.innerHTML = '@keyframes swingBar'+i+'{0%{transform:scaleX(0.1);opacity:0.25;}50%{transform:scaleX(1);opacity:1;}100%{transform:scaleX(0.1);opacity:0.25;}}';
    document.head.appendChild(styleTag);

    var bar = document.createElement('div');
    bar.className = 'dna-bar';
    bar.style.left = '-40px'; bar.style.width = '80px';
    bar.style.animation = 'swingBar'+i+' '+duration+'s linear infinite';
    bar.style.animationDelay = delay+'s';

    var dotL = document.createElement('div');
    dotL.className = 'dna-dot l'; dotL.style.left = '-6px';
    dotL.style.animation = 'swingL '+duration+'s linear infinite';
    dotL.style.animationDelay = delay+'s';

    var dotR = document.createElement('div');
    dotR.className = 'dna-dot r'; dotR.style.left = '-6px';
    dotR.style.animation = 'swingR '+duration+'s linear infinite';
    dotR.style.animationDelay = delay+'s';

    rung.appendChild(bar); rung.appendChild(dotL); rung.appendChild(dotR);
    helix.appendChild(rung);
  }

  var start = null, target = new Date().getFullYear() - 1989, dur = 1500;
  function step(ts) {
    if (!start) start = ts;
    var p = Math.min((ts - start) / dur, 1);
    document.getElementById('banner_cnt1').textContent = Math.floor(p * target);
    if (p < 1) requestAnimationFrame(step);
    else document.getElementById('banner_cnt1').textContent = target;
  }
  requestAnimationFrame(step);

  var targets = {A: 29, T: 24, G: 27, C: 20};
  document.querySelectorAll('.banner-bb').forEach(function(el){
    var t = targets[el.getAttribute('data-base')];
    setTimeout(function(){ el.style.transition = 'height 1.2s ease-out'; el.style.height = (t*2.3)+'px'; }, 200);
    setInterval(function(){
      var jitter = t + (Math.random()*6-3);
      el.style.height = Math.max(8,jitter*2.3)+'px';
    }, 2200);
  });

});
"))
		)
	),
	
	# ══ OBJECTIF GLOBAL DE L'APPLICATION ═══════════════════════
	fluidRow(
		column(12,
			div(
				style = "background: linear-gradient(135deg, #001f54 0%, #0F4C81 100%);
                 border-radius: 12px; padding: 20px 26px; margin-bottom: 22px;
                 box-shadow: 0 4px 14px rgba(0,0,0,0.15);
                 display:flex; align-items:center; gap:18px;",
				div(style = "font-size:38px; flex-shrink:0;", "🎯"),
				div(
					h4(style = "color:#FFD700; margin:0 0 8px 0; font-weight:800;
                       letter-spacing:0.5px; font-size:15px; text-transform:uppercase;",
						"Objective"),
					p(style = "color:#F0F4FA; margin:0; font-size:15px; line-height:1.6;",
						"IMGT/PROMISE provides an integrated environment for exploring and ",
						"analyzing regulatory features of 5\u2032UTR regions in IG and TR V-genes available in IMGT/LIGM-DB, ",
						"combining statistical analyses, visualization, comparative approaches, ",
						"and automatic motif discovery for species lacking curated motif annotations.")
				)
			)
		)
	),
	
	# ══ PHASE 1 ═════════════════════════════════════════════════
	fluidRow(
		column(12,
			div(class = "phase-band",
				
				# Titre
				div(class = "phase-band-title",
					style = "background: linear-gradient(135deg, #001f54 0%, #28a745 100%);",
					div(class = "phase-icon", icon("database")),
					div(class = "phase-num", "1"),
					div(class = "phase-name", "Dataset Construction", br(), "& Descriptive Statistics")
				),
				
				# Analyses
				div(class = "phase-band-items",
					style = "background: linear-gradient(135deg, rgba(0,31,84,0.90) 0%, rgba(40,167,69,0.90) 100%);",
					div(class = "phase-item-row", div(class="dot"), "Import 5UTR files"),
					div(class = "phase-item-row", div(class="dot"), "Import motif files "),
					div(class = "phase-item-row", div(class="dot"), "Join data & Final Dataset"),
					div(class = "phase-item-row", div(class="dot"), "Descriptive table by species"),
					div(class = "phase-item-row", div(class="dot"), "Sequence & Functionality Distribution")
				),
				
				# 3 figures
				div(class = "phase-band-imgs",
					style = "background: linear-gradient(135deg, rgba(40,167,69,0.30) 0%, rgba(0,31,84,0.50) 100%);",
					div(class = "phase-fig",
						tags$img(src = "monj.PNG", alt = "Final Dataset"),
						div(class = "phase-fig-label", "Final Dataset")),
					div(class = "phase-fig",
						tags$img(src = "sousgroupe.png", alt = "Functionality"),
						div(class = "phase-fig-label", "Functionality")),
					div(class = "phase-fig",
						tags$img(src = "tableau.png", alt = "Descriptive Table"),
						div(class = "phase-fig-label", "Descriptive Table"))
				)
			)
		)
	),
	
	# ══ PHASE 2 ═════════════════════════════════════════════════
	fluidRow(
		column(12,
			div(class = "phase-band",
				
				div(class = "phase-band-title",
					style = "background: linear-gradient(135deg, #28a745 0%, #001f54 100%);",
					div(class = "phase-icon", icon("chart-bar")),
					div(class = "phase-num", "2"),
					div(class = "phase-name", "Conservation", br(), "& Diversity Motifs ")
				),
				
				div(class = "phase-band-items",
					style = "background: linear-gradient(135deg, rgba(40,167,69,0.90) 0%, rgba(0,31,84,0.90) 100%);",
					div(class = "phase-item-row", div(class="dot"), "Motif presence frequencies"),
					div(class = "phase-item-row", div(class="dot"), "Logistic Regression"),
					div(class = "phase-item-row", div(class="dot"), "Sequence logos"),
					div(class = "phase-item-row", div(class="dot"), "Clustering & Sequence Identity"),
					div(class = "phase-item-row", div(class="dot"), "Diversity indices")
				),
				
				div(class = "phase-band-imgs",
					style = "background: linear-gradient(135deg, rgba(0,31,84,0.30) 0%, rgba(40,167,69,0.50) 100%);",
					div(class = "phase-fig",
						tags$img(src = "ROC.png", alt = "ROC and AUC curves"),
						div(class = "phase-fig-label", "ROC and AUC curves")),
					div(class = "phase-fig",
						tags$img(src = "logoseq.png", alt = "Sequence logos"),
						div(class = "phase-fig-label", "Sequence logos")),
					div(class = "phase-fig",
						tags$img(src = "diversite.png", alt = "Sequence Identity"),
						div(class = "phase-fig-label", "Sequence Identity"))
				)
			)
		)
	),
	
	# ══ PHASE 3 ═════════════════════════════════════════════════
	fluidRow(
		column(12,
			div(class = "phase-band",
				
				div(class = "phase-band-title",
					style = "background: linear-gradient(135deg, #001f54 0%, #28a745 60%, #ffd700 100%);",
					div(class = "phase-icon", icon("microscope")),
					div(class = "phase-num", "3"),
					div(class = "phase-name", "Advanced", br(), "Analyses")
				),
				
				div(class = "phase-band-items",
					style = "background: linear-gradient(135deg, rgba(0,31,84,0.90) 0%, rgba(40,167,69,0.90) 60%, rgba(255,215,0,0.90) 100%);",
					div(class = "phase-item-row", div(class="dot"), "Genomic Map"),
					div(class = "phase-item-row", div(class="dot"), "Co-occurrence & Network"),
					div(class = "phase-item-row", div(class="dot"), "Motif Phylogeny"),
					div(class = "phase-item-row", div(class="dot"), "Alluvial Flow"),
					div(class = "phase-item-row", div(class="dot"), "Distance Distribution"),
					div(class = "phase-item-row", div(class="dot"), "PERMANOVA")
				),
				
				div(class = "phase-band-imgs",
					style = "background: linear-gradient(135deg, rgba(184,134,11,0.30) 0%, rgba(0,31,84,0.50) 100%);",
					div(class = "phase-fig",
						tags$img(src = "motifss.png", alt = "Genomic Map"),
						div(class = "phase-fig-label", "Genomic Map")),
					div(class = "phase-fig",
						tags$img(src = "phy.png", alt = "Phylogeny"),
						div(class = "phase-fig-label", "Phylogeny")),
					div(class = "phase-fig",
						tags$img(src = "Distances.png", alt = "Distance Distribution"),
						div(class = "phase-fig-label", "Distance Distribution"))
				)
			)
		)
	),
	
	# ══ PHASE 4 ═════════════════════════════════════════════════
	fluidRow(
		column(12,
			div(class = "phase-band",
				
				div(class = "phase-band-title",
					style = "background: linear-gradient(135deg, #ffd700 0%, #28a745 50%, #001f54 100%);",
					div(class = "phase-icon", icon("brain")),
					div(class = "phase-num", "4"),
					div(class = "phase-name", "5UTR", br(), "Motif Scanner")
				),
				
				div(class = "phase-band-items",
					style = "background: linear-gradient(135deg, rgba(255,215,0,0.90) 0%, rgba(40,167,69,0.90) 50%, rgba(0,31,84,0.90) 100%);",
					div(class = "phase-item-row", div(class="dot"), "Import 5UTR files"),
					div(class = "phase-item-row", div(class="dot"), "Import L_PART1 files"),
					div(class = "phase-item-row", div(class="dot"), "Visualisation"),
					div(class = "phase-item-row", div(class="dot"), "Join Data & compute distances"),
					div(class = "phase-item-row", div(class="dot"), "Final Dataset & Export")
				),
				
				div(class = "phase-band-imgs",
					style = "background: linear-gradient(135deg, rgba(255,215,0,0.30) 0%, rgba(40,167,69,0.40) 50%, rgba(0,31,84,0.50) 100%);",
					div(class = "phase-fig",
						tags$img(src = "motif3.PNG", alt = "Consensus-based detection"),
						div(class = "phase-fig-label", "Consensus-based detection")),
					div(class = "phase-fig",
						tags$img(src = "motif3.PNG", alt = "PWM-based probabilistic detection"),
						div(class = "phase-fig-label", "PWM-based probabilistic detection")),
					div(class = "phase-fig",
						tags$img(src = "compa.PNG", alt =  "comparison and performance"),
						div(class = "phase-fig-label", "comparison and performance"))
				)
			)
		)
	),
	
	# ══ NOTE  ════════════════════════════════════════════════
	fluidRow(
		column(12,
			div(class = "info-note",
				icon("info-circle"),
				strong(" Welcome to IMGT/PROMISE —"),
				" Navigate in the left menu to access the different analysis phases.",
				tags$span(
					style = "float: right; color: #ffd700; font-size: 11px;",
					icon("code"), " IMGT/PROMISE — M2 EDSB "
				)
			)
		)
	)
	
) # end tabItem home
