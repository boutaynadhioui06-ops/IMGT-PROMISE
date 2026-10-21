# IMGT/PROMISE

**PROmoter Motif Identification and Sequence Explorer**

Statistical analysis of regulatory regions (5′UTR) of IG and TR V-genes.
R Shiny application developed within the IMGT® team (IGH, CNRS – Université de Montpellier).

## Features

- Dataset construction from 5′UTR and motif FASTA files
- Descriptive statistics (sequences, functionality)
- Motif conservation (frequencies, logos, logistic regression)
- Motif positioning (genomic map, distance distributions)
- Motif similarity (clustering, identity, phylogeny)
- Motif diversity, PERMANOVA and co-occurrence network
- Automatic motif discovery (consensus and PWM-based scanners)
- Dynamic analysis report (HTML / Word / PDF)

## Installation

```bash
git clone https://github.com/<ton-compte>/IMGT-PROMISE.git
cd IMGT-PROMISE
```

Then, in R:

```r
install.packages("renv")
renv::restore()
shiny::runApp()
```

PDF reports additionally require a LaTeX distribution (e.g. TinyTeX).

## Data

The application expects FASTA files exported from IMGT/LIGM-DB
(5′UTR, L-PART1 and motif files). See *Resources → Documentation*
inside the app for the expected format.

## Intended use

Research use only. Results depend on the quality and completeness of the
imported annotations and must be interpreted by a domain expert.

## Contact

Boutayna Dhioui, IMGT® Team, IGH – boutayna.dhioui@igh.cnrs.fr
https://www.imgt.org

## License

To be defined before public release.

## Author

Boutayna DHIOUI — IMGT® team.