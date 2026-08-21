annotables
================

[![DOI](https://zenodo.org/badge/3882/stephenturner/annotables.svg)](https://zenodo.org/badge/latestdoi/3882/stephenturner/annotables)

Provides tables for converting and annotating Ensembl Gene IDs.

## Installation

``` r
install.packages("devtools")
devtools::install_github("stephenturner/annotables")
```

## Rationale

Many bioinformatics tasks require converting gene identifiers from one
convention to another, or annotating gene identifiers with gene symbol,
description, position, etc. Sure,
[biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)
does this for you, but I got tired of remembering biomaRt syntax and
hammering Ensembl’s servers every time I needed to do this.

This package has basic annotation information from **Ensembl Genes 116**
for the species below. Gene counts are computed from the bundled tables;
assemblies, accessions, release dates and taxonomy IDs come from the
[Ensembl REST API](https://rest.ensembl.org).

| Dataset | Species | Taxon | Assembly | Accession | Assembly released | Genes | Coding | Non-coding | tx2gene |
|:---|:---|---:|:---|:---|:---|---:|---:|---:|---:|
| `grch38` | Human | 9606 | GRCh38 | GCA_000001405.29 | 2013-12 | 91,743 | 24,105 | 67,638 | 670,670 |
| `grch37` | Human | 9606 | GRCh37 | GCA_000001405.14 | 2009-02 | 66,978 | 25,166 | 41,812 | 215,170 |
| `grcm38` | Mouse | 10090 | GRCm39 | GCA_000001635.9 | 2020-06 | 78,718 | 22,087 | 56,631 | 481,956 |
| `rnor6` | Rat | 10116 | GRCr8 | GCA_036323735.1 | 2024-01 | 57,760 | 22,384 | 35,376 | 95,472 |
| `galgal5` | Chicken | 9031 | bGalGal1.mat.broiler.GRCg7b | GCA_016699485.1 | 2021-01 | 34,332 | 18,464 | 15,868 | 72,689 |
| `wbcel235` | Worm | 6239 | WBcel235 | GCA_000002985.3 | 2012-12 | 46,926 | 19,985 | 26,941 | 60,000 |
| `bdgp6` | Fly | 7227 | BDGP6.54 | GCA_000001215.4 | — | 28,759 | 18,484 | 10,275 | 41,600 |
| `drerio` | Zebrafish | 7955 | GRCz11 | GCA_000002035.4 | 2017-05 | 92,137 | 34,779 | 57,358 | 65,905 |
| `mmul801` | Macaque | 9544 | Mmul_10 | GCA_003339765.3 | 2019-02 | 37,169 | 22,236 | 14,933 | 64,228 |
| `cfamiliaris` | Dog | 9615 | ROS_Cfam_1.0 | GCA_014441545.1 | 2020-09 | 34,012 | 21,501 | 12,511 | 55,335 |
| `sscrofa` | Pig | 9823 | Sscrofa11.1 | GCA_000003025.6 | 2017-02 | 35,819 | 22,150 | 13,669 | 60,440 |

Coding versus non-coding is `biotype == "protein_coding"` against
everything else; totals count gene rows, not unique symbols.

**Several dataset names no longer match the assembly they contain.**
Recipes select by species, and BioMart returns whatever assembly is
current for that species, so `rnor6` now holds GRCr8, `grcm38` holds
GRCm39, and `galgal5`, `mmul801` and `cfamiliaris` likewise carry newer
assemblies than their names suggest. Use the Assembly column above, not
the dataset name, when coordinates matter.

New datasets are named after the genome assembly rather than the
species, so this cannot recur. Existing names are unchanged because
renaming them would break existing code. See `NEWS.md`.

Where each table contains:

- `ensgene`: Ensembl gene ID
- `entrez`: Entrez gene ID
- `symbol`: Gene symbol
- `chr`: Chromosome
- `start`: Start
- `end`: End
- `strand`: Strand
- `biotype`: Protein coding, pseudogene, mitochondrial tRNA, etc.
- `description`: Full gene name/description

Additionally, there are `tx2gene` tables that link Ensembl gene IDs to
Ensembl transcript IDs.

## Usage

``` r
library(annotables)
```

Look at the human genes table (note the description column gets cut off
because the table becomes too wide to print nicely):

``` r
grch38
```

    ## # A tibble: 91,743 × 9
    ##    ensgene         entrez symbol chr     start    end strand biotype description
    ##    <chr>            <int> <chr>  <chr>   <int>  <int>  <int> <chr>   <chr>      
    ##  1 ENSG00000000003   7105 TSPAN6 X      1.01e8 1.01e8     -1 protei… tetraspani…
    ##  2 ENSG00000000005  64102 TNMD   X      1.01e8 1.01e8      1 protei… tenomodulin
    ##  3 ENSG00000000419   8813 DPM1   20     5.09e7 5.10e7     -1 protei… dolichyl-p…
    ##  4 ENSG00000000457  57147 SCYL3  1      1.70e8 1.70e8     -1 protei… SCY1 like …
    ##  5 ENSG00000000460  55732 FIRRM  1      1.70e8 1.70e8      1 protei… FIGNL1 int…
    ##  6 ENSG00000000938   2268 FGR    1      2.76e7 2.76e7     -1 protei… FGR proto-…
    ##  7 ENSG00000000971   3075 CFH    1      1.97e8 1.97e8      1 protei… complement…
    ##  8 ENSG00000001036   2519 FUCA2  6      1.43e8 1.44e8     -1 protei… alpha-L-fu…
    ##  9 ENSG00000001084   2729 GCLC   6      5.35e7 5.36e7     -1 protei… glutamate-…
    ## 10 ENSG00000001167   4800 NFYA   6      4.11e7 4.11e7      1 protei… nuclear tr…
    ## # ℹ 91,733 more rows

Look at the human genes-to-transcripts table:

``` r
grch38_tx2gene
```

    ## # A tibble: 670,670 × 2
    ##    enstxp          ensgene        
    ##    <chr>           <chr>          
    ##  1 ENST00000373020 ENSG00000000003
    ##  2 ENST00000612152 ENSG00000000003
    ##  3 ENST00000496771 ENSG00000000003
    ##  4 ENST00000494424 ENSG00000000003
    ##  5 ENST00000867886 ENSG00000000003
    ##  6 ENST00000867887 ENSG00000000003
    ##  7 ENST00000867888 ENSG00000000003
    ##  8 ENST00000867889 ENSG00000000003
    ##  9 ENST00000867890 ENSG00000000003
    ## 10 ENST00000867891 ENSG00000000003
    ## # ℹ 670,660 more rows

Tables are saved in [tibble](http://tibble.tidyverse.org) format,
pipe-able with [dplyr](http://dplyr.tidyverse.org):

``` r
grch38 %>% 
    dplyr::filter(biotype == "protein_coding" & chr == "1") %>% 
    dplyr::select(ensgene, symbol, chr, start, end, description) %>% 
    head %>% 
    knitr::kable(.)
```

| ensgene | symbol | chr | start | end | description |
|:---|:---|:---|---:|---:|:---|
| ENSG00000000457 | SCYL3 | 1 | 169816039 | 169894267 | SCY1 like pseudokinase 3 |
| ENSG00000000460 | FIRRM | 1 | 169662007 | 169873258 | FIGNL1 interacting regulator of recombination and mitosis |
| ENSG00000000938 | FGR | 1 | 27611388 | 27635565 | FGR proto-oncogene, Src family tyrosine kinase |
| ENSG00000000971 | CFH | 1 | 196651754 | 196753075 | complement factor H |
| ENSG00000001460 | STPG1 | 1 | 24355830 | 24416934 | sperm tail PG-rich repeat containing 1 |
| ENSG00000001461 | NIPAL3 | 1 | 24413515 | 24476735 | NIPA like domain containing 3 |

Example with
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
results from the
[airway](https://bioconductor.org/packages/release/data/experiment/html/airway.html)
package:

``` r
library(DESeq2)
library(airway)

data(airway)
airway <- DESeqDataSet(airway, design = ~cell + dex)
airway <- DESeq(airway)
res <- results(airway)

# tidy the results
res_tidy <- res %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tibble::as_tibble() %>%
    dplyr::rename(estimate = log2FoldChange,
                  stderror = lfcSE,
                  statistic = stat,
                  p.value = pvalue,
                  p.adjusted = padj)
head(res_tidy)
```

    ## # A tibble: 6 × 7
    ##   gene            baseMean estimate stderror statistic   p.value p.adjusted
    ##   <chr>              <dbl>    <dbl>    <dbl>     <dbl>     <dbl>      <dbl>
    ## 1 ENSG00000000003  709.      0.381     0.101     3.79   0.000152    0.00128
    ## 2 ENSG00000000005    0      NA        NA        NA     NA          NA      
    ## 3 ENSG00000000419  520.     -0.207     0.112    -1.84   0.0653      0.196  
    ## 4 ENSG00000000457  237.     -0.0379    0.143    -0.264  0.792       0.911  
    ## 5 ENSG00000000460   57.9     0.0882    0.287     0.307  0.759       0.895  
    ## 6 ENSG00000000938    0.318   1.38      3.50      0.394  0.694      NA

``` r
res_tidy %>% 
    dplyr::arrange(p.adjusted) %>% 
    head(20) %>% 
    dplyr::inner_join(grch38, by = c("gene" = "ensgene")) %>% 
    dplyr::select(gene, estimate, p.adjusted, symbol, description) %>% 
    knitr::kable(.)
```

| gene | estimate | p.adjusted | symbol | description |
|:---|---:|---:|:---|:---|
| ENSG00000152583 | -4.574919 | 0 | SPARCL1 | SPARC like 1 |
| ENSG00000165995 | -3.291062 | 0 | CACNB2 | calcium voltage-gated channel auxiliary subunit beta 2 |
| ENSG00000120129 | -2.947810 | 0 | DUSP1 | dual specificity phosphatase 1 |
| ENSG00000101347 | -3.766995 | 0 | SAMHD1 | SAM and HD domain containing deoxynucleoside triphosphate triphosphohydrolase 1 |
| ENSG00000189221 | -3.353580 | 0 | MAOA | monoamine oxidase A |
| ENSG00000211445 | -3.730403 | 0 | GPX3 | glutathione peroxidase 3 |
| ENSG00000157214 | -1.976773 | 0 | STEAP2 | STEAP2 metalloreductase |
| ENSG00000162614 | -2.035665 | 0 | NEXN | nexilin F-actin binding protein |
| ENSG00000125148 | -2.210979 | 0 | MT2A | metallothionein 2A |
| ENSG00000154734 | -2.345604 | 0 | ADAMTS1 | ADAM metallopeptidase with thrombospondin type 1 motif 1 |
| ENSG00000139132 | -2.228903 | 0 | FGD4 | FYVE, RhoGEF and PH domain containing 4 |
| ENSG00000162493 | -1.891217 | 0 | PDPN | podoplanin |
| ENSG00000134243 | -2.195712 | 0 | SORT1 | sortilin 1 |
| ENSG00000179094 | -3.191750 | 0 | PER1 | period circadian regulator 1 |
| ENSG00000162692 | 3.692661 | 0 | VCAM1 | vascular cell adhesion molecule 1 |
| ENSG00000163884 | -4.459128 | 0 | KLF15 | KLF transcription factor 15 |
| ENSG00000178695 | 2.528175 | 0 | KCTD12 | potassium channel tetramerization domain containing 12 |
| ENSG00000198624 | -2.918436 | 0 | CCDC69 | coiled-coil domain containing 69 |
| ENSG00000107562 | 1.911670 | 0 | CXCL12 | C-X-C motif chemokine ligand 12 |
| ENSG00000148848 | 1.814543 | 0 | ADAM12 | ADAM metallopeptidase domain 12 |
