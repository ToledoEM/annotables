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

This package has basic annotation information from **Ensembl Genes 114**
for:

- Human build 38 (`grch38`)
- Human build 37 (`grch37`)
- Mouse (`grcm38`)
- Rat (`rnor6`)
- Chicken (`galgal5`)
- Worm (`wbcel235`)
- Fly (`bdgp6`)
- Macaque (`mmul801`)
- Pig (`Sscrofa11.1`)
- Dog (`ROS_Cfam_1.0`)
- Zebrafish (`GRCz11`)

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

    ## # A tibble: 91,673 × 9
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
    ## # ℹ 91,663 more rows

Look at the human genes-to-transcripts table:

``` r
grch38_tx2gene
```

    ## # A tibble: 412,044 × 2
    ##    enstxp          ensgene        
    ##    <chr>           <chr>          
    ##  1 ENST00000373020 ENSG00000000003
    ##  2 ENST00000612152 ENSG00000000003
    ##  3 ENST00000496771 ENSG00000000003
    ##  4 ENST00000494424 ENSG00000000003
    ##  5 ENST00000373031 ENSG00000000005
    ##  6 ENST00000485971 ENSG00000000005
    ##  7 ENST00000466152 ENSG00000000419
    ##  8 ENST00000371582 ENSG00000000419
    ##  9 ENST00000683048 ENSG00000000419
    ## 10 ENST00000371588 ENSG00000000419
    ## # ℹ 412,034 more rows

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
| ENSG00000000457 | SCYL3 | 1 | 169849631 | 169894267 | SCY1 like pseudokinase 3 |
| ENSG00000000460 | FIRRM | 1 | 169662007 | 169854080 | FIGNL1 interacting regulator of recombination and mitosis |
| ENSG00000000938 | FGR | 1 | 27612064 | 27635185 | FGR proto-oncogene, Src family tyrosine kinase |
| ENSG00000000971 | CFH | 1 | 196651754 | 196752476 | complement factor H |
| ENSG00000001460 | STPG1 | 1 | 24356999 | 24416934 | sperm tail PG-rich repeat containing 1 |
| ENSG00000001461 | NIPAL3 | 1 | 24415802 | 24476735 | NIPA like domain containing 3 |

Example with
[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
results from the
[airway](https://bioconductor.org/packages/release/data/experiment/html/airway.html)
package, made tidy with
[biobroom](http://www.bioconductor.org/packages/devel/bioc/html/biobroom.html):

``` r
library(DESeq2)
```

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     explain

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     count

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## The following object is masked from 'package:rlang':
    ## 
    ##     exprs

``` r
library(airway)

data(airway)
airway <- DESeqDataSet(airway, design = ~cell + dex)
airway <- DESeq(airway)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
res <- results(airway)

# tidy results with biobroom
library(biobroom)
```

    ## Loading required package: broom

    ## Registered S3 methods overwritten by 'biobroom':
    ##   method      from 
    ##   glance.list broom
    ##   tidy.list   broom

``` r
res_tidy <- tidy.DESeqResults(res)
```

    ## Warning: `tbl_df()` was deprecated in dplyr 1.0.0.
    ## ℹ Please use `tibble::as_tibble()` instead.
    ## ℹ The deprecated feature was likely used in the biobroom package.
    ##   Please report the issue at <https://github.com/StoreyLab/biobroom/issues>.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
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
