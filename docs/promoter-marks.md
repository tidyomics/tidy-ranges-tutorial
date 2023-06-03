# Promoter marks

Objective: determine if tissue-specific promoter marks (e.g. H3K27ac)
are often near genes that are expressed in a tissue-specific manner.

We will load expression data from the GTEx project [@gtex], which gives
median expression in TPM for each tissue. We will use H3K27ac ChIP-seq
data from the ENCODE project [@encode].


```r
library(tidyr)
file <- "data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
gtex <- read.delim(file, skip=2)
```

We select two tissues, bladder and kidney, and convert the data from a
wide format into a tidy format.


```r
tissues <- gtex %>% dplyr::select(Name, Bladder, Kidney...Cortex) %>%
  dplyr::rename(gene = Name, Kidney = Kidney...Cortex) %>%
  dplyr::mutate(gene = sub("\\..*","",gene)) %>%
  pivot_longer(!gene, names_to="tissue", values_to="tpm")
```

Now define two vectors of genes that are specific to bladder and kidney:


```r
bladder_expr <- tissues %>%
  dplyr::filter(tissue == "Bladder" & tpm > 10) %>%
  dplyr::pull(gene)
kidney_expr <- tissues %>%
  dplyr::filter(tissue == "Kidney" & tpm > 10) %>%
  dplyr::pull(gene)
int <- intersect(bladder_expr, kidney_expr)
bladder_expr <- setdiff(bladder_expr, int)
kidney_expr <- setdiff(kidney_expr, int)
# save(bladder_expr, kidney_expr, file="data/bladder_kidney_expr.rda")
```

Next, use an existing TxDb to locate these genes in the genomes. While
we usually recommend to use GENCODE genes for human analysis, because
the ENCODE chromatin modification peak files on AnnotationHub are in
hg19, we use the UCSC hg19 genes here for simplicity of the code:


```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
```

Add the ENSEMBL ID and pull out the two tissue-specific sets. 


```r
g <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
```

```
##   403 genes were dropped because they have exons located on both strands of the same
##   reference sequence or on more than one reference sequence, so cannot be represented by a
##   single genomic range.
##   Use 'single.strand.genes.only=FALSE' to get all the genes in a GRangesList object, or use
##   suppressMessages() to suppress this message.
```


```r
library(plyranges)
g <- g %>% mutate(ensembl = mapIds(org.Hs.eg.db, gene_id, "ENSEMBL", "ENTREZID"))
bladder_g <- g %>% filter(ensembl %in% bladder_expr)
kidney_g <- g %>% filter(ensembl %in% kidney_expr)
```

Finally we combine the two sets with `bind_ranges`, and we change
the feature size from the whole gene extent (the range from the
leftmost exon to rightmost exon), to just the TSS, using `anchor_5p`
and `mutate`.


```r
tss <- bind_ranges(bladder=bladder_g,
                   kidney=kidney_g,
                   .id="gtissue") %>%
  anchor_5p() %>%
  mutate(width=1)
```

Now we will obtain the H3K27ac peak sets:


```r
library(AnnotationHub)
ah <- AnnotationHub()
# query(ah, c("Homo sapiens", "bladder", "H3K27ac", "narrowPeak"))
bladder_pks <- ah[["AH44180"]]
# query(ah, c("Homo sapiens", "kidney", "H3K27ac", "narrowPeak"))
kidney_pks <- ah[["AH43443"]]
save(bladder_pks, kidney_pks, file="data/peaks.rda")
```



We download these and scale so they have the same 90% quantile of `signalValue`.


```r
ninety <- function(x) quantile(x, .9, names=FALSE)
bladder_pks <- bladder_pks %>%
  mutate(signal = signalValue / ninety(signalValue))
kidney_pks <- kidney_pks %>%
  mutate(signal = signalValue / ninety(signalValue))
```

Combine the peaks from bladder and kidney, filter to those with < 0.1%
FDR, and center the peak on the summit (the `peak` column gives the
shift from the left side to the summit).


```r
pks <- bind_ranges(bladder=bladder_pks,
                   kidney=kidney_pks,
                   .id="ptissue") %>%
  filter(qValue > 3, width <= 1000) %>%
  mutate(start = start + peak) %>%
  select(-peak) %>%
  mutate(width = 1)
```

Finally, once we have two tidy range sets, we can perform the analysis
by a join, followed by two lines that take care of multiple overlaps,
followed by two lines that give us our tallies of interest.

It appears that tissue-specific peaks are enriched near the
tissue-specific genes for both bladder and kidney. 


```r
tss %>%
  join_overlap_left(pks, maxgap=500) %>%
  group_by(ptissue) %>% # within peak tissue...
  filter(!duplicated(gene_id)) %>% # ...just take the first overlap per gene
  group_by(gtissue, ptissue) %>%
  summarize(count = n())
```

```
## DataFrame with 6 rows and 3 columns
##   gtissue     ptissue     count
##     <Rle> <character> <integer>
## 1 bladder     bladder      2163
## 2 bladder      kidney      1568
## 3 bladder          NA       859
## 4  kidney     bladder       103
## 5  kidney      kidney       208
## 6  kidney          NA       242
```

The above number could also be found with four `countOverlaps` calls,
by considering all four pairs of overlaps of the two sets of genes and
peaks. 

Another way to avoid counting overlaps more than once per gene is to
use the *plyranges* function, `n_distinct()`:


```r
tss %>%
  join_overlap_left(pks, maxgap=500) %>%
  group_by(gtissue, ptissue) %>%
  summarize(count = n_distinct(gene_id))
```

```
## DataFrame with 6 rows and 3 columns
##   gtissue     ptissue     count
##     <Rle> <character> <integer>
## 1 bladder     bladder      2163
## 2 bladder      kidney      1568
## 3 bladder          NA       859
## 4  kidney     bladder       103
## 5  kidney      kidney       208
## 6  kidney          NA       242
```

If we want more information per gene, e.g. suppose we want to compute
the average signal per gene of peaks nearby, we need to group twice,
once also by gene ID, and the second time integrating over gene
ID. While here we add a few more lines of code, performing such an
operation with base Bioconductor functions would require adding code
to perform the loops, adding many intermediate variables to store
results, etc.


```r
tss %>%
  join_overlap_left(pks, maxgap=500) %>%
  group_by(gtissue, ptissue, gene_id) %>% # need per gene stats
  summarize(num_overlaps = n(), signal = mean(signal)) %>%
  as_tibble() %>% # DataFrame to tibble for further processing
  group_by(gtissue, ptissue) %>%
  summarize(sum_any_overlaps = sum(num_overlaps > 0),
            mean_signal=mean(signal))
```

```
## `summarise()` has grouped output by 'gtissue'. You can override using the `.groups` argument.
```

```
## # A tibble: 6 × 4
## # Groups:   gtissue [2]
##   gtissue ptissue sum_any_overlaps mean_signal
##   <chr>   <chr>              <int>       <dbl>
## 1 bladder bladder             2163       1.18 
## 2 bladder kidney              1568       0.929
## 3 bladder <NA>                 859      NA    
## 4 kidney  bladder              103       0.991
## 5 kidney  kidney               208       0.714
## 6 kidney  <NA>                 242      NA
```

What's wrong with this analysis?

* We didn't figure out the expressed promoter, we just looked at the
  left or rightmost isoform (for + or - strand genes, respectively). 
