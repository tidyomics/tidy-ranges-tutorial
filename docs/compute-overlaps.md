# Compute overlaps

Objective: compute overlaps and summary statistics between two sets of
genomic ranges. In particular, suppose we want to compute the mean
genomic extent (distance from left-most to right-most basepair) of 
genes overlapping a set of query ranges.

We move on from the "classroom example" by seeing how we compute
overlaps when the features are in genomic space. We will use *GRanges*
in the Bioconductor package *GenomicRanges* to represent the features
and *plyranges* to compute the overlaps, similarly to how we used
*dplyr* to compute the overlaps in the previous analysis. So
data.frame is to *dplyr* as GRanges is to *plyranges*.


```r
library(plyranges)
```

Note the structure of the *GRanges* object. We can create a *GRanges*
from a data.frame by specifying two of: `start`, `end`, or `width`.


```r
df <- data.frame(
  seqnames="chr1",
  start=1 + c(34e6,36e6,36.6e6),
  width=c(2e5,2e5,1e5),
  strand=c("+","-","-"),
  range_id=factor(c("foo","bar","boo")))
r <- as_granges(df)
r
```

```
## GRanges object with 3 ranges and 1 metadata column:
##       seqnames            ranges strand | range_id
##          <Rle>         <IRanges>  <Rle> | <factor>
##   [1]     chr1 34000001-34200000      + |      foo
##   [2]     chr1 36000001-36200000      - |      bar
##   [3]     chr1 36600001-36700000      - |      boo
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

In case you haven't seen this before, *GRanges* objects have specific
functions to pull out information. See `?GRanges` for details.


```r
length(r)
```

```
## [1] 3
```

```r
seqnames(r)
```

```
## factor-Rle of length 3 with 1 run
##   Lengths:    3
##   Values : chr1
## Levels(1): chr1
```

```r
strand(r)
```

```
## factor-Rle of length 3 with 2 runs
##   Lengths: 1 2
##   Values : + -
## Levels(3): + - *
```

Let's find which genes overlap a region of interest. We will load the
Ensembl genes, here from a Bioconductor data package (usually we would
obtain these from AnnotationHub). The last lines use convenience
functions to convert to, e.g. `"chr1"`, ..., `"chrY"`.


```r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
```


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
g <- keepStandardChromosomes(g, pruning.mode="coarse")
```

Now we are ready to test for overlaps. A left join gives us all the
overlaps for ranges on the left side (here `r`). If a range on the
left has no overlaps it appears with `NA` for the metadata columns of
the right side ranges. If a range on the left side has multiple
overlaps with the right side, it will appear multiple times in the
output. 

Keeping reading below on how to deal with this, if it is desired to
have statistics on per-range overlaps.


```r
r %>% join_overlap_left(g)
```

```
## GRanges object with 10 ranges and 2 metadata columns:
##        seqnames            ranges strand | range_id     gene_id
##           <Rle>         <IRanges>  <Rle> | <factor> <character>
##    [1]     chr1 34000001-34200000      + |      foo      114784
##    [2]     chr1 36000001-36200000      - |      bar      127703
##    [3]     chr1 36000001-36200000      - |      bar       23154
##    [4]     chr1 36000001-36200000      - |      bar      339488
##    [5]     chr1 36000001-36200000      - |      bar        5690
##    [6]     chr1 36000001-36200000      - |      bar       63967
##    [7]     chr1 36000001-36200000      - |      bar       79932
##    [8]     chr1 36600001-36700000      - |      boo       27095
##    [9]     chr1 36600001-36700000      - |      boo       55700
##   [10]     chr1 36600001-36700000      - |      boo        9967
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

If we want to exclude the zero matches cases, we can use an inner
join:


```r
r %>% join_overlap_inner(g)
```

```
## GRanges object with 10 ranges and 2 metadata columns:
##        seqnames            ranges strand | range_id     gene_id
##           <Rle>         <IRanges>  <Rle> | <factor> <character>
##    [1]     chr1 34000001-34200000      + |      foo      114784
##    [2]     chr1 36000001-36200000      - |      bar      127703
##    [3]     chr1 36000001-36200000      - |      bar       23154
##    [4]     chr1 36000001-36200000      - |      bar      339488
##    [5]     chr1 36000001-36200000      - |      bar        5690
##    [6]     chr1 36000001-36200000      - |      bar       63967
##    [7]     chr1 36000001-36200000      - |      bar       79932
##    [8]     chr1 36600001-36700000      - |      boo       27095
##    [9]     chr1 36600001-36700000      - |      boo       55700
##   [10]     chr1 36600001-36700000      - |      boo        9967
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

**What about other types of overlap?**
Note that we can specify a distance, called a "gap", if we want to
also find when ranges are *near* each other, up to a maximum allowed
distance. We can also specify the minimum amount of overlapping
basepairs. Every overlap function in plyranges has arguments `maxgap`
and `minoverlap`. For example, to see if ranges are 50kb from each
other, we would specify `maxgap=5e4`. If we want to know if ranges are
50kb from a particular endpoint of another set of ranges, for example
TSS, we could perform the operations `anchor_5p()` followed by
`mutate(width=1)`, before overlapping the sets.

We can also perform summarization by columns either in the `r` or the
`g` object:


```r
r %>% join_overlap_inner(g) %>%
  group_by(range_id) %>%
  summarize(count=n())
```

```
## DataFrame with 3 rows and 2 columns
##   range_id     count
##   <factor> <integer>
## 1      bar         6
## 2      boo         3
## 3      foo         1
```

This is giving us the same information as the following:


```r
r %>% count_overlaps(g)
```

```
## [1] 1 6 3
```

Which can be added to the range data with a `mutate` call:


```r
r %>% mutate(overlaps = count_overlaps(., g))
```

```
## GRanges object with 3 ranges and 2 metadata columns:
##       seqnames            ranges strand | range_id  overlaps
##          <Rle>         <IRanges>  <Rle> | <factor> <integer>
##   [1]     chr1 34000001-34200000      + |      foo         1
##   [2]     chr1 36000001-36200000      - |      bar         6
##   [3]     chr1 36600001-36700000      - |      boo         3
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

If we don't care about multiple overlaps, but just want a binary 
variable that records if there was one or more overlaps or not,
we can ask if the count of overlaps is greater than 0:


```r
r %>% mutate(overlaps_any = count_overlaps(., g) > 0)
```

```
## GRanges object with 3 ranges and 2 metadata columns:
##       seqnames            ranges strand | range_id overlaps_any
##          <Rle>         <IRanges>  <Rle> | <factor>    <logical>
##   [1]     chr1 34000001-34200000      + |      foo         TRUE
##   [2]     chr1 36000001-36200000      - |      bar         TRUE
##   [3]     chr1 36600001-36700000      - |      boo         TRUE
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

If we want to keep the information about the gene ranges, we swap the
order of the ranges in the command:


```r
g %>% join_overlap_inner(r)
```

```
## GRanges object with 10 ranges and 2 metadata columns:
##          seqnames            ranges strand |     gene_id range_id
##             <Rle>         <IRanges>  <Rle> | <character> <factor>
##   114784     chr1 33979609-34631443      - |      114784      foo
##   127703     chr1 36179477-36184790      - |      127703      bar
##    23154     chr1 36023393-36032380      + |       23154      bar
##    27095     chr1 36602170-36621654      - |       27095      boo
##   339488     chr1 36038971-36060927      + |      339488      bar
##    55700     chr1 36621803-36646441      + |       55700      boo
##     5690     chr1 36035413-36107445      - |        5690      bar
##    63967     chr1 36197713-36235551      - |       63967      bar
##    79932     chr1 35899091-36023551      - |       79932      bar
##     9967     chr1 36690017-36770957      + |        9967      boo
##   -------
##   seqinfo: 25 sequences (1 circular) from hg19 genome
```

If we want strand specific overlaps, we can add `_directed`:


```r
g %>% join_overlap_inner_directed(r)
```

```
## GRanges object with 5 ranges and 2 metadata columns:
##          seqnames            ranges strand |     gene_id range_id
##             <Rle>         <IRanges>  <Rle> | <character> <factor>
##   127703     chr1 36179477-36184790      - |      127703      bar
##    27095     chr1 36602170-36621654      - |       27095      boo
##     5690     chr1 36035413-36107445      - |        5690      bar
##    63967     chr1 36197713-36235551      - |       63967      bar
##    79932     chr1 35899091-36023551      - |       79932      bar
##   -------
##   seqinfo: 25 sequences (1 circular) from hg19 genome
```

By turning the join around, we have access to the genomic range 
information about the genes. Now we can compute, e.g. the average 
genomic extent of the genes (first base to last base), per overlapping
range.


```r
g %>% join_overlap_inner_directed(r) %>%
  group_by(range_id) %>%
  summarize(count=n(), mean_width=mean(width))
```

```
## DataFrame with 2 rows and 3 columns
##   range_id     count mean_width
##   <factor> <integer>  <numeric>
## 1      bar         4    59911.8
## 2      boo         1    19485.0
```

What about `"boo"`? We need to add a `complete()` call to account for
the fact that we are missing those overlaps after the join. We need to
call the function explicitly from the *tidyr* package but by not
loading the package we can avoid some function name conflicts with
*plyranges*. Also we need to convert to tibble (explanation
follows).


```r
library(tibble)
g %>% join_overlap_inner_directed(r) %>%
  group_by(range_id) %>%
  summarize(count=n(), mean_width=mean(width)) %>%
  as_tibble() %>%
  tidyr::complete(range_id, fill=list(count=0))
```

```
## # A tibble: 3 × 3
##   range_id count mean_width
##   <fct>    <int>      <dbl>
## 1 bar          4     59912.
## 2 boo          1     19485 
## 3 foo          0        NA
```

Why did we have to convert to tibble before running `complete()`?
This is because metadata columns of *GRanges* objects are in a format
called *DataFrame* which the *tidyr* / *dplyr* functions don't know
how to operate on. To access these metadata columns, you can use any
of these types of calls:


```r
mcols(r)
```

```
## DataFrame with 3 rows and 1 column
##   range_id
##   <factor>
## 1      foo
## 2      bar
## 3      boo
```

```r
mcols(r)$range_id
```

```
## [1] foo bar boo
## Levels: bar boo foo
```

```r
r$range_id # this works also
```

```
## [1] foo bar boo
## Levels: bar boo foo
```

```r
mcols(r)[["range_id"]] # for programmatic access
```

```
## [1] foo bar boo
## Levels: bar boo foo
```

But if you want to work on them in with *tidyr* / *dplyr*, you need to
first convert to tibble (or data.frame):


```r
mcols(r) %>% as_tibble()
```

```
## # A tibble: 3 × 1
##   range_id
##   <fct>   
## 1 foo     
## 2 bar     
## 3 boo
```

**Reduce instead of summarize:**
Above when we used `group_by` and `summarize` we lost the original
range data. Another option, to preserve the range data, is to use
the function `reduce_ranges` within groups that we define. If we want
to preserve the range information for the `r` object, we can start
with `r` and proceed to join, group by, and reduce within groups. In
order to compute on the gene widths, we have to add that as a metadata
column within the join. To keep the no-gene-overlapping ranges in `r`,
we can count when the gene ID is not `NA`.


```r
r %>%
  join_overlap_left_directed(g %>% mutate(gene_width=width)) %>%
  group_by(range_id) %>%
  reduce_ranges(count=sum(!is.na(gene_id)),
                mean_width=mean(gene_width))
```

```
## GRanges object with 3 ranges and 3 metadata columns:
##       seqnames            ranges strand | range_id     count mean_width
##          <Rle>         <IRanges>  <Rle> | <factor> <integer>  <numeric>
##   [1]     chr1 34000001-34200000      * |      foo         0         NA
##   [2]     chr1 36000001-36200000      * |      bar         4    59911.8
##   [3]     chr1 36600001-36700000      * |      boo         1    19485.0
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

Hopefully, you've seen that there are many routes to compute the types
of statistics of interest. The best way to decide which one to use is
to think first: what do I want the final output to look like, and do I
need to keep track of non-overlapping ranges? This will help dictate
the way you set up the call, whether a `join` or a `mutate` to just
tally a column of overlaps, etc.
