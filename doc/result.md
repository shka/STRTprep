# Interpretation of the results

## Table of contents

1. [Quality check](#quality-check)
2. [Differential expression tests](#differential-expression-tests)
3. [For Genome Browsers](#for-genome-browsers) to understand the results

## Quality check

There are four important files about the quality check report, (i) `out/byGene/samples_all.csv`, (ii) `out/byGene/samples.xls` (or `out/byGene/samples.csv`, which is identical content but different format), (iii) `out/byGene/plugin_heatmap_global.pdf` and (iv) `out/byGene/plugin_pca_global.pdf`.

> `out/byGene/fig_heatmap_global.pdf` or `out/byGene/fig_pca.pdf` instead of `out/byGene/plugin_*_global.pdf` in case of processing by version 2-beta

The former table is the statistical report about the qualities by samples on the rows, and the the latter one is the report only about the qualified samples. You must check whether there are enough number of qualified samples for your study.

The two figures are for overview of differential expression in all the qualified samples. According to the results, you can find (i) major variation in your samples, (ii) library bias, (iii) sample batch bias, (iv) outlier samples in the expression profile, due to sampling/experimental errors in the sample preparation steps (rather than technical errors during STRT library preparation/sequencing), and so on. You can design the studies again by re-editing of `src/samples.csv` based on this figure.

### `out/byGene/samples_all.csv`

This table is a quality check report for all samples. The following is the table legend.

- Columns `LIBRARY`, `WELL`, and `NAME`: Library name, well name and sample name, which were given by `src/samples.csv`.
- `QUALIFIED_READS`: Read count just before the redundant-read exclusion
- `TOTAL_READS`: Read count just after the redundant-read exclusion
- `REDUNDANCY`:	= QUALIFIED_READ / TOTAL_READS
- `MAPPED_READS`: Read count, which are aligned at only one location
- `MAPPED_RATE`: = MAPPED_READS / TOTAL_READS
- `SPIKEIN_READS`: Read count, which are aligned on spike-in sense strand
- `MAPPED/SPIKEIN`: ~ Relative abundance of endogenous polyA+ transcripts
- `SPIKEIN_5END_READS`: Read count, which are aligned within 5'-end 50 nt region of spike-in sense strand
- `SPIKEIN_5END_RATE`: = SPIKEIN_5END_READS / SPIKEIN_READS
- `CODING_READS`: Read count, which are aligned within any exon or the proximal upstream (500bp) of coding gene sense strand
- `CODING_5END_READS`: Read count, which are aligned within 5'-UTR or the proximal upstream of coding gene sense strand
- `CODING_5END_RATE` = CODING_5END_READS / CODING_READS
- `SPIKEIN_READS.OUTLIER`: TRUE if it is outlier of SPIKEIN_READS in the library
- `MAPPED/SPIKEIN.OUTLIER`: TRUE if it is outlier of MAPPED/SPIKEIN in the library
- `SPIKEIN_5END_RATE.OUTLIER`: TRUE if it is outlier of SPIKEIN_5END_RATE in the library
- `CODING_5END_RATE.OUTLIER`: TRUE if it is outlier of CODING_5END_RATE in the library
- ... and the other columns given by `src/samples.csv`; see a [section in "Preprocessing and analysis protocol"](protocol.md#design-of-experiments).

The outlier wells are not used for further statistics, for example differential expression tests, because it suspects any experimental errors. You can find the distribution and the outliers by figures `out/byGene/fig_outliers_*.pdf`.

### `out/byGene/samples.xls`

This is also quality check report, but it contains only the qualified samples in `out/byGene/samples_all.xls` - all samples must be that (i) `NAME` is not `NA`, and (ii) all '*.OUTLIER' values are FALSE.

### `out/byGene/plugin_heatmap_global.pdf`

This is heatmap and clustering of fluctuated gene expression in the qualified samples. The pipeline selects fluctuated genes automatically by significance of coefficient of variance between the qualified samples; the significance is estimated by comparison with the expected technical variation in the spike-in RNAs.

> The figure file name might be `out/byGene/fig_heatmap_global.pdf` if processed by version 2-beta.

### `out/byGene/plugin_pca_global.pdf`

This is PCA plot in the qualified samples with the fluctuated gene expression.

> The figure file name might be `out/byGene/fig_pca.pdf` if processed by version 2-beta.

## Differential expression tests

The table `out/byGene/diffexp.xls` (or `out/byGene/diffexp.csv`, which is identical content but different format) contains statistical results by protein coding genes on the rows, and the columns are as below. There are three important purposes for this table.

1. You can find the differentially regulated genes from this table. Reasonable thresholding is by both of `qvalue` and `fluctuation`, but sometimes `pvalue` and `fluctuation` only work when you prepared a few biological replicates.
2. You can find the normalized expression levels in your genes of interest, and draw the figure by any your favorite tool, for example R or Microsoft Excel.
3. You can apply the other statistics using the other third-party tools with the raw read counts or the normalized expression levels (and we are grateful if you could provide the analysis code!).

Also, there are several heatmaps `out/byGene/plugin_heatmap_*.pdf` automatically (where `*` is the test number given in `src/samples.csv`) if there are significantly regulated genes. The pipeline selects significantly regulated genes automatically by both of `qvalue` and `fluctuation` thresholds given in `conf.yaml`.

> The figure file name might be `out/byGene/fig_heatmap_diffexp*.pdf` if processed by version 2-beta.

### `out/byGene/diffexp.xls`

This is very big table containing expression levels and statistical results of the differential expression in the qualified samples.

- Column `Gene`: Gene name
- `fluctuationScore.global`: Statistical score in degree of variation of the normalized expression levels between the qualified wells towards expected technical noise levels
- `fluctuation.global`: Adjusted p-value of the degree of variation between the qualified wells
- `diffexpScore.n` (n=0, 1, ...): Statistical score of differential expression in comparison n; n is the test number given in `src/samples.csv`. Positive score in two class comparison is up-regulation of the class 2 samples than the class 1, and negative is down-regulation. In case of multiclass comparison, the score is all positive. Larger absolute score is more significant.
- `pvalue.n` (n=0, 1, ...): P-value of differential expression in comparison n
- `qvalue.n` (n=0, 1, ...): Q-value (a.k.a. FDR) of differential expression in comparison n
- `fluctuationScore.n` (n=0, 1, ...): Statistical score in degree of variation of the normalized expression levels between targets of the comparison n towards expected technical noise levels
- `fluctuation.n` (n=0, 1, ...): Adjusted p-value of the degree of variation between the targets of the comparison n
- `N|*library.well*|*name*`: Spike-in based normalized expression level
- `R|*library.well*|*name*`: Raw aligned read count within 5'-UTR of the gene, or the proximal (500 nt) upstream; such transcripts starting these regions would be intact mRNAs as template for translation into proteins. You can find the regions [[by Genome Browsers|Results#for-genome-browsers]] with file `out/web/regions_byGene.bed.gz`.

> No `fluctuationScore.*` columns, and `Score.n` instead of `diffexpScore.n`, if processed by version 2-beta.

> `diffexpScore.n`, `pvalue.n` and `qvalue.n` are calculated by SAMstrt [[Katayama et al. 2013](http://www.ncbi.nlm.nih.gov/pubmed/?term=23995393)]; the `diffexpScore.n` by STRTprep is `Score(d)` by SAMstrt. Briefly, in case of two class comparison, the score is average of Wilcoxon statistics with multiple Poisson resampling. In contrast, the score in multiple class comparison is average of Kruskal-Wallis statistic. See also SAMseq [[Li and Tibshirani 2013](http://www.ncbi.nlm.nih.gov/pubmed?term=22127579)] and the manual in the [official web page](http://statweb.stanford.edu/~tibs/SAM/) to understand the statistical background and the outputs.

> `fluctuationScore.*` and `fluctuation.*` are calculated by this pipeline [Krjutškov and Katayama et al., submitted].

> [AutoFilter](https://support.office.com/en-ca/article/Quick-start-Filter-data-by-using-an-AutoFilter-08647e19-11d1-42f6-b376-27b932e186e0?ui=en-US&rs=en-CA&ad=CA) by Microsoft Excel is very useful to select the significantly fluctuated genes. After the filtering, copy the gene names, then paste them into your favorite analysis tools!

## For Genome Browsers

Files in a folder `out/web` are appendix for understanding of the results on Genome Browsers; for example, you can [upload to UCSC Genome Browser](http://genome-euro.ucsc.edu/FAQ/FAQcustom.html#custom1).

### `out/web/regions_byGene.bed.gz`

The bands in this track are 5'-UTR of protein coding genes, or the proximal (500 nt) upstream. Aligned reads within the bands were used for the `byGene` analysis.
