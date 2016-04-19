# Utilities in this package

There were several scripts of utility for the analysis.

- [`utils/intersect_simple_diffexps.R`](#utils_intersect_simple_diffexps_R)

## `utils/intersect_simple_diffexps.R`

This script outputs intersection of differentially expressed genes (or TFEs) between two comparisons. Thresholds for the two comparisons must be configured by parameters of [`simple_diffexp`](plugin.md#plugin-simple_diffexp) plugin. The output contains five statistic values (`diffexpScore.n`, `pvalue.n`, `qvalue.n`, `fluctuationScore.n` and `fluctuation.n`; see also [`out/byGene/diffexp.xls`](result.md#outbygenediffexpxls)).

### Usage

`utils/intersect_simple_diffexps.R` _quantificationType_ _class1_ _class2_


Options|Value
-------|-----
_quantificationType_ | Either `byGene` or `byTFE`
_class1_ | Test ID 1
_class2_ | Test ID 2

### Example

To get intersection of differentially expressed genes in the tests 1 and 2,

```
utils/intersect_simple_diffexps.R byGene 1 2 > out/160419-intersect_simple_diffexps_1vs2.csv
```
