# CFunStats: Compact Functional Statistics

Statistics to preserve continuity, capture discontinuity, and identify locality of discrete functional patterns. The input is discrete data of two variables. A corresponding contingency table is used to encode the discrete pattern. The output is functional statistics on a compact table representing the original table. Rows of the compact table are the trimmed and compressed version of the input table while columns are the same with the input table. The compact table is obtained by dynamic programming. Let the row variable be the independent variable and the column variable be the dependent variable. We characterize the statistical significance and the strength (effect size) of the pattern being a function by the model-free functional chi-squared test. The compact functional statistics can be applied to time-course signals to recognize both temporal continuity and discontinuity while reducing model biases.

# Installation

`devtools` is required to install CFunStats:

```r
library(devtools)
install_github("https://github.com/sajal49/CFunStats")
```


# Usage

The function `c_fun_stats` in CFunStats takes two discrete variables as input, whether or not trimming to the most active window is desired and whether log of p-value should be returned. It return functional statistics for the compacted joint pattern between the two discrete variables.

```r
library(CFunStats)
x = c(1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,3, 6,6,6,6,6, 1,1,1,1,1, 2,2,2,2,2, 4,4,4,4,4, 6,6,6,6,6, 1,1,1,1,1, 2,2,2,2,2, 5,5,5,5,5, 6,6,6,6,6)
y = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)

# statistics on trimmed and compressed table
trimmed_stats = c_fun_stats(x, y)
print(trimmed_stats)

# statistics on full table with compression but no trimming
full_stats = c_fun_stats(x, y, full=TRUE)
print(full_stats)
```

To see the compressed / trimmed table `construct_compact_table` can be used.

```r
library(CFunStats)
x = c(1,1,1,1,1, 2,2,2,2,2, 3,3,3,3,3, 6,6,6,6,6, 1,1,1,1,1, 2,2,2,2,2, 4,4,4,4,4, 6,6,6,6,6, 1,1,1,1,1, 2,2,2,2,2, 5,5,5,5,5, 6,6,6,6,6)
y = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3)

# orignal table
print(table(x, y))

# trimmed and compressed table
trimmed_table = construct_compact_table(x, y)
print(trimmed_table)

# only compressed table
full_table = construct_compact_table(x, y, full=TRUE)
print(full_table)
```


