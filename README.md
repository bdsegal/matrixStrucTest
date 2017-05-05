# matrixTest

This is an R package that tests for block-diagonal structure in symmetric matrices. The null hypothesis is that the off-diagonal elements are exchangeable. Monte Carlo methods are used to approximate the permutation p-value with Hubert's Gamma (Hubert, 1976) and a t-statistic with unequal variance. This package also implements a Chi-squared statistic described by Steiger (1980)

Please see Segal, et. al for more information.

## Installation

```{r}
library("devtools")
install_github("bdsegal/matrixTest")
```

## Examples

```{r}
library(matrixTest)

# prepare data for matrixTest -------------------------------------------------
data("big5")

# get column numbers for questionnaire items
items <- grep("[0-9]", colnames(big5))

# compute Spearman's correlation matrix
A <- cor(big5[, items], use = "complete.obs", method = "spearman")

# get column numbers for items within each block/group
extrovert <- grep("E", colnames(A))
neurotic <- grep("N", colnames(A))
agreeable <- grep("A", colnames(A))
conscientious <- grep("C", colnames(A))
open <- grep("O", colnames(A))

# put blocks/groups in list for matrixTest function
group_list <- list(extrovert = extrovert, 
                   neurotic = neurotic, 
                   agreeable = agreeable, 
                   conscientious = conscientious,
                   open = open)

# compute permutation p-values ------------------------------------------------
matrixTest(A = A, group_list = group_list, B = 1000, absolute = TRUE)

```

## References

Lawrence Hubert and James Schultz. Quadratic assignment as a general data analysis
strategy. British journal of mathematical and statistical psychology, 29(2):190–241, 1976.

James H. Steiger. Tests for comparing elements of a correlation matrix. Psychological
bulletin, 87(2):245–251, 1980.

Brian D. Segal, Thomas Braun, Richard Gonzalez, and Michael Elliott. Tests of matrix structure for construct validation. Submitted.
<!-- Working draft available online at [https://arxiv.org/abs/[id]](https://arxiv.org/abs/[id]). -->