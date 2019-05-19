# matrixTest

Tests for block-diagonal structure in symmetric matrices (e.g. correlation matrices) under the null hypothesis of exchangeable off-diagonal elements. As described in Segal et al. (2019), these tests can be useful for construct validation either by themselves or as a complement to confirmatory factor analysis. Monte Carlo methods are used to approximate the permutation p-value with Hubert's Gamma (Hubert, 1976) and a t-statistic with unequal variance. This package also implements the Chi-squared statistic described by Steiger (1980).

## Installation

```{r}
library("devtools")
install_github("bdsegal/matrixTest")
```

## Examples

```{r}
library(matrixTest)
library(ggplot2)
library(reshape2)
     
# prepare data for matrixTest -------------------------------------------------
data("big5")

# get column numbers for questionnaire items
items <- grep("[0-9]", colnames(big5))

# compute Spearman's correlation matrix
A <- cor(big5[, items], use = "complete.obs", method = "spearman")

# Specify groups
groups <- "extrovert ~ E1 + E2 + E3 + E4 + E5 + E6 + E7 + E8 + E9 + E10
           neurotic ~ N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9 + N10
           agreeable ~ A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9 + A10
           conscientious ~ C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10
           open ~ O1 + O2 + O3 + O4 + O5 + O6 + O7 + O8 + O9 + O10"
     
# compute permutation p-values ------------------------------------------------
result <- matrixTest(A = A, groups = groups, B = 1000, absolute = TRUE)
result

# Visualize groups ------------------------------------------------------------
ord <- unlist(result$group_list)
diag(A) <- NA # remove diagonals from color scale
Am <- melt(A[ord, ord])
names(Am) <- c("x", "y", "value")
Am$y <- factor(Am$y, levels = rev(levels(Am$y)))

ggplot(aes(x = x, y = y, fill = abs(value)), data = Am)+
  geom_tile()+
  theme_bw(18)+
  scale_fill_gradient2(space="Lab", name="abs(Cor)", lim = c(0, 1))+
  labs(x = "", y = "")+
  theme(axis.text.x = element_text(angle = 90, vjust = .35,hjust=1))
```

## References

Hubert, L., Schultz, J. (1976) Quadratic assignment as a general data analysis strategy. British Journal of Mathematical and Statistical Psychology, 29(2):190–241.

Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. Psychological Bulletin, 87(2):245–251.

Segal, B. D., Braun, T., Gonzalez, R., and Elliott, M. R. (2019). Tests of matrix structure for construct validation. Psychometrika, 84(1), 65-83.
