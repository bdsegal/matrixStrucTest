# MC estimates of permutation p-values with Hubert's Gamma and
# t-statistic with unequal variance

# functions for making Delta matrix -------------------------------------------
deltaSub <- function(i,j, group_list){
  #' Sub-routine to create Delta matrix
  #'
  #' This sub-routine outputs 1 if i and j are in at least one group together, 
  #' and 0 otherwise, and is called by \code{matrixStrucTest} and 
  #' \code{prepBoxPlots}.
  #' @param i First index
  #' @param j Second index
  #' @param group_list List of indices for each block

  i_group <- which(sapply(group_list, function(group){i %in% group}))
  j_group <- which(sapply(group_list, function(group){j %in% group}))
  
  # allows for multiple membership
  as.numeric(length(intersect(i_group, j_group)) >= 1)
}

multiSub <- function(i, j, group){
  #' Sub-routine to create Delta matrix for block-specific tests
  #'
  #' This sub-routine outputs TRUE if either \code{i} or \code{j} are in 
  #' \code{group}, FALSE otherwise, and is called by \code{matrixStrucTest} 
  #' and \code{prepBoxPlots}.
  #' @param i First index
  #' @param j Second index
  #' @param group Indices for items in group

  i_in <- i %in% group
  j_in <- j %in% group
  
  (i_in + j_in) >= 1
}

makeGroupList <- function(groups, A) {
  #' Convert character string in lavaan syntax into a list of indices
  #'
  #' This sub-routine is called by \code{matrixStrucTest} and \code{prepBoxPlots}.
  #' @param groups Character string in lavaan syntax specifying groups
  #' @param A A Distance or similarity matrix. Must have column names
  #'
  #' @return group_list List of column indices of A corresponding to each group
  #' @export

  split <- gsub(" ", "", unlist(strsplit(groups, "\\n")))
  group_list <- vector(length = length(split), mode = "list")
  for (i in 1:length(split)) {
    sub_split <- unlist(strsplit(split[i], "~"))
    names(group_list)[i] <- sub_split[1]

    labels <- unlist(strsplit(sub_split[2], "\\+"))
    group_list[[i]] <- which(colnames(A) %in% labels)
  }

  rm_group <- which(is.na(names(group_list)))
  if (length(rm_group) >= 1) {
    group_list <- group_list[-rm_group]
  }

  return(group_list)
}

# t-statistic and Hubert's gamma ----------------------------------------------
matrixStrucTestSub <- function(A, group_list_ord, Delta, multi_group_ind,
                         A_upper_ind, K){
  #' Compute Gamma and t-statistics for a single permutation
  #'
  #' This sub-routine is called by \code{matrixStrucTest} and \code{prepBoxPlots}.
  #' @param A Distance or similarity matrix, e.g. correlation
  #' @param group_list_ord List of groupings for ordered matrix A
  #' @param Delta Delta matrix
  #' @param multi_group_ind List of indicator matrices for membership in block k test
  #' @param A_upper_ind indicator matrix for upper triangular elements
  #' @param K Total number of hypothesized blocks
  #' @return Gamma_overall: Overall Hubert's Gamma
  #' @return Gamma_multi: Block-specific Hubert's Gamma
  #' @return t_overall: Overall t-statistic (unequal variance)
  #' @return t_multi: Block-specific t-statistics
  #' @return Ak_list: List of values in A used for each block-specific test
  #' @return Deltak_list: List of values in Delta used for each block-specific test
  #' @return A_upper: Upper triangular elements of A (used for box plot function)
  #' @return Delta_upper: Upper triangular elements of Delta (used for box plot function)

  # Calculate Hubert's Gamma with upper triangular elements
  A_upper <- A[upper.tri(A)]
  Delta_upper <- Delta[upper.tri(Delta)]

  Gamma_overall <- cor(A_upper, Delta_upper)

  g1 <- A_upper[which(Delta_upper == 1)]
  g0 <- A_upper[which(Delta_upper == 0)]
  t_overall <- t.test(g1, g0, equal.var = FALSE)$statistic

  # multiple hypothesis testing ---------------------------
  Ak_list <- list()
  Deltak_list <- list()

  Gamma_multi <- rep(NA, K)
  t_multi <- rep(NA, K)
  for (k in 1:K){

   Ak_list[[k]] <- A[A_upper_ind & multi_group_ind[[k]]]
   Deltak_list[[k]] <- Delta[A_upper_ind & multi_group_ind[[k]]]
   Gamma_multi[k] <- cor(Ak_list[[k]], Deltak_list[[k]])

   g0 <- Ak_list[[k]][which(Deltak_list[[k]] == 0)]
   g1 <- Ak_list[[k]][which(Deltak_list[[k]] == 1)]
   t_multi[k] <- t.test(g1, g0, equal.var = FALSE)$statistic

  }

  return(list(Gamma_overall = Gamma_overall,
              Gamma_multi = Gamma_multi,
              t_overall = t_overall,
              t_multi = t_multi,
              Ak_list = Ak_list,
              Deltak_list = Deltak_list,
              A_upper = A_upper,
              Delta_upper = Delta_upper))
}


matrixStrucTest <- function(A, group_list = NULL, groups = NULL, B = 1000, absolute = TRUE){
  #' Permutation p-values for Gamma and t-statistics
  #'
  #' This function computes permutation p-values for Hubert's Gamma and t-statistics
  #' for both overall and block-specific tests.
  #' @param A Distance or similarity matrix, e.g. correlation
  #' @param groups CFA model in lavaan syntax. Either \code{groups} or \code{group_list} but not both must be supplied.
  #' @param group_list List of column indices of A for each group. Either \code{groups} or \code{group_list} but not both must be supplied.
  #' @param B Number of Monte Carlo resamples (defaults to B=1000)
  #' @param absolute Use the absolute values of A (defaults to TRUE)

  #' @return pt_overall_one_sided: Overall one-sided p-value using t statistic
  #' @return pt_overall_two_sided: Overall two-sided p-value using t statistic
  #' @return pt_multi_one_sided: Block-specific one-sided p-values using t statistic
  #' @return pt_multi_two_sided: Block-specific two-sided p-values using t statistic
  #' @return t0 Observed overall: t statistic
  #' @return t0k: Observed block-specific t statistic
  #' @return t_overall: Vector of overall t statistics from permuted A
  #' @return t_max_one_sided: Vector of max t statistics from permuted A (one-sided)
  #' @return t_max_two_sided: Vector of max t statistics from permuted A (two-sided)
  #' @return pG_overall_one_sided: Overall one-sided p-value using Hubert's Gamma
  #' @return pG_overall_two_sided: Overall two-sided p-value using Hubert's Gamma
  #' @return pG_multi_one_sided: Block-specific one-sided p-values using Hubert's Gamma
  #' @return pG_multi_two_sided: Block-specific two-sided p-values using Hubert's Gamma
  #' @return Gamma0: Observed overall Hubert's Gamma
  #' @return Gamma0k: Observed block-specific Hubert's Gamma
  #' @return Gamma_overall: Vector of Hubert's Gamma statistics from permuted A
  #' @return Gamma_max_one_sided: Vector of max Hubert's Gamma statistics from permuted A (one-sided)
  #' @return Gamma_max_two_sided: Vector of max Hubert's Gamma statistics from permuted A (two-sided)
  #' @return B: number of Monte Carlo resamples
  #' @return group_list: List of column/row indices corresponding to each group
  #' @keywords Hubert's Gamma, t-statistic, permutation, MC resampling
  #' @export
  #' @examples
  #' # example for matrixStrucTest package
  #' library(matrixStrucTest)
  #' data("big5")
  #'
  #' # get column numbers for questionnaire items
  #' items <- grep("[0-9]", colnames(big5))
  #' 
  #' # compute Spearman's correlation matrix
  #' A <- cor(big5[, items], use = "complete.obs", method = "spearman")
  #' 
  #' # specify the groups
  #' groups <- "extrovert ~ E1 + E2 + E3 + E4 + E5 + E6 + E7 + E8 + E9 + E10
  #'            neurotic ~ N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9 + N10
  #'            agreeable ~ A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9 + A10
  #'            conscientious ~ C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10
  #'            open ~ O1 + O2 + O3 + O4 + O5 + O6 + O7 + O8 + O9 + O10"
  #'
  #' # compute permutation p-values
  #' # Note: Using small B for fast checking on CRAN. Set B >= 1000 in practice.
  #' result <- matrixStrucTest(A = A, groups = groups, B = 100, absolute = TRUE)
  #'
  #' # Note: two-sided p-values from Hubert's Gamma printed by default
  #' #       other results available by directing accessing them from the
  #' #       returned object
  #' result
  #'
  #' # Alternative approach for specifying the groups as a list of column/row indices
  #' extrovert <- grep("E", colnames(A))
  #' neurotic <- grep("N", colnames(A))
  #' agreeable <- grep("A", colnames(A))
  #' conscientious <- grep("C", colnames(A))
  #' open <- grep("O", colnames(A))
  #' 
  #' # put blocks/groups in list
  #' group_list <- list(extrovert = extrovert, 
  #'                    neurotic = neurotic, 
  #'                    agreeable = agreeable, 
  #'                    conscientious = conscientious,
  #'                    open = open)
  #' 
  #' # Note: Using small B for fast checking on CRAN. Set B >= 1000 in practice.
  #' result <- matrixStrucTest(A = A, group_list = group_list, B = 100, absolute = TRUE)
  #'
  #' # Note: two-sided p-values from Hubert's Gamma printed by default
  #' #       other results available by directing accessing them from the
  #' #       returned object
  #' result
  #' 
  #' # Visualize groups
  #' library(ggplot2)
  #' library(reshape2)
  #' 
  #' ord <- unlist(result$group_list)
  #' diag(A) <- NA # remove diagonals from color scale
  #' Am <- melt(A[ord, ord])
  #' names(Am) <- c("x", "y", "value")
  #' Am$y <- factor(Am$y, levels = rev(levels(Am$y)))
  #' 
  #' ggplot(aes(x = x, y = y, fill = abs(value)), data = Am)+
  #'   geom_tile()+
  #'   theme_bw(18)+
  #'   scale_fill_gradient2(space="Lab", name="abs(Cor)", lim = c(0, 1))+
  #'   labs(x = "", y = "")+
  #'   theme(axis.text.x = element_text(angle = 90, vjust = .35,hjust=1))

  if(!is.null(groups) & !is.null(group_list)) {
    stop("Both a 'groups' and 'group_list' were supplied as arguments. Only one of the two can be given.")
  } else if(is.null(groups) & is.null(group_list)) {
    stop("Either 'groups' or 'group_list' must be supplied as an argument")
  }

  # if groups given, convert to group_list
  if(!is.null(groups)) {
    if(is.null(colnames(A))) {
      stop("When 'groups' is supplied, the matrix A must have column names")
    }

  group_list <- makeGroupList(groups, A)
  }

  if(absolute){
    A <- abs(A)
  }

  p <- ncol(A)
  K <- length(group_list)

  if(sum(sapply(group_list, length)) > p){
    warning("Some variables may be assigned to multiple clusters")
  } else if(sum(sapply(group_list, length)) < p){
    warning("Some variables not assigned to a cluster")
  }

  # clean up names of groups
  if(is.null(names(group_list))) {
    names(group_list) <- paste("block", 1:K)
  } else {
    no_name <- which(sapply(names(group_list), nchar) == 0)
    names(group_list)[no_name] <- paste("block", no_name)
  }

  # create delta matrix
  Delta <- matrix(nrow=p,ncol=p)
  for (i in 1:p){
    for (j in 1:p){
       Delta[i,j] <- deltaSub(i, j, group_list)
    }
  }

  # order matrices by hypothesized groups
  ord <- unlist(group_list)
  A_ord <- A[ord, ord]
  Delta_ord <- Delta[ord, ord]

  # get group list for reordered matrices
  nk <- sapply(group_list, length)
  firstIndex <- 1
  group_list_ord <- list()
  for (k in 1:K) {
    group_list_ord[[k]] <- firstIndex:(firstIndex+nk[k] -1)
    firstIndex <- firstIndex + nk[k]
  }

  # get indicator matrices for block-specific tests
  newOrd <- unlist(group_list_ord)
  A_upper_ind <- upper.tri(A)
  multi_group_ind <- list()
  for (k in 1:K){
    multi_group_ind[[k]] <- outer(newOrd, newOrd, multiSub, 
                                  group = group_list_ord[[k]])
  }

  # observed test statistic
  out0 <- matrixStrucTestSub(A = A_ord,
  	                  group_list_ord = group_list_ord,
  	                  Delta = Delta_ord, 
   	                  multi_group_ind = multi_group_ind,
   	                  A_upper_ind = A_upper_ind,
   	                  K = K)

  Gamma0 <- out0$Gamma_overall
  Gamma0k <- out0$Gamma_multi
  t0 <- out0$t_overall
  t0k <- out0$t_multi

  # permutations
  Gamma_overall <- rep(NA, B)
  Gamma_k_mat <- matrix(nrow=B, ncol=K) 
  t_overall <- rep(NA, B)
  t_k_mat <- matrix(nrow=B, ncol=K)

  for (b in 1:B){
    perm <- sample(1:p)
    Ab <- A_ord[perm, perm]
    out <- matrixStrucTestSub(A = Ab,
    	                 group_list_ord = group_list_ord,
    	                 Delta = Delta_ord,
     	                 multi_group_ind = multi_group_ind,
    	                 A_upper_ind = A_upper_ind,
     	                 K = K)
    Gamma_overall[b] <- out$Gamma_overall
    Gamma_k_mat[b,] <- out$Gamma_multi
    t_overall[b] <- out$t_overall
    t_k_mat[b,] <- out$t_multi
  }
  
  # t-statistic 
  pt_overall_one_sided <- mean(c(t0, t_overall) >= t0)
  pt_overall_two_sided <- mean(abs(c(t0, t_overall)) >= abs(t0))
  t_max_one_sided <- apply(t_k_mat, 1, max)
  t_max_two_sided <- apply(t_k_mat, 1, function(x){max(abs(x))})
  pt_multi_one_sided <- sapply(t0k, function(x){
                               mean(c(t_max_one_sided, x) >= x)})
  pt_multi_two_sided <- sapply(t0k, function(x){
                               mean(c(t_max_two_sided, abs(x)) >= abs(x))})
  names(pt_multi_one_sided) <- names(pt_multi_two_sided) <- names(group_list)
  names(t0k) <- colnames(t_k_mat) <- names(group_list)
   
  # Hubert's Gamma
  pG_overall_one_sided <- mean(c(Gamma0, Gamma_overall) >= Gamma0)
  pG_overall_two_sided <- mean(abs(c(Gamma0, Gamma_overall)) >= abs(Gamma0))
  Gamma_max_one_sided <- apply(Gamma_k_mat, 1, max)
  Gamma_max_two_sided <- apply(Gamma_k_mat, 1, function(x){max(abs(x))})
  pG_multi_one_sided <- sapply(Gamma0k, function(x){
                               mean(c(Gamma_max_one_sided, x) >= x)})
  pG_multi_two_sided <- sapply(Gamma0k, function(x){
                               mean(c(Gamma_max_two_sided, abs(x)) >= abs(x))})
  names(pG_multi_one_sided) <- names(group_list)
  names(Gamma0k) <- names(group_list)
   
  ret <- list(pt_overall_one_sided = pt_overall_one_sided,
              pt_overall_two_sided = pt_overall_two_sided,
              pt_multi_one_sided = pt_multi_one_sided,
              pt_multi_two_sided = pt_multi_two_sided,
              t0 = t0,
              t0k = t0k,
              t_overall = t_overall,
              t_max_one_sided = t_max_one_sided,
              t_max_two_sided = t_max_two_sided,
              pG_overall_one_sided = pG_overall_one_sided,
              pG_overall_two_sided = pG_overall_two_sided,
              pG_multi_one_sided = pG_multi_one_sided,
              pG_multi_two_sided = pG_multi_two_sided,
              Gamma0 = Gamma0,
              Gamma0k = Gamma0k,
              Gamma_overall = Gamma_overall,
              Gamma_max_one_sided = Gamma_max_one_sided,
              Gamma_max_two_sided = Gamma_max_two_sided,
              group_list = group_list,
              B = B)
  class(ret) <- "mst"
  return(ret)
}

print.mst <- function(x, ...) {
  #' Print results from \code{matrixStrucTest}
  #'
  #' This function prints results from an object returned by \code{matrixStrucTest}.
  #' @param x Output from matrixStrucTest
  #' @param ... Further arguments passed to print
  #' @keywords matrixStrucTest print
  #' @export

  if(x$pG_overall_two_sided > 1/(x$B + 1)) {
    result <- paste("    Test of matrix structure   \n\n",
      prettyNum(x$B, big.mark = ","), " Monte Carlo resamples, two-sided p-values",
      "\n\nOverall Hubert's Gamma",
      "\nGamma0 = ", signif(x$Gamma0, 3), ", p-val = ", signif(x$pG_overall_two_sided, 3),
      "\n\nBlock-specific Hubert's Gamma", sep = "")
  } else {
    result <- paste("    Test of matrix structure   \n\n",
      prettyNum(x$B, big.mark = ","), " Monte Carlo resamples, two-sided p-values",
      "\n\nOverall Hubert's Gamma",
      "\nGamma0 = ", signif(x$Gamma0, 3), ", p-val < ", 1/x$B,
      "\n\nBlock-specific Hubert's Gamma", sep = "")
  }

  p_multi <- x$pG_multi_two_sided
  hit_min <- p_multi == 1/(x$B + 1)
  p_multi <- as.character(signif(p_multi, 3))
  p_multi[hit_min] <- paste("<", 1/(x$B))
  writeLines(result)
  print(data.frame(Gamma0 = signif(x$Gamma0k, 3),
                   pval = p_multi),
        ...)
}

prepBoxPlots <- function(A, groups = NULL, group_list = NULL, absolute = TRUE){
  #' Prepare data for box plots
  #'
  #' This function prepares the data for making box plots.
  #' @param A Distance or similarity matrix, e.g. correlation
  #' @param groups CFA model in lavaan syntax. Either \code{groups} or \code{group_list} but not both must be supplied.
  #' @param group_list List of groupings. Either \code{groups} or \code{group_list} but not both must be supplied.
  #' @param absolute Use the absolute values of A (defaults to TRUE)
  #' @return multi: data frame for making box plots for block-specific tests
  #' @return overall: data frame for making box plots for overall test
  #' @export
  #' @examples
  #' library(matrixStrucTest)
  #' library(ggplot2)
  #' 
  #' data("big5")
  #' 
  #' # get column numbers for questionnaire items
  #' items <- grep("[0-9]", colnames(big5))
  #' 
  #' # compute Spearman's correlation matrix
  #' A <- cor(big5[, items], use = "complete.obs", method = "spearman")
  #' 
  #' groups <- "extrovert ~ E1 + E2 + E3 + E4 + E5 + E6 + E7 + E8 + E9 + E10
  #'            neurotic ~ N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9 + N10
  #'            agreeable ~ A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9 + A10
  #'            conscientious ~ C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10
  #'            open ~ O1 + O2 + O3 + O4 + O5 + O6 + O7 + O8 + O9 + O10"
  #' 
  #' # Make box plots contrasting within and between group correlations
  #' box <- prepBoxPlots(A = A, groups = groups, absolute = TRUE)
  #' 
  #' ggplot(aes(x = as.factor(delta), y = a), data = box$overall)+
  #'   geom_boxplot()+
  #'   theme_bw(22)+
  #'   labs(x = expression(Delta), y="|a|")
  #' 
  #' dev.new(width = 12, height = 5)
  #' ggplot(aes(x = as.factor(delta), y = a), data = box$multi)+
  #'   geom_boxplot()+
  #'   facet_grid(~block)+
  #'   theme_bw(22)+
  #' labs(x = expression(Delta), y = "|a|")
   
  if(!is.null(groups) & !is.null(group_list)) {
    stop("Both a 'groups' and 'group_list' were supplied as arguments. Only one of the two can be given.")
  } else if(is.null(groups) & is.null(group_list)) {
    stop("Either 'groups' or 'group_list' must be supplied as an argument")
  }

  # if groups given, convert to group_list
  if(!is.null(groups)) {
    if(is.null(colnames(A))) {
      stop("When 'groups' is supplied, the matrix A must have column names")
    }

    group_list <- makeGroupList(groups, A)
  }

   if(absolute){
      A <- abs(A)
   }
   
   p <- ncol(A)
   K <- length(group_list)
   
   if(sum(sapply(group_list,length)) > p){
      warning("Some variables may be assigned to multiple clusters")
   } else if(sum(sapply(group_list,length)) < p){
      warning("Some variables not assigned to a cluster")
   }
   
  # clean up names of groups
  if(is.null(names(group_list))) {
    names(group_list) <- paste("block", 1:K)
  } else {
    no_name <- which(sapply(names(group_list), nchar) == 0)
    names(group_list)[no_name] <- paste("block", no_name)
  }

   # create delta matrix
   Delta <- matrix(nrow = p, ncol = p)
   for (i in 1:p){
      for (j in 1:p){
         Delta[i,j] <- deltaSub(i, j, group_list)
      }
   }

  # order matrices by hypothesized groups
  ord <- unlist(group_list)
  A_ord <- A[ord, ord]
  Delta_ord <- Delta[ord, ord]

  # get group list for reordered matrices
  nk <- sapply(group_list, length)
  firstIndex <- 1
  group_list_ord <- list()
  for (k in 1:K) {
    group_list_ord[[k]] <- firstIndex:(firstIndex + nk[k] -1)
    firstIndex <- firstIndex + nk[k]
  }

  # get indicator matrices for block-specific tests
  newOrd <- unlist(group_list_ord)
  A_upper_ind <- upper.tri(A)
  multi_group_ind <- list()
  for (k in 1:K){
    multi_group_ind[[k]] <- outer(newOrd, newOrd, multiSub, 
                                  group = group_list_ord[[k]])
  }
 
  # observed test statistic
  out <- matrixStrucTestSub(A = A_ord,
                     group_list_ord = group_list_ord,
                     Delta = Delta_ord, 
                     multi_group_ind = multi_group_ind,
                     A_upper_ind = A_upper_ind,
                     K = K)
   
  box_overall <- data.frame(a = out$A_upper, delta = out$Delta_upper)
   
  box_list <- list()
  for (k in 1:K){
    box_list[[k]] <- data.frame(a = out$Ak_list[[k]], 
                                delta = out$Deltak_list[[k]],
                                block = names(group_list)[k])
  }
   
  box_multi <- do.call(rbind, box_list)
  
  return(list(multi = box_multi, overall = box_overall))
}
