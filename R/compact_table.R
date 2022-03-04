# Employs dynamic programming to trim and compress a pattern to its most compact cross-tabular representation
# Input : x (an ordinal factor), y (a factor), full (logical; if no trimming is desired)
# Output: A compacted contingency table of X (rows) and Y (columns)
# Created by: Sajal Kumar and Dr. Mingzhou (Joe) Song
# Date Created: 18th January, 2020

#' construct_compact_table
#'
#' Employs dynamic programming to trim and compress a pattern to its most compact cross-tabular representation
#'
#' @importFrom DescTools Untable
#' @importFrom FunChisq fun.chisq.test
#' @importFrom stats pchisq
#' @import Rcpp
#' @useDynLib CFunStats
#'
#' @param x (factor, integer or numeric) An ordinal categorical independent variable. An ordered factor is recommended,
#' with pre-defined ordering scheme as y and x are sorted w.r.t x.
#' @param y (factor, integer or numeric) A categorical dependent variable.
#' @param full (logical) Disables trimming and returns the compressed global pattern. Default: FALSE
#'
#' @return
#' \code{ctable}: most compact cross cross-tabular representation of \code{x} and \code{y}.
#'
#' @keywords Compact representation, Dynamic programming, Trimming, Local pattern.
#' @export
#' @examples
#' library(DescTools)
#' xy_table = matrix(c( 20,20,20,
#'                      20,20,20,
#'                      20,20,20,
#'                      10,0,0,
#'                      0,10,0,
#'                      0,0,10,
#'                      0,10,0,
#'                      10,0,0,
#'                      20,20,20,
#'                      20,20,20,
#'                      20,20,20 ), nrow=11, ncol=3, byrow = TRUE)
#'  xy = Untable(xy_table)
#'  x = xy$Var1
#'  y = xy$Var2
#'  print(construct_compact_table(x, y))
construct_compact_table = function(x, y, full=FALSE){

  # X and Y should be a factor, integer or numeric
  class_x = class(x)
  class_y = class(y)

  if(!("numeric" %in% class_x  || "factor" %in% class_x || "integer" %in% class_x)){
    stop("x must be an integer, numeric or factor!!")
  }

  if(!("numeric" %in% class_y  || "factor" %in% class_y || "integer" %in% class_y)){
    stop("y must be an integer, numeric or factor!!")
  }

  # save original X and Y
  org_x = as.character(x)
  org_y = as.character(y)

  # convert to consecutive numbers
  x = as.numeric(as.factor(x))
  y = as.numeric(as.factor(y))

  # sort x and y according to x
  ord_x = order(x, decreasing = FALSE)
  y = y[ord_x]
  org_x = org_x[ord_x]
  org_y = org_y[ord_x]
  x = x[ord_x]

  # save original labels for X and Y
  org_unx = unique(org_x)
  org_uny = unique(org_y)

  # construct original table
  org_table = table(x, y)

  if(!full){ # if trimming is desired

    # compute funchisq for  all subtables
    subtabs = as.data.frame(matrix(0, ncol = 5))

    count = 1
    for(i in 3:(nrow(org_table))){
      for(j in 1:(nrow(org_table)-i+1)){
        q = fun.chisq.test(org_table[j:(j+i-1),], log.p = TRUE)
        if(q$statistic == 0){
          subtabs[count,] = c(beg = j,
                              end = (j+i-1),
                              statistic = 0,
                              pvalue = 0,
                              estimate = 0)
        } else {
          subtabs[count,] = c(beg = j,
                              end = (j+i-1),
                              statistic = q$statistic,
                              pvalue = q$p.value,
                              estimate = q$estimate)
        }

        count = count + 1
      }
    }

    # adjust p-values by Benjamini Hochberg
    pvalues = subtabs[,4]
    pvalues_r = rank(as.numeric(pvalues))
    pvalues = as.numeric(pvalues) + log(nrow(subtabs)) - log(pvalues_r)
    pvalues[pvalues > 0] = 0
    subtabs[,4] = pvalues

    # filter
    subtabs = subtabs[subtabs[,4] < (-3),] # equivalent to 0.05

    # min pvalue per starting point
    min_pval_perx = sapply(c(1:nrow(org_table)), function(i){
      perx_pval = subtabs[which(subtabs[,1]==i),4]
      return(ifelse((length(perx_pval)==0 || all(is.na(min(perx_pval)))), 0, min(perx_pval, na.rm = TRUE)))
    })

    # # [conservative]
    # min_pval_perx = round(min_pval_perx)

    # get the smaller index [conservative]
    ind = min(which(min_pval_perx == min(min_pval_perx)))
  } else { # if trimming is not desired
    ind = 1
  }

  # subset the original table and optimally compress it
  xy = Untable(org_table[c(ind:nrow(org_table)),])
  x = as.numeric(as.factor(xy$x))
  y = as.numeric(as.factor(xy$y))

  # sort x and y
  ord_x = order(x, decreasing = FALSE)
  y = y[ord_x]
  x = x[ord_x]

  # prepare for compression
  slices_x = which(!duplicated(x))
  slices_x = c(slices_x[2:length(slices_x)],length(x)+1) - 1

  out_dyn = DynTable(X = x-1, Y = y-1, sorted_x = c(1:max(x)),
                     sorted_y = c(1:max(y)), slices_x = slices_x,
                     method = "fchisq")

  # get mar dyn table
  mar_dyn = out_dyn$mar_dyn_table

  # get cluster info
  clust_dyn = out_dyn$clust_table

  # get dyn table
  out_dyn = out_dyn$dyn_table

  if(!full){

    # we now find the best compressed ending
    pval_table = vector("list", length = length(out_dyn))
    for(i in 1:length(out_dyn)){
      pvals = rep(0, length(out_dyn))
      for(j in 1:length(out_dyn)){

        fc = (out_dyn[[i]][j] - mar_dyn[[i]][j])
        pvals[j] = pchisq(fc, df = ((j-1)*(ncol(org_table)-1)),
                          lower.tail = FALSE, log.p = TRUE)
      }
      pval_table[[i]] = pvals
    }

    min_pval = which.min(unlist(lapply(pval_table, min)))
    comp_lvls = lapply(pval_table, which.min)[[min_pval]]

    # prepare the trimmed and compressed table
    xy = Untable(org_table[c(ind:(ind+min_pval-1)),])
    x = as.numeric(as.factor(xy$x))
    y = as.numeric(as.factor(xy$y))

    # sort x and y
    ord_x = order(x, decreasing = FALSE)
    y = y[ord_x]
    x = x[ord_x]

    # get sub-slices
    slices_x = which(!duplicated(x))
    slices_x = c(slices_x[2:length(slices_x)],length(x))

    # Get the transformed x from clust_table
    temp = 1
    rec_k = clust_dyn[[min_pval]][comp_lvls] + 1
    rec_j = comp_lvls - 1
    x_new = rep(0, slices_x[length(slices_x)])
    for(itr in length(x_new) : 1)
    {
      if(rec_k!=0 && itr<(slices_x[rec_k]))
      {
        temp = temp + 1
        rec_k = clust_dyn[[rec_k]][rec_j] + 1
        if(rec_j>1)
          rec_j = rec_j - 1
      }
      x_new[itr] = temp
    }

    # adjust x_new
    x_new = abs(x_new - (temp+1))

    # adjust old_x and old_y
    new_org_x = c()
    for(i in ind:(ind+min_pval-1)){
      new_org_x = c(new_org_x, subset(org_x, org_x==org_unx[i]))
    }

    # names of x_new
    unx_new = sort(unique(x_new))
    rname_x_new = rep("",length(unx_new))
    for(i in 1:length(unx_new))
      rname_x_new[i] = paste(sort(unique(new_org_x[x_new == unx_new[i]])),collapse = "#")

  } else {

    # Adjust slices_x
    slices_x[1:(length(out_dyn) - 1)] = slices_x[1:(length(out_dyn) - 1)] + 1

    # Compute p-value for the last row
    last_row_stat = out_dyn[[length(out_dyn)]] - mar_dyn[[length(out_dyn)]]
    last_row_pvalue = rep(1, length(last_row_stat))
    for(i in 2:length(last_row_stat))
    {
      stat = last_row_stat[i]
      df = (ncol(org_table) - 1) * (i-1)
      last_row_pvalue[i] = pchisq(q=stat, df=df, lower.tail=FALSE, log.p = TRUE)
    }

    # optimal cluster from the dyn_table
    opt_clust = which.min(last_row_pvalue)

    # Get the transformed x from clust_table
    temp = 1
    rec_k = clust_dyn[[length(out_dyn)]][opt_clust] + 1
    rec_j = opt_clust - 1
    x_new = rep(0, length(x))
    for(itr in length(x) : 1)
    {
      if(rec_k!=0 && itr<(slices_x[rec_k]))
      {
        temp = temp + 1
        rec_k = clust_dyn[[rec_k]][rec_j] + 1
        if(rec_j>1)
          rec_j = rec_j - 1
      }
      x_new[itr] = temp
    }

    # adjust x_new
    x_new = abs(x_new - (temp+1))

    # names of x_new
    unx_new = sort(unique(x_new))
    rname_x_new = rep("",length(unx_new))
    for(i in 1:length(unx_new))
      rname_x_new[i] = paste(sort(unique(org_x[x_new == unx_new[i]])),collapse = "#")
  }

  # compact table
  ctable = tableRcpp(x_new-1, y-1, ylevels = length(org_uny))
  rownames(ctable) = rname_x_new
  colnames(ctable) = org_uny
  return(ctable)

}
