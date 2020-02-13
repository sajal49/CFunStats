# Computes the most compact contingency table and calculates the compact functional statistics
# Input : x (an ordinal factor), y (a factor), full (logical; if no trimming is desired),
# log.p (logical; if log of pvalue is desired)
# Output: Compact functional statistics for the compacted contingency table
# Created by: Sajal Kumar and Dr. Mingzhou (Joe) Song
# Date Created: 18th January, 2020

#' c_fun_stats
#'
#' Computes the most compact contingency table and calculates the compact functional statistics
#'
#' @importFrom FunChisq fun.chisq.test
#'
#' @param x (factor, integer or numeric) An ordinal categorical independent variable. An ordered factor is recommended,
#' with pre-defined ordering scheme as Y and X are sorted w.r.t X.
#' @param y (factor, integer or numeric) A categorical dependent variable.
#' @param full (logical) Disables trimming and returns the compressed global pattern. Default: FALSE
#' @param log.p (logical) If natural log of the p-value should be returned. Default: FALSE
#'
#' @return Compact functional statistics for the most compact pattern.
#' A list with class "\code{htest}" containing the following components:
#' \describe{
#'  \item{\code{statistic}}{the functional chi-squared statistic}
#'
#'  \item{\code{parameter}}{degrees of freedom for the functional chi-squared statistic}
#'
#'  \item{\code{p.value}}{p-value of the functional test computed by an asymptotic chi-squared distribution}
#'
#'  \item{\code{estimate}}{an estimate of function index between 0 and 1.
#' The value of 1 indicates a strictly mathematical function. It is
#' asymmetrical with respect to transpose of the input contingency table,
#' different from the symmetrical Cramer's V based on the Pearson's
#' chi-squared test statistic.}
#' }
#' @keywords Compact representation, Local pattern, Functional dependency
#'
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
#'  print(c_fun_stats(x, y))
c_fun_stats = function(x, y, full=FALSE, log.p=FALSE) {

  # get the compact pattern table
  ctable = construct_compact_table(x, y, full)

  # caluclate the funchisq stats
  stats = fun.chisq.test(ctable, log.p = log.p)

  # organize and return result
  return(structure(list(statistic = stats$statistic, p.value = stats$p.value,
                        estimate = stats$estimate, data.name = "Compacted x->y",
                        method = "Compact Functional Statistics"),
                   class = "htest"))

}
