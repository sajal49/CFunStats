# tests c_fun_stats and construct_compact_table functions
# Created by: Sajal Kumar and Dr. Mingzhou (Joe) Song
# Date Created: 18th January, 2020

library(testthat)
library(CFunStats)
library(DescTools)

context("Testing CFunStats")

test_that("Testing construct_compact_table and c_fun_stats", {

  # test 1

  x <- c(1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4, 5,5,5,5,5)
  y <- c(1,1,1,1, 1,1,2,2, 2,2,2,2, 2,2, 3,3,3,3,3)

  ctable <- construct_compact_table(x, y)

  expect_equivalent(ctable, matrix(c(6, 2, 0,
                                     0, 6, 0,
                                     0, 0, 5), nrow=3, ncol=3, byrow=TRUE))
  cfunstats <- c_fun_stats(x, y)

  expect_equivalent(round(cfunstats$statistic, digits=3),
                    28.263)
  expect_equivalent(signif(cfunstats$p.value, digits=3),
                    1.1e-05)
  expect_equivalent(round(cfunstats$estimate, digits=3),
                    0.871)

  # test 2
  xy_table = matrix(c( 5,  5,  5,
                      10, 10, 10,
                      20, 20, 20,
                      10,  0,  0,
                       0, 10,  0,
                       0,  0, 10,
                      20, 20, 20,
                      10, 10, 10), nrow=8, ncol = 3, byrow = TRUE)
  xy = Untable(xy_table)
  x = as.factor(xy$Var1)
  y = as.factor(xy$Var2)

  ctable = construct_compact_table(x, y)

  expect_equivalent(ctable, matrix(c(10,  0,  0,
                                      0, 10,  0,
                                      0,  0, 10), nrow=3, ncol=3, byrow=TRUE))

  cfunstats = c_fun_stats(x, y)

  expect_equivalent(round(cfunstats$statistic, digits=3),
                    60)
  expect_equivalent(signif(cfunstats$p.value, digits=3),
                    2.9e-12)
  expect_equivalent(round(cfunstats$estimate, digits=3),
                    1)


  # test 3
  xy_table = matrix(c(  5,  5,  5,
                       10, 10, 10,
                       20, 20, 20,
                       10,  0,  0,
                        0, 10,  0,
                        0,  0, 10,
                       20, 20, 20,
                       10, 10, 10), nrow=8, ncol = 3, byrow = TRUE)
  xy = Untable(xy_table)
  x = as.factor(xy$Var1)
  y = as.factor(xy$Var2)

  ctable = construct_compact_table(x, y, full = TRUE)

  expect_equivalent(ctable, matrix(c(35, 35, 35,
                                     10,  0,  0,
                                      0, 10,  0,
                                      0,  0, 10,
                                     30, 30, 30), nrow=5, ncol=3, byrow=TRUE))

  cfunstats = c_fun_stats(x, y, full = TRUE)

  expect_equivalent(round(cfunstats$statistic, digits=3),
                    60)
  expect_equivalent(signif(cfunstats$p.value, digits=3),
                    4.6e-10)
  expect_equivalent(round(cfunstats$estimate, digits=3),
                    0.365)


  # test 4
  xy_table = matrix(c( 50,  0,  0,
                       10, 10, 10,
                       20, 20, 20,
                       10,  0,  0,
                        0, 10,  0,
                        0,  0, 10,
                       20, 20, 20,
                       10, 10, 10,
                        0,  0, 50), nrow=9, ncol = 3, byrow = TRUE)
  xy = Untable(xy_table)
  x = as.factor(xy$Var1)
  y = as.factor(xy$Var2)

  ctable = construct_compact_table(x, y)

  expect_equivalent(ctable, matrix(c(50,  0,  0,
                                     30, 30, 30,
                                     10,  0,  0,
                                      0, 10,  0,
                                      0,  0, 10,
                                     30, 30, 30,
                                      0,  0, 50), nrow=7, ncol=3, byrow=TRUE))

  cfunstats = c_fun_stats(x, y)

  expect_equivalent(round(cfunstats$statistic, digits=3),
                    243.871)
  expect_equivalent(signif(cfunstats$p.value, digits=3),
                    2.59e-45)
  expect_equivalent(round(cfunstats$estimate, digits=3),
                    0.635)


  # test 5
  xy_table = matrix(c( 50,  0,  0,
                       10, 10, 10,
                       20, 20, 20,
                       10,  0,  0,
                        0, 10,  0,
                        0,  0, 10,
                       20, 20, 20,
                       10, 10, 10,
                        0,  0, 50), nrow=9, ncol = 3, byrow = TRUE)
  xy = Untable(xy_table)
  x = as.factor(xy$Var1)
  y = as.factor(xy$Var2)

  ctable = construct_compact_table(x, y, full = TRUE)

  expect_equivalent(ctable, matrix(c(50,  0,  0,
                                     30, 30, 30,
                                     10,  0,  0,
                                      0, 10,  0,
                                      0,  0, 10,
                                     30, 30, 30,
                                      0,  0, 50), nrow=7, ncol=3, byrow=TRUE))

  cfunstats = c_fun_stats(x, y, full=TRUE)

  expect_equivalent(round(cfunstats$statistic, digits=3),
                    243.871)
  expect_equivalent(signif(cfunstats$p.value, digits=3),
                    2.59e-45)
  expect_equivalent(round(cfunstats$estimate, digits=3),
                    0.635)


})
