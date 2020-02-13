#include <Rcpp.h>
#include "define.h"
using namespace Rcpp;

// [[Rcpp::export]]
List DynTable(NumericVector X, NumericVector Y, NumericVector sorted_x, NumericVector sorted_y,
              NumericVector slices_x, std::string method){

  // Conversions :
  std::vector<int> x = as<std::vector<int>>(X);
  std::vector<int> y = as<std::vector<int>>(Y);

  std::vector<int> sx = as<std::vector<int>>(sorted_x);
  std::vector<int> sy = as<std::vector<int>>(sorted_y);

  std::vector<int> slc_x = as<std::vector<int>>(slices_x);

  // clusttable is a 2d matrix with the same dimension as dyntable. It keeps track of
  // the floor of the last cluster, to aid in backtracking.
  frame<int> clusttable(sx.size(), std::vector<int>(sx.size(),-1));

  // Get the dyntable
  frame<double> tb = makeDynTable(x,y,sx,sy,slc_x,method,clusttable);

  // Return as list
  return List::create(Named("dyn_table") = tb, Named("clust_table") = clusttable);
}

// [[Rcpp::export]]
NumericMatrix tableRcpp(NumericVector x, NumericVector y, int xlevels=-1, int ylevels=-1) {

  //Convert numeric vector to standard vector
  std::vector<int> x_vec = as< std::vector<int> >(x);
  std::vector<int> y_vec = as< std::vector<int> >(y);

  // Get table from CPP
  frame<int> tablecpp = tableCpp(x_vec, y_vec, xlevels, ylevels);

  // Convert to numeric matrix
  NumericMatrix table(tablecpp.size(), tablecpp[0].size());

  for(int i=0; i<tablecpp.size(); i++){
    for(int j=0; j<tablecpp[0].size(); j++){
      table(i,j) = tablecpp[i][j];
    }
  }

  return table;
}

