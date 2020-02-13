//
//  define.h
//  eft-mhg
//
//  Created by Hua Zhong on 6/28/15.
//  Modified by Sajal Kumar
//  Copyright (c) 2015 New Mexico State University. All rights reserved.
//

#if defined _WIN32
typedef double mydouble;
//typedef long double mydouble;
#else
typedef long double mydouble;
#endif

#include <vector>
#include <string>
#include <algorithm>
template <typename T>
using frame = std::vector<std::vector<T> >;

template <typename T>
using vec = std::vector<T>;


mydouble funchisq(const std::vector<std::vector<int> > & O, mydouble & estimate,
                  const std::string index_kind, const std::string method);

mydouble funchisq(const std::vector<std::vector<int> > & O, const std::vector<int> & rowsums,
                  const std::vector<int> & colsums, int n);

mydouble getFCstats(const frame<int> & O, const std::string method);


frame<double> makeDynTable(const vec<int> & x, const vec<int> & y,
                           const vec<int> & sorted_unx, const vec<int> & sorted_uny,
                           const vec<int> & slices_x, const std::string method,
                           frame<int> & clusttable);

frame<int> tableCpp(std::vector<int> x_vec, std::vector<int> y_vec, int xlevels, int ylevels);

void print_tableCpp(frame<int> table);
