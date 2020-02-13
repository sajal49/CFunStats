#include<iostream>
#include<climits>
#include<math.h>
#include "define.h"

mydouble getFCstats(const std::vector<std::vector<int> > & O, const std::string method){
  // For internal use

  mydouble estimate, fc;

  // Given a table, get its FunChisq statistics
  fc = funchisq(O,estimate,"conditional","non-constant");
  int df = (O.size() - 1)*(O[0].size() - 1);

  if(method == "nfchisq") {
    fc = (fc - df)/(sqrt(2*df));
  }

  return fc;
}


frame<double> makeDynTable(const vec<int> & x, const vec<int> & y,
                   const vec<int> & sorted_unx, const vec<int> & sorted_uny,
                   const vec<int> & slices_x, const std::string method,
                   frame<int> & clusttable)
{
  // dyntable is a 2d matrix (frame) with number of points on rows and number of clusters
  // on columns.
  frame<double> dyntable(sorted_unx.size(), std::vector<double>(sorted_unx.size(),0));

  // helper variables
  std::vector<int> tempx, tempy;
  frame<int> temp_table;
  vec<int> bestclust;
  vec<double> bestfc;
  mydouble fc;
  int temp, rec_k, rec_j;

  // Base cases
  // if j = 0, dyntable[i,0] = 0; If there is only one cluster, then stats is 0
  for(int i=0; i<dyntable.size(); i++)
    dyntable[i][0] = 0;


  // if i==j, dyntable[i,i] = stats for each all x in their own cluster
  for(int i=1; i<dyntable.size(); i++) {

    // clear the temp vectors
    tempx.clear();
    tempy.clear();

    // count the subset of data denoted by slices
    for(int j=0; j<slices_x[i]; j++){
      tempx.push_back(x[j]);
      tempy.push_back(y[j]);
    }

    // make table of the subset and compute the funchisq statistics
    temp_table = tableCpp(tempx, tempy, -1, sorted_uny.size());
    fc = getFCstats(temp_table, method);

    dyntable[i][i] = fc;
    clusttable[i][i] = i-1;
  }

  // Recursive Case
  // For all other cases we look for the best solution in [i-1 : j - 1, j-1], neighboring cells.
  for(int i=1; i<dyntable.size(); i++){
    for(int j=1; j<i; j++){

      bestclust.clear(); // keep track of all possible cluster floors.
      bestfc.clear(); // keep track of best fchisq score so far.

      for(int k=i-1; k>=j-1; k--){

        // prepare helper variables
        tempx.clear();
        tempy.clear();
        temp = 0;
        rec_k = k;
        rec_j = j-1;

        for(int itr=slices_x[i]-1; itr>=0; itr--){

          // Only updates itself when a cluster floor is being hit in clusttable.
          // -1 means we have reached the first column (also denoted by rec_j reaching 0).
          if(rec_k!=-1 && itr<slices_x[rec_k]){
            temp++;
            rec_k = clusttable[rec_k][rec_j];
            if(rec_j>0)
              rec_j--;
          }

          tempx.push_back(temp);
          tempy.push_back(y[itr]);
        }

        // Fix the reverse order of clusters.
        temp = *std::max_element(tempx.begin(), tempx.end());
        for(int i=0; i<tempx.size(); i++){
          tempx[i] = (int)std::abs((tempx[i] - temp));
        }

        // make table of the subset and compute the funchisq statistics
        temp_table = tableCpp(tempx, tempy, -1, sorted_uny.size());
        fc = getFCstats(temp_table, method);

        bestfc.push_back(fc);
        bestclust.push_back(k);
      }

      // Find the best fc and update dyntable , clusttable accordingly.
      fc = bestfc[0];
      temp = bestclust[0];

      for(int i=1; i<bestclust.size(); i++){
        if(fc<bestfc[i]){
          fc = bestfc[i];
          temp = bestclust[i];
        }
      }

      dyntable[i][j] = fc;
      clusttable[i][j] = temp;
    }
  }

  return dyntable;
}

frame<int> tableCpp(std::vector<int> x_vec, std::vector<int> y_vec, int xlevels, int ylevels) {

  int min_val_x = INT_MAX, min_val_y = INT_MAX, max_val_x = INT_MIN, max_val_y = INT_MIN;

  // Find min and max for both vectors
  for(int i=0; i<x_vec.size(); i++){
    if(x_vec[i] < min_val_x)
      min_val_x = x_vec[i];

    if(x_vec[i] > max_val_x)
      max_val_x = x_vec[i];

    if(y_vec[i] < min_val_y)
      min_val_y = y_vec[i];

    if(y_vec[i] > max_val_y)
      max_val_y = y_vec[i];
  }

  //If smallest element is greater than 1 then subtract min_val to offset all elements to 0,
  if(min_val_x >= 1)
  {
    for(int i=0; i<x_vec.size(); i++){
      x_vec[i]-=min_val_x;
    }
    max_val_x-=min_val_x;
  }

  if(min_val_y >= 1)
  {
    for(int i=0; i<y_vec.size(); i++){
      y_vec[i]-=min_val_y;
    }
    max_val_y-=min_val_y;
  }

  if(xlevels == -1 && ylevels == -1){
    xlevels = max_val_x+1;
    ylevels = max_val_y+1;
  } else if(xlevels == -1){
    xlevels = max_val_x+1;
  } else if(ylevels == -1){
    ylevels = max_val_y+1;
  }

  std::vector<std::vector<int > > table(xlevels, std::vector<int>(ylevels,0));

  for(int i=0; i<x_vec.size(); i++)
    table[x_vec[i]][y_vec[i]]++;

  return table;
}

/*void print_tableCpp(frame<int> table){
  // For internal use

  // column names
  std::cout<<"   ";
  for(int j=0; j<table[0].size(); j++){
    std::cout<<j<<" ";
  }
  std::cout<<std::endl;

  // table
  for(int i=0; i<table.size(); i++){
    // rownames
    std::cout<<i<<": ";
    for(int j=0; j<table[0].size(); j++){
      std::cout<<table[i][j]<<" ";
    }
    std::cout<<std::endl;
  }

  std::cout<<std::endl;

}*/

