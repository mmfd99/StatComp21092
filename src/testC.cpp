#include <Rcpp.h>
using namespace Rcpp;
//' @title A test function
//' @description A test function
//' @param x a number
//' @export
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param thin the number of between-sample random numbers
//' @return random sample list
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int thin) {
  NumericMatrix mat(N, 2);
  double x=0, y=0;
  int n=20;
  double a=2,b=1;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rbinom(1, n, y)[0];
      y = rbeta(1, x+a, n-x+b)[0];
    }
    mat(i, 0)=x;
    mat(i, 1)=y;
  }
  return(mat);
}
