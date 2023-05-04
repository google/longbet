#include <ctime>
#include <RcppArmadillo.h>
#include "common.h"

using namespace std;


void rcpp_to_std2(arma::mat y, arma::mat X, arma::mat Xtest, Rcpp::NumericMatrix &y_std, double &y_mean, Rcpp::NumericMatrix &X_std, Rcpp::NumericMatrix &Xtest_std, matrix<size_t> &Xorder_std);

void arma_to_std(const arma::mat &matrix_in, std::vector<double> &vector_out);

// transfers data from an armadillo matrix object (column 0) to an std vector object
void arma_to_std(const arma::mat &matrix_in, std::vector<size_t> &vector_out);

void arma_to_rcpp(const arma::mat &matrix_in, Rcpp::NumericMatrix &matrix_out);

void arma_to_std_ordered(const arma::mat &matrix_in, matrix<size_t> &matrix_ordered_std);

void std_to_rcpp(const matrix<double> &matrix_in, Rcpp::NumericMatrix &matrix_out);

double compute_vec_mean(const std::vector<double> &vec);

double compute_mat_mean(const Rcpp::NumericMatrix &matrix);

void rcpp_to_std2(arma::mat X, arma::mat Xtest, Rcpp::NumericMatrix &X_std, Rcpp::NumericMatrix &Xtest_std, matrix<size_t> &Xorder_std);

void std_to_arma(matrix<double> &mat_std, arma::mat &mat_arma);
