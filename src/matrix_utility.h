#ifndef THIRD_PARTY_R_PACKAGES_LONGBET_LONGBET_SRC_MATRIX_UTILITY_H_
#define THIRD_PARTY_R_PACKAGES_LONGBET_LONGBET_SRC_MATRIX_UTILITY_H_

# include "common.h"

// all code belongs to https://www.codeproject.com/Articles/1268576/Singular-Values-Decomposition-SVD-In-Cplusplus11-B#:~:text=Singular%20value%20decomposition%20(Singular%20Value,to%20visualize%20the%20available%20data.

void jordan_gaussian_transform(matrix<double> mat,
std::vector<double>& eigenvector);

void get_hermitian_matrix(std::vector<double> eigenvector,
matrix<double>& h_matrix);

void get_hermitian_matrix_inverse(std::vector<double> eigenvector,
matrix<double>& ih_matrix);

void matrix_by_matrix(matrix<double> matrix1, matrix<double>& matrix2,
matrix<double>& matrix3);

void get_reduced_matrix(matrix<double> mat, matrix<double>& r_matrix,
std::size_t new_size);

void compute_evd(matrix<double> mat, std::vector<double> &eigenvalues,
matrix<double> &eigenvectors, std::size_t eig_count);

void matrix_transpose(matrix<double> matrix1, matrix<double>& matrix2);

void get_inverse_diagonal_matrix(matrix<double> mat,
matrix<double>& inv_matrix);

void svd(matrix<double> mat, matrix<double>& s, matrix<double>& u,
matrix<double>& v);

void ginv(matrix<double> &mat, matrix<double> &mat_inv);

#endif  // THIRD_PARTY_R_PACKAGES_LONGBET_LONGBET_SRC_MATRIX_UTILITY_H_
