
#include "matrix_utility.h"
#include <cstddef>
#include <iostream>


void jordan_gaussian_transform(matrix<double> mat,
std::vector<double>& eigenvector)
{
  const double eps = 0.000001; bool eigenv_found = false;
  for (std::uint32_t s = 0; s < mat.size() - 1 && !eigenv_found; s++)
  {
    std::uint32_t col = s; double alpha = mat[s][s];
    while ((col < mat[s].size()) && (alpha != 0) && (alpha != 1)){;
      mat[s][col++] /= alpha;
    }

    for (std::uint32_t col = s; col < mat[s].size() && !alpha; col++){
      std::swap(mat[s][col], mat[s + 1][col]);
    }

    for (std::uint32_t row = 0; row < mat.size(); row++)
    {
      double gamma = mat[row][s];
      for (std::uint32_t col = s; col < mat[row].size() && (row != s); col++)
        mat[row][col] = mat[row][col] - mat[s][col] * gamma;
    }

    if ((s == (mat.size() - 2)) || (std::fabs(mat[s + 1][s + 1]) < eps))
    {
      if (s == (mat.size() - 2)) {std::cout << "last column?" << endl;}
      for (size_t row = 0; row < mat.size(); row++){
        eigenvector.push_back(-mat[row][s + 1]);
      }
    }

    if (eigenvector.size() == mat.size())
    {
      eigenv_found = true; eigenvector[s + 1] = 1;
      for (std::uint32_t index = s + 1; index < eigenvector.size(); index++){
        eigenvector[index] = (std::fabs(eigenvector[index]) >= eps) ?
        eigenvector[index] : 0;
      }
    }
  }
}


void get_hermitian_matrix(std::vector<double> eigenvector,
matrix<double>& h_matrix)
{
  std::cout << "eigenvector = " << eigenvector << endl;
  h_matrix.resize(eigenvector.size());
  for (std::uint32_t row = 0; row < eigenvector.size(); row++)
    h_matrix[row].resize(eigenvector.size());

  h_matrix[0][0] = 1 / eigenvector[0];
  for (std::uint32_t row = 1; row < eigenvector.size(); row++)
    h_matrix[row][0] = -eigenvector[row] / eigenvector[0];

  for (std::uint32_t row = 1; row < eigenvector.size(); row++)
    h_matrix[row][row] = 1;
}

void get_hermitian_matrix_inverse(std::vector<double> eigenvector,
matrix<double>& ih_matrix)
{
  ih_matrix.resize(eigenvector.size());
  for (std::uint32_t row = 0; row < eigenvector.size(); row++)
    ih_matrix[row].resize(eigenvector.size());

  ih_matrix[0][0] = eigenvector[0];
  for (std::uint32_t row = 1; row < eigenvector.size(); row++)
    ih_matrix[row][0] = -eigenvector[row];

  for (std::uint32_t row = 1; row < eigenvector.size(); row++)
    ih_matrix[row][row] = 1;
}

void matrix_by_matrix(matrix<double> matrix1, matrix<double>& matrix2,
matrix<double>& matrix3)
{
  matrix3.resize(matrix1.size());
  for (std::uint32_t row = 0; row < matrix1.size(); row++)
  {
    matrix3[row].resize(matrix1[row].size());
    for (std::uint32_t col = 0; col < matrix1[row].size(); col++)
    {
      matrix3[row][col] = 0.00;
      for (std::uint32_t k = 0; k < matrix1[row].size(); k++)
        matrix3[row][col] += matrix1[row][k] * matrix2[k][col];
    }
  }
}

void get_reduced_matrix(matrix<double> mat, matrix<double>& r_matrix,
std::size_t new_size)
{
  r_matrix.resize(new_size);
  std::size_t index_d = mat.size() - new_size;
  std::uint32_t row = index_d, row_n = 0;
  while (row < mat.size())
  {
    r_matrix[row_n].resize(new_size);
    std::uint32_t col = index_d, col_n = 0;
    while (col < mat.size())
      r_matrix[row_n][col_n++] = mat[row][col++];

    row++; row_n++;
  }
}

void compute_evd(matrix<double> mat, std::vector<double> &eigenvalues,
matrix<double> &eigenvectors, std::size_t eig_count)
{
  std::size_t m_size = mat.size();
  std::vector<double> vec(m_size, 1);

  static matrix<double> matrix_i;
  if (eigenvalues.size() == 0 && eigenvectors.size() == 0)
  {
    eigenvalues.resize(m_size);
    eigenvectors.resize(eigenvalues.size());
    matrix_i = mat;
  }

  matrix<double> m;
  m.resize(m_size);
  for (std::uint32_t row = 0; row < m_size; row++) m[row].resize(100);

  double lambda_old = 0;
  std::uint32_t index = 0; bool is_eval = false;
  while (is_eval == false)
  {
    for (std::uint32_t row = 0; row < m_size && (index % 100) == 0; row++)
      m[row].resize(m[row].size() + 100);

    for (std::uint32_t row = 0; row < m_size; row++)
    {
      m[row][index] = 0;
      for (std::uint32_t col = 0; col < m_size; col++)
        m[row][index] += mat[row][col] * vec[col];
    }

    for (std::uint32_t col = 0; col < m_size; col++)
      vec[col] = m[col][index];
    if (index > 0)
    {
      double lambda = (m[0][index - 1] != 0) ? (m[0][index] / m[0][index - 1]) : m[0][index];
      is_eval = (std::fabs(lambda - lambda_old) < 10e-7) ? true : false;
      lambda = (std::fabs(lambda) >= 10e-6) ? lambda : 0;
      eigenvalues[eig_count] = lambda; lambda_old = lambda;
    }

    index++;
    if (index > 10000) {
      std::cout << "Eigenvalue decomposition can not converge" << endl;
      abort();
    }
  }

  matrix<double> matrix_new;

  if (m_size > 1)
  {
    matrix<double> matrix_tTet;
    matrix_tTet.resize(m_size);
    for (std::uint32_t row = 0; row < m_size; row++)
    {
      matrix_tTet[row].resize(m_size);
      for (std::uint32_t col = 0; col < m_size; col++)
        matrix_tTet[row][col] = (row == col) ? (mat[row][col] - eigenvalues[eig_count]) : mat[row][col];
    }

    std::vector<double> eigenvector;
    jordan_gaussian_transform(matrix_tTet, eigenvector);

    matrix<double> hermitian_matrix;
    get_hermitian_matrix(eigenvector, hermitian_matrix);

    matrix<double> ha_matrix_product;
    matrix_by_matrix(hermitian_matrix, mat, ha_matrix_product);

    matrix<double> inverse_hermitian_matrix;
    get_hermitian_matrix_inverse(eigenvector, inverse_hermitian_matrix);
   
    matrix<double> iha_matrix_product;
    matrix_by_matrix(ha_matrix_product, inverse_hermitian_matrix,
    iha_matrix_product);

    get_reduced_matrix(iha_matrix_product, matrix_new, m_size - 1);
  }

  if (m_size <= 1)
  {
    std::cout << "msize = 1" << endl;
    for (std::uint32_t index = 0; index < eigenvalues.size(); index++)
    {
      double lambda = eigenvalues[index];
      matrix<double> matrix_tTet;
      matrix_tTet.resize(matrix_i.size());

      for (std::uint32_t row = 0; row < matrix_i.size(); row++)
      {
        matrix_tTet[row].resize(matrix_i.size());
        for (std::uint32_t col = 0; col < matrix_i.size(); col++)
          matrix_tTet[row][col] = (row == col) ? (matrix_i[row][col] - lambda) : matrix_i[row][col];
      }

      eigenvectors.resize(matrix_i.size());
      jordan_gaussian_transform(matrix_tTet, eigenvectors[index]);

      double eigsum_sq = 0;
      for (std::uint32_t v = 0; v < eigenvectors[index].size(); v++)
        eigsum_sq += std::pow(eigenvectors[index][v], 2.0);

      for (std::uint32_t v = 0; v < eigenvectors[index].size(); v++)
        eigenvectors[index][v] /= sqrt(eigsum_sq);

      eigenvalues[index] = std::sqrt(eigenvalues[index]);
    }
    return;
  }
  compute_evd(matrix_new, eigenvalues, eigenvectors, eig_count + 1);
}

void matrix_transpose(matrix<double> matrix1, matrix<double>& matrix2)
{
  matrix2.resize(matrix1.size());
  for (std::uint32_t row = 0; row < matrix1.size(); row++)
  {
    matrix2[row].resize(matrix1[row].size());
    for (std::uint32_t col = 0; col < matrix1[row].size(); col++)
      matrix2[row][col] = matrix1[col][row];
  }
}

void get_inverse_diagonal_matrix(matrix<double> mat, matrix<double>& inv_matrix)
{
  inv_matrix.resize(mat.size());
  for (std::uint32_t index = 0; index < mat.size(); index++)
  {
    inv_matrix[index].resize(mat[index].size());
    inv_matrix[index][index] = 1.0 / mat[index][index];
  }
}

void svd(matrix<double> mat, matrix<double>& s, matrix<double>& u,
matrix<double>& v)
{
  // this code belongs to https://www.codeproject.com/Articles/1268576/Singular-Values-Decomposition-SVD-In-Cplusplus11-B#:~:text=Singular%20value%20decomposition%20(Singular%20Value,to%20visualize%20the%20available%20data.
  matrix<double> matrix_t;
  matrix_transpose(mat, matrix_t);

  matrix<double> matrix_product1;
  matrix_by_matrix(mat, matrix_t, matrix_product1);

  matrix<double> matrix_product2;
  matrix_by_matrix(matrix_t, mat, matrix_product2);

  matrix<double> u_1;
  matrix<double> v_1;

  std::vector<double> eigenvalues;
  compute_evd(matrix_product2, eigenvalues, v_1, 0);
  matrix_transpose(v_1, v);
  std::cout << "finish compute_evd" << endl;

  s.resize(mat.size());
  for (std::uint32_t index = 0; index < eigenvalues.size(); index++)
  {
    s[index].resize(eigenvalues.size());
    s[index][index] = eigenvalues[index];
  }

  matrix<double> s_inverse;
  get_inverse_diagonal_matrix(s, s_inverse);

  matrix<double> av_matrix;
  matrix_by_matrix(mat, v, av_matrix);
  matrix_by_matrix(av_matrix, s_inverse, u);
}


void ginv(matrix<double> &mat, matrix<double> &mat_inv)
{
  // get Moore-Penrose generalized inverse of a matrix
  matrix<double> u, v, s;
  svd(mat, s, u, v);

  matrix<double> s_inverse;
  get_inverse_diagonal_matrix(s, s_inverse);

  matrix<double> u_t;
  matrix_transpose(u, u_t);

  matrix<double> vs;
  matrix_by_matrix(v, s_inverse, vs);
  matrix_by_matrix(vs, u_t, mat_inv);

  // check ginv is correct
  matrix<double> identity;
  matrix_by_matrix(mat, mat_inv, identity);
  std::cout << "mat * ginv(mat) = " << identity << endl;
}