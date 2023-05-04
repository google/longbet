#include <rcpp_utility.h>

void rcpp_to_std2(arma::mat y, arma::mat X, arma::mat Xtest, Rcpp::NumericMatrix &y_std, double &y_mean, Rcpp::NumericMatrix &X_std, Rcpp::NumericMatrix &Xtest_std, matrix<size_t> &Xorder_std)
{
    // The goal of this function is to convert RCPP object to std objects

    // TODO: Refactor code so for loops are self contained functions
    // TODO: Why RCPP and not std?
    // TODO: inefficient Need Replacement?

    size_t N = X.n_rows;
    size_t p = X.n_cols;
    size_t N_test = Xtest.n_rows;
    size_t p_y = y.n_cols;

    // Create y_std
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < p_y; j++){
            y_std(i, j) = y(i, j);
            y_mean = y_mean + y_std(i, j);
        }
    }
    y_mean = y_mean / (double)N / (double) p_y;

    // X_std
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < p; j++)
        {
            X_std(i, j) = X(i, j);
        }
    }

    //X_std_test
    for (size_t i = 0; i < N_test; i++)
    {
        for (size_t j = 0; j < p; j++)
        {
            Xtest_std(i, j) = Xtest(i, j);
        }
    }

    // Create Xorder
    // Order
    arma::umat Xorder(X.n_rows, X.n_cols);
    for (size_t i = 0; i < X.n_cols; i++)
    {
        Xorder.col(i) = arma::sort_index(X.col(i));
    }
    // Create
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < p; j++)
        {
            Xorder_std[j][i] = Xorder(i, j);
        }
    }

    return;
}

// FUNCTION arma_to_std (instance 1)
// transfers data from an armadillo matrix object (column 0) to an std vector object
void arma_to_std(const arma::mat &matrix_in, std::vector<double> &vector_out)
{
    size_t dim = matrix_in.n_rows;

    for (size_t i = 0; i < dim; i++)
    {
        vector_out[i] = matrix_in(i, 0);
    }

    return;
}

// transfers data from an armadillo matrix object (column 0) to an std vector object
void arma_to_std(const arma::mat &matrix_in, std::vector<size_t> &vector_out)
{
    size_t dim = matrix_in.n_rows;

    for (size_t i = 0; i < dim; i++)
    {
        vector_out[i] = (size_t)matrix_in(i, 0);
    }

    return;
}

// FUNCTION arma_to_rcpp (instance 1)                    ?? Rcpp matrix or std matrix ??
// transfers data from an armadillo matrix object to an Rcpp matrix object
void arma_to_rcpp(const arma::mat &matrix_in, Rcpp::NumericMatrix &matrix_out)
{
    size_t dim_x = matrix_in.n_rows;
    size_t dim_y = matrix_in.n_cols;
    // cout << "mat dim " << dim_x << " " << matrix_in.n_cols << " first value " << matrix_in(0,0) <<endl;

    for (size_t i = 0; i < dim_x; i++)
    {
        for (size_t j = 0; j < dim_y; j++)
        {
            matrix_out(i, j) = matrix_in(i, j);
        }
    }

    return;
}

// FUNCTION arma_to_std_ordered
// transfers data from an armadillo matrix object to an std matrix object with indeces [carries the pre-sorted features]
void arma_to_std_ordered(const arma::mat &matrix_in, matrix<size_t> &matrix_ordered_std)
{
    size_t dim_x = matrix_in.n_rows;
    size_t dim_y = matrix_in.n_cols;

    arma::umat matrix_ordered(dim_x, dim_y);
    for (size_t i = 0; i < dim_y; i++)
    {
        matrix_ordered.col(i) = arma::sort_index(matrix_in.col(i));
    }

    for (size_t i = 0; i < dim_x; i++)
    {
        for (size_t j = 0; j < dim_y; j++)
        {
            matrix_ordered_std[j][i] = matrix_ordered(i, j);
        }
    }

    return;
}

// FUNCTION std_to_Rcpp
// transfers data from an std matrix object to an Rcpp NumericMatrix object
void std_to_rcpp(const matrix<double> &matrix_in, Rcpp::NumericMatrix &matrix_out)
{
    size_t dim_x = matrix_in.size();
    size_t dim_y = matrix_in[0].size();
    for (size_t i = 0; i < dim_y; i++)
    {
        for (size_t j = 0; j < dim_x; j++)
        {
            matrix_out(i, j) = matrix_in[j][i];
        }
    }

    return;
}

double compute_vec_mean(const std::vector<double> &vec)
{
    double mean = 0;
    int length = vec.size();

    for (size_t i = 0; i < length; i++)
    {
        mean = mean + vec[i];
    }
    mean = mean / (double)length;
    return mean;
}

double compute_mat_mean(const Rcpp::NumericMatrix &matrix)
{
    double mean = 0;

    for (size_t i = 0; i < matrix.nrow(); i++)
    {
        for (size_t j = 0; j < matrix.ncol(); j++){
            mean = mean + matrix(i, j);
        }
    }
    mean = mean / (double) matrix.nrow() / (double) matrix.ncol();
    return mean;
}

void rcpp_to_std2(arma::mat X, arma::mat Xtest, Rcpp::NumericMatrix &X_std, Rcpp::NumericMatrix &Xtest_std, matrix<size_t> &Xorder_std)
{
    // The goal of this function is to convert RCPP object to std objects

    // TODO: Refactor code so for loops are self contained functions
    // TODO: Why RCPP and not std?
    // TODO: inefficient Need Replacement?

    size_t N = X.n_rows;
    size_t p = X.n_cols;
    size_t N_test = Xtest.n_rows;

    // X_std
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < p; j++)
        {
            X_std(i, j) = X(i, j);
        }
    }

    //X_std_test
    for (size_t i = 0; i < N_test; i++)
    {
        for (size_t j = 0; j < p; j++)
        {
            Xtest_std(i, j) = Xtest(i, j);
        }
    }

    // Create Xorder
    // Order
    arma::umat Xorder(X.n_rows, X.n_cols);
    for (size_t i = 0; i < X.n_cols; i++)
    {
        Xorder.col(i) = arma::sort_index(X.col(i));
    }
    // Create
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < p; j++)
        {
            Xorder_std[j][i] = Xorder(i, j);
        }
    }

    return;
}

void std_to_arma(matrix<double> &mat_std, arma::mat &mat_arma)
{
    mat_arma.resize(mat_std.size(), mat_std[0].size());
    for (size_t i = 0; i < mat_std.size(); i++)
    {
        for (size_t j = 0; j < mat_std[0].size(); j++)
        {
            mat_arma(i, j) = mat_std[i][j];
        }
    }
}