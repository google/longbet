#include <ctime>
// #include <Rcpp.h>
#include "tree.h"
#include "forest.h"
#include <chrono>
#include "mcmc_loop.h"
#include "utility.h"
#include "json_io.h"
#include "model.h"
#include <RcppArmadillo.h>
#include "rcpp_utility.h"


using namespace arma;
// [[Rcpp::export]]
Rcpp::List predict(arma::mat X, arma::mat t, 
                        Rcpp::XPtr<std::vector<std::vector<tree>>> tree_pnt)
{
    // Process the prognostic input
    size_t N = X.n_rows;
    size_t p = X.n_cols;
    size_t p_y = t.n_rows;

    if (N * p * p_y == 0)
    {
        std::cout << "Wrong dimensions " << " N " << N << " p " << p << " p_y " << p_y << endl;
        abort();
    }

    // Init X_std matrix
    Rcpp::NumericMatrix X_std(N, p);
    Rcpp::NumericMatrix t_std(N, p_y);

    arma_to_rcpp(X, X_std);
    arma_to_rcpp(t, t_std);

    double *Xpointer = &X_std[0];
    double *tpointer = &t_std[0];

    // Trees
    std::vector<std::vector<tree>> *trees = tree_pnt;

    // Result Container
    size_t num_sweeps = (*trees).size();
    std::vector<matrix<double>> pred_xinfo(num_sweeps);
    for (size_t i = 0; i < num_sweeps; i++)
    { 
        ini_matrix(pred_xinfo[i], p_y, N);
    }

    longBetModel *model = new longBetModel();

    // Predict
    model->predict_std(Xpointer, tpointer, N, p_y, num_sweeps, pred_xinfo, *trees);

    // Convert back to Rcpp
    Rcpp::NumericMatrix preds(N * p_y, num_sweeps);

    for (size_t sweep_ind = 0; sweep_ind < num_sweeps; sweep_ind++)
    {
        for (size_t col = 0; col < p_y; col ++)
        {
            for (size_t row = 0; row < N; row ++)
            {
                preds(col * N +  row, sweep_ind) = pred_xinfo[sweep_ind][row][col];
            }
        }
    }


    return Rcpp::List::create(Rcpp::Named("preds") = preds);
}

// [[Rcpp::export]]
Rcpp::List predict_beta(arma::mat t_test, arma::mat t_train,
    arma::mat res, arma::mat A_diag, arma::mat Sig_diag,
    double sig_knl, double lambda_knl)
{
    size_t tr_size = t_train.n_rows;
    size_t te_size = t_test.n_rows;
    size_t num_sweeps = res.n_cols;

    if (tr_size * te_size * num_sweeps == 0)
    {
        std::cout << "Wrong dimensions " << " tr_size " << tr_size << " te_size " << te_size << " num_sweeps " << num_sweeps << endl;
        abort();
    }
    // std::cout << "tr_size " << tr_size << endl;
    // std::cout << "t_train " << t_train << endl;
    // std::cout << "t_test " << t_test << endl;

    // Init matrix
    std::vector<double> tr_std(tr_size);
    std::vector<double> te_std(te_size);
    // Rcpp::NumericMatrix res_std(t_size, num_sweeps);
    // Rcpp::NumericMatrix A_diag_std(t_size, num_sweeps);
    // Rcpp::NumericMatrix Sig_diag_std(t_size, num_sweeps);
    
    arma_to_std(t_train, tr_std);
    arma_to_std(t_test, te_std);
    
    // std::cout << "tr_std " << tr_std << endl;
    // std::cout << "te_std " << te_std << endl;
    
    // Calculate Sigma matrix
    matrix<double> Sigma_tr_std;
    matrix<double> Sigma_te_std;
    matrix<double> Sigma_tt_std;
    ini_matrix(Sigma_tr_std, tr_size, tr_size);
    ini_matrix(Sigma_te_std, te_size, te_size);
    ini_matrix(Sigma_tt_std, te_size, tr_size);

    cov_kernel(tr_std, tr_std, sig_knl, lambda_knl, Sigma_tr_std);
    cov_kernel(te_std, te_std, sig_knl, lambda_knl, Sigma_te_std);
    cov_kernel(tr_std, te_std, sig_knl, lambda_knl, Sigma_tt_std);

    longBetModel *model = new longBetModel();

    // cout << "tr_size" << tr_size << endl;

    // mat A(tr_size, tr_size, fill::zeros);
    // mat Sig(tr_size, tr_size, fill::zeros);
    // mat res_vec(tr_size, 1, fill::zeros);
    // cout << "A " << A << endl;

    std::vector<double> a_vec(tr_size);
    std::vector<double> sig_vec(tr_size);
    std::vector<double> res_vec(tr_size);

    // output
    std::random_device rd;
    std::mt19937 gen = std::mt19937(rd());;
    std::vector<double> beta_std(te_size);
    Rcpp::NumericMatrix beta(te_size, num_sweeps);

    for (size_t sweeps = 0; sweeps < num_sweeps; sweeps++)
    {
        for (size_t i = 0; i < tr_size; i++)
        {
            a_vec[i] = A_diag(i, sweeps);
            sig_vec[i] = Sig_diag(i, sweeps);
            res_vec[i] = res(i, sweeps);
        }

        model->predict_beta(beta_std, res_vec, a_vec, sig_vec, Sigma_tr_std, Sigma_te_std, Sigma_tt_std, gen);

        for (size_t i = 0; i < te_size; i++)
        {
            beta(i, sweeps) = beta_std[i];
        }
    }
    return Rcpp::List::create(Rcpp::Named("beta") = beta);
}

// // [[Rcpp::export]]
// Rcpp::List predict_beta(arma::mat t, arma::mat train_t, arma::mat beta,
//     arma::mat time_residuals, arma::mat time_diag_A, arma::mat time_diag_Sig,
//     double sig_knl, double lambda_knl)
// {
//     cout << "predict beta time t = " << t << endl;
//     size_t N_t = t.n_rows;
//     size_t t_size = beta.n_rows;
//     size_t num_sweeps = beta.n_cols;

//     std::vector<double> t_std(N_t);
//     std::vector<double> train_t_std(t_size);

//     for (size_t i = 0; i < N_t; i++) t_std[i] = t(i, 0);
//     for (size_t i = 0; i < t_size; i++) train_t_std[i] = train_t(i, 0);


//     matrix<double> beta_std;
//     matrix<double> residuals;
//     matrix<double> diag_A;
//     matrix<double> diag_Sig;

//     ini_matrix(beta_std, t_size, num_sweeps);
//     ini_matrix(residuals, t_size, num_sweeps);
//     ini_matrix(diag_A, t_size, num_sweeps);
//     ini_matrix(diag_Sig, t_size, num_sweeps);

//     for (size_t i = 0; i < num_sweeps; i++)
//     {
//         for (size_t j = 0; j < t_size; j++)
//         {
//             beta_std[i][j] = beta(j, i);
//             residuals[i][j] = time_residuals(j, i);
//             diag_A[i][j] = time_diag_A(j, i);
//             diag_Sig[i][j] = time_diag_Sig(j, i);
//         }
//     }

//     // ini output
//     matrix<double> beta_draws;
//     ini_matrix(beta_draws, N_t, num_sweeps);

//     // check if extrapolation is needed
//     // initialize inverse matrix
//     size_t idx;
//     bool extrapolate = false;
//     for (size_t i = 0; i < N_t; i++){
//         idx = std::find(train_t_std.begin(), train_t_std.end(), t_std[i]) -
//             train_t_std.begin();
//         if ((idx == t_size) & (t_std[i] != train_t_std[t_size - 1]))
//         {
//             extrapolate = true;
//             break;
//         }
//     }

//     matrix<double> cov_kernel;
//     arma::mat res_cov(t_size, t_size, arma::fill::zeros);
//     arma::mat res_cov_inv(t_size, t_size, arma::fill::zeros);
//     arma::mat t_cov(t_size, 1);
//     arma::mat res(t_size, 1);
//     double mu;
//     double var;
//     // arma::mat mu(1, 1);
//     // arma::mat var(1, 1);

//     if (extrapolate)
//     {
//         std::random_device rd;
//         std::mt19937 gen = std::mt19937(rd());
//         std::normal_distribution<double> normal_samp(0.0, 1.0);
//         // initialize covariance kernel
//         double sigma2 = pow(sig_knl, 2);
//         double lambda2 = pow(lambda_knl, 2);

//         ini_matrix(cov_kernel, t_size, t_size);
//         double cov_kernel_diag = squared_exponential(train_t_std[0], train_t_std[0],
//         sigma2, lambda2);

//         for (size_t i = 0; i < t_size; i++)
//         {
//             // calculate diagonal element
//             cov_kernel[i][i] = cov_kernel_diag;
//             for (size_t j = 0; j < i; j++)
//             {
//                 cov_kernel[i][j] = squared_exponential(train_t_std[i],
//                     train_t_std[j], sigma2, lambda2);
//                 cov_kernel[j][i] = cov_kernel[i][j];
//             }
//         }

//         for (size_t sweeps = 0; sweeps < num_sweeps; sweeps++)
//         {
//             // get residuals
//             for (size_t i = 0; i < t_size; i++)
//             {
//                 res(i, 0) = residuals[sweeps][i];
//             }
//             cout << "res = " << res << endl;

//             // calculate residual covariance
//             res_cov.zeros();
//             res_cov_inv.zeros();
//             for (size_t i = 0; i < t_size; i++)
//             {
//                 res_cov(i, i) += diag_Sig[sweeps][i];
//                 for (size_t j = 0; j < t_size; j++)
//                 {
//                     res_cov(i, j) += diag_A[sweeps][i] * cov_kernel[i][j] *
//                         diag_A[sweeps][j];
//                 }
//             }
//             cout << "res_cov = " << res_cov << endl;
//             // get inverse of residual covaraince
//             res_cov_inv = inv_sympd(res_cov);
//             cout << "res_cov_inv = " << res_cov_inv << endl;

//             for (size_t i = 0; i < N_t; i++)
//             {
//                 idx = std::find(train_t_std.begin(), train_t_std.end(), t_std[i]) - train_t_std.begin();
//                 if ((idx < t_size) | (t_std[i] == train_t_std[t_size - 1])){
//                     // t exist in training set
//                     beta_draws[sweeps][i] = beta_std[sweeps][idx];
//                 } else {
//                     cout << "extrapolate t = " << t_std[i] << endl;
//                     // t needs to be extrapolated
//                     // get covariance kernel for t v.s. train_t
//                     for (size_t j = 0; j < t_size; j++)
//                     {
//                         t_cov(j, 0) = diag_A[sweeps][j] * squared_exponential(
//                             t_std[i], train_t_std[j], sigma2, lambda2);

//                     }
//                     cout << "t_cov = " << t_cov << endl;

//                     // get mean and variance
//                     mu = as_scalar(trans(t_cov) * res_cov_inv * res);
//                     cout << "mu = " << mu << endl;
//                     var = cov_kernel_diag - as_scalar(trans(t_cov) * res_cov_inv * t_cov);

//                     // TODO: add check on dimension, both mu and var should be 1 by 1 matrix (scalar).
//                     cout << "mu = " << mu << ", var = " << var << endl;
//                     beta_draws[sweeps][i] = mu + pow(var, 0.5) * normal_samp(gen);
//                 }
//             }
//         }
//     } else {
//         for (size_t i = 0; i < N_t; i++)
//         {
//             idx = std::find(train_t_std.begin(), train_t_std.end(), t_std[i]) -
//                 train_t_std.begin();
//             for (size_t sweeps = 0; sweeps < num_sweeps; sweeps++)
//             {
//                 beta_draws[sweeps][i] = beta_std[sweeps][idx];
//             }
//         }
//     }

//     Rcpp::NumericMatrix beta_draws_rcpp(N_t, num_sweeps);
//     for (size_t i = 0; i < num_sweeps; i++)
//     {
//         for (size_t j = 0; j < N_t; j++)
//         {
//             beta_draws_rcpp(j, i) = beta_draws[i][j];
//         }
//     }

//     return Rcpp::List::create(Rcpp::Named("beta_draws") = beta_draws_rcpp);
// }