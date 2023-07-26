#include <ctime>
// #include <Rcpp.h>
#include "tree.h"
#include "forest.h"
#include <chrono>
// #include "mcmc_loop.h"
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

    longbetModel *model = new longbetModel();

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

    longbetModel *model = new longbetModel();

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

// [[Rcpp::export]]
Rcpp::StringVector r_to_json(double y_mean, Rcpp::XPtr<std::vector<std::vector<tree>>> tree_pnt)
{

    Rcpp::StringVector result(1);
    std::vector<std::vector<tree>> *trees = tree_pnt;
    json j = get_forest_json(*trees, y_mean);
    result[0] = j.dump(4);
    return result;
}

// [[Rcpp::export]]
Rcpp::List json_to_r(Rcpp::StringVector json_string_r)
{

    std::vector<std::string> json_string(json_string_r.size());
    //std::string json_string = json_string_r(0);
    json_string[0] = json_string_r(0);
    double y_mean;

    // Create trees
    vector<vector<tree>> *trees2 = new std::vector<vector<tree>>();

    // Load
    from_json_to_forest(json_string[0], *trees2, y_mean);

    // Define External Pointer
    Rcpp::XPtr<std::vector<std::vector<tree>>> tree_pnt(trees2, true);

    return Rcpp::List::create(Rcpp::Named("model_list") = Rcpp::List::create(Rcpp::Named("tree_pnt") = tree_pnt,
                                                                             Rcpp::Named("y_mean") = y_mean));
}