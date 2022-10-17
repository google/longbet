#include <ctime>
#include <RcppArmadillo.h>
// #include <Rcpp.h>
// #include <armadillo>
#include "tree.h"
#include "forest.h"
#include <chrono>
// #include "mcmc_loop.h"
#include "utility.h"
#include "json_io.h"
#include "model.h"
#include "rcpp_utility.h"


// [[Rcpp::export]]
Rcpp::List predict(arma::mat X, arma::mat t, 
                        Rcpp::XPtr<std::vector<std::vector<tree>>> tree_pnt)
{

    // Process the prognostic input
    size_t N = X.n_rows;
    size_t p = X.n_cols;
    size_t p_y = t.n_rows;

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

    // Init matrix
    std::vector<double> tr_std(tr_size);
    std::vector<double> te_std(te_size);
    // Rcpp::NumericMatrix res_std(t_size, num_sweeps);
    // Rcpp::NumericMatrix A_diag_std(t_size, num_sweeps);
    // Rcpp::NumericMatrix Sig_diag_std(t_size, num_sweeps);
    
    arma_to_std(t_train, tr_std);
    arma_to_std(t_test, te_std);
    // arma_to_rcpp(res, res_std);
    // arma_to_rcpp(A_diag, A_diag_std);
    // arma_to_rcpp(Sig_diag, Sig_diag_std);
    
    // Calculate Sigma matrix
    matrix<double> Sigma_tr;
    matrix<double> Sigma_te;
    matrix<double> Sigma_tt;
    ini_matrix(Sigma_tr, tr_size, tr_size);
    ini_matrix(Sigma_te, te_size, te_size);
    ini_matrix(Sigma_tt, te_size, tr_size);

    cov_kernel(tr_std, tr_std, sig_knl, lambda_knl, Sigma_tr);
    cov_kernel(te_std, te_std, sig_knl, lambda_knl, Sigma_te);
    cov_kernel(tr_std, te_std, sig_knl, lambda_knl, Sigma_tt);


    // // output
    // Rcpp::NumericMatrix preds(pred_t_size, num_sweeps);

    // for (size_t sweeps = 0; sweeps < num_sweeps; sweeps++)
    // {

    // }

    return Rcpp::List::create();
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