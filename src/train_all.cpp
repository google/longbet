#include <ctime>
#include <RcppArmadillo.h>
#include "tree.h"
#include <chrono>
// #include "mcmc_loop.h"
#include "X_struct.h"
#include "mcmc_loop.h"
#include "common.h"

using namespace std;
using namespace chrono;
////////////////////////////////////////////////////////////////////////
//                                                                    //
//                                                                    //
//  Full function, support both continuous and categorical variables  //
//                                                                    //
//                                                                    //
////////////////////////////////////////////////////////////////////////

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

// FUNCTION longBet
// preprocesses input received from R
// feeds data into main loop function 'mcmc_loop_longBet'
// returns the list of objects to later become the output in R
// general attributes: y (vector of responses), X (matrix of covariates),
//                     z (vector of treatment assignments),
//                     time (tmpry input for time trend),
//                     num_sweeps, burnin (# of burn-in sweeps),
//                     max_depth (of a tree), n_min (minimum node size),
//                     num_cutpoints (# of cutpoints considered at each split),
//                     no_split_penalty,
//                     mtry (# of variables considered at each split),
//                     p_categorical (# of categorical regressors)
// per forest:         alpha, beta, tau, (BART prior parameters)
//                     num_trees,
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
Rcpp::List longBet_cpp(arma::mat y, arma::mat X, arma::mat X_tau, arma::mat z,
                    arma::mat t_con, arma::mat t_mod,
                    size_t num_sweeps, size_t burnin = 1,
                    size_t max_depth = 1, size_t n_min = 5,
                    size_t num_cutpoints = 1,
                    double no_split_penality = 0.001,
                    size_t mtry_pr = 0, size_t mtry_trt = 0,
                    size_t p_categorical_pr = 0,
                    size_t p_categorical_trt = 0,
                    size_t num_trees_pr = 200,
                    double alpha_pr = 0.95, double beta_pr = 2,
                    double tau_pr = 0.5, double kap_pr = 16, double s_pr = 4,
                    bool pr_scale = false,
                    size_t num_trees_trt = 50,
                    double alpha_trt = 0.25, double beta_trt = 3,
                    double tau_trt = 0.5, double kap_trt = 16, double s_trt = 4,
                    bool trt_scale = false,
                    bool verbose = false, bool parallel = true,
                    bool set_random_seed = false, size_t random_seed = 0,
                    bool sample_weights_flag = true,
                    bool a_scaling = true, bool b_scaling = true,
                    bool split_t_mod = true, bool split_t_con = false,
                    double sig_knl = 1, double lambda_knl = 2)
{
    auto start = system_clock::now();

    size_t N = X.n_rows;

    // number of total variables
    size_t p_pr = X.n_cols;
    size_t p_trt = X_tau.n_cols;
    size_t p_y = y.n_cols;  // dimension of response variables

    // number of continuous variables
    size_t p_continuous_pr = p_pr - p_categorical_pr;
    size_t p_continuous_trt = p_trt - p_categorical_trt;

    // suppose first p_continuous variables are continuous, then categorical

    assert(mtry_pr <= p_pr);
    assert(mtry_trt <= p_trt);

    assert(burnin <= num_sweeps);

    if (mtry_pr == 0)
    {
        mtry_pr = p_pr;
    }

    if (mtry_pr != p_pr)
    {
        COUT << "Sample " << mtry_pr << " out of " << p_pr <<
         " variables when grow each tree." << endl;
    }

    if (mtry_trt == 0)
    {
        mtry_trt = p_trt;
    }

    if (mtry_trt != p_trt)
    {
        COUT << "Sample " << mtry_trt << " out of " << p_trt <<
        " variables when grow each tree." << endl;
    }

    arma::umat Xorder(X.n_rows, X.n_cols);
    matrix<size_t> Xorder_std;
    ini_matrix(Xorder_std, N, p_pr);

    matrix<size_t> torder_mu_std;
    ini_matrix(torder_mu_std, t_con.n_rows, t_con.n_cols);

    matrix<size_t> torder_tau_std;  // post treatment time;
    ini_matrix(torder_tau_std, t_mod.n_rows, t_mod.n_cols);

    arma::umat Xorder_tau(X_tau.n_rows, X_tau.n_cols);
    matrix<size_t> Xorder_tau_std;
    ini_matrix(Xorder_tau_std, N, p_trt);

    // std::vector<double> y_std(N);
    // std::vector<size_t> z_std(N);
    // std::vector<double> time_std(N);
    double y_mean = 0.0;

    Rcpp::NumericMatrix z_std(N, p_y);
    Rcpp::NumericMatrix y_std(N, p_y);
    Rcpp::NumericMatrix X_std(N, p_pr);
    Rcpp::NumericMatrix X_tau_std(N, p_trt);
    Rcpp::NumericMatrix tcon_std(N, t_con.n_cols);
    Rcpp::NumericMatrix tmod_std(N, t_mod.n_cols);

    arma_to_rcpp(y, y_std);
    arma_to_rcpp(z, z_std);
    arma_to_rcpp(t_con, tcon_std);
    arma_to_rcpp(t_mod, tmod_std);
    arma_to_rcpp(X, X_std);
    arma_to_std_ordered(X, Xorder_std);
    arma_to_std_ordered(t_con, torder_mu_std);
    arma_to_std_ordered(t_mod, torder_tau_std);
    arma_to_rcpp(X_tau, X_tau_std);
    arma_to_std_ordered(X_tau, Xorder_tau_std);
    y_mean = compute_mat_mean(y_std);

    ///////////////////////////////////////////////////////////////////
    std::vector<double> sigma_vec(2);  // vector of sigma0, sigma1
    sigma_vec[0] = 1.0;
    sigma_vec[1] = 1.0;

    double bscale0 = -0.5;
    double bscale1 = 0.5;

    std::vector<double> b_vec(2);  // vector of sigma0, sigma1
    b_vec[0] = bscale0;
    b_vec[1] = bscale1;

    std::vector<size_t> num_trees(2); // vector of tree number for each of mu and tau
    num_trees[0] = num_trees_pr;
    num_trees[1] = num_trees_trt;

    std::vector<size_t> n_trt(p_y, 0); // number of treated individuals TODO: remove from here and from constructor as well

    // assuming we have presorted data (treated individuals first, then control group)

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < p_y; j++){
            if (z_std(i, j)== 1) n_trt[j]++;
        }
    }

    double *Xpointer = &X_std[0];
    double *Xpointer_tau = &X_tau_std[0];
    double *ypointer = &y_std[0];
    double *zpointer = &z_std[0];
    double *tpointer_mu = &tcon_std[0];
    double *tpointer_tau = &tmod_std[0];
    // double *Xtestpointer = &Xtest_std[0];

    std::vector<matrix<double>> tauhats_xinfo(num_sweeps);
    std::vector<matrix<double>> muhats_xinfo(num_sweeps);

    for (size_t i = 0; i < num_sweeps; i++)
    { 
        ini_matrix(tauhats_xinfo[i], p_y, N);
        ini_matrix(muhats_xinfo[i], p_y, N);
    }
    // matrix<double> total_fit;
    // ini_matrix(total_fit, N, num_sweeps);

    matrix<double> sigma0_draw_xinfo;
    ini_matrix(sigma0_draw_xinfo, num_trees_trt + num_trees_pr, num_sweeps);

    matrix<double> sigma1_draw_xinfo;
    ini_matrix(sigma1_draw_xinfo, num_trees_trt + num_trees_pr, num_sweeps);

    matrix<double> a_xinfo;
    ini_matrix(a_xinfo, num_sweeps, 1);

    matrix<double> b_xinfo;
    ini_matrix(b_xinfo, num_sweeps, 2);

    matrix<double> beta_xinfo;
    ini_matrix(beta_xinfo, p_y, num_sweeps);

    // // Create trees
    vector<vector<tree>> *trees_pr = new vector<vector<tree>>(num_sweeps);
    for (size_t i = 0; i < num_sweeps; i++)
    {
        (*trees_pr)[i] = vector<tree>(num_trees_pr);
    }

    // // Create trees
    vector<vector<tree>> *trees_trt = new vector<vector<tree>>(num_sweeps);
    for (size_t i = 0; i < num_sweeps; i++)
    {
        (*trees_trt)[i] = vector<tree>(num_trees_trt);
    }
    // define the model for the prognostic term
    longBetModel *model_pr = new longBetModel(kap_pr, s_pr, tau_pr, alpha_pr, beta_pr);
    model_pr->setNoSplitPenality(no_split_penality);

    // define the model for the treatment term
    longBetModel *model_trt = new longBetModel(kap_trt, s_trt, tau_trt, alpha_trt, beta_trt);
    model_trt->setNoSplitPenality(no_split_penality);

    // State settings for the prognostic term
    std::unique_ptr<State> state(new longBetState(Xpointer, Xorder_std, N,
    n_trt, p_pr, p_trt, p_y, num_trees, p_categorical_pr, p_categorical_trt,
    p_continuous_pr, p_continuous_trt, set_random_seed, random_seed, n_min,
    num_cutpoints, parallel, mtry_pr, mtry_trt, Xpointer, num_sweeps,
    sample_weights_flag, ypointer, zpointer, sigma_vec, b_vec, max_depth, y_mean,
    burnin, model_trt->dim_suffstat));

    // initialize X_struct for the prognostic term
    std::vector<double> initial_theta_pr(1, y_mean / (double)num_trees_pr);
    std::unique_ptr<X_struct> x_struct_pr(new X_struct(Xpointer, ypointer, tpointer_mu, N, p_y, Xorder_std, torder_mu_std, p_categorical_pr, p_continuous_pr, &initial_theta_pr, num_trees_pr, sig_knl, lambda_knl));

    // initialize X_struct for the treatment term
    std::vector<double> initial_theta_trt(1, 0);
    std::unique_ptr<X_struct> x_struct_trt(new X_struct(Xpointer_tau, ypointer, tpointer_tau, N, p_y, Xorder_tau_std, torder_tau_std, p_categorical_trt, p_continuous_trt, &initial_theta_trt, num_trees_trt, sig_knl, lambda_knl));

    // mcmc_loop returns tauhat [N x sweeps] matrix
    mcmc_loop_longBet(Xorder_std, Xorder_tau_std, Xpointer, Xpointer_tau, torder_mu_std, torder_tau_std, verbose, sigma0_draw_xinfo, sigma1_draw_xinfo, b_xinfo, a_xinfo, beta_xinfo, *trees_pr, *trees_trt, no_split_penality,
                   state, model_pr, model_trt, x_struct_pr, x_struct_trt, a_scaling, b_scaling, split_t_mod, split_t_con);

    // predict tauhats and muhats
    // cout << "predict " << endl;
    model_pr->predict_std(Xpointer, tpointer_mu, N, p_y, num_sweeps, muhats_xinfo, *trees_pr);
    model_trt->predict_std(Xpointer_tau, tpointer_tau, N, p_y, num_sweeps, tauhats_xinfo, *trees_trt);
    // cout << "finish predict " << endl;

    // R Objects to Return
    Rcpp::NumericMatrix tauhats(N * p_y, num_sweeps);
    Rcpp::NumericMatrix muhats(N * p_y, num_sweeps);
    Rcpp::NumericMatrix sigma0_draws(num_trees_trt + num_trees_pr, num_sweeps);
    Rcpp::NumericMatrix sigma1_draws(num_trees_trt + num_trees_pr, num_sweeps);
    Rcpp::NumericMatrix b_draws(num_sweeps, 2);
    Rcpp::NumericMatrix a_draws(num_sweeps, 1);
    Rcpp::NumericMatrix beta_draws(p_y, num_sweeps);
    Rcpp::XPtr<std::vector<std::vector<tree>>> tree_pnt_pr(trees_pr, true);
    Rcpp::XPtr<std::vector<std::vector<tree>>> tree_pnt_trt(trees_trt, true);

    for (size_t sweep_ind = 0; sweep_ind < num_sweeps; sweep_ind++)
    {
        for (size_t col = 0; col < p_y; col ++)
        {
            for (size_t row = 0; row < N; row ++)
            {
                tauhats(col * N +  row, sweep_ind) = tauhats_xinfo[sweep_ind][row][col];
                muhats(col * N + row, sweep_ind) = muhats_xinfo[sweep_ind][row][col];
            }
        }
    }

    std_to_rcpp(sigma0_draw_xinfo, sigma0_draws);
    std_to_rcpp(sigma1_draw_xinfo, sigma1_draws);
    std_to_rcpp(b_xinfo, b_draws);
    std_to_rcpp(a_xinfo, a_draws);
    std_to_rcpp(beta_xinfo, beta_draws);


    auto end = system_clock::now();

    auto duration = duration_cast<microseconds>(end - start);

    // print out tree structure, for usage of BART package

    std::stringstream treess_pr;
    std::stringstream treess_trt;

    Rcpp::StringVector output_tree_pr(num_sweeps);
    Rcpp::StringVector output_tree_trt(num_sweeps);

    for (size_t i = 0; i < num_sweeps; i++)
    {
        treess_pr.precision(10);
        treess_trt.precision(10);

        treess_pr.str(std::string());
        treess_pr << num_trees_pr << " " << p_pr << endl;

        treess_trt.str(std::string());
        treess_trt << num_trees_trt << " " << p_trt << endl;

        for (size_t t = 0; t < num_trees_pr; t++)
        {
            treess_pr << (*trees_pr)[i][t];
        }

        for (size_t t = 0; t < num_trees_trt; t++)
        {
            treess_trt << (*trees_trt)[i][t];
        }

        output_tree_pr(i) = treess_pr.str();
        output_tree_trt(i) = treess_trt.str();
    }

    // clean memory
    delete model_pr;
    delete model_trt;
    state.reset();
    x_struct_pr.reset();
    x_struct_trt.reset();
    // R Objects to Return
    return Rcpp::List::create(
        Rcpp::Named("tauhats") = tauhats,
        Rcpp::Named("muhats") = muhats,
        Rcpp::Named("sigma0_draws") = sigma0_draws,
        Rcpp::Named("sigma1_draws") = sigma1_draws,
        Rcpp::Named("b_draws") = b_draws,
        Rcpp::Named("a_draws") = a_draws,
        Rcpp::Named("beta_draws") = beta_draws,
        Rcpp::Named("model_list") = Rcpp::List::create(Rcpp::Named("tree_pnt_pr") = tree_pnt_pr,
                                                       Rcpp::Named("tree_pnt_trt") = tree_pnt_trt,
                                                       Rcpp::Named("y_mean") = y_mean),
        Rcpp::Named("treedraws_pr") = output_tree_pr,
        Rcpp::Named("treedraws_trt") = output_tree_trt,
        Rcpp::Named("sdy_use") = NULL,
        Rcpp::Named("sdy") = NULL,
        Rcpp::Named("meany") = NULL,
        Rcpp::Named("tauhats.adjusted") = NULL,
        Rcpp::Named("muhats.adjusted") = NULL,
        Rcpp::Named("model_params") = Rcpp::List::create(Rcpp::Named("num_sweeps") = num_sweeps,
                                                         Rcpp::Named("burnin") = burnin,
                                                         Rcpp::Named("max_depth") = max_depth,
                                                         Rcpp::Named("Nmin") = n_min,
                                                         Rcpp::Named("num_cutpoints") = num_cutpoints,
                                                         Rcpp::Named("alpha_pr") = alpha_pr,
                                                         Rcpp::Named("beta_pr") = beta_pr,
                                                         Rcpp::Named("tau_pr") = tau_pr,
                                                         Rcpp::Named("p_categorical_pr") = p_categorical_pr,
                                                         Rcpp::Named("num_trees_pr") = num_trees_pr,
                                                         Rcpp::Named("alpha_trt") = alpha_trt,
                                                         Rcpp::Named("beta_trt") = beta_trt,
                                                         Rcpp::Named("tau_trt") = tau_trt,
                                                         Rcpp::Named("p_categorical_trt") = p_categorical_trt,
                                                         Rcpp::Named("num_trees_trt") = num_trees_trt),
        Rcpp::Named("input_var_count") = Rcpp::List::create(Rcpp::Named("x_con") = p_pr-1,
                                                            Rcpp::Named("x_mod") = p_trt)

    );
}