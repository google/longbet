#include <ctime>
#include "tree.h"
#include "forest.h"
#include <chrono>
#include "model.h"
#include "state.h"
#include "cdf.h"
#include "X_struct.h"
//#include "MH.h"

void mcmc_loop_longBet(matrix<size_t> &Xorder_std, matrix<size_t> &Xorder_tau_std,
                    const double *X_std, const double *X_tau_std,
                    matrix<size_t> &torder_mu_std,
                    matrix<size_t> &torder_tau_std,
                    bool verbose,
                    matrix<double> &sigma0_draw_xinfo,
                    matrix<double> &sigma1_draw_xinfo,
                    matrix<double> &b_xinfo,
                    matrix<double> &a_xinfo,
                    matrix<double> &beta_xinfo,
                    vector<vector<tree>> &trees_ps,
                    vector<vector<tree>> &trees_trt,
                    double no_split_penality,
                    std::unique_ptr<State> &state,
                    //std::unique_ptr<State> &state_trt,
                    longBetModel *model_ps,
                    longBetModel *model_trt,
                    std::unique_ptr<X_struct> &x_struct_ps,
                    std::unique_ptr<X_struct> &x_struct_trt,
                    bool a_scaling,
                    bool b_scaling,
                    bool split_t_mod,
                    bool split_t_con,
                    matrix<double> &resid_info,
                    matrix<double> &A_diag_info,
                    matrix<double> &Sig_diag_info
                    );