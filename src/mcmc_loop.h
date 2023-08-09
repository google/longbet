#include <ctime>
#include "tree.h"
#include "forest.h"
#include <chrono>
#include "model.h"
#include "state.h"
#include "cdf.h"
#include "X_struct.h"
//#include "MH.h"

void mcmc_loop_longbet( 
  std::unique_ptr<split_info> &split_ps,  
  std::unique_ptr<split_info> &split_trt,
  bool verbose,
  matrix<double> &sigma0_draw_xinfo,
  matrix<double> &sigma1_draw_xinfo,
  matrix<double> &b_xinfo,
  matrix<double> &a_xinfo,
  matrix<double> &beta_info,
  matrix<double> &beta_xinfo,
  vector<vector<tree>> &trees_ps,
  vector<vector<tree>> &trees_trt,
  double no_split_penality,
  std::unique_ptr<State> &state,
  longbetModel *model_ps,
  longbetModel *model_trt,
  std::unique_ptr<X_struct> &x_struct_ps,
  std::unique_ptr<X_struct> &x_struct_trt,
  bool a_scaling,
  bool b_scaling,
  bool split_time_ps,
  bool split_time_trt,
  matrix<double> &resid_info,
  matrix<double> &A_diag_info,
  matrix<double> &Sig_diag_info
  );