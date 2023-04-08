#include "tree.h"
#include "model.h"
// #include <armadillo>
#include <cstddef>
#include <memory>
#include <numeric>

using namespace arma;

//////////////////////////////////////////////////////////////////////////////////////
//
//
//  XBCF Model
//
//
//////////////////////////////////////////////////////////////////////////////////////

// adds residual to suff stats
// called from calcSuffStat_categorical, calcSuffStat_continuous in tree.cpp
// void longBetModel::incSuffStat(std::unique_ptr<State> &state,
// size_t index_next_obs, matrix<double> &suffstats)
// {
//   // TODO: reconstruct this for n*t y matrix
//   double gp;
//   double resid;
//   for (size_t j = 0; j < state->p_y; j++){
//     gp = state->z[index_next_obs][j];  // treatment group
//     resid = *(state->y_std + state->n_y * j + index_next_obs) -
//     state->a * state->mu_fit[index_next_obs][j] - state->b_vec[gp] *
//     state->beta_t[j] * state->tau_fit[index_next_obs][j];

//     if (state->fl == 0)  // suff stat for prognostic trees
//     {
//       if (state->z[index_next_obs][j] == 1)
//       {
//         suffstats[j][1] += resid / state->a;
//         suffstats[j][3] += 1;
//       } else {
//         suffstats[j][0] += resid / state->a;
//         suffstats[j][2] += 1;
//       }
//     } else {  // suff stat for treatment trees
//       if (state->z[index_next_obs][j] == 1)
//       {
//         // beta_t^2 * r / b / beta_t = beta_t * r / b
//         suffstats[j][1] += state->beta_t[j] * resid / state->b_vec[1];
//         suffstats[j][3] += pow(state->beta_t[j], 2);
//       } else {
//         suffstats[j][0] += state->beta_t[j] * resid / state->b_vec[0];
//         suffstats[j][2] += pow(state->beta_t[j], 2);
//       }
//     }
//   }
// }

void longBetModel::incSuffStat(std::unique_ptr<State> &state,
size_t index_next_obs, size_t index_next_t, std::vector<double> &suffstats)
{
  double gp = *(state->z + index_next_t * state->n_y + index_next_obs);
  double resid = *(state->y_std + state->n_y * index_next_t + index_next_obs) -
  state->a * state->mu_fit[index_next_obs][index_next_t] - state->b_vec[gp] *
  state->beta_t[index_next_t] * state->tau_fit[index_next_obs][index_next_t];

  if (state->fl == 0)  // suff stat for prognostic trees
  {
    if (gp == 1)
    {
      suffstats[1] += resid / state->a;
      suffstats[3] += 1;
    } else {
      suffstats[0] += resid / state->a;
      suffstats[2] += 1;
    }
  } else {  // suff stat for treatment trees
  // cout << "resid = " << resid << endl;
    if (gp == 1)
    {
      // beta_t^2 * r / b / beta_t = beta_t * r / b
      suffstats[1] += state->beta_t[index_next_t] * resid / state->b_vec[1];
      suffstats[3] += pow(state->beta_t[index_next_t], 2);
    } else {
      suffstats[0] += state->beta_t[index_next_t] * resid / state->b_vec[0];
      suffstats[2] += pow(state->beta_t[index_next_t], 2);
    }
  }

}

// samples leaf parameter
// called from GFR in tree.cpp
void longBetModel::samplePars(std::unique_ptr<State> &state,
std::vector<double> &suff_stat, std::vector<double> &theta_vector,
double &prob_leaf)
{
  std::normal_distribution<double> normal_samp(0.0, 1.0);
  double s0 = 0;
  double s1 = 0;

  if (state->fl == 0)  // no sum of beta_t sufficient for prognostic trees
  {
    s0 = pow(state->a / state->sigma_vec[0], 2);
    s1 = pow(state->a / state->sigma_vec[1], 2);
  } else {
    s0 = pow(state->b_vec[0] / state->sigma_vec[0], 2);
    s1 = pow(state->b_vec[1] / state->sigma_vec[1], 2);
  }

  double denominator = suff_stat[2] * s0 + suff_stat[3] * s1  + 1 / tau;
  double m1 = (suff_stat[0] * s0 + suff_stat[1] * s1) / denominator;
  double v1 = 1 / denominator;

  // sample leaf parameter
  theta_vector[0] = m1 + sqrt(v1) * normal_samp(state->gen);
  // cout << "suff_stat = " << suff_stat << ", s0 = " << s0 << ", s1 = " << s1 << endl;
  // cout << "m1 = " << m1 << ", v1 = " << v1 << endl;
  // cout << "theta = " << theta_vector[0] << endl;

  // also update probability of leaf parameters
  prob_leaf = 1.0;
  if (isnan(theta_vector[0])) {
    cout << "theta is nan, m1 = " << m1 << ", v1 = " << v1 << endl;
    cout << "suff_stat = " << suff_stat << ", s0 = " << s0 << ", s1 = " << s1 << endl;
    abort();
  }
}

// updates sigmas (new)
void longBetModel::draw_sigma(std::unique_ptr<State> &state, size_t ind)
{
  double m1 = 0;
  double v1 = 0;
  double sigma;
  double squared_resid = 0;
  double N_trt = std::accumulate(state->n_trt.begin(), state->n_trt.end(), 0.0);

  if (ind == 0){  // update sigma0
    for (size_t i = 0; i < state->p_y; i++){
      squared_resid += sum_squared(state->full_residual_ctrl[i]);
    }
    m1 = (state->n_y * state->p_y - N_trt + kap) / 2.0;
    v1 =  2.0 / (squared_resid + s);
  } else {
    for (size_t i = 0; i < state->p_y; i++){
      squared_resid += sum_squared(state->full_residual_trt[i]);
    }
    m1 = (N_trt + kap) / 2.0;
    v1 = 2.0 / (squared_resid + s);
  }

  // computing both sigmas here due to structural complexity of splitting them
  std::gamma_distribution<double> gamma_samp(m1, v1);
  sigma = 1.0 / sqrt(gamma_samp(state->gen));

  // update the corresponding value in the state object
  state->update_sigma(sigma, ind);
}

// initializes root suffstats
// called from mcmc_loop_longBet in mcmc_loop.cpp
void longBetModel::initialize_root_suffstat(std::unique_ptr<State> &state,
std::vector<double> &suff_stat)
{
  // ini_matrix(suff_stat, 4, state->p_y);
  // for (size_t i = 0; i < state->p_y; i++)
  //   std::fill(suff_stat[i].begin(), suff_stat[i].end(), 0.0);
  suff_stat.resize(4);
  std::fill(suff_stat.begin(), suff_stat.end(), 0.0);

  for (size_t i = 0; i < state->n_y; i++)
  {
    for (size_t j = 0; j < state->p_y; j++)
    {
      incSuffStat(state, i, j, suff_stat);
    }
  }
}

// updates node suffstats for the split
// called from split_xorder_std_continuous, split_xorder_std_categorical in tree.cpp
// it is executed after suffstats for the node has been initialized by suff_stats_ini [defined in tree.h]
void longBetModel::updateNodeSuffStat(std::vector<double> &suff_stat, std::unique_ptr<State> &state, matrix<size_t> &Xorder_std, matrix<size_t> &torder_std, size_t &split_var, size_t row_ind)
{
  for (size_t i = 0; i < torder_std[0].size(); i++){
    incSuffStat(state, Xorder_std[split_var][row_ind], torder_std[0][i], suff_stat);
  }
}

// updates the other side node's side suffstats for the split
// called from split_xorder_std_continuous, split_xorder_std_categorical in tree.cpp
void longBetModel::calculateOtherSideSuffStat(std::vector<double> &parent_suff_stat, std::vector<double> &lchild_suff_stat, std::vector<double> &rchild_suff_stat, bool &compute_left_side)
{

  // in function split_xorder_std_categorical, for efficiency, the function only calculates suff stat of ONE child
  // this function calculate the other side based on parent and the other child

  if (compute_left_side)
  {
    rchild_suff_stat = parent_suff_stat - lchild_suff_stat;
  }
  else
  {
    lchild_suff_stat = parent_suff_stat - rchild_suff_stat;
  }
}

// updates partial residual for the next tree to fit
// called from mcmc_loop_xbcf in xbcf_mcmc_loop.cpp
void longBetModel::state_sweep(size_t tree_ind, matrix<double> &fit, std::unique_ptr<X_struct> &x_struct) const
{
  matrix<double> mu_ft;
  ini_matrix(mu_ft, fit[0].size(), fit.size());
  for (size_t i = 0; i < fit.size(); i++)
  {
    for (size_t j = 0; j < fit[i].size(); j++){
      // fit[i][j] += (*(x_struct->data_pointers[tree_ind][i * fit[i].size() + j]))[0];
      mu_ft[i][j] = (*(x_struct->data_pointers[tree_ind][i * fit[i].size() + j]))[0];
      fit[i][j] += mu_ft[i][j];
    }
  }

  for (size_t i = 0; i < fit.size(); i++){
    // cout << "mu " << i << " = " << mu_ft[i] << endl;
  }
}

// computes likelihood of a split
// called from GFR in tree.cpp
double longBetModel::likelihood(std::vector<double> &temp_suff_stat,
std::vector<double> &suff_stat_all, bool left_side,
bool no_split, std::unique_ptr<State> &state) const
{
  // helper variables
  double s0 = 0;
  double s1 = 0;
  double denominator = 1;   // (1 + tau * precision_squared)
  double s_psi_squared = 0;  // (residual * precision_squared)^2

  if (state->fl == 0)  // no sum of beta_t sufficient for prognostic trees
  {
    s0 = pow(state->a, 2) / pow(state->sigma_vec[0], 2);
    s1 = pow(state->a, 2) / pow(state->sigma_vec[1], 2);
  } else {
    s0 = pow(state->b_vec[0], 2) / pow(state->sigma_vec[0], 2);
    s1 = pow(state->b_vec[1], 2) / pow(state->sigma_vec[1], 2);
  }

  if (no_split)
  {
    denominator = 1 + tau  * (suff_stat_all[2] * s0 + suff_stat_all[3] * s1);
    s_psi_squared = suff_stat_all[0] * s0 + suff_stat_all[1] * s1;
  } else {
    if (left_side)
    {
      denominator = 1 + tau * (temp_suff_stat[2] * s0 + temp_suff_stat[3]*s1);
      s_psi_squared = temp_suff_stat[0] * s0 + temp_suff_stat[1] * s1;
    } else {
      denominator = 1 + tau * ((suff_stat_all[2] - temp_suff_stat[2])
      * s0 + (suff_stat_all[3] - temp_suff_stat[3]) * s1);
      s_psi_squared = (suff_stat_all[0] - temp_suff_stat[0]) * s0 +
      (suff_stat_all[1] - temp_suff_stat[1]) * s1;
    }
  }

  return 0.5 * log(1 / denominator) +
  0.5 * pow(s_psi_squared, 2) * tau / denominator;
}

// makes a prediction for treatment effect on the given Xtestpointer data
void longBetModel::predict_std(const double *Xtestpointer, const double *tpointer, size_t N_test, size_t p, size_t num_sweeps, std::vector<matrix<double>> &yhats_test_xinfo, vector<vector<tree>> &trees)
{
  std::vector<double> output(this->dim_theta, 0.0);
  for (size_t sweeps = 0; sweeps < num_sweeps; sweeps++)
  {
    for (size_t data_ind = 0; data_ind < N_test; data_ind++)
    {
      for (size_t time_ind = 0; time_ind < p; time_ind++)
      {
        for (size_t tree_ind  = 0; tree_ind < trees[0].size(); tree_ind++){
          // cout << "data = " << data_ind << ", time = " << time_ind << ", tree = " << tree_ind << endl;
          getThetaForObs_Outsample(output, trees[sweeps][tree_ind],
          data_ind, time_ind, Xtestpointer, tpointer, N_test, p);

          yhats_test_xinfo[sweeps][data_ind][time_ind] += output[0];
        }
      }
    }
  }
}

// updates parameter a
// called from mcmc_loop_xbcf in xbcf_mcmc_loop.cpp
void longBetModel::update_a_value(std::unique_ptr<State> &state)
{
  std::normal_distribution<double> normal_samp(0.0, 1.0);

  double mu2sum_ctrl = 0;
  double mu2sum_trt = 0;
  double muressum_ctrl = 0;
  double muressum_trt = 0;
  double s0 = pow(state->sigma_vec[0], 2);
  double s1 = pow(state->sigma_vec[1], 2);

  // compute the residual y-b*beta_t*tau(x)
  for (size_t i = 0; i < state->n_y; i++)
  {
    for (size_t j = 0; j < state->p_y; j++){
      if ((*(state->z + j * state->n_y + i)) == 1)
      {
        state->residual[i][j] = *(state->y_std + state->n_y * j + i) -
        state->b_vec[1] * state->beta_t[j] * state->tau_fit[i][j];

        mu2sum_trt += state->mu_fit[i][j] * state->mu_fit[i][j];
        muressum_trt += state->mu_fit[i][j] * state->residual[i][j];
      } else {
        state->residual[i][j] = *(state->y_std + state->n_y * j + i) -
        state->b_vec[0] * state->beta_t[j] * state->tau_fit[i][j];

        mu2sum_ctrl += state->mu_fit[i][j] * state->mu_fit[i][j];
        muressum_ctrl += state->mu_fit[i][j] * state->residual[i][j];
      }
    }
  }

  // update parameters
  double denominator = mu2sum_ctrl / s0 + mu2sum_trt / s1 + 1;
  double m1 = (muressum_ctrl / s0 + muressum_trt / s1) / denominator;
  double v1 = 1 / denominator;

  // sample a
  state->a = m1 + sqrt(v1) * normal_samp(state->gen);
}

// updates parameters b0, b1
// called from mcmc_loop_xbcf in xbcf_mcmc_loop.cpp
void longBetModel::update_b_values(std::unique_ptr<State> &state)
{
  std::normal_distribution<double> normal_samp(0.0, 1.0);

  double tau2sum_ctrl = 0;
  double tau2sum_trt = 0;
  double tauressum_ctrl = 0;
  double tauressum_trt = 0;
  double s0 = pow(state->sigma_vec[0], 2);
  double s1 = pow(state->sigma_vec[1], 2);


  for (size_t i = 0; i < state->n_y; i++)
  {
    for (size_t j = 0; j < state->p_y; j++){
      state->residual[i][j] = *(state->y_std + state->n_y * j + i) -
      state->a * state->mu_fit[i][j];

      if (*(state->z + j * state->n_y + i) == 1)
      {
        tau2sum_trt += pow(state->tau_fit[i][j] * state->beta_t[j], 2);
        tauressum_trt += state->beta_t[j] * state->tau_fit[i][j] *
        state->residual[i][j];
      } else {
        tau2sum_ctrl += pow(state->tau_fit[i][j] * state->beta_t[j], 2);
        tauressum_ctrl += state->beta_t[j] * state->tau_fit[i][j] *
        state->residual[i][j];
      }
    }
  }

  // update parameters
  double v0 = 1 / (2 + tau2sum_ctrl / s0);
  double v1 = 1 / (2 + tau2sum_trt / s1);

  double m0 = v0 * (tauressum_ctrl) / s0;
  double m1 = v1 * (tauressum_trt) / s1;

  // sample b0, b1
  double b0 = m0 + sqrt(v0) * normal_samp(state->gen);
  double b1 = m1 + sqrt(v1) * normal_samp(state->gen);

  state->b_vec[1] = b1;
  state->b_vec[0] = b0;
}

void longBetModel::update_time_coef(std::unique_ptr<State> &state, std::unique_ptr<X_struct> &x_struct,
  matrix<size_t> &torder_std, std::vector<double> &resid, std::vector<double> &diag, std::vector<double> &sig, std::vector<double> &beta)
{
  // get total number of time
  double t_size = x_struct->t_values.size();
  double n = state->n_y;  // n obs per period. TODO: need update
  std::vector<double> res_ctrl(t_size, 0);  // residuals
  std::vector<double> res_trt(t_size, 0);

  // diagonal element of matrix A: sigma_{z_i}^{-1} * b_{z_i} * tau_i
  std::vector<double> diag_ctrl(t_size, 0);
  std::vector<double> diag_trt(t_size, 0);
  // std::vector<double> diag(t_size, 0);

  double sig02 = pow(state->sigma_vec[0], 2);
  double sig12 = pow(state->sigma_vec[1], 2);
  // vec sig(t_size, fill::zeros);


  std::vector<size_t> idx(state->p_y);  // keep track of t-values
  size_t t_idx;
  size_t counts = 0;
  const double *z_pointer;
  const double *y_pointer;

  // if(t_size <= 1){
  //   cout << "unique t values need to be greater than 1" << endl;
  //   throw;
  // }

  for (size_t i = 0; i < t_size; i++)
  {
    for (size_t j = 0; j < x_struct->t_counts[i]; j++)
    {
      t_idx = torder_std[0][counts];
      counts++;
      idx[t_idx] = i;
      z_pointer = state->z + state->n_y * t_idx;
      y_pointer = state->y_std + state->n_y * t_idx;
      for (size_t k = 0; k < state->n_y; k++)
      {
        if (*(z_pointer + k) == 0)
        {
          res_ctrl[i] += *(y_pointer + k) - state->a * state->mu_fit[k][t_idx];
          diag_ctrl[i] += state->tau_fit[k][t_idx];
          sig[i] += sig02;
        } else {
          res_trt[i] += *(y_pointer + k) - state->a * state->mu_fit[k][t_idx];
          diag_trt[i] += state->tau_fit[k][t_idx];
          sig[i] += sig12;
        }
      }
    }
  }

  for (size_t i = 0; i < t_size; i++){
    resid[i] = (res_trt[i] + res_ctrl[i]) / n / x_struct->t_counts[i];
    diag[i] = (state->b_vec[1] * diag_trt[i] + state->b_vec[0] * diag_ctrl[i])/n;
    sig[i] = 1 / (sig[i] / pow(n, 2) / x_struct->t_counts[i]);
  }
  // solve by var = (Sigma0^-1 + Sigma^-1)^-1
  // Sigma0 = A*cov_kernel*A'
  // Sigma = diag(sig)
  // mu = var * (Sigma0^-1 * mu0 + res)
  arma::mat Sigma0(t_size, t_size);
  // arma::mat Sigma_inv(t_size, t_size);
  arma::mat Sigma_inv = diagmat(conv_to<mat>::from(sig));
  for (size_t i = 0; i < t_size; i++){
    for (size_t j = 0; j < t_size; j++){
      Sigma0(i, j) = diag[i] * x_struct->cov_kernel[i][j] * diag[j];
    }
    // Sigma_inv(i, i) = 1 / sig[i];
  }

  arma::mat Sigma0_inv = pinv(Sigma0);
  arma::mat var_inv = Sigma0_inv + Sigma_inv;
  arma::mat var = pinv(var_inv);

  arma::mat U, V;
  arma::vec s;
  svd(U, s, V, var);

  arma::mat L = U * diagmat(s);
  // mean
  arma::mat res_vec(t_size, 1);
  for (size_t i = 0; i < t_size; i++){
    res_vec(i, 0) = resid[i];
  }
  arma::mat mu = var * Sigma_inv * res_vec;

  std::normal_distribution<double> normal_samp(0.0, 1.0);
  arma::mat draws(t_size, 1);
  for (size_t i = 0; i < t_size; i++){ draws(i, 0) = normal_samp(state->gen); }

  arma::mat beta_tilde = mu + L * draws;

  // beta = diag^-1 * beta_tilde
  // arma::mat beta(t_size, 1);
  for (size_t i = 0; i < t_size; i++){
    beta[i] = beta_tilde(i, 0) / diag[i];
  }

  // match beta to beta_t
  for (size_t i = 0; i < state->p_y; i++){
    state->beta_t[i] = beta[idx[i]];
  }
}


// subtracts old tree contribution from the fit
// called from mcmc_loop_xbcf in xbcf_mcmc_loop.cpp
void longBetModel::subtract_old_tree_fit(size_t tree_ind, matrix<double> &fit, std::unique_ptr<X_struct> &x_struct)
{
  for (size_t i = 0; i < fit.size(); i++)  // N
  {
    for (size_t j = 0; j < fit[i].size(); j++){  // p_y
      fit[i][j] -= (*(x_struct->data_pointers[tree_ind][i * fit[i].size() + j]))[0];
    }
  }
}

// sets unique term parameters in the state object depending on the term being updated
// called from mcmc_loop_xbcf in xbcf_mcmc_loop.cpp
void longBetModel::set_state_status(std::unique_ptr<State> &state, size_t value, const double *X, matrix<size_t> &Xorder, const double *t_std)
{
  state->fl = value; // value can only be 0 or 1 (to alternate between arms)
  state->iniSplitStorage(state->fl);
  state->adjustMtry(state->fl);
  state->X_std = X;
  state->Xorder_std = Xorder;
  state->t_std = t_std;
  if(value == 0)
  {
    state->p = state->p_pr;
    state->p_categorical = state->p_categorical_pr;
    state->p_continuous = state->p_continuous_pr;
  } else {
    state->p = state->p_trt;
    state->p_categorical = state->p_categorical_trt;
    state->p_continuous = state->p_continuous_trt;
  }

}

void longBetModel::predict_beta(std::vector<double> &beta,
  std::vector<double> &res_vec, std::vector<double> &a_vec, std::vector<double> &sig_vec, 
  matrix<double> &Sigma_tr_std, matrix<double> &Sigma_te_std, matrix<double> &Sigma_tt_std,
  std::mt19937 &gen)
{
  vec a_diag = conv_to<vec>::from(a_vec);
  vec sig_diag = conv_to<vec>::from(sig_vec);
  vec res = conv_to<vec>::from(res_vec);

  size_t tr_size = Sigma_tr_std.size();
  size_t te_size = Sigma_te_std.size();
  mat Sigma_tr(tr_size, tr_size);
  mat Sigma_te(te_size, te_size);
  mat Sigma_tt(tr_size, te_size);

  std_to_arma(Sigma_tr_std, Sigma_tr);
  std_to_arma(Sigma_te_std, Sigma_te);
  std_to_arma(Sigma_tt_std, Sigma_tt);

  mat A = diagmat(a_diag);
  mat Sig = diagmat(sig_diag);
  mat Sig_inv = pinv(Sig + A * Sigma_tr * A.t());
  mat common_mat = Sigma_tt.t() * A.t() * Sig_inv;
  
  mat mu = common_mat * res;
  mat var = Sigma_te - common_mat * A * Sigma_tt;

  arma::mat U, V;
  arma::vec s;
  svd(U, s, V, var);

  arma::mat L = U * diagmat(s);

  std::normal_distribution<double> normal_samp(0.0, 1.0);
  arma::mat draws(te_size, 1);
  for (size_t i = 0; i < te_size; i++){ draws(i, 0) = normal_samp(gen); }

  arma::mat beta_tilde = mu + L * draws;

  for (size_t i = 0; i < te_size; i++){
    beta[i] = beta_tilde(i, 0);
  }
  
}
