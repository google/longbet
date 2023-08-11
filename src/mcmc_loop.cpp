#include "mcmc_loop.h"
#include <memory>
#include <ostream>
#include "split_info.h"
// BCF main loop
// input includes information about two sets of trees 
// (one for prognostic term, the other for treatment term)

void mcmc_loop_longbet( 
  std::unique_ptr<split_info> &split_pr,  
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
  )
{

  if (state->parallel)
    thread_pool.start();

  verbose = true;

  for (size_t sweeps = 0; sweeps < state->num_sweeps; sweeps++)
  {
    if (verbose == true)
    {
      COUT << "--------------------------------" << endl;
      COUT << "number of sweeps " << sweeps << endl;
      COUT << "--------------------------------" << endl;
    }

    model_ps->set_state_status(state, 0, x_struct_ps->X_std, split_pr->Xorder_std, x_struct_ps->t_std);

    ////////////// Prognostic term loop
    for (size_t tree_ind = 0; tree_ind < state->num_trees_vec[0]; tree_ind++)
    {
      if (verbose == true)
      {
        COUT << "--------------------------------" << endl;
        COUT << "number of prognostic trees " << tree_ind << endl;
        COUT << "--------------------------------" << endl;
      }
      state->update_residuals();  // update residuals
      model_ps->draw_sigma(state, 0);  // draw sigmas
      // store sigma draws
      sigma0_draw_xinfo[sweeps][tree_ind] = state->sigma_vec[0];
      sigma1_draw_xinfo[sweeps][tree_ind] = state->sigma_vec[1];
      // cout << "sigma = " << state->sigma_vec << endl;

      if (state->use_all && (sweeps > state->burnin) &&
      (state->mtry_pr != state->p_pr))
      {
        state->use_all = false;
      }

      // clear counts of splits for one tre
      std::fill(state->split_count_current_tree.begin(),
      state->split_count_current_tree.end(), 0.0);

      if (state->sample_weights_flag)
      {
        // subtract old tree for sampling case
        state->mtry_weight_current_tree = state->mtry_weight_current_tree -
        state->split_count_all_tree_pr[tree_ind];
      }

      // get partial mu_fit -- thus take out the old fitted values
      model_ps->subtract_old_tree_fit(tree_ind, state->mu_fit, x_struct_ps);
      // initialize suff stat using partial fit
      model_ps->initialize_root_suffstat(state,
      trees_ps[sweeps][tree_ind].suff_stat);
      // cout << "root suffstat = " << trees_ps[sweeps][tree_ind].suff_stat << endl;
      // GFR
      trees_ps[sweeps][tree_ind].grow_from_root(state, split_pr, model_ps,
      x_struct_ps, sweeps, tree_ind, split_time_ps);
      model_ps->state_sweep(tree_ind, state->mu_fit, x_struct_ps);  // update total mu_fit by adding just fitted values

      state->update_split_counts(tree_ind, 0);  // update split counts for mu 
    }

    model_ps->set_state_status(state, 1, x_struct_trt->X_std, split_trt->Xorder_std, x_struct_trt->t_std);
    
    ////////////// Treatment term loop
    for (size_t tree_ind = 0; tree_ind < state->num_trees_vec[1]; tree_ind++)
    {
       if (verbose == true)
      {
        COUT << "--------------------------------" << endl;
        COUT << "number of treatment trees " << tree_ind << endl;
        COUT << "--------------------------------" << endl;
      }
      // cout << "beta_t = " << state->beta_t << endl;
      // cout << "b_vec = " << state->b_vec << endl;
      state->update_residuals(); // update residuals
      model_trt->draw_sigma(state, 1); // draw sigmas (and update them in the state obj)

      // store sigma draws
      sigma0_draw_xinfo[sweeps][state->num_trees_vec[0]+tree_ind] = state->sigma_vec[0]; // storing sigmas
      sigma1_draw_xinfo[sweeps][state->num_trees_vec[0]+tree_ind] = state->sigma_vec[1]; // storing sigmas

      if (state->use_all && (sweeps > state->burnin) && (state->mtry_trt != state->p_trt))
      {
        state->use_all = false;
      }

      std::fill(state->split_count_current_tree.begin(), state->split_count_current_tree.end(), 0.0); // clear counts of splits for one tree

      if (state->sample_weights_flag)
      {
        state->mtry_weight_current_tree = state->mtry_weight_current_tree - state->split_count_all_tree_trt[tree_ind]; // subtract old tree for sampling case
      }

      model_trt->subtract_old_tree_fit(tree_ind, state->tau_fit, x_struct_trt); // for GFR we will need partial tau_fit -- thus take out the old fitted values

      model_trt->initialize_root_suffstat(state, trees_trt[sweeps][tree_ind].suff_stat); // initialize suff stat using partial fit
      // cout << "root suffstat = " << trees_trt[sweeps][tree_ind].suff_stat << endl;
      // GFR
      trees_trt[sweeps][tree_ind].grow_from_root(state, split_trt, model_trt, x_struct_trt, sweeps, tree_ind, split_time_trt);
      // cout << "finish treatment " << tree_ind << endl;

      model_trt->state_sweep(tree_ind, state->tau_fit, x_struct_trt); // update total tau_fit by adding just fitted values

      state->update_split_counts(tree_ind, 1); // update split counts for tau
      
    }

    if (sweeps != 0)
    {
      if (a_scaling) // in case b_scaling on, we update b0 and b1
      {
        model_ps->update_a_value(state);
      }
      if (b_scaling) // in case b_scaling on, we update b0 and b1
      {
        model_trt->update_b_values(state);
      }
    }

    // TODO: replace torder with sorder
    model_ps->update_time_coef(state, x_struct_trt, split_trt->torder_std, 
      resid_info[sweeps], A_diag_info[sweeps], Sig_diag_info[sweeps], beta_info[sweeps]); 

    std::copy(state->beta_t.begin(), state->beta_t.end(),
    beta_xinfo[sweeps].begin());
    // store draws for b0, b1 and a, although they are updated per tree, we save results per forest (sweep)
    b_xinfo[0][sweeps] = state->b_vec[0];
    b_xinfo[1][sweeps] = state->b_vec[1];
    a_xinfo[0][sweeps] = state->a;
    // cout << "finish " << endl;
  }

  thread_pool.stop();
}