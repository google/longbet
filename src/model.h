
#ifndef model_h
#define model_h

#include "common.h"
#include "utility.h"
#include "matrix_utility.h"
#include <memory>
#include "state.h"
#include "X_struct.h"
#include "cdf.h"

using namespace std;

class tree;

class Model
{

public:
    size_t dim_theta;

    size_t dim_suffstat;

    size_t dim_residual;

    /////////////////////////////////////
    //
    //  suff_stat_model and suff_stat_total
    //  are useless for NormalModel now
    //  They are still here because CLT class depends on them
    //  Delelte them later
    //
    /////////////////////////////////////
    std::vector<double> suff_stat_model;

    std::vector<double> suff_stat_total;

    double no_split_penality;

    // tree prior
    double alpha;

    double beta;

    Model(size_t dim_theta, size_t dim_suff)
    {
        this->dim_theta = dim_theta;
        this->dim_suffstat = dim_suff;
    };

    // Abstract functions
    virtual void incSuffStat(std::unique_ptr<State> &state, size_t index_next_obs, size_t index_next_t, std::vector<double> &suffstats) { };

    virtual void samplePars(std::unique_ptr<State> &state, std::vector<double> &suff_stat, std::vector<double> &theta_vector, double &prob_leaf) { };

    virtual void update_state(std::unique_ptr<State> &state, size_t tree_ind, std::unique_ptr<X_struct> &x_struct) { };

    virtual void initialize_root_suffstat(std::unique_ptr<State> &state, std::vector<double> &suff_stat) { };

    virtual void updateNodeSuffStat(std::vector<double> &suff_stat, std::unique_ptr<State> &state, matrix<size_t> &Xorder_std, matrix<size_t> &torder_std, size_t &split_var, size_t row_ind) { };

    virtual void calculateOtherSideSuffStat(std::vector<double> &parent_suff_stat, std::vector<double> &lchild_suff_stat, std::vector<double> &rchild_suff_stat, bool &compute_left_side) { };
    
    virtual void state_sweep(size_t tree_ind, size_t M, matrix<double> &residual_std, std::unique_ptr<X_struct> &x_struct) const { };
    
    virtual double likelihood(std::vector<double> &temp_suff_stat, std::vector<double> &suff_stat_all, bool left_side, bool no_split, std::unique_ptr<State> &state) const { return 0.0; };

    // virtual double likelihood_no_split(std::vector<double> &suff_stat, std::unique_ptr<State> &state) const { return 0.0; };

    virtual void ini_residual_std(std::unique_ptr<State> &state) { };

    // virtual double predictFromTheta(const std::vector<double> &theta_vector) const { return 0.0; };

    virtual void predict_std(const double *Xtestpointer, const double *tpointer, size_t N_test, size_t p, size_t num_trees, size_t num_sweeps, std::vector<matrix<double>> &yhats_test_xinfo, vector<vector<tree>> &trees) { };

    virtual void subtract_old_tree_fit(size_t tree_ind, matrix<double> &fit,
    std::unique_ptr<X_struct> &x_struct) { };

    virtual Model *clone() { return nullptr; };

    // Getters and Setters
    // num classes
    size_t getNumClasses() const { return dim_theta; };

    void setNumClasses(size_t n_class) { dim_theta = n_class; };

    // dim suff stat
    size_t getDimSuffstat() const { return dim_suffstat; };

    void setDimSuffStat(size_t dim_suff) { dim_suffstat = dim_suff; };

    //penality
    double getNoSplitPenality()
    {
        return no_split_penality;
        ;
    };
    void setNoSplitPenality(double pen) { this->no_split_penality = pen; };
};

//////////////////////////////////////////////////////////////////////////////////////
//
//
//  longBet Model
//
//
//////////////////////////////////////////////////////////////////////////////////////

class longBetModel : public Model
{
public:
  size_t dim_suffstat = 4;

  // model prior
  // prior on sigma
  double kap;
  double s;
  // prior on leaf parameter
  double tau;

  longBetModel(double kap, double s, double tau, double alpha, double beta) : Model(1, 4)
  {
    this->kap = kap;
    this->s = s;
    this->tau = tau;
    this->alpha = alpha;
    this->beta = beta;
    this->dim_residual = 1;
  }

  longBetModel() : Model(1, 4) {}

  Model *clone() { return new longBetModel(*this); }

  void incSuffStat(std::unique_ptr<State> &state, size_t index_next_obs, size_t index_next_t, std::vector<double> &suffstats);

  void samplePars(std::unique_ptr<State> &state, std::vector<double> &suff_stat_vec, std::vector<double> &theta_vector, double &prob_leaf);

  void draw_sigma(std::unique_ptr<State> &state, size_t ind);

  void initialize_root_suffstat(std::unique_ptr<State> &state, std::vector<double> &suff_stat);

  void updateNodeSuffStat(std::vector<double> &suff_stat, std::unique_ptr<State> &state, matrix<size_t> &Xorder_std, matrix<size_t> &torder_std, size_t &split_var, size_t row_ind);

  void calculateOtherSideSuffStat(std::vector<double> &parent_suff_stat, std::vector<double> &lchild_suff_stat, std::vector<double> &rchild_suff_stat, bool &compute_left_side);

  void state_sweep(size_t tree_ind, matrix<double> &fit, std::unique_ptr<X_struct> &x_struct) const;

  void update_xinfo(matrix<double> &yhats_xinfo, size_t sweep_num, size_t num_trees, size_t N, std::unique_ptr<X_struct> &x_struct);

  double likelihood(std::vector<double> &temp_suff_stat, std::vector<double> &suff_stat_all, bool left_side, bool no_split, std::unique_ptr<State> &state) const;

  void predict_std(const double *Xtestpointer, const double *tpointer, size_t N_test, size_t p, size_t num_sweeps, std::vector<matrix<double>> &yhats_test_xinfo, vector<vector<tree>> &trees);

  void update_a_value(std::unique_ptr<State> &state);

  void update_b_values(std::unique_ptr<State> &state);

  void update_time_coef(std::unique_ptr<State> &state, std::unique_ptr<X_struct> &x_struct,
    matrix<size_t> &torder_std, std::vector<double> &resid);
  
  void subtract_old_tree_fit(size_t tree_ind, matrix<double> &fit,
  std::unique_ptr<X_struct> &x_struct);

  void set_state_status(std::unique_ptr<State> &state, size_t value,
  const double *X, matrix<size_t> &Xorder, const double *t_std);
};


#endif
