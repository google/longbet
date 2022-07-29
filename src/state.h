#ifndef GUARD_fit_info_h
#define GUARD_fit_info_h

#include <algorithm>
#include <cstddef>
#include <ctime>
#include "common.h"
#include "utility.h"
#include <chrono>
#include <vector>

class State
{
public:
    size_t dim_suffstat;


    // residual vectors
    // total residual vector
    matrix<double> residual;
    // residual for treated group, length n_trt
    matrix<double> full_residual_trt;
    matrix<double> full_residual_ctrl;

    std::vector<double> beta_t;
    std::vector<double> t;
    std::vector<double> unique_t;

    // Random
    std::vector<double> prob;
    std::random_device rd;
    std::mt19937 gen;
    std::discrete_distribution<> d;

    // Splits
    matrix<double> split_count_all_tree;
    matrix<double> split_count_all_tree_pr;  // TODO: move to xbcfState
    matrix<double> split_count_all_tree_trt; // TODO: move to xbcfState
    std::vector<double> split_count_current_tree;
    std::vector<double> mtry_weight_current_tree;

    // mtry
    bool use_all = true;

    // fitinfo
    size_t n_min;
    size_t n_cutpoints;
    bool parallel;
    size_t p_categorical;
    size_t p_continuous;
    size_t p; // total number of variables = p_categorical + p_continuous
    size_t mtry;
    size_t n_y;  // number of total data points in root node
    size_t p_y;  // dimension of response variables

    const double *X_std;  // pointer to original data
    const double *y_std;  // pointer to y data
    const double *z;            // the scaled treatment vector            TODO: move to xbcfState
    
    std::vector<size_t> n_trt;                     // the number of treated individuals      TODO: check if it's used anywhere after restructuring
    matrix<double> mu_fit;       // total mu_fit                           TODO: move to xbcfState
    matrix<double> tau_fit;      // total tau_fit                          TODO: move to xbcfState
    std::vector<double> b_vec;        // scaling parameters for tau (b0,b1)     TODO: move to xbcfState
    std::vector<double> sigma_vec;    // residual standard deviations           TODO: move to xbcfState
    double a;                         // scaling parameter for mu               TODO: move to xbcfState
    size_t p_categorical_pr;          // TODO: move to xbcfState
    size_t p_continuous_pr;           // TODO: move to xbcfState
    size_t p_categorical_trt;         // TODO: move to xbcfState
    size_t p_continuous_trt;          // TODO: move to xbcfState
    size_t p_pr;                      // total number of variables for mu          TODO: move to xbcfState
    size_t p_trt;                     // total number of variables for tau          TODO: move to xbcfState
    size_t mtry_pr;                   // TODO: move to xbcfState
    size_t mtry_trt;                  // TODO: move to xbcfState

    size_t max_depth;
    size_t num_trees;
    std::vector<size_t> num_trees_vec;
    size_t num_sweeps;
    size_t burnin;
    bool sample_weights_flag;
    double ini_var_yhat;
    size_t fl; // flag for likelihood function to alternate between mu loop and tau loop calculations  TODO: move to xbcfState

    matrix<size_t> Xorder_std;

    // residual standard deviation
    double sigma;
    double sigma2; // sigma squared

    // time information
    const double *t_std;
    size_t n_t;
    size_t p_t;
    //std::vector<double> precision_squared;

    void update_residuals()
    {
        size_t index_trt;
        size_t index_ctrl;
        const double *temp_pointer;

        for (size_t j = 0; j < this->p_y; j++){
            index_trt = 0;
            index_ctrl = 0;
            temp_pointer = this->z + j * this->n_y;

            for (size_t i = 0; i < this->n_y; i++)
            {
                if (*(temp_pointer + i) == 1)
                {
                    this->full_residual_trt[j][index_trt] = *(this->y_std + this->n_y * j + i) - this->a * this->mu_fit[i][j] - this->b_vec[1] * this->beta_t[j] * this->tau_fit[i][j];
                    index_trt++;
                }
                else
                {
                    this->full_residual_ctrl[j][index_ctrl] = *(this->y_std + this->n_y * j + i) - this->a * this->mu_fit[i][j] - this->b_vec[0] * this->beta_t[j] * this->tau_fit[i][j];
                    index_ctrl++;
                }
            }
        }
    }

    void update_sigma(double sigma)
    {
        this->sigma = sigma;
        this->sigma2 = pow(sigma, 2);
   }

    // sigma update for longBetModel       TODO: move to xbcfClass
    void update_sigma(double sigma, size_t ind)
    {
        this->sigma_vec[ind] = sigma;    }

    // sigma update for longBetModel       TODO: move to xbcfClass
    void update_bscales(double b0, double b1)
    {
        this->b_vec[0] = b0;  // sigma for the control group
        this->b_vec[1] = b1;
    }

    //  TODO: update the constructor / get rid of it
    State(const double *Xpointer, matrix<size_t> &Xorder_std, size_t N,
    size_t p_pr, size_t p_trt, size_t p_y, std::vector<size_t> num_trees_vec,
    size_t p_categorical_pr, size_t p_categorical_trt, size_t p_continuous_pr,
    size_t p_continuous_trt, bool set_random_seed, size_t random_seed,
    size_t n_min, size_t n_cutpoints, bool parallel, size_t mtry_pr,
    size_t mtry_trt, const double *X_std, size_t num_sweeps, bool
    sample_weights_flag, const double *y_std,
    const double *z, std::vector<double> sigma_vec,
    std::vector<double> b_vec, size_t max_depth, double ini_var_yhat,
    size_t burnin)
    {
        // Random
        this->prob = std::vector<double>(2, 0.5);
        this->gen = std::mt19937(rd());
        if (set_random_seed)
        {
            gen.seed(random_seed);
        }
        this->d = std::discrete_distribution<>(prob.begin(), prob.end());

        // Splits
        ini_xinfo(this->split_count_all_tree_pr, p_pr, num_trees_vec[0]);
        ini_xinfo(this->split_count_all_tree_trt, p_trt, num_trees_vec[1]);

        this->n_min = n_min;
        this->n_cutpoints = n_cutpoints;
        this->parallel = parallel;
        this->p_categorical_pr = p_categorical_pr;
        this->p_continuous_pr = p_continuous_pr;
        this->p_categorical_trt = p_categorical_trt;
        this->p_continuous_trt = p_continuous_trt;
        this->mtry_pr = mtry_pr;
        this->mtry_trt = mtry_trt;
        this->X_std = X_std;
        this->p_pr = p_categorical_pr + p_continuous_pr;
        this->p_trt = p_categorical_trt + p_continuous_trt;
        this->n_y = N;
        this->p_y = p_y;
        this->num_trees_vec = num_trees_vec;  // stays the same even for vector
        this->num_sweeps = num_sweeps;
        this->sample_weights_flag = sample_weights_flag;
        this->y_std = y_std;
        this->max_depth = max_depth;
        this->burnin = burnin;
        this->ini_var_yhat = ini_var_yhat;
        this->Xorder_std = Xorder_std;

   }

    void update_split_counts(size_t tree_ind)
    {
        mtry_weight_current_tree = mtry_weight_current_tree + split_count_current_tree;
        split_count_all_tree[tree_ind] = split_count_current_tree;
    }

    void update_split_counts(size_t tree_ind, size_t flag)
    {
        mtry_weight_current_tree = mtry_weight_current_tree + split_count_current_tree;
        if (flag == 0)
        {
            split_count_all_tree_pr[tree_ind] = split_count_current_tree;
        }
        else
        {
            split_count_all_tree_trt[tree_ind] = split_count_current_tree;
        }
    }

    void iniSplitStorage(size_t flag)
    {
        if (flag == 0)
        {
            this->split_count_current_tree = std::vector<double>(this->p_pr, 0);
            this->mtry_weight_current_tree = std::vector<double>(this->p_pr, 0);
        }
        else if (flag == 1)
        {
            this->split_count_current_tree = std::vector<double>(this->p_trt, 0);
            this->mtry_weight_current_tree = std::vector<double>(this->p_trt, 0);
        }
    }

    void adjustMtry(size_t flag)
    {
        if (flag == 0)
        {
            this->mtry = this->mtry_pr;
        }
        else if (flag == 1)
        {
            this->mtry = this->mtry_trt;
        }
    }

    void set_t_info(const double *t_std, size_t n_t)
    {
        this->t_std = t_std;
        this->n_t = n_t;
    }
};

class longBetState : public State
{
 public:
    longBetState(const double *Xpointer, matrix<size_t> &Xorder_std, size_t N,
    std::vector<size_t> n_trt, size_t p, size_t p_tau, size_t p_y,
    std::vector<size_t> num_trees_vec,
    size_t p_categorical_pr, size_t p_categorical_trt, size_t p_continuous_pr,
    size_t p_continuous_trt, bool set_random_seed, size_t random_seed,
    size_t n_min, size_t n_cutpoints, bool parallel, size_t mtry_pr,
    size_t mtry_trt, const double *X_std, size_t num_sweeps,
    bool sample_weights_flag, const double *y_std,
    const double *z,
    std::vector<double> sigma_vec, std::vector<double> b_vec, size_t max_depth,
    double ini_var_yhat, size_t burnin, size_t dim_suffstat) :
    State(Xpointer, Xorder_std, N, p, p_tau, p_y, num_trees_vec,
    p_categorical_pr, p_categorical_trt, p_continuous_pr, p_continuous_trt,
    set_random_seed, random_seed, n_min, n_cutpoints, parallel, mtry_pr,
    mtry_trt, X_std, num_sweeps, sample_weights_flag, y_std, z,
    sigma_vec, b_vec, max_depth, ini_var_yhat, burnin)
    {
        this->sigma_vec = sigma_vec;
        this->b_vec = b_vec;
        this->n_trt = n_trt;
        this->num_trees_vec = num_trees_vec;
        this->z = z;
        this->a = 1;  // initialize a at 1 for now

        this->dim_suffstat = dim_suffstat;

        ini_matrix(this->mu_fit, p_y, N);
        ini_matrix(this->tau_fit, p_y, N);
        this->beta_t = std::vector<double>(p_y, 1);

        // those are for XBCF, initialize at a length 1 vector
        // this->residual = std::vector<double>(N, 0);
        // this->full_residual_ctrl = std::vector<double>(N - n_trt, 0);
        // this->full_residual_trt = std::vector<double>(n_trt, 0);
        ini_matrix(this->residual, p_y, N);
        ini_residuals(N, p_y, n_trt);
    }

    void ini_residuals(size_t N, size_t p_y, std::vector<size_t> n_trt){
        this->full_residual_trt.resize(p_y);
        this->full_residual_ctrl.resize(p_y);
        for (size_t i = 0; i < p_y; i++){
            this->full_residual_trt[i].resize(n_trt[i]);
            this->full_residual_ctrl[i].resize(N - n_trt[i]);
        }
    }

};


#endif