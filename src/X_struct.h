#ifndef GUARD_X_struct_h
#define GUARD_X_struct_h
#include <iostream>
#include <vector>
#include "common.h"
#include "utility.h"

struct X_struct
{
public:
    // Vector pointers
    // std::vector<matrix<std::vector<double> *>> data_pointers;
    matrix<std::vector<double> *> data_pointers;

    std::vector<double> X_values;
    std::vector<size_t> X_counts;
    std::vector<size_t> variable_ind;
    std::vector<size_t> X_num_unique;

    std::vector<double> t_values;
    std::vector<size_t> t_counts;
    std::vector<size_t> t_variable_ind;
    std::vector<size_t> t_num_unique;
    matrix<double> cov_kernel;

    std::vector<double> s_values;
    std::vector<size_t> s_counts;
    std::vector<size_t> s_variable_ind;
    std::vector<size_t> s_num_unique;

    const double *X_std;  // pointer to original data
    const double *y_std;  // pointer to y data
    const double *t_std;  // pointer to t data
    const double *s_std;  // pointer to s (cumulative treatment time)
    size_t n_y;  // number of total data points in root node
    size_t p_y;
    size_t p_continuous;
    size_t p_x;

    size_t n_t;
    size_t p_t;

    X_struct(const double *X_std, const double *y_std, const double *t_std, const double *s_std,
    size_t n_y, size_t p_y, std::vector<std::vector<size_t>> &Xorder_std,
    std::vector<std::vector<size_t>> &torder_std, std::vector<std::vector<size_t>> &sorder_std,
    size_t p_categorical, size_t p_continuous,
    std::vector<double> *initial_theta, size_t num_trees,
    double &sig_knl, double &lambda_knl)
    {

        this->variable_ind = std::vector<size_t>(p_categorical + 1);
        this->X_num_unique = std::vector<size_t>(p_categorical);

        init_tree_pointers(initial_theta, num_trees, n_y, p_y);

        // std::cout << "ini dp size = " << this->data_pointers[0].size() << endl;

        unique_value_count2(X_std, Xorder_std, X_values, X_counts, variable_ind, n_y, X_num_unique, p_categorical, p_continuous);

        size_t t_categorical = 1;
        size_t t_continuous = 0;
        this->t_variable_ind = std::vector<size_t>(t_categorical + 1);
        this->t_num_unique = std::vector<size_t>(t_categorical);

        unique_value_count2(t_std, torder_std, t_values, t_counts, t_variable_ind, p_y, t_num_unique, t_categorical, t_continuous);

        size_t s_categorical = n_y;
        size_t s_continuous = 0;
        this->s_variable_ind = std::vector<size_t>(s_categorical + 1);
        this->s_num_unique = std::vector<size_t>(s_categorical);

        unique_value_count2(s_std, sorder_std, s_values, s_counts, s_variable_ind, p_y, s_num_unique, s_categorical, s_continuous);

        this->X_std = X_std;
        this->y_std = y_std;
        this->t_std = t_std;
        this->s_std = s_std;
        this->n_y = n_y;
        this->p_y = p_y;
        this->p_continuous = p_continuous;
        this->p_x = Xorder_std.size();
        this->n_t = torder_std[0].size();
        this->p_t = torder_std.size();

        ini_cov_kernel(sig_knl, lambda_knl);
    }

    void init_tree_pointers(std::vector<double> *initial_theta, size_t num_trees, size_t N, size_t p_y)
    {
        ini_matrix(this->data_pointers, p_y * N, num_trees);

        for (size_t i = 0; i < num_trees; i++)
        {
            for (size_t j = 0; j < N; j++)
            {
                std::vector<std::vector<double> *> &pointer_vec = this->data_pointers[i];
                for (size_t k = 0; k < p_y; k++){
                    pointer_vec[j * p_y + k] = initial_theta;
                    // data_pointers[i][j * p_y + k] = initial_theta;
                }
            }
        }
    }

    void ini_cov_kernel(double &sig_knl, double &lambda_knl){
        double sigma2 = pow(sig_knl, 2);
        double lambda2 = pow(lambda_knl, 2);

        size_t t_size = t_values.size();
        ini_matrix(this->cov_kernel, t_size, t_size);
        double diag = squared_exponential(t_values[0], t_values[0], sigma2, lambda2);

        for (size_t i = 0; i < t_size; i++){
            // calculate diagonal element
            cov_kernel[i][i] = diag;
            for (size_t j = 0; j < i; j++){
                cov_kernel[i][j] = squared_exponential(t_values[i], t_values[j],
                sigma2, lambda2);
                cov_kernel[j][i] = cov_kernel[i][j];
            }
        }
        // std::cout << "cov_kernel = " << cov_kernel << endl;
    }
};

#endif