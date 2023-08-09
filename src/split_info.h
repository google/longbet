#ifndef THIRD_PARTY_R_PACKAGES_LONGBET_LONGBET_SRC_SPLIT_INFO_H_
#define THIRD_PARTY_R_PACKAGES_LONGBET_LONGBET_SRC_SPLIT_INFO_H_
#include <iostream>
#include <memory>
#include <vector>
#include "common.h"
#include "state.h"
#include "model.h"
#include "X_struct.h"

struct split_info
{
public:
    matrix<size_t> Xorder_std;
    matrix<size_t> torder_std;

    matrix<size_t> sorder_std;
    std::vector<double> s_values;

    std::vector<size_t> X_counts;
    std::vector<size_t> X_num_unique;
    std::vector<size_t> t_counts;
    std::vector<size_t> t_num_unique;

    size_t N_Xorder;
    size_t p_Xorder;
    size_t N_torder;
    size_t p_torder;

    split_info(std::unique_ptr<X_struct> &x_struct, matrix<size_t> &Xorder_std, matrix<size_t> &torder_std, matrix<size_t> &sorder_std, std::vector<double> s_values){
      this->Xorder_std = Xorder_std;
      this->torder_std = torder_std;
      this->X_counts  = x_struct->X_counts;
      this->X_num_unique = x_struct->X_num_unique;
      this->t_counts  = x_struct->t_counts;
      this->t_num_unique = x_struct->t_num_unique;

      this->sorder_std = sorder_std;
      this->s_values = s_values;

      this->N_Xorder = Xorder_std[0].size();
      this->p_Xorder = Xorder_std.size();
      this->N_torder = torder_std[0].size();
      this->p_torder = torder_std.size();
    }

    split_info(std::unique_ptr<split_info> &parent, size_t &split_var,
    size_t &split_point, bool split_t, bool left)
    {
      this->X_num_unique.resize(parent->X_num_unique.size());
      this->X_counts.resize(parent->X_counts.size());

      this->t_num_unique.resize(parent->t_num_unique.size());
      this->t_counts.resize(parent->t_counts.size());

      if (split_t)
      {
        ini_xinfo_sizet(this->Xorder_std, parent->N_Xorder, parent->p_Xorder);
        // copy Xorder
        for (size_t i = 0; i < Xorder_std.size(); i++){
          std::copy(parent->Xorder_std[i].begin(), parent->Xorder_std[i].end(), this->Xorder_std[i].begin());
          std::copy(parent->Xorder_std[i].begin(), parent->Xorder_std[i].end(), this->Xorder_std[i].begin());
        }

        if (left)
        {
          ini_xinfo_sizet(this->torder_std, split_point + 1, parent->p_torder);
        } else {
          ini_xinfo_sizet(this->torder_std, parent->N_torder - split_point - 1, parent->p_torder);
        }
      } else {
        ini_xinfo_sizet(this->torder_std, parent->N_torder, parent->p_torder);
        // copy torder
        for (size_t i = 0; i < torder_std.size(); i++){
          std::copy(parent->torder_std[i].begin(), parent->torder_std[i].end(), this->torder_std[i].begin());
          std::copy(parent->torder_std[i].begin(), parent->torder_std[i].end(), this->torder_std[i].begin());
        }
        if (left) 
        {
          ini_xinfo_sizet(this->Xorder_std, split_point + 1, parent->p_Xorder);
        } else {
          ini_xinfo_sizet(this->Xorder_std, parent->N_Xorder - split_point - 1, parent->p_Xorder);
        }
      }

      this->N_Xorder = Xorder_std[0].size();
      this->p_Xorder = Xorder_std.size();
      this->N_torder = torder_std[0].size();
      this->p_torder = torder_std.size();
    }

    // split(std::unique_ptr<split_info> &parent,)

  void split_this(std::unique_ptr<split_info> &split_left,
  std::unique_ptr<split_info> &split_right,
  size_t split_var, size_t split_point, bool split_t,
  Model *model, std::unique_ptr<X_struct> &x_struct,
  std::unique_ptr<State> &state, std::vector<double> &suff_stat,
  std::vector<double> &left_suff_stat, std::vector<double> &right_suff_stat)
  {
    if (split_t){

    // cout << "split_t, split_point = " << split_point << endl;
    split_torder_std(split_left, split_right, split_var, split_point,
    model, x_struct, state, suff_stat, left_suff_stat, right_suff_stat);
    
    // cout << "left_suff_stat = " << left_suff_stat << endl;
    // cout << "right_suff_stat = " << right_suff_stat << endl;
    } else {
      
      if (state->p_categorical > 0)
      {
        split_xorder_std_categorical(split_left, split_right, split_var, split_point, model, x_struct, state,
        suff_stat, left_suff_stat, right_suff_stat);
      }

      if (state->p_continuous > 0)
      {
        split_xorder_std_continuous(split_left, split_right, split_var, split_point, model, x_struct, state,
        suff_stat, left_suff_stat, right_suff_stat);
      }
    }

    double tol = -0.1;
    if ((left_suff_stat[2] < tol) | (right_suff_stat[2] < tol) | (left_suff_stat[3] < tol) | (right_suff_stat[3] < tol) ){
      cout << "suff stat < 0" << endl;
      cout << "parent suff stat = " << suff_stat << endl;
      cout << "left suff stat = " << left_suff_stat << endl;
      cout << "right suff stat = " << right_suff_stat << endl;
      abort();
    }
  }

  void split_xorder_std_continuous(std::unique_ptr<split_info> &split_left,
  std::unique_ptr<split_info> &split_right, size_t split_var, size_t split_point,
  Model *model, std::unique_ptr<X_struct> &x_struct,
  std::unique_ptr<State> &state, std::vector<double> &suff_stat,
  std::vector<double> &left_suff_stat, std::vector<double> &right_suff_stat)
  {

    // when find the split point, split Xorder matrix to two sub matrices for both subnodes

    // preserve order of other variables
    size_t N_Xorder = Xorder_std[0].size();
    size_t left_ix = 0;
    size_t right_ix = 0;
    size_t N_Xorder_left = split_left->Xorder_std[0].size();
    size_t N_Xorder_right = split_right->Xorder_std[0].size();

    // if the left side is smaller, we only compute sum of it
    bool compute_left_side = N_Xorder_left < N_Xorder_right;

    double cutvalue = *(state->X_std + state->n_y * split_var + Xorder_std[split_var][split_point]);

    const double *temp_pointer = state->X_std + state->n_y * split_var;

    // ini suff stat
    std::fill(left_suff_stat.begin(), left_suff_stat.end(), 0.0);
    std::fill(right_suff_stat.begin(), right_suff_stat.end(), 0.0);

    for (size_t j = 0; j < N_Xorder; j++)
    {
        if (compute_left_side)
        {
            if (*(temp_pointer + Xorder_std[split_var][j]) <= cutvalue)
            {
                model->updateNodeSuffStat(left_suff_stat, state, Xorder_std, torder_std, split_var, j);
            }
        }
        else
        {
            if (*(temp_pointer + Xorder_std[split_var][j]) > cutvalue)
            {
                model->updateNodeSuffStat(right_suff_stat, state, Xorder_std, torder_std, split_var, j);
            }
        }
    }

    const double *split_var_x_pointer = state->X_std + state->n_y * split_var;

    for (size_t i = 0; i < x_struct->p_continuous; i++) // loop over variables
    {
        // lambda callback for multithreading
        auto split_i = [&, i]() {
            size_t left_ix = 0;
            size_t right_ix = 0;

            std::vector<size_t> &xo = Xorder_std[i];
            std::vector<size_t> &xo_left = split_left->Xorder_std[i];
            std::vector<size_t> &xo_right = split_right->Xorder_std[i];

            for (size_t j = 0; j < N_Xorder; j++)
            {
                if (*(split_var_x_pointer + xo[j]) <= cutvalue)
                {
                    xo_left[left_ix] = xo[j];
                    left_ix = left_ix + 1;
                }
                else
                {
                    xo_right[right_ix] = xo[j];
                    right_ix = right_ix + 1;
                }
            }
        };
        if (thread_pool.is_active())
            thread_pool.add_task(split_i);
        else
            split_i();
    }
    if (thread_pool.is_active())
        thread_pool.wait();

    model->calculateOtherSideSuffStat(suff_stat, left_suff_stat, right_suff_stat, compute_left_side);

  }

  void split_xorder_std_categorical(std::unique_ptr<split_info> &split_left,
  std::unique_ptr<split_info> &split_right,
  size_t split_var, size_t split_point,
  Model *model, std::unique_ptr<X_struct> &x_struct,
  std::unique_ptr<State> &state, std::vector<double> &suff_stat,
  std::vector<double> &left_suff_stat, std::vector<double> &right_suff_stat)
  {
    // preserve order of other variables
    size_t N_Xorder = Xorder_std[0].size();
    size_t left_ix = 0;
    size_t right_ix = 0;
    size_t N_Xorder_left = split_left->Xorder_std[0].size();
    size_t N_Xorder_right = split_right->Xorder_std[0].size();

    size_t X_counts_index = 0;

    // if the left side is smaller, we only compute sum of it
    bool compute_left_side = N_Xorder_left < N_Xorder_right;

    size_t start;
    size_t end;

    const double *temp_pointer = state->X_std + state->n_y * split_var;
    double cutvalue = *(state->X_std + state->n_y * split_var +
    Xorder_std[split_var][split_point]);

    for (size_t i = x_struct->p_continuous; i < x_struct->p_x; i++)
    {
      // loop over variables
      left_ix = 0;
      right_ix = 0;

      // start <= i <= end;
      start = x_struct->variable_ind[i - x_struct->p_continuous];
      end = x_struct->variable_ind[i + 1 - x_struct->p_continuous];

      if (i == split_var)
      {
        // split the split_variable, only need to find row of cutvalue
        // COUT << "compute left side " << compute_left_side << endl;
        ///////////////////////////////////////////////////////////
        //
        // We should be able to run this part in parallel
        //
        //  just like split_xorder_std_continuous
        //
        ///////////////////////////////////////////////////////////

        if (compute_left_side)
        {
          for (size_t j = 0; j < N_Xorder; j++)
          {
            if (*(temp_pointer + Xorder_std[i][j]) <= cutvalue)
            {
              model->updateNodeSuffStat(left_suff_stat, state,
              Xorder_std, torder_std, split_var, j);
              split_left->Xorder_std[i][left_ix] = Xorder_std[i][j];
              left_ix = left_ix + 1;
            } else {
              // go to right side
              split_right->Xorder_std[i][right_ix] = Xorder_std[i][j];
              right_ix = right_ix + 1;
            }
          }
        } else {
          for (size_t j = 0; j < N_Xorder; j++)
          {
            if (*(temp_pointer + Xorder_std[i][j]) <= cutvalue)
            {
              split_left->Xorder_std[i][left_ix] = Xorder_std[i][j];
              left_ix = left_ix + 1;
            } else {
              // cout << "right, update node suff stat" << endl;
              model->updateNodeSuffStat(right_suff_stat, state,
              Xorder_std, torder_std, split_var, j);
              split_right->Xorder_std[i][right_ix] = Xorder_std[i][j];
              right_ix = right_ix + 1;
            }
          }
        }

        // for the cut variable, it's easy to counts X_counts_left and 
        // X_counts_right, simply cut X_counts to two pieces.
        for (size_t k = start; k < end; k++)
        {
          // loop from start to end!

          if (x_struct->X_values[k] <= cutvalue)
          {
            // smaller than cutvalue, go left
            split_left->X_counts[k] = X_counts[k];
          } else {
            // otherwise go right
            split_right->X_counts[k] = X_counts[k];
          }
        }

      } else {
        X_counts_index = start;
        for (size_t j = 0; j < N_Xorder; j++)
        {
          while (*(state->X_std + state->n_y * i + Xorder_std[i][j]) != x_struct->X_values[X_counts_index])
          {
            // find location of corresponding unique values
            X_counts_index++;
          }
          // cout << "j = "<< j << ", value =  " << *(temp_pointer + Xorder_std[i][j]) <<  endl;
          if (*(temp_pointer + Xorder_std[i][j]) <= cutvalue)
          {
            // go to left side
            split_left->Xorder_std[i][left_ix] = Xorder_std[i][j];
            left_ix = left_ix + 1;
            split_left->X_counts[X_counts_index]++;
          } else {
            // go to right side
            split_right->Xorder_std[i][right_ix] = Xorder_std[i][j];
            right_ix = right_ix + 1;

            split_right->X_counts[X_counts_index]++;
          }
        }
      }
    }

    model->calculateOtherSideSuffStat(suff_stat, left_suff_stat, right_suff_stat, compute_left_side);

    // update X_num_unique
    std::fill(split_left->X_num_unique.begin(), split_left->X_num_unique.end(), 0.0);
    std::fill(split_right->X_num_unique.begin(), split_right->X_num_unique.end(), 0.0);

    for (size_t i = x_struct->p_continuous; i < x_struct->p_x; i++)
    {
      start = x_struct->variable_ind[i - x_struct->p_continuous];
      end = x_struct->variable_ind[i + 1 - x_struct->p_continuous];

      // COUT << "start " << start << " end " << end << " size " << X_counts_left.size() << endl;
      for (size_t j = start; j < end; j++)
      {
        if (split_left->X_counts[j] > 0)
        {
          split_left->X_num_unique[i - x_struct->p_continuous] = split_left->X_num_unique[i - x_struct->p_continuous] + 1;
        }
        if (split_right->X_counts[j] > 0)
        {
          split_right->X_num_unique[i - x_struct->p_continuous] = split_right->X_num_unique[i - x_struct->p_continuous] + 1;
        }
      }
    }

  }

  void split_torder_std(std::unique_ptr<split_info> &split_left,
  std::unique_ptr<split_info> &split_right,
  size_t split_var, size_t split_point,
  Model *model, std::unique_ptr<X_struct> &x_struct,
  std::unique_ptr<State> &state, std::vector<double> &suff_stat,
  std::vector<double> &left_suff_stat, std::vector<double> &right_suff_stat)
  {
    // TODO: 
    //      Restructure sindex as vector of vector with varying size.
    //      All the splitting code need to be
    // initialize X_struct for the treatment term

    // cout << "split point " << split_point << " value " << *(state->t_std + state->n_t * split_var + torder_std[split_var][split_point]) << endl;
    // for (size_t i = 0; i < )
    // split t as categorical variable
    // preserve order of other variables
    size_t N_torder = torder_std[0].size();
    size_t left_ix = 0;
    size_t right_ix = 0;
    size_t N_torder_left = split_left->torder_std[0].size();
    size_t N_torder_right = split_right->torder_std[0].size();

    size_t X_counts_index = 0;

    // if the left side is smaller, we only compute sum of it
    bool compute_left_side = N_torder_left < N_torder_right;
    // cout << "Nleft " << N_torder_left << " Nright " << N_torder_right << endl; 

    size_t start;
    size_t end;

    const double *temp_pointer = state->t_std + state->n_t * split_var;
    double cutvalue = *(state->t_std + state->n_t * split_var +
    torder_std[split_var][split_point]);
    // cout << "cutvalue = " << cutvalue << endl;

    for (size_t i = 0; i < x_struct->p_t; i++)
    {
      // loop over variables
      left_ix = 0;
      right_ix = 0;

      start = x_struct->t_variable_ind[i];
      end = x_struct->t_variable_ind[i + 1];
      // cout << "i = " << i << " start = " << start << " end = " << end << endl;
      if (i == split_var)
      {
        ///////////////////////////////////////////////////////////
        //
        // We should be able to run this part in parallel
        //
        //  just like split_xorder_std_continuous
        //
        ///////////////////////////////////////////////////////////
        
        if (compute_left_side)
        {
          for (size_t j = 0; j < N_torder; j++)
          {

            if (*(temp_pointer + torder_std[i][j]) <= cutvalue)
            {
              for (size_t k = 0; k < Xorder_std[0].size(); k++){
                model->incSuffStat(state, Xorder_std[0][k],
                torder_std[split_var][j], left_suff_stat);
              }
              split_left->torder_std[i][left_ix] = torder_std[i][j];
              left_ix = left_ix + 1;
            } else {
              // go to right side
              split_right->torder_std[i][right_ix] = torder_std[i][j];
              right_ix = right_ix + 1;
            }
          }
        } else {
          for (size_t j = 0; j < N_torder; j++)
          {
            if (*(temp_pointer + torder_std[i][j]) <= cutvalue)
            {
              split_left->torder_std[i][left_ix] = torder_std[i][j];
              left_ix = left_ix + 1;
            } else {
              for (size_t k = 0; k < Xorder_std[0].size(); k++){
                model->incSuffStat(state, Xorder_std[0][k],
                torder_std[split_var][j], right_suff_stat);
              }
              split_right->torder_std[i][right_ix] = torder_std[i][j];
              right_ix = right_ix + 1;
            }
          }
        }

        for (size_t k = start; k < end; k++)
        {
          // loop from start to end!
          if (x_struct->t_values[k] <= cutvalue)
          {
            // smaller than cutvalue, go left
            split_left->t_counts[k] = t_counts[k];
          } else {
            // otherwise go right
            split_right->t_counts[k] = t_counts[k];
          }
        }
      } else {
        X_counts_index = start;

        // split other variables, need to compare each row
        for (size_t j = 0; j < N_torder; j++)
        {
          while (*(state->t_std + state->n_t * i + torder_std[i][j]) != x_struct->t_values[X_counts_index])
          {
            // for the current observation, find location of corresponding unique values
            X_counts_index++;
          }

          if (*(temp_pointer + torder_std[i][j]) <= cutvalue)
          {
            // go to left side
            split_left->torder_std[i][left_ix] = torder_std[i][j];
            left_ix = left_ix + 1;
            split_left->t_counts[X_counts_index]++;
          } else {
            // go to right side

            split_right->torder_std[i][right_ix] = torder_std[i][j];
            right_ix = right_ix + 1;

            split_right->t_counts[X_counts_index]++;
          }
        }
      }
    }

    model->calculateOtherSideSuffStat(suff_stat, left_suff_stat, right_suff_stat, compute_left_side);

    // update X_num_unique

    std::fill(split_left->t_num_unique.begin(), split_left->t_num_unique.end(), 0.0);
    std::fill(split_right->t_num_unique.begin(), split_right->t_num_unique.end(), 0.0);

    for (size_t i = 0; i < x_struct->p_t; i++)
    {
      start = x_struct->t_variable_ind[i];
      end = x_struct->t_variable_ind[i + 1];

      // COUT << "start " << start << " end " << end << " size " << X_counts_left.size() << endl;
      for (size_t j = start; j < end; j++)
      {
        if (split_left->t_counts[j] > 0)
        {
          split_left->t_num_unique[i] = split_left->t_num_unique[i] + 1;
        }
        if (split_right->t_counts[j] > 0)
        {
          split_right->t_num_unique[i] = split_right->t_num_unique[i] + 1;
        }
      }
    }

  }

};

#endif  // THIRD_PARTY_R_PACKAGES_LONGBET_LONGBET_SRC_SPLIT_INFO_H_
