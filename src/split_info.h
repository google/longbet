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
          
          this->s_values.resize(split_point + 1);
          std::copy(parent->s_values.begin(), parent->s_values.begin() + split_point, this->s_values.begin());

        } else {
          ini_xinfo_sizet(this->torder_std, parent->N_torder - split_point - 1, parent->p_torder);
          
          this->s_values.resize(parent->s_values.size() - split_point - 1);
          std::copy(parent->s_values.begin() + split_point + 1, parent->s_values.end(), this->s_values.begin());
        }
      } else {
        ini_xinfo_sizet(this->torder_std, parent->N_torder, parent->p_torder);
        // copy torder
        for (size_t i = 0; i < torder_std.size(); i++){
          std::copy(parent->torder_std[i].begin(), parent->torder_std[i].end(), this->torder_std[i].begin());
          std::copy(parent->torder_std[i].begin(), parent->torder_std[i].end(), this->torder_std[i].begin());
        }

        this->s_values = parent->s_values;
        
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

    void split_this(std::unique_ptr<split_info> &split_left,
      std::unique_ptr<split_info> &split_right,
      size_t split_var, size_t split_point, bool split_t,
      Model *model, std::unique_ptr<X_struct> &x_struct,
      std::unique_ptr<State> &state, std::vector<double> &suff_stat,
      std::vector<double> &left_suff_stat, std::vector<double> &right_suff_stat);

    void split_xorder_std_continuous(std::unique_ptr<split_info> &split_left,
      std::unique_ptr<split_info> &split_right, size_t split_var, size_t split_point,
      Model *model, std::unique_ptr<X_struct> &x_struct,
      std::unique_ptr<State> &state, std::vector<double> &suff_stat,
      std::vector<double> &left_suff_stat, std::vector<double> &right_suff_stat);

    void split_xorder_std_categorical(std::unique_ptr<split_info> &split_left,
      std::unique_ptr<split_info> &split_right,
      size_t split_var, size_t split_point,
      Model *model, std::unique_ptr<X_struct> &x_struct,
      std::unique_ptr<State> &state, std::vector<double> &suff_stat,
      std::vector<double> &left_suff_stat, std::vector<double> &right_suff_stat);

    void split_torder_std(std::unique_ptr<split_info> &split_left,
      std::unique_ptr<split_info> &split_right,
      size_t split_var, size_t split_point,
      Model *model, std::unique_ptr<X_struct> &x_struct,
      std::unique_ptr<State> &state, std::vector<double> &suff_stat,
      std::vector<double> &left_suff_stat, std::vector<double> &right_suff_stat);

};

#endif  // THIRD_PARTY_R_PACKAGES_LONGBET_LONGBET_SRC_SPLIT_INFO_H_
