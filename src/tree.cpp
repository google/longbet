#include "tree.h"
// #include <RcppArmadilloExtensions/sample.h>
#include <chrono>
#include <cstddef>
#include <memory>
#include <ostream>

using namespace std;
using namespace chrono;

//--------------------
// node id
size_t tree::nid() const
{
    if (!p)
        return 1; //if you don't have a parent, you are the top
    if (this == p->l)
        return 2 * (p->nid()); //if you are a left child
    else
        return 2 * (p->nid()) + 1; //else you are a right child
}
//--------------------
tree::tree_p tree::getptr(size_t nid)
{
    if (this->nid() == nid)
        return this; //found it
    if (l == 0)
        return 0; //no children, did not find it
    tree_p lp = l->getptr(nid);
    if (lp)
        return lp; //found on left
    tree_p rp = r->getptr(nid);
    if (rp)
        return rp; //found on right
    return 0;      //never found it
}
//--------------------

//--------------------
//depth of node
// size_t tree::depth()
// {
//     if (!p)
//         return 0; //no parents
//     else
//         return (1 + p->depth());
// }
//--------------------
//tree size
size_t tree::treesize()
{
    if (l == 0)
        return 1; //if bottom node, tree size is 1
    else
        return (1 + l->treesize() + r->treesize());
}
//--------------------
//node type
char tree::ntype()
{
    //t:top, b:bottom, n:no grandchildren, i:internal
    if (!p)
        return 't';
    if (!l)
        return 'b';
    if (!(l->l) && !(r->l))
        return 'n';
    return 'i';
}
//--------------------
//print out tree(pc=true) or node(pc=false) information
void tree::pr(bool pc)
{
    size_t d = this->depth;
    size_t id = nid();

    size_t pid;
    if (!p)
        pid = 0; //parent of top node
    else
        pid = p->nid();

    std::string pad(2 * d, ' ');
    std::string sp(", ");
    if (pc && (ntype() == 't'))
        COUT << "tree size: " << treesize() << std::endl;
    COUT << pad << "(id,parent): " << id << sp << pid;
    COUT << sp << "(v,c): " << v << sp << c;
    COUT << sp << "theta: " << theta_vector;
    COUT << sp << "type: " << ntype();
    COUT << sp << "depth: " << this->depth;
    COUT << sp << "pointer: " << this << std::endl;

    if (pc)
    {
        if (l)
        {
            l->pr(pc);
            r->pr(pc);
        }
    }
}

//--------------------
//is the node a nog node
bool tree::isnog()
{
    bool isnog = true;
    if (l)
    {
        if (l->l || r->l)
            isnog = false; //one of the children has children.
    }
    else
    {
        isnog = false; //no children
    }
    return isnog;
}
//--------------------
size_t tree::nnogs()
{
    if (!l)
        return 0; //bottom node
    if (l->l || r->l)
    { //not a nog
        return (l->nnogs() + r->nnogs());
    }
    else
    { //is a nog
        return 1;
    }
}
//--------------------
size_t tree::nbots()
{
    if (l == 0)
    { //if a bottom node
        return 1;
    }
    else
    {
        return l->nbots() + r->nbots();
    }
}
//--------------------
//get bottom nodes
void tree::getbots(npv &bv)
{
    if (l)
    { //have children
        l->getbots(bv);
        r->getbots(bv);
    }
    else
    {
        bv.push_back(this);
    }
}
//--------------------
//get nog nodes
void tree::getnogs(npv &nv)
{
    if (l)
    { //have children
        if ((l->l) || (r->l))
        { //have grandchildren
            if (l->l)
                l->getnogs(nv);
            if (r->l)
                r->getnogs(nv);
        }
        else
        {
            nv.push_back(this);
        }
    }
}
//--------------------
//get pointer to the top tree
tree::tree_p tree::gettop()
{
    if (!p)
    {
        return this;
    }
    else
    {
        return p->gettop();
    }
}
//--------------------
//get all nodes
void tree::getnodes(npv &v)
{
    v.push_back(this);
    if (l)
    {
        l->getnodes(v);
        r->getnodes(v);
    }
}
void tree::getnodes(cnpv &v) const
{
    v.push_back(this);
    if (l)
    {
        l->getnodes(v);
        r->getnodes(v);
    }
}
//--------------------
tree::tree_p tree::bn(double *x, matrix<double> &xi)
{

    // original BART function, v and c are index of split point in matrix<double>& xi

    if (l == 0)
        return this; //no children
    if (x[v] <= xi[v][c])
    {
        // if smaller than or equals to the cutpoint, go to left child

        return l->bn(x, xi);
    }
    else
    {
        // if greater than cutpoint, go to right child
        return r->bn(x, xi);
    }
}

tree::tree_p tree::bn_std(double *x)
{
    // v is variable to split, c is raw value
    // not index in matrix<double>, so compare x[v] with c directly

    if (l == 0)
        return this;
    if (x[v] <= c)
    {
        return l->bn_std(x);
    }
    else
    {
        return r->bn_std(x);
    }
}

tree::tree_p tree::search_bottom_std(const double *X, const double *t, const size_t &i, const size_t &j, const size_t &N, const size_t &N_t)
{
    // X is a matrix, std vector of vectors, stack by column, N rows and p columns
    // // i is index of row in X to predict
    // cout << "v = " << v << ", c = " << c << ", no split = " << this->no_split <<  ", split_t = " << this->split_t << endl;
    // cout << "value = " << *(X + N * v + i) << ", l = " << l <<   endl;
    
    if (this->l)
    {
        if (this->split_t)
        {
            if (*(t + N_t * v + j) <= c)
            {
                return l->search_bottom_std(X, t, i, j, N, N_t);
            } else {
                return r->search_bottom_std(X, t, i, j, N, N_t);
            }
        } else {
            // X[v][i], v-th column and i-th row
            // if(X[v][i] <= c){
            if (*(X + N * v + i) <= c)
            {
                return l->search_bottom_std(X, t, i, j, N, N_t);
            } else {
                return r->search_bottom_std(X, t, i, j, N, N_t);
            }
        }
    } else {
        return this;
    }
}

//--------------------
//find region for a given variable
void tree::rg(size_t v, size_t *L, size_t *U)
{
    if (this->p == 0)
    {
        return;
    }
    if ((this->p)->v == v)
    { //does my parent use v?
        if (this == p->l)
        { //am I left or right child
            if ((size_t)(p->c) <= (*U))
                *U = (p->c) - 1;
            p->rg(v, L, U);
        }
        else
        {
            if ((size_t)(p->c) >= *L)
                *L = (p->c) + 1;
            p->rg(v, L, U);
        }
    }
    else
    {
        p->rg(v, L, U);
    }
}
//--------------------
//cut back to one node
void tree::tonull()
{
    size_t ts = treesize();
    //loop invariant: ts>=1
    while (ts > 1)
    { //if false ts=1
        npv nv;
        getnogs(nv);
        for (size_t i = 0; i < nv.size(); i++)
        {
            delete nv[i]->l;
            delete nv[i]->r;
            nv[i]->l = 0;
            nv[i]->r = 0;
        }
        ts = treesize(); //make invariant true
    }
    v = 0;
    c = 0;
    p = 0;
    l = 0;
    r = 0;
}
//--------------------
//copy tree tree o to tree n
void tree::cp(tree_p n, tree_cp o)
//assume n has no children (so we don't have to kill them)
//recursion down
// create a new copy of tree in NEW memory space
{
    if (n->l)
    {
        COUT << "cp:error node has children\n";
        return;
    }

    n->v = o->v;
    n->c = o->c;
    n->prob_split = o->prob_split;
    n->prob_leaf = o->prob_leaf;
    n->drawn_ind = o->drawn_ind;
    n->loglike_node = o->loglike_node;
    n->tree_like = o->tree_like;
    n->theta_vector = o->theta_vector;
    n->no_split = o->no_split;
    n->split_t = o->split_t;

    if (o->l)
    { //if o has children
        n->l = new tree;
        (n->l)->p = n;
        cp(n->l, o->l);
        n->r = new tree;
        (n->r)->p = n;
        cp(n->r, o->r);
    }
}

void tree::copy_only_root(tree_p o)
//assume n has no children (so we don't have to kill them)
//NOT LIKE cp() function
//this function pointer new root to the OLD structure
{
    this->v = o->v;
    this->c = o->c;
    this->prob_split = o->prob_split;
    this->prob_leaf = o->prob_leaf;
    this->drawn_ind = o->drawn_ind;
    this->loglike_node = o->loglike_node;
    this->tree_like = o->tree_like;
    this->theta_vector = o->theta_vector;

    if (o->l)
    {
        // keep the following structure, rather than create a new tree in memory
        this->l = o->l;
        this->r = o->r;
        // also update pointers to parents
        this->l->p = this;
        this->r->p = this;
    }
    else
    {
        this->l = 0;
        this->r = 0;
    }
}

json tree::to_json()
{
    json j;
    if (l == 0)
    {
        j = this->theta_vector;
    }
    else
    {
        j["variable"] = this->v;
        j["cutpoint"] = this->c;
        j["nodeid"] = this->nid();
        j["left"] = this->l->to_json();
        j["right"] = this->r->to_json();
    }
    return j;
}

void tree::from_json(json &j3, size_t dim_theta)
{
    if (j3.is_array())
    {
        std::vector<double> temp;
        j3.get_to(temp);
        if (temp.size() > 1)
        {
            this->theta_vector = temp;
        }
        else
        {
            this->theta_vector[0] = temp[0];
        }
    }
    else
    {
        j3.at("variable").get_to(this->v);
        j3.at("cutpoint").get_to(this->c);

        tree *lchild = new tree(dim_theta);
        lchild->from_json(j3["left"], dim_theta);
        tree *rchild = new tree(dim_theta);
        rchild->from_json(j3["right"], dim_theta);

        lchild->p = this;
        rchild->p = this;
        this->l = lchild;
        this->r = rchild;
    }
}

//--------------------------------------------------
//operators
tree &tree::operator=(const tree &rhs)
{
    if (&rhs != this)
    {
        tonull();       //kill left hand side (this)
        cp(this, &rhs); //copy right hand side to left hand side
    }
    return *this;
}
//--------------------------------------------------
std::ostream& operator<<(std::ostream& os, const tree& t)
{
   tree::cnpv nds;
   t.getnodes(nds);
   os << nds.size() << std::endl;
   for(size_t i=0;i<nds.size();i++) {
      os << nds[i]->nid() << " ";
      os << nds[i]->getv() << " ";
      os << nds[i]->getc_index() << " ";
      os << nds[i]->getc() << " ";
      os << nds[i]->theta_vector[0] << std::endl;
   }
   return os;
}

std::istream &operator>>(std::istream &is, tree &t)
{
    size_t tid, pid;                    //tid: id of current node, pid: parent's id
    std::map<size_t, tree::tree_p> pts; //pointers to nodes indexed by node id
    size_t nn;                          //number of nodes

    t.tonull(); // obliterate old tree (if there)

    //read number of nodes----------
    is >> nn;
    if (!is)
    {
        return is;
    }

    // The idea is to dump string to a lot of node_info structure first, then link them as a tree, by nid

    //read in vector of node information----------
    std::vector<node_info> nv(nn);
    for (size_t i = 0; i != nn; i++)
    {
        is >> nv[i].id >> nv[i].v >> nv[i].c_index >> nv[i].c >> nv[i].theta_vector[0]; // Only works on first theta for now, fix latex if needed
        if (!is)
        {
            return is;
        }
    }

    //first node has to be the top one
    pts[1] = &t; //be careful! this is not the first pts, it is pointer of id 1.
    t.setv(nv[0].v);
    t.setc(nv[0].c);
    t.setc_index(nv[0].c_index);
    t.settheta(nv[0].theta_vector);
    t.p = 0;

    //now loop through the rest of the nodes knowing parent is already there.
    for (size_t i = 1; i != nv.size(); i++)
    {
        tree::tree_p np = new tree;
        np->v = nv[i].v;
        np->c_index = nv[i].c_index;
        np->c = nv[i].c;
        np->theta_vector = nv[i].theta_vector;
        tid = nv[i].id;
        pts[tid] = np;
        pid = tid / 2;
        if (tid % 2 == 0)
        { //left child has even id
            pts[pid]->l = np;
        }
        else
        {
            pts[pid]->r = np;
        }
        np->p = pts[pid];
    }
    return is;
}


void tree::grow_from_root(std::unique_ptr<State> &state,
std::unique_ptr<split_info> &split, Model *model,
std::unique_ptr<X_struct> &x_struct, const size_t &sweeps,
const size_t &tree_ind, bool control_split_t)
{
    // grow a tree, users can control number of split points
    size_t N_Xorder = split->Xorder_std[0].size();
    size_t N_t = split->s_values.size();
    size_t p = split->Xorder_std.size();
    size_t ind;
    size_t split_var;
    size_t split_point;
    // cout << "grow_from_root" << ", N_x = " << N_Xorder << ", N_t =  " << N_t << ", p = " << p << endl;

    this->N = N_Xorder;

    // tau is prior VARIANCE, do not take squares

    // update_theta
    // cout << "depth = " << this->depth <<  " suff_stat " << this->suff_stat << endl;
    model->samplePars(state, this->suff_stat, this->theta_vector, this->prob_leaf);

    if (N_Xorder <= state->n_min)
    {
        return;
    }

    if (this->depth >= state->max_depth - 1)
    {
        return;
    }

    std::vector<size_t> subset_vars(p);

    if (state->use_all)
    {
        std::iota(subset_vars.begin(), subset_vars.end(), 0);
    } else {
        // TODO: add time variable for sampling
        if (state->sample_weights_flag)
        {
            std::vector<double> weight_samp(p);
            double weight_sum;

            // Sample Weights Dirchelet
            for (size_t i = 0; i < p; i++)
            {
                std::gamma_distribution<double> temp_dist(state->mtry_weight_current_tree[i], 1.0);
                weight_samp[i] = temp_dist(state->gen);
            }
            weight_sum = accumulate(weight_samp.begin(), weight_samp.end(), 0.0);
            for (size_t i = 0; i < p; i++)
            {
                weight_samp[i] = weight_samp[i] / weight_sum;
            }

            subset_vars = sample_int_ccrank(p, state->mtry, weight_samp, state->gen);
        }
        else
        {
            subset_vars = sample_int_ccrank(p, state->mtry, state->mtry_weight_current_tree, state->gen);
        }
    }

    // cout << "start BART" << endl;
    BART_likelihood_all(split, this->no_split, split_var, split_point, subset_vars, model, x_struct, state, this, this->split_t, control_split_t);
    // cout << "no split = " << this->no_split << ", split_t = " << this->split_t << ", var = " << split_var << ", split_point  = " << split_point << endl;

    this->loglike_node = model->likelihood(this->suff_stat, this->suff_stat,
    false, true, state);
    // cout << "finish BART" << endl;

    if (!this->no_split)
    {
        this->v = split_var;

        if (this->split_t)
        { 
            int place_holder = 1;
            // cout << "split_t, split_point = " << split_point << " value = " << *(x_struct->t_std + x_struct->n_t * split_var + split->torder_std[split_var][split_point]) << endl;
            // this->c = *(x_struct->t_std + x_struct->n_t * split_var + split->torder_std[split_var][split_point]);
            // while ((split_point < N_torder - 1) && (*(x_struct->t_std + x_struct->n_t * split_var + split->torder_std[split_var][split_point + 1]) == this->c))
            // {
            //     split_point = split_point + 1;
            // }

            // // If our current split is same as parent, exit
            // if ( (this->p) && ((this->p)->split_t) && (this->v == (this->p)->v) && (this->c == (this->p)->c))
            // {
            //     split_t = false;
            //     no_split = true;
            // }
            // cout << "split point = " << split_point << endl;
        } else {
            // cout << "split on X" << endl;
            this->c = *(state->X_std + state->n_y * split_var + split->Xorder_std[split_var][split_point]);
            // Update Cutpoint to be a true seperating point
            // Increase split_point (index) until it is no longer equal to cutpoint value
            while ((split_point < N_Xorder - 1) && (*(state->X_std + state->n_y * split_var + split->Xorder_std[split_var][split_point + 1]) == this->c))
            {
                split_point = split_point + 1;
            }
            // If our current split is same as parent, exit
            if ( (this->p) && ((this->p)->split_t == false) && (this->v == (this->p)->v) && (this->c == (this->p)->c))
            {
                no_split = true;
            }
            state->split_count_current_tree[split_var] += 1;
        }
        
    }

    if (this->no_split == true)
    {
        // cout << "no split " << endl;
        for (size_t i = 0; i < N_Xorder; i++)
        {
            for (auto j: split->sorder_std[i])
            {
                x_struct->data_pointers[tree_ind][split->Xorder_std[0][i] * state->p_y + j] = &this->theta_vector;
            }
        }

        this->l = 0;
        this->r = 0;

        return;
    }

    // If do not update split prob ONLY
    // grow from root, initialize new nodes

    tree::tree_p lchild = new tree(model->getNumClasses(), this, model->dim_suffstat);
    tree::tree_p rchild = new tree(model->getNumClasses(), this, model->dim_suffstat);

    this->l = lchild;
    this->r = rchild;

    lchild->depth = this->depth + 1;
    rchild->depth = this->depth + 1;

    lchild->ID = 2 * (this->ID);
    rchild->ID = lchild->ID + 1;

    this->l->ini_suff_stat();
    this->r->ini_suff_stat();

    // cout << "create split info " << endl;
    std::unique_ptr<split_info> split_left( new split_info(split, split_var,
    split_point, split_t, true));
    std::unique_ptr<split_info> split_right( new split_info(split, split_var,
    split_point, split_t, false));

    std::vector<double> left_suff_stat(this->l->suff_stat.size(), 0.0);
    std::vector<double> right_suff_stat(this->r->suff_stat.size(), 0.0);

    // cout << "split this" << endl;
    // TODO: potential bug using split_info structure without categorical variable
    split->split_this(split_left, split_right, split_var, split_point, split_t,
    model, x_struct, state, this->suff_stat, left_suff_stat, right_suff_stat);
    // cout << "end split" << endl;

    // cout << "split_t = " << split_t << ", split_var = " << split_var << ", split_point = " << split_point << ", depth = " << this->depth + 1 << endl;
    // cout << "left suff stat = " << left_suff_stat << endl;
    // cout << "right suff stat = " << right_suff_stat << endl;

    std::copy(left_suff_stat.begin(), left_suff_stat.end(), this->l->suff_stat.begin());
    std::copy(right_suff_stat.begin(), right_suff_stat.end(), this->r->suff_stat.begin());

    this->l->grow_from_root(state, split_left, model, x_struct, sweeps, tree_ind, control_split_t);

    this->r->grow_from_root(state, split_right, model, x_struct, sweeps, tree_ind, control_split_t);
}

void BART_likelihood_all(std::unique_ptr<split_info> &split_info,
 bool &no_split, size_t &split_var, size_t &split_point,
 const std::vector<size_t> &subset_vars, Model *model,
 std::unique_ptr<X_struct> &x_struct, std::unique_ptr<State> &state,
 tree *tree_pointer, bool &split_t, bool control_split_t)
{
    // cout << "start BART likelihood" << endl;
    // compute BART posterior (loglikelihood + logprior penalty)

    // subset_vars: a vector of indexes of varibles to consider (like random forest)

    // use stacked vector loglike instead of a matrix, stacked by column
    // length of loglike is p * (N - 1) + 1
    // N - 1 has to be greater than 2 * Nmin

    size_t N = split_info->Xorder_std[0].size();
    size_t p = split_info->Xorder_std.size();
    size_t ind;
    size_t N_Xorder = N;
    size_t total_categorical_split_candidates = 0;

    // double sigma2 = pow(sigma, 2);

    double loglike_max = -INFINITY;

    std::vector<double> loglike;

    size_t loglike_start;
    size_t loglike_time_start;
    size_t loglike_t_size;

    if (control_split_t){
        loglike_t_size = split_info->s_values.size();
    } else {
        loglike_t_size = 0;
    }

    // decide lenght of loglike vector
    if (N <= state->n_cutpoints + 1 + 2 * state->n_min)
    {
        loglike.resize((N_Xorder - 1) * state->p_continuous + x_struct->X_values.size() + loglike_t_size + 1, -INFINITY);
        loglike_start = (N_Xorder - 1) * state->p_continuous;
    }
    else
    {
        loglike.resize(state->n_cutpoints * state->p_continuous + x_struct->X_values.size() + loglike_t_size + 1, -INFINITY);
        loglike_start = state->n_cutpoints * state->p_continuous;
    }
    loglike_time_start = loglike_start + x_struct->X_values.size();

    // cout << "loglike_start = " << loglike_start << ", x_values = " << x_struct->X_values.size() << ", loglike_t_size = " << loglike_t_size << ", loglike size = " << loglike.size() <<endl;

    // calculate for each cases
    if (state->p_continuous > 0)
    {
        // cout << "calculate_loglikelihood_continuous" << endl;
        calculate_loglikelihood_continuous(loglike, subset_vars, N_Xorder, split_info->Xorder_std, split_info->sorder_std, loglike_max, model, x_struct, state, tree_pointer);
        // cout << "finish " << endl;
    }

    if (state->p_categorical > 0)
    {
        // cout << "calculate_loglikelihood_categorical" << endl;
        calculate_loglikelihood_categorical(loglike, loglike_start, subset_vars, N_Xorder, split_info->Xorder_std, split_info->sorder_std, loglike_max, split_info->X_counts, split_info->X_num_unique, model, x_struct, total_categorical_split_candidates, state, tree_pointer);
        // cout << "finish" << endl;
    }

    if (control_split_t)
    {
        // cout << "calc time like "  << endl;
        calculate_loglikelihood_time(loglike, loglike_time_start, loglike_max,
        model, x_struct, split_info, state, tree_pointer);
        // cout << "finish " << endl;
    }

    // calculate likelihood of no-split option
    calculate_likelihood_no_split(loglike, N_Xorder, loglike_max, model, x_struct, total_categorical_split_candidates, state, tree_pointer);

    // transfer loglikelihood to likelihood
    for (size_t ii = 0; ii < loglike.size(); ii++)
    {
        // if a variable is not selected, take exp will becomes 0
        loglike[ii] = exp(loglike[ii] - loglike_max);
    }
    
    // cout << "loglike = " << loglike << endl;

    // sampling cutpoints
    if (N <= state->n_cutpoints + 1 + 2 * state->n_min)
    {
        // N - 1 - 2 * Nmin <= Ncutpoints, consider all data points

        // if number of observations is smaller than Ncutpoints, all data are splitpoint candidates
        // note that the first Nmin and last Nmin cannot be splitpoint candidate

        if ((N - 1) > 2 * state->n_min)
        {
            // for(size_t i = 0; i < p; i ++ ){
            for (auto &&i : subset_vars)
            {
                if (i < state->p_continuous)
                {
                    // delete some candidates, otherwise size of the new node can be smaller than Nmin
                    std::fill(loglike.begin() + i * (N - 1), loglike.begin() + i * (N - 1) + state->n_min + 1, 0.0);
                    std::fill(loglike.begin() + i * (N - 1) + N - 2 - state->n_min, loglike.begin() + i * (N - 1) + N - 2 + 1, 0.0);
                }
            }
        }
        else
        {
            // do not use all continuous variables
            std::fill(loglike.begin(), loglike.begin() + (N_Xorder - 1) * state->p_continuous - 1, 0.0);
        }
        std::discrete_distribution<> d(loglike.begin(), loglike.end());

        // for MH update usage only
        tree_pointer->num_cutpoint_candidates = count_non_zero(loglike);

        // sample one index of split point
        ind = d(state->gen);
        tree_pointer->drawn_ind = ind;
        // cout << "ind = " << ind << endl;

        // save the posterior of the chosen split point
        vec_sum(loglike, tree_pointer->prob_split);
        tree_pointer->prob_split = loglike[ind] / tree_pointer->prob_split;
        if (ind == loglike.size() - 1)
        {
            // no split
            no_split = true;
            split_var = 0;
            split_point = 0;
        }
        else if ((N - 1) <= 2 * state->n_min)
        {
            // np split

            /////////////////////////////////
            //
            // Need optimization, move before calculating likelihood
            //
            /////////////////////////////////

            no_split = true;
            split_var = 0;
            split_point = 0;
        }
        // else if (ind < loglike_start - x_struct->X_values.size())
        else if (ind < loglike_start)
        {
            // split at continuous variable
            split_var = ind / (N - 1);
            split_point = ind % (N - 1);
        }
        // else if (ind < loglike_start)
        else if (ind < loglike_time_start)
        {
            // split at categorical variable
            size_t start;
            // ind = ind - loglike_start - x_struct->X_values.size();
            ind = ind - loglike_start;
            for (size_t i = 0; i < (x_struct->variable_ind.size() - 1); i++)
            {
                if (x_struct->variable_ind[i] <= ind && x_struct->variable_ind[i + 1] > ind)
                {
                    split_var = i;
                }
            }
            start = x_struct->variable_ind[split_var];
            // count how many
            split_point = std::accumulate(split_info->X_counts.begin() + start, split_info->X_counts.begin() + ind + 1, 0);
            // minus one for correct index (start from 0)
            split_point = split_point - 1;
            split_var = split_var + state->p_continuous;
        } else {
            // split at time variable
            split_t = true;
            size_t start;
            // cout << "ind = " << ind << ", loglike_time_start = " << loglike_time_start << endl;
            ind = ind - loglike_time_start;
            split_var = 0;
            split_point = ind;
            // for (size_t i = 0; i < (x_struct->t_variable_ind.size() - 1); i++)
            // {
            //     if (x_struct->t_variable_ind[i] <= ind && x_struct->t_variable_ind[i + 1] > ind)
            //     {
            //         split_var = i;
            //     }
            // }
            // start = x_struct->t_variable_ind[split_var];
            // // count how many
            // split_point = std::accumulate(split_info->t_counts.begin() + start, split_info->t_counts.begin() + ind + 1, 0);
            // cout << "split_point = " << split_point << endl;
            // // minus one for correct index (start from 0)
            // split_point = split_point - 1;
        }
    }
    else
    {
        // use adaptive number of cutpoints

        std::vector<size_t> candidate_index(state->n_cutpoints);

        seq_gen_std(state->n_min, N - state->n_min, state->n_cutpoints, candidate_index);

        std::discrete_distribution<size_t> d(loglike.begin(), loglike.end());

        // For MH update usage only
        tree_pointer->num_cutpoint_candidates = count_non_zero(loglike);

        // // sample one index of split point
        ind = d(state->gen);
        tree_pointer->drawn_ind = ind;
        // cout << "loglike size = " << loglike.size() << ", no_split = " << loglike[loglike.size() - 1] << ", ind = " << ind << endl;

        // save the posterior of the chosen split point
        vec_sum(loglike, tree_pointer->prob_split);
        tree_pointer->prob_split = loglike[ind] / tree_pointer->prob_split;

        if (ind == loglike.size() - 1)
        {
            // no split
            no_split = true;
            split_var = 0;
            split_point = 0;
        } else if (ind < loglike_start) {
            // split at continuous variable
            split_var = ind / state->n_cutpoints;
            split_point = candidate_index[ind % state->n_cutpoints];
        } else if (ind < loglike_time_start) {
            // split at categorical variable
            size_t start;
            // ind = ind - loglike_start - x_struct->X_values.size();
            ind = ind - loglike_start;
            for (size_t i = 0; i < (x_struct->variable_ind.size() - 1); i++)
            {
                if (x_struct->variable_ind[i] <= ind && x_struct->variable_ind[i + 1] > ind)
                {
                    split_var = i;
                }
            }
            start = x_struct->variable_ind[split_var];
            // count how many
            split_point = std::accumulate(split_info->X_counts.begin() + start, split_info->X_counts.begin() + ind + 1, 0);
            // minus one for correct index (start from 0)
            split_point = split_point - 1;
            split_var = split_var + state->p_continuous;
        } else {
            // split at time variable
            split_t = true;
            size_t start;
            ind = ind - loglike_time_start;
            split_var = 0;
            split_point = ind;
            // for (size_t i = 0; i < (x_struct->t_variable_ind.size() - 1); i++)
            // {
            //     if (x_struct->t_variable_ind[i] <= ind && x_struct->t_variable_ind[i + 1] > ind)
            //     {
            //         split_var = i;
            //     }
            // }
            // start = x_struct->t_variable_ind[split_var];
            // // count how many
            // split_point = std::accumulate(split_info->t_counts.begin() + start, split_info->t_counts.begin() + ind + 1, 0);
            // // minus one for correct index (start from 0)
            // split_point = split_point - 1;
        }
    }
}

void calculate_loglikelihood_continuous(std::vector<double> &loglike,
const std::vector<size_t> &subset_vars, size_t &N_Xorder,
matrix<size_t> &Xorder_std, matrix<size_t> &torder_std,
double &loglike_max, Model *model,
std::unique_ptr<X_struct> &x_struct, std::unique_ptr<State> &state,
tree *tree_pointer)
{
    size_t N = N_Xorder;

    std::vector<double> temp_suff_stat(model->dim_suffstat);
    std::vector<double> temp_suff_stat2(model->dim_suffstat);

    if (N_Xorder <= state->n_cutpoints + 1 + 2 * state->n_min)
    {
        // if we only have a few data observations in current node
        // use all of them as cutpoint candidates

        double n1tau;
        double n2tau;
        // double Ntau = N_Xorder * model->tau;

        // to have a generalized function, have to pass an empty candidate_index object for this case
        // is there any smarter way to do it?
        std::vector<size_t> candidate_index(1);

        for (auto &&i : subset_vars)
        {
            if (i < state->p_continuous)
            {
                std::vector<size_t> &xorder = Xorder_std[i];

                // initialize sufficient statistics
                std::fill(temp_suff_stat.begin(), temp_suff_stat.end(), 0.0);

                ////////////////////////////////////////////////////////////////
                //
                //  This part can be run in parallel, just like continuous case below, Ncutpoint case
                //
                //  If run in parallel, need to redefine model class for each thread
                //
                ////////////////////////////////////////////////////////////////

                for (size_t j = 0; j < N_Xorder - 1; j++)
                {
                    calcSuffStat_continuous(temp_suff_stat, xorder, torder_std, candidate_index, j, false, model, state);

                    loglike[(N_Xorder - 1) * i + j] = model->likelihood(temp_suff_stat, tree_pointer->suff_stat, true, false, state) + model->likelihood(temp_suff_stat, tree_pointer->suff_stat, false, false, state);

                    if (loglike[(N_Xorder - 1) * i + j] > loglike_max)
                    {
                        loglike_max = loglike[(N_Xorder - 1) * i + j];
                    }
                }
            }
        }
    }
    else
    {

        // otherwise, adaptive number of cutpoints
        // use Ncutpoints

        std::vector<size_t> candidate_index2(state->n_cutpoints + 1);
        seq_gen_std2(state->n_min, N - state->n_min, state->n_cutpoints, candidate_index2);

        // double Ntau = N_Xorder * model->tau;

        std::mutex llmax_mutex;

        for (auto &&i : subset_vars)
        {
            if (i < state->p_continuous)
            {

                // Lambda callback to perform the calculation
                auto calcllc_i = [i, &loglike, &loglike_max, &Xorder_std, &torder_std, &state, &candidate_index2, &model, &llmax_mutex, N_Xorder, &tree_pointer]() {
                    std::vector<size_t> &xorder = Xorder_std[i];
                    double llmax = -INFINITY;

                    std::vector<double> temp_suff_stat(model->dim_suffstat);

                    std::fill(temp_suff_stat.begin(), temp_suff_stat.end(), 0.0);

                    for (size_t j = 0; j < state->n_cutpoints; j++)
                    {

                        calcSuffStat_continuous(temp_suff_stat, xorder, torder_std, candidate_index2, j, true, model, state);

                        loglike[(state->n_cutpoints) * i + j] = model->likelihood(temp_suff_stat, tree_pointer->suff_stat, true, false, state) + model->likelihood(temp_suff_stat, tree_pointer->suff_stat, false, false, state);

                        if (loglike[(state->n_cutpoints) * i + j] > llmax)
                        {
                            llmax = loglike[(state->n_cutpoints) * i + j];
                        }
                    }
                    llmax_mutex.lock();
                    if (llmax > loglike_max)
                        loglike_max = llmax;
                    llmax_mutex.unlock();
                };

                if (thread_pool.is_active())
                    thread_pool.add_task(calcllc_i);
                else
                    calcllc_i();
            }
        }
        if (thread_pool.is_active())
            thread_pool.wait();
    }
}

void calculate_loglikelihood_time(std::vector<double> &loglike,
size_t &loglike_start, double &loglike_max, Model *model,
std::unique_ptr<X_struct> &x_struct, std::unique_ptr<split_info> &split_info,
std::unique_ptr<State> &state, tree *tree_pointer)
{

//     std::vector<double> loglike_copy = loglike;

//     size_t N = split_info->Xorder_std[0].size();
//     size_t start;
//     size_t end;
//     size_t end2;
//     size_t n1;
//     size_t temp;
//     size_t effective_cutpoints = 0;

//     std::vector<double> temp_suff_stat(tree_pointer->suff_stat.size());
//     std::vector<double> total_suff_stat(model->dim_suffstat);
//     std::vector<double> left_suff_stat(model->dim_suffstat);

//     for (size_t i = 0; i < split_info->torder_std.size(); i++)
//     {
//         if (split_info->t_num_unique[i] > 1)
//         {
//             // more than one unique values
//             start = x_struct->t_variable_ind[i];
//             end = x_struct->t_variable_ind[i + 1] - 1;
//             end2 = end;

//             while (split_info->t_counts[end2] == 0)
//             {
//                 end2 = end2 - 1;
//             }
//             end2 = end2 - 1;

//             std::fill(temp_suff_stat.begin(), temp_suff_stat.end(), 0.0);

//             ////////////////////////////////////////////////////////////////
//             //
//             //  This part can be run in parallel, just like continuous case
//             //
//             //  If run in parallel, need to redefine model class for each thread
//             //
//             ////////////////////////////////////////////////////////////////

//             n1 = 0;

//             for (size_t j = start; j <= end2; j++)
//             {

//                 if (split_info->t_counts[j] != 0)
//                 {

//                     temp = n1 + split_info->t_counts[j] - 1;
//                     calcSuffStat_time(temp_suff_stat, split_info->Xorder_std[0],
//                     split_info->torder_std[i], n1, temp, model, state);
//                     n1 = n1 + split_info->t_counts[j];

//                     loglike_copy[loglike_start + j] = model->likelihood(temp_suff_stat, tree_pointer->suff_stat, true, false, state)
//                     + model->likelihood(temp_suff_stat, tree_pointer->suff_stat, false, false, state);              
                    
//                     // count total number of cutpoint candidates
//                     effective_cutpoints++;

//                     // if (loglike[loglike_start + j] > loglike_max)
//                     // {
//                     //     loglike_max = loglike[loglike_start + j];
//                     // }
//                 }
//             }
//         }
//     } 

    std::vector<double> temp_suff_stat(tree_pointer->suff_stat.size());
    std::fill(temp_suff_stat.begin(), temp_suff_stat.end(), 0.0);

    for (size_t s = 0; s <split_info->s_values.size() - 1; s++){
        
        for (auto i: split_info->Xorder_std[0])
        {
            if (split_info->sorder_std[i].size() == 0) {continue;}

            for (auto j: split_info->sorder_std[i]){ 

                // this can be improved.
                if (x_struct->Tpt[i + j*state->n_y] == split_info->s_values[s]){ // ==? <=
                    model->incSuffStat(state, i, j, temp_suff_stat);
                }
            }
        }

        loglike[loglike_start + s] = model->likelihood(temp_suff_stat, tree_pointer->suff_stat, true, false, state)
                    + model->likelihood(temp_suff_stat, tree_pointer->suff_stat, false, false, state);       
        if (loglike[loglike_start + s] > loglike_max)
        {
            loglike_max = loglike[loglike_start + s];
        }       
    }

}


void calculate_loglikelihood_categorical(std::vector<double> &loglike,
size_t &loglike_start, const std::vector<size_t> &subset_vars, size_t &N_Xorder,
matrix<size_t> &Xorder_std, matrix<size_t> &torder_std, double &loglike_max, std::vector<size_t> &X_counts,
std::vector<size_t> &X_num_unique, Model *model, std::unique_ptr<X_struct> &x_struct, size_t &total_categorical_split_candidates, std::unique_ptr<State> &state, tree *tree_pointer)
{

    // loglike_start is an index to offset
    // consider loglikelihood start from loglike_start

    size_t start;
    size_t end;
    size_t end2;
    double y_cumsum = 0.0;
    size_t n1;
    size_t n2;
    size_t temp;
    size_t N = N_Xorder;

    size_t effective_cutpoints = 0;

    std::vector<double> temp_suff_stat(tree_pointer->suff_stat.size());
    std::vector<double> total_suff_stat(model->dim_suffstat);
    std::vector<double> left_suff_stat(model->dim_suffstat);

    for (auto &&i : subset_vars)
    {
        if ((i >= state->p_continuous) && (X_num_unique[i - state->p_continuous] > 1))
        {
            // more than one unique values
            start = x_struct->variable_ind[i - state->p_continuous];
            end = x_struct->variable_ind[i + 1 - state->p_continuous] - 1; // minus one for indexing starting at 0
            end2 = end;

            while (X_counts[end2] == 0)
            {
                // move backward if the last unique value has zero counts
                end2 = end2 - 1;
            }
            // move backward again, do not consider the last unique value as cutpoint
            end2 = end2 - 1;

            y_cumsum = 0.0;
            //model -> suff_stat_fill(0.0); // initialize sufficient statistics
            std::fill(temp_suff_stat.begin(), temp_suff_stat.end(), 0.0);

            ////////////////////////////////////////////////////////////////
            //
            //  This part can be run in parallel, just like continuous case
            //
            //  If run in parallel, need to redefine model class for each thread
            //
            ////////////////////////////////////////////////////////////////

            n1 = 0;

            for (size_t j = start; j <= end2; j++)
            {

                if (X_counts[j] != 0)
                {

                    temp = n1 + X_counts[j] - 1;
                    calcSuffStat_categorical(temp_suff_stat, Xorder_std[i], torder_std, n1, temp, model, state);

                    n1 = n1 + X_counts[j];
                    // n1tau = (double)n1 * model->tau;
                    // n2tau = ntau - n1tau;

                    loglike[loglike_start + j] = model->likelihood(temp_suff_stat, tree_pointer->suff_stat, true, false, state) + model->likelihood(temp_suff_stat, tree_pointer->suff_stat, false, false, state);

                    // count total number of cutpoint candidates
                    effective_cutpoints++;

                    if (loglike[loglike_start + j] > loglike_max)
                    {
                        loglike_max = loglike[loglike_start + j];
                    }
                }
            }
        }
    }
}

void calculate_likelihood_no_split(std::vector<double> &loglike, size_t &N_Xorder, double &loglike_max, Model *model, std::unique_ptr<X_struct> &x_struct, size_t &total_categorical_split_candidates, std::unique_ptr<State> &state, tree *tree_pointer)
{

    loglike[loglike.size() - 1] = model->likelihood(tree_pointer->suff_stat, tree_pointer->suff_stat, false, true, state) + log(pow(1.0 + tree_pointer->getdepth(), model->beta) / model->alpha - 1.0) + log((double)loglike.size() - 1.0) + log(model->getNoSplitPenality());

    // this is important, update maximum of loglike vector
    if (loglike[loglike.size() - 1] > loglike_max)
    {
        loglike_max = loglike[loglike.size() - 1];
    }
}

// void predict_from_tree(tree &tree, const double *X_std, size_t N, size_t p, std::vector<double> &output, Model *model)
// {
//     tree::tree_p bn;
//     for (size_t i = 0; i < N; i++)
//     {
//         bn = tree.search_bottom_std(X_std, i, p, N);
//         output[i] = model->predictFromTheta(bn->theta_vector);
//     }
//     return;
// }

// void predict_from_datapointers(size_t tree_ind, Model *model, std::unique_ptr<State> &state, std::unique_ptr<X_struct> &x_struct)
// {
//     // // tree search, but read from the matrix of pointers to end node directly
//     // // easier to get fitted value of training set
//     // for (size_t i = 0; i < state->n_y; i++)
//     // {
//     //     state->predictions_std[tree_ind][i] = model->predictFromTheta(*(x_struct->data_pointers[tree_ind][i]));
//     // }
//     // return;
// }

void calcSuffStat_categorical(std::vector<double> &temp_suff_stat,
std::vector<size_t> &xorder, matrix<size_t> &torder, size_t &start,
size_t &end, Model *model, std::unique_ptr<State> &state)
{
    // calculate sufficient statistics for categorical variables

    // compute sum of y[Xorder[start:end, var]]
    for (size_t i = start; i <= end; i++)
    {
        for (auto j: torder[xorder[i]]){
            model->incSuffStat(state, xorder[i], j, temp_suff_stat);
        }
    }
}

void calcSuffStat_time(std::vector<double> &temp_suff_stat,
std::vector<size_t> &xorder, std::vector<size_t> &torder, size_t &start,
size_t &end, Model *model, std::unique_ptr<State> &state)
{
    // calculate sufficient statistics for time variable
    if (state->fl == 0){
        for (size_t i = start; i <= end; i++)
        {
            for (size_t j = 0; j < xorder.size(); j++){
                model->incSuffStat(state, xorder[j], torder[i],temp_suff_stat);
            }
        }
    } else {
        for (size_t i = start; i <= end; i++)
        {
            // TODO: split based on s values for each observation
            for (size_t j = 0; j < xorder.size(); j++){
                // check how many s for the i-th individual.
                model->incSuffStat(state, xorder[j], torder[i],temp_suff_stat);
            }
        }

    }
   
}

void calcSuffStat_continuous(std::vector<double> &temp_suff_stat,
std::vector<size_t> &xorder, matrix<size_t> &torder,
std::vector<size_t> &candidate_index,
size_t index, bool adaptive_cutpoint, Model *model,
std::unique_ptr<State> &state)
{
    // calculate sufficient statistics for continuous variables

    if (adaptive_cutpoint)
    {

        if (index == 0)
        {
            // initialize, only for the first cutpoint candidate, thus index == 0
            for (auto j: torder[xorder[0]]){
                model->incSuffStat(state, xorder[0], j, temp_suff_stat);

            }
        }

        // if use adaptive number of cutpoints, calculated based on vector candidate_index
        for (size_t q = candidate_index[index] + 1; q <= candidate_index[index + 1]; q++)
        {
            for (auto j: torder[xorder[q]]){
                model->incSuffStat(state, xorder[q], j, temp_suff_stat);
            }
        }
    }
    else
    {
        // use all data points as candidates
        for (auto j: torder[xorder[index]]){
            model->incSuffStat(state, xorder[index], j, temp_suff_stat);
        }
    }
}

void getThetaForObs_Outsample(std::vector<double> &output, tree &tree, size_t x_index, size_t t_index, const double *Xtest, const double *tpointer, size_t N_Xtest, size_t p)
{
    // get theta of ONE observation of ALL trees, out sample fit
    // input is a pointer to testing set matrix because it is out of sample
    // tree is a vector of all trees

    // output should have dimension (dim_theta, num_trees)

    tree::tree_p bn; // pointer to bottom node

    bn = tree.search_bottom_std(Xtest, tpointer, x_index, t_index, N_Xtest, p);
    output = bn->theta_vector;
}

void getThetaForObs_Outsample_ave(matrix<double> &output, std::vector<tree> &tree, size_t x_index, const double *Xtest, size_t N_Xtest, size_t p)
{
    // This function takes AVERAGE of ALL thetas on the PATH to leaf node

    // get theta of ONE observation of ALL trees, out sample fit
    // input is a pointer to testing set matrix because it is out of sample
    // tree is a vector of all trees

    // output should have dimension (dim_theta, num_trees)

    tree::tree_p bn; // pointer to bottom node
    size_t count = 1;

    for (size_t i = 0; i < tree.size(); i++)
    {

        // loop over trees
        // tree search
        bn = &tree[i]; // start from root node

        std::fill(output[i].begin(), output[i].end(), 0.0);
        count = 0;

        while (bn->getl())
        {
            // while bn has child (not bottom node)

            output[i] = output[i] + bn->theta_vector;
            count++;

            // move to the next level
            if (*(Xtest + N_Xtest * bn->getv() + x_index) <= bn->getc())
            {
                bn = bn->getl();
            }
            else
            {
                bn = bn->getr();
            }
        }

        // bn is the bottom node

        output[i] = output[i] + bn->theta_vector;
        count ++ ;

        // take average of the path
        for (size_t j = 0; j < output[i].size(); j++)
        {
            output[i][j] = output[i][j] / (double)count;
        }

    }

    return;
}

#ifndef NoRcpp
#endif
