#' Full longBet model with all parameters available for tuning
#'
#' @param y An n by t matrix of outcome variables
#' @param z An n by t matrix of treatment assignments
#' @param x_con An input matrix for the prognostic term of size n by p1. Column order matters: continuos features should all bgo before of categorical.
#' @param x_mod An input matrix for the treatment term of size n by p2 (default x_mod = x_con). Column order matters: continuos features should all go beforeof categorical.
#' @param t_con A t by 1 vector of time variables (potentially can be t by x matrix)
#' @param t_mod A t by 1 vector of post-treatment time (TODO: change to n by t to accomodate time-varying treatments)
#' @param pihat An array of propensity score estimates (default is NULL). (Currently deprecated)
#' @param pcat_con The number of categorical inputs in the prognostic term input matrix x_con.
#' @param pcat_mod The number of categorical inputs in the treatment term input matrix x_mod.
#' @param num_sweeps The total number of sweeps over data (default is 60).
#' @param num_burnin The number of burn-in sweeps (default is 20).
#' @param n_trees_con The number of trees in the prognostic forest (default is 30).
#' @param n_trees_mod The number of trees in the treatment forest (default is 10).
#' @param max_depth The maximum possible depth of a tree. (default is 50)
#' @param Nmin The minimum node size. (default is 1)
#' @param num_cutpoints The number of adaptive cutpoints considered at each split for continuous variables (default is 20).
#' @param no_split_penality As the name states..
#' @param sample_weights_flag bool variable, enable weight sampling for split candiate sampling.
#' @param mtry_con number of variables to be sampled as split candidate in prognostic term
#' @param mtry_mod number of variables to be sampled as split candidate in treatment term
#' @param alpha_con Base parameter for tree prior on trees in prognostic forest (default is 0.95).
#' @param beta_con Power parameter for tree prior on trees in prognostic forest (default is 1.25).
#' @param tau_con Prior variance over the mean on on trees in prognostic forest (default is 0.6*var(y)/n_trees_con).
#' @param kap_con
#' @param s_con
#' @param pr_scale
#' @param alpha_mod Base parameter for tree prior on trees in treatment forest (default is 0.25).
#' @param beta_mod Power parameter for tree prior on trees in treatment forest (default is 3).
#' @param tau_mod Prior variance over the mean on on trees in treatment forest (default is 0.1*var(y)/n_trees_mod).
#' @param kap_mod
#' @param s_mod
#' @param trt_scale
#' @param verbose bool variable to print 
#' @param parallel bool variable for parallel (not sure if it works)
#' @param random_seed random seed for sampling
#' @param a_scaling
#' @param b_scaling TODO: does not mix will with GP. Set to FALSE by default
#' @param sig_knl variance parameter for squared exponential kernel
#' @param lambda_knl lengthscale parameter for squared exponential kernel
#'
#' @return A fit file, which contains the draws from the model as well as parameter draws at each sweep.
#' @export

longBet_full <- function(y, z, x_con, x_mod = x_con, t_con = NULL, t_mod = NULL, 
                pihat = NULL,
                pcat_con = NULL,
                pcat_mod = pcat_con,
                num_sweeps = 60, burnin = 20,
                n_trees_con = 30L, 
                n_trees_mod = 10L,
                max_depth = 50, Nmin = 1L,
                num_cutpoints = 20,
                no_split_penality = "Auto", 
                sample_weights_flag = TRUE,
                mtry_con = 0L, mtry_mod = 0L,
                alpha_con = 0.95, beta_con = 1.25, tau_con = NULL,
                kap_con = 16, s_con = 4,
                pr_scale = FALSE,
                alpha_mod = 0.25, beta_mod = 3, tau_mod = NULL,
                kap_mod = 16, s_mod = 4,
                trt_scale = FALSE,
                verbose = FALSE, parallel = TRUE,
                random_seed = NULL, 
                a_scaling = TRUE, b_scaling = FALSE,
                sig_knl = 1, lambda_knl = 1) {

    if(!("matrix" %in% class(x_con))){
        cat("Msg: input x_con is not a matrix, try to convert type.\n")
        x_con = as.matrix(x_con)
    }
    if(!("matrix" %in% class(x_mod))){
        cat("Msg: input x_mod is not a matrix, try to convert type.\n")
        x_mod = as.matrix(x_mod)
    }
    if(!("matrix" %in% class(z))){
        cat("Msg: input z is not a matrix, try to convert type.\n")
        z = as.matrix(z)
    }
    if(!("matrix" %in% class(y))){
        cat("Msg: input y is not a matrix, try to convert type.\n")
        y = as.matrix(y)
    }

    if (is.null(t_con)){
        t_con = rep(1, ncol(y))
    }
    if(!("matrix" %in% class(t_con))){
        t_con = as.matrix(t_con)
    }
    if (is.null(t_mod)){
        if (ncol(z) > 1){
            # create post treatment time
            t_mod <- matrix(0, nrow(z), ncol(z))
            for (i in 1:nrow(z)){
                for (j in 1:ncol(z)){
                    if (z[i, j] > 0){
                        t_mod[i, j:ncol(z)] = 1:(ncol(z) - j + 1)
                        break;
                    }
                }
            }
        } else {
            t_mod <- rep(1, ncol(z))
        }
       
    }
    if(!("matrix" %in% class(t_mod))){
        t_mod = as.matrix(t_mod)
    }

    if (!is.null(pihat)) {
      x_con <- cbind(pihat, x_con)
    }
    # # compute pihat if it wasn't provided with the call
    # if(is.null(pihat)) {
    #     sink("/dev/null") # silence output
    #     fitz = nnet::nnet(z~.,data = x_con, size = 3,rang = 0.1, maxit = 1000, 
    #     abstol = 1.0e-8, decay = 5e-2)
    #     sink() # close the stream
    #     pihat = fitz$fitted.values
    # }
    # if(!("matrix" %in% class(pihat))){
    #     cat("Msg: input pihat is not a matrix, try to convert type.\n")
    #     pihat = as.matrix(pihat)
    # }
    # x_con <- cbind(pihat, x_con)
    
    p_X <- ncol(x_con)
    p_Xt <- ncol(x_mod)


    if(nrow(x_con) != nrow(x_mod)) {
        stop('row number mismatch for the two input matrices')
    }
    if (nrow(x_con) != nrow(y)) {
        stop(paste0('row number mismatch between X (', nrow(x_con), ') and y (', nrow(y), ')'))
    }
    if (nrow(x_con) != nrow(z)) {
        stop(paste0('row number mismatch between X (', nrow(x_con), ') and z (', nrow(z), ')'))
    }

    # check if p_categorical was not provided
    if(is.null(pcat_con)) {
        stop('number of categorical variables pcat_con is not specified')
    }
    if(is.null(pcat_mod)) {
        stop('number of categorical variables pcat_mod is not specified')
    }

    # check if p_categorical exceeds the number of columns
    if(pcat_con > p_X) {
        stop('number of categorical variables (pcat_con) cannot exceed number of columns')
    }
    if(pcat_mod > p_Xt) {
        stop('number of categorical variables (pcat_mod) cannot exceed number of columns')
    }

    # check if p_categorical is negative
    if(pcat_con < 0 || pcat_mod < 0) {
        stop('number of categorical values can not be negative: check pcat_con and pcat_mod')
    }

    # check if mtry exceeds the number of columns
    if(mtry_con > p_X) {
        cat('Msg: mtry value cannot exceed number of columns; set to default.\n')
        mtry_con <- 0
    }
    if(mtry_mod > p_Xt) {
        cat('Msg: mtry value cannot exceed number of columns; set to default.\n')
        mtry_mod <- 0
    }

    # check if mtry is negative
    if(mtry_con < 0) {
        cat('Msg: mtry value cannot exceed number of columns; set to default.\n')
        mtry_con <- 0
    }
    if(mtry_mod < 0) {
        cat('Msg: mtry value cannot exceed number of columns; set to default.\n')
        mtry_mod <- 0
    }

    # set defaults for taus if it wasn't provided with the call
    if(is.null(tau_con)) {
        tau_con = 0.6 * var(as.vector(y)) / n_trees_con
    }
    # set defaults for taus if it wasn't provided with the call
    if(is.null(tau_mod)) {
        tau_mod = 0.1 * var(as.vector(y)) / n_trees_mod
    }

    # meany = mean(y) # disable meany temporarily
    meany = 0 
    y = y - meany
    sdy = sd(y)

    if(sdy == 0) {
        stop('y is a constant variable; sdy = 0')
    } else {
        y = y / sdy
    }

    # compute default values for taus if none provided
    if(is.null(tau_con)) {
        tau_con <- 0.6*var(y)/n_trees_con
    }

    if(is.null(tau_mod)) {
        tau_mod <- 0.1*var(y)/n_trees_mod
    }

    if(is.null(random_seed)){
        set_random_seed = FALSE
        random_seed = 0;
    }else{
        cat("Set random seed as ", random_seed, "\n")
        set_random_seed = TRUE
    }

    if(burnin >= num_sweeps){
        stop(paste0('burnin (',burnin,') cannot exceed or match the total number of sweeps (',sweeps,')'))
    }
    if(no_split_penality == "Auto"){
        no_split_penality = log(num_cutpoints)
    }
    
    obj = longBet_cpp(y, x_con, x_mod, z, t_con, t_mod,
                         num_sweeps, burnin,
                         max_depth, Nmin,
                         num_cutpoints,
                         no_split_penality, mtry_con, mtry_mod,
                         pcat_con,
                         pcat_mod,
                         n_trees_con,
                         alpha_con, beta_con, tau_con,
                         kap_con, s_con,
                         pr_scale,
                         n_trees_mod,
                         alpha_mod, beta_mod, tau_mod,
                         kap_mod, s_mod,
                         trt_scale,
                         verbose, parallel, set_random_seed,
                         random_seed, sample_weights_flag,
                         a_scaling, b_scaling,
                         sig_knl, lambda_knl)
    class(obj) = "longBet"

    #obj$sdy_use = sdy_use
    obj$sdy = sdy
    obj$meany = meany
    obj$tauhats = obj$tauhats * sdy
    obj$muhats = obj$muhats * sdy

    obj$tauhats.adjusted <- array(NA, dim = c(nrow(y), ncol(y), num_sweeps-burnin))
    obj$muhats.adjusted <- array(NA, dim = c(nrow(y), ncol(y), num_sweeps-burnin))
    seq <- (burnin+1):num_sweeps
    for (i in seq) {
        obj$tauhats.adjusted[,, i - burnin] = matrix(obj$tauhats[,i], nrow(y), ncol(y)) * (obj$b_draws[i,2] - obj$b_draws[i,1])
        obj$tauhats.adjusted[,,i - burnin] = obj$tauhats.adjusted[,,i - burnin] * t(matrix(rep(obj$beta_draws[,i], nrow(y)), ncol(y), nrow(y)))
        obj$muhats.adjusted[,, i - burnin] = matrix(obj$muhats[,i], nrow(y), ncol(y)) * (obj$a_draws[i]) + meany
    }
    
    obj$beta_draws = obj$beta_draws[, (burnin+1):num_sweeps]
    return(obj)
}
