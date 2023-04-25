#' longBet model
#'
#' @param y An n by t matrix of outcome variables.
#' @param x n by p input matrix of covariates. (If the covariates matrix is different for the prognostic and treatment term, please use longBet_full).
#' @param z An n by t matrix of treatment assignments.
#' @param t time variable (post-treatment time for treatment term will be infered based on input t and z).
#' @param pcat The number of categorical inputs in matrix x.
#' @param num_sweeps The total number of sweeps over data (default is 60).
#' @param num_burnin The number of burn-in sweeps (default is 20).
#' @param num_trees_pr The number of trees in the prognostic forest (default is 50).
#' @param num_trees_trt The number of trees in the treatment forest (default is 20).
#' @param mtry number of variables to be sampled as split candidate per tree.
#' @param n_min The minimum node size. (default is 1)
#' @param sig_knl variance parameter for squared exponential kernel (default is 1).
#' @param lambda_knl lengthscale parameter for squared exponential kernel (default is 5).
#'
#' @return A fit file, which contains the draws from the model as well as parameter draws at each sweep.
#' @export
longbet <- function(y, x, z, t, pcat, 
                    num_sweeps = 60, num_burnin = 20,
                    num_trees_pr = 50, num_trees_trt = 20,
                    mtry = 0L, n_min = 10,
                    sig_knl = 1, lambda_knl = 5,
                    a_scaling = TRUE, b_scaling = FALSE) {

    if(!("matrix" %in% class(x))){
        cat("Msg: input x is not a matrix, try to convert type.\n")
        x = as.matrix(x)
    }
    if(!("matrix" %in% class(z))){
        cat("Msg: input z is not a matrix, try to convert type.\n")
        z = as.matrix(z)
    }
    if(!("matrix" %in% class(y))){
        cat("Msg: input y is not a matrix, try to convert type.\n")
        y = as.matrix(y)
    }

    if(any(dim(z) != dim(y))) {
        stop("Dimensions of response y and treatment z do not match. \n")
    }

    if (length(t) != ncol(y)){
        stop("Lenght of input t should match the columns of y. \n")
    }

    # check if treatment all start at the same time
    # number of treated periods per unit should only be 0 or t1 - t0 + 1
    unique_z_sum <- unique(rowSums(z))
    if (length(unique_z_sum) != 2) {
        stop("Current version can only handle treamtments occured at the same time. \n")
    }

    # get post-treatment time variable
    if (is.null(t)){
        t_con = 1:ncol(y)
    } else {
        t_con = t
    }

    if (ncol(y) > 1) {
        post_t <- sort(unique_z_sum)[2]
        t0 <- t_con[ncol(y) - post_t]
        t_mod <- sapply(t_con, function(x) max(x - t0, 0))
        print("Adjusted treatment time:")
        print(t_mod)
        # t_mod <- c(rep(0, ncol(y) - post_t), 1:post_t)
    } else {
        t_mod <- c(1)
        t0 <- NULL
    }

    if(!("matrix" %in% class(t_con))){
        t_con = as.matrix(t_con)
    }
    if(!("matrix" %in% class(t_mod))){
        t_mod = as.matrix(t_mod)
    }


    if (nrow(x) != nrow(y)) {
        stop(paste0('row number mismatch between X (', nrow(x), ') and y (', nrow(y), ')'))
    }

    # check if p_categorical was not provided
    if(is.null(pcat)) {
        stop('number of categorical variables pcat_con is not specified')
    }

    # check if p_categorical exceeds the number of columns
    if(pcat > ncol(x)) {
        stop('number of categorical variables (pcat_con) cannot exceed number of columns')
    }

    # check if p_categorical is negative
    if(pcat < 0) {
        stop('number of categorical values can not be negative: check pcat_con and pcat_mod')
    }

    # check if mtry exceeds the number of columns
    if(mtry > ncol(x)) {
        cat('Msg: mtry value cannot exceed number of columns; set to default.\n')
        mtry <- 0
    }
    # check if mtry is negative
    if(mtry < 0) {
        cat('Msg: mtry value cannot exceed number of columns; set to default.\n')
        mtry <- 0
    }

    # set defaults for taus if it wasn't provided with the call
    tau_con = 0.6 * var(as.vector(y)) / num_trees_pr
    tau_mod = 0.1 * var(as.vector(y)) / num_trees_trt
    
    meany = mean(y) # disable meany temporarily
    y = y - meany
    sdy = sd(y)

    if(sdy == 0) {
        stop('y is a constant variable; sdy = 0')
    } else {
        y = y / sdy
    }

    if(num_burnin >= num_sweeps){
        stop(paste0('num_burnin (',num_burnin,') cannot exceed or match the total number of sweeps (',num_sweeps,')'))
    }

    # deprecated hyperparameters for user experience
    max_depth = 50
    num_cutpoints = 20
    no_split_penality = log(num_cutpoints)
    alpha_con = 0.95; beta_con = 1.25
    kap_con = 16; s_con = 4
    pr_scale = FALSE
    alpha_mod = 0.95; beta_mod = 1.25
    kap_mod = 16; s_mod = 4
    trt_scale = FALSE
    verbose = FALSE; parallel = TRUE
    set_random_seed = FALSE; random_seed = 0
    sample_weights_flag = TRUE
    split_t_mod = TRUE; split_t_con = TRUE
    print( "call longbet function" )
    
    obj = longBet_cpp(y = y,
                    X = x, 
                    X_tau = x, 
                    z = z, 
                    t_con = t_con, 
                    t_mod = t_mod,
                    num_sweeps = num_sweeps, 
                    burnin = num_burnin,
                    max_depth = max_depth, 
                    n_min = n_min,
                    num_cutpoints = num_cutpoints,
                    no_split_penality = no_split_penality,
                    mtry_pr = mtry, 
                    mtry_trt = mtry,
                    p_categorical_pr = pcat, 
                    p_categorical_trt = pcat,
                    num_trees_pr = num_trees_pr,
                    alpha_pr = alpha_con, 
                    beta_pr = beta_con, 
                    tau_pr = tau_con,
                    kap_pr = kap_con, 
                    s_pr = s_con,
                    pr_scale = pr_scale,
                    num_trees_trt = num_trees_trt,
                    alpha_trt = alpha_mod, 
                    beta_trt = beta_mod,
                    tau_trt = tau_mod,
                    kap_trt = kap_mod, 
                    s_trt = s_mod,
                    trt_scale = trt_scale,
                    verbose = verbose, 
                    parallel = parallel, 
                    set_random_seed = set_random_seed,
                    random_seed = random_seed, 
                    sample_weights_flag = sample_weights_flag,
                    a_scaling = a_scaling, 
                    b_scaling = b_scaling,
                    split_t_mod = split_t_mod, 
                    split_t_con = split_t_con,
                    sig_knl = sig_knl, 
                    lambda_knl = lambda_knl)
    class(obj) = "longBet"

    obj$t0 = t0
    obj$sdy = sdy
    obj$meany = meany
    obj$tauhats = obj$tauhats * sdy
    obj$muhats = obj$muhats * sdy

    obj$tauhats.adjusted <- array(NA, dim = c(nrow(y), ncol(y), num_sweeps - num_burnin))
    obj$muhats.adjusted <- array(NA, dim = c(nrow(y), ncol(y), num_sweeps - num_burnin))
    seq <- (num_burnin+1):num_sweeps
    for (i in seq) {
        obj$tauhats.adjusted[,, i - num_burnin] = matrix(obj$tauhats[,i], nrow(y), ncol(y)) * (obj$b_draws[i,2] - obj$b_draws[i,1])
        obj$tauhats.adjusted[,,i - num_burnin] = obj$tauhats.adjusted[,,i - num_burnin] * t(matrix(rep(obj$beta_draws[,i], nrow(y)), ncol(y), nrow(y)))
        obj$muhats.adjusted[,, i - num_burnin] = matrix(obj$muhats[,i], nrow(y), ncol(y)) * (obj$a_draws[i]) + meany
    }
    
    # obj$beta_draws = obj$beta_draws[, (num_burnin+1):num_sweeps]
    return(obj)
}
