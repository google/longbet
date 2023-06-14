#' Get post-burnin draws from longbet model
#'
#' @param model A trained longbet model.
#' @param x An input matrix for size n by p1. Column order matters: continuos features should all bgo before of categorical.
#' @param z n by p_y treatment matrix indicating whether each unit get treated at each step, should match the training period
#' @param gp bool, predict time coefficient beta using gaussian process
#'
#' @return A matrix for predicted prognostic effect and a matrix for predicted treatment effect. 
#' @export
predict.longBet <- function(model, x, z, sigma = NULL, lambda = NULL, ...) {

    print(dim(x))
    if(!("matrix" %in% class(x))) {
        cat("Msg: input x is not a matrix, try to convert type.\n")
        x = as.matrix(x)
    }

    if(ncol(x) != model$input_var_count$x_con) {
        stop(paste0('Check dimensions of input matrices. The model was trained on
        x with ', model$input_var_count$x_con,
        ' columns; trying to predict on x with ', ncol(x),' columns.'))
    }

    if(!("matrix" %in% class(z))) {
        cat("Msg: input z is not a matrix, try to convert type.\n")
        z = as.matrix(z)
    }

    t_con <- as.matrix(1:ncol(z))
    # t_mod <- as.matrix(rep(0, ncol(z)))
    t_mod <- as.matrix(sapply(t_con, function(x) max(x - model$t0, 0)))
    # print("Adjusted treatment time to predict:")
    # print(t_mod)
    
    obj_mu = .Call(`_longBet_predict`, x, t_con, model$model_list$tree_pnt_pr)

    obj_tau = .Call(`_longBet_predict`, x, t_mod, model$model_list$tree_pnt_trt)

    # Match post treatment periods
    n <- nrow(z)
    p <- ncol(z)

    num_sweeps <- ncol(model$tauhats)
    num_burnin <- model$model_params$burnin

    if(num_burnin >= num_sweeps) {
        stop(paste0('burnin (',num_burnin,') cannot exceed or match the total number of sweeps (',num_sweeps,')'))
    }


    post_trt <- t(apply(z, 1, cumsum)) + 1
    beta_preds <- array(NA, dim = c(n, p, num_sweeps))

    max_post_trt <- max(post_trt)
    if (max_post_trt > nrow(longbet.fit$beta_draws)){
        # predict beta
        stop("TODO: update extrapolation code for staggered adoption, \n")
        if (is.null(sigma)) { 
            if (nrow(model$beta_draws) > t0 + 1){
                # if the training has more than 1 treatment period
                sigma = sigma_knl = mean( sqrt( apply(longbet.fit$beta_draws[model$t0:t1,], 2, var) ))
            } else {
                sigma = 1
            }
        }
        if (is.null(lambda)) { lambda = (nrow(model$beta_draws) - t0) / 2}
        print(paste("predict beta with GP, sigma = ", sigma, ", lambda = ", lambda, sep = ""))

        # obj_beta = .Call(`_longBet_predict_beta`, as.matrix(nrow(model$beta_draws)), 
        #     model$gp_info$t_values, model$gp_info$resid, model$gp_info$A_diag, model$gp_info$Sig_diag,
        #     sigma, lambda)
        # beta[is.na(idx), ] <- obj_beta$beta
    }
    for (i in 1:num_sweeps){
        beta_preds[,,i] <- t(apply(post_trt, 1, function(x, beta) beta[x], beta = model$beta_draws[,i]))
    }

    obj_mu$preds <- obj_mu$preds * model$sdy
    obj_tau$preds <- obj_tau$preds * model$sdy


    obj <- list()
    class(obj) = "longBet.pred"
    
    obj$muhats <- array(NA, dim = c(n, p, num_sweeps - num_burnin))
    obj$tauhats <- array(NA, dim = c(n, p, num_sweeps - num_burnin))
    seq <- (num_burnin+1):num_sweeps
    for (i in seq) {
        obj$muhats[,, i - num_burnin] = matrix(obj_mu$preds[,i], n, p) * (model$a_draws[i]) + model$meany +  matrix(obj_tau$preds[,i], n, p) *  model$b_draws[i,1] * model$beta_draws[1, i]
        # obj$tauhats[,, i - num_burnin] = matrix(obj_tau$preds[,i], n, p) * (model$b_draws[i,2] - model$b_draws[i,1]) * beta_preds[,,i]
        obj$tauhats[,, i - num_burnin] = matrix(obj_tau$preds[,i], n, p) * (model$b_draws[i,2] * beta_preds[,,i] - model$b_draws[i,1] * model$beta_draws[1, i]) # * beta_preds[,,i]
        # TODO: change tauhat to b1 * beta_s * tau_s - b0 * beta_0 * tau_0 when tau can split on post-treatment time
    }
    obj$beta_preds <- beta_preds
    # obj$beta_draws = beta
    return(obj)
}

get_ate <- function(object, alpha = 0.05, ...) {
    if(class(object) != "longBet.pred"){
        stop("Input object should be output from predict.longBet function")    
    }
    ate_full <- apply(object$tauhats, c(2, 3), mean)
    obj <- list()
    obj$ate <- rowMeans(ate_full)
    obj$interval <- apply(ate_full, 1, quantile, probs = c(alpha / 2, 1- alpha / 2))
    obj$ate_full <- ate_full
    return(obj)
}

get_att <- function(object, z, alpha = 0.05, ...){
    if(class(object) != "longBet.pred"){
        stop("Input object should be output from predict.longBet function")    
    }
    att_full <- apply(object$tauhats[z,,], c(2, 3), mean)
    obj <- list()
    obj$att <- rowMeans(att_full)
    obj$interval <- apply(att_full, 1, quantile, probs = c(alpha / 2, 1- alpha / 2))
    obj$att_full <- att_full
    return(obj)
}

get_cate <- function(object, alpha = 0.05, ...){
    if(class(object) != "longBet.pred"){
        stop("Input object should be output from predict.longBet function")    
    }
    obj <- list()
    obj$cate <- apply(object$tauhats, c(1, 2), mean)
    obj$interval <- apply(object$tauhats, c(1, 2), quantile, probs = c(alpha / 2, 1 - alpha / 2))
    return(obj)
}
