#' Get post-burnin draws from longbet model
#'
#' @param model A trained longbet model.
#' @param x An input matrix for size n by p1. Column order matters: continuos features should all bgo before of categorical.
#' @param z n by p_y treatment matrix indicating whether each unit get treated at each step, should match the training period
#' @param gp bool, predict time coefficient beta using gaussian process
#'
#' @return A matrix for predicted prognostic effect and a matrix for predicted treatment effect. 
#' @export
predict.longbet <- function(model, x, x_trt, z, t = NULL, sigma = NULL, lambda = NULL, ...) {

    if(!("matrix" %in% class(x))) {
        cat("Msg: input x is not a matrix, try to convert type.\n")
        x = as.matrix(x)
    }

    if(ncol(x) != model$input_var_count$x_con) {
        stop(paste0('Check dimensions of input matrices. The model was trained on
        x with ', model$input_var_count$x_con,
        ' columns; trying to predict on x with ', ncol(x),' columns.'))
    }

    if(!("matrix" %in% class(x_trt))) {
        cat("Msg: input x is not a matrix, try to convert type.\n")
        x_trt = as.matrix(x_trt)
    }

    if(ncol(x_trt) != model$input_var_count$x_mod) {
        stop(paste0('Check dimensions of input matrices. The model was trained on
        x with ', model$input_var_count$x_mod,
        ' columns; trying to predict on x with ', ncol(x_trt),' columns.'))
    }

    if(!("matrix" %in% class(z))) {
        cat("Msg: input z is not a matrix, try to convert type.\n")
        z = as.matrix(z)
    }

    if (nrow(z) != nrow(x)){
        stop("X and Z should have the same number of rows. \n")
    }

    
    if (is.null(t)){
        print(paste(c("Predicting from time", longbet.fit$time), collapse = " "))
        t_con <-  matrix(rep(model$time, nrow(x)), nrow = nrow(x), byrow = T)
    } else {
        if (length(t) != ncol(z)){
            stop("Msg: lenght of t should match the size of z. \n")
        }
        print(paste(c("Predicting from time", t), collapse = " "))
        t_con <-  matrix(rep(t, nrow(x)), nrow = nrow(x), byrow = T)
    }

    t_mod <- t( apply(z, 1, cumsum) )
    
    obj_mu = .Call(`_longbet_predict_longbet`, x, t_con, model$model_list$tree_pnt_pr)

    obj_tau = .Call(`_longbet_predict_longbet`, x_trt, t_mod, model$model_list$tree_pnt_trt)
    
    obj_tau0 = .Call(`_longbet_predict_longbet`, x_trt, matrix(rep(0, nrow(x)), ncol = 1), model$model_list$tree_pnt_trt)

    # Match post treatment periods
    n <- nrow(z)
    p <- ncol(z)

    num_sweeps <- ncol(model$tauhats)
    num_burnin <- model$model_params$burnin

    if(num_burnin >= num_sweeps) {
        stop(paste0('burnin (',num_burnin,') cannot exceed or match the total number of sweeps (',num_sweeps,')'))
    }


    post_trt <- t(apply(z, 1, cumsum))
    beta_preds <- array(NA, dim = c(n, p, num_sweeps))

    max_post_trt <- max(post_trt)
    S <- nrow(model$beta_values) - 1 # max S observed
    if (max_post_trt > S){
        # predict beta
        # stop("TODO: update extrapolation code for staggered adoption, \n")
        if (is.null(sigma)) {  sigma = 1 }
        if (is.null(lambda)) { lambda = nrow(model$beta_values) / 2}
        print(paste("predict beta with GP, sigma = ", sigma, ", lambda = ", lambda, sep = ""))

        # beta to be predicted?
        beta_test <- as.matrix((S + 1) : max_post_trt)
        obj_beta = .Call(`_longbet_predict_beta`, beta_test,
            as.matrix(model$gp_info$t_values), model$gp_info$resid, model$gp_info$A_diag, model$gp_info$Sig_diag,
            sigma, lambda)
        model$beta_values <- rbind(model$beta_values, obj_beta$beta)
    }
    for (i in 1:num_sweeps){
        beta_preds[,,i] <- t(apply(post_trt, 1, function(x, beta) beta[x + 1], beta = model$beta_values[,i]))
    }

    obj_mu$preds <- obj_mu$preds * model$sdy
    obj_tau$preds <- obj_tau$preds * model$sdy
    obj_tau0$preds <- obj_tau0$preds * model$sdy

    obj <- list()
    class(obj) = "longbet.pred"
    
    obj$muhats <- array(NA, dim = c(n, p, num_sweeps - num_burnin))
    obj$tauhats <- array(NA, dim = c(n, p, num_sweeps - num_burnin))
    seq <- (num_burnin+1):num_sweeps
    for (i in seq) {
        obj$muhats[,, i - num_burnin] = matrix(obj_mu$preds[,i], n, p) * (model$a_draws[i]) + model$meany +  matrix(obj_tau$preds[,i], n, p) *  model$b_draws[i,1] * model$beta_draws[1, i]
        # obj$tauhats[,, i - num_burnin] = matrix(obj_tau$preds[,i], n, p) * (model$b_draws[i,2] * beta_preds[,,i] - model$b_draws[i,1] * model$beta_draws[1, i]) # * beta_preds[,,i]
        obj$tauhats[,, i - num_burnin] = model$b_draws[i,2] * beta_preds[,,i] * matrix(obj_tau$preds[,i], n, p)  - matrix(rep(model$b_draws[i,1] * model$beta_values[1, i] * obj_tau0$preds[,i], p), ncol = p)
        # TODO: change tauhat to b1 * beta_s * tau_s - b0 * beta_0 * tau_0 when tau can split on post-treatment time
    }
    obj$beta_values <- model$beta_values
    obj$beta_preds <- beta_preds
    obj$z <- z
    return(obj)
}

get_att <- function(object, alpha = 0.05, ...){
    if(class(object) != "longbet.pred"){
        stop("Input object should be output from predict.longbet function")    
    }
    # att_full <- apply(object$tauhats[z,,], c(2, 3), mean)

    n <- dim(object$tauhats)[1]
    treatment_period <- nrow(object$beta_values) - 1
    num_sweeps <- dim(object$tauhats)[3]

    # Align treatment effect 
    align_catt <- array(NA, dim = c(n, treatment_period, num_sweeps))

    for (i in 1:n){
        if (sum(object$z[i,]) == 0) {next}
        post_t <- 1:sum(object$z[i,])
        align_catt[i,post_t,] = object$tauhats[i, object$z[i,]== 1,]
    }

    att_hat <- apply(align_catt, c(2, 3), mean, na.rm = T)

    obj <- list()
    obj$att <- rowMeans(att_hat)
    obj$intervals <- apply(att_hat, 1, quantile, probs = c(alpha / 2, 1- alpha / 2))
    obj$att_full <- att_hat
    return(obj)
}

get_catt <- function(object, alpha = 0.05, ...){
    if(class(object) != "longbet.pred"){
        stop("Input object should be output from predict.longbet function")    
    }
    obj <- list()
    obj$catt <- apply(object$tauhats, c(1, 2), mean)
    obj$intervals <- apply(object$tauhats, c(1, 2), quantile, probs = c(alpha / 2, 1 - alpha / 2))
    return(obj)
}
