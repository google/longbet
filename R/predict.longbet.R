#' Get post-burnin draws from longbet model
#'
#' @param model A trained longbet model.
#' @param x An input matrix for size n by p1. Column order matters: continuos features should all bgo before of categorical.
#' @param t time variable (post-treatment time for treatment term will be infered based on input t and z).
#' @param gp bool, predict time coefficient beta using gaussian process
#'
#' @return A list with two matrices. Each matrix corresponds to a set of draws of predicted values; rows are datapoints, columns are iterations.
#' @export
predict.longBet <- function(model, x, t, gp = FALSE, ...) {

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

    t_con <- as.matrix(t)
    t_mod <- as.matrix(sapply(t_con, function(x) max(x - model$t0, 0)))
    # print("Adjusted treatment time to predict:")
    # print(t_mod)
    if (any(t_mod != model$t_mod)){
        print('Adjusted treatment time does not trained treatment time:')
        print(model$t_mod)
        stop()
    }
    print("Predict con")
    obj_mu = .Call(`_longBet_predict`, x, t_con, model$model_list$tree_pnt_pr)
    print("Predict mod")
    obj_tau = .Call(`_longBet_predict`, x, t_mod, model$model_list$tree_pnt_trt)
    print("Predict beta")
    obj_beta = .Call(`_longBet_predict_beta`, t_mod, 
        model$t_mod, model$gp_info$resid, model$gp_info$A_diag, model$gp_info$Sig_diag,
        model$model_params$sig_knl, model$model_params$lambda_knl)

    num_sweeps <- ncol(model$tauhats)
    num_burnin <- model$model_params$burnin

    if(num_burnin >= num_sweeps) {
        stop(paste0('burnin (',num_burnin,') cannot exceed or match the total number of sweeps (',num_sweeps,')'))
    }

    n <- nrow(x)
    p <- length(t)

    obj_mu$preds <- obj_mu$preds * model$sdy
    obj_tau$preds <- obj_tau$preds * model$sdy


    obj <- list()
    obj$muhats <- obj_mu$preds
    obj$tauhats <- obj_tau$preds

    obj$tauhats.adjusted <- array(NA, dim = c(n, p, num_sweeps - num_burnin))
    obj$muhats.adjusted <- array(NA, dim = c(n, p, num_sweeps - num_burnin))
    seq <- (num_burnin+1):num_sweeps
    for (i in seq) {
        obj$muhats.adjusted[,, i - num_burnin] = matrix(obj_mu$preds[,i], n, p) * (model$a_draws[i]) + model$meany
        obj$tauhats.adjusted[,, i - num_burnin] = matrix(obj_tau$preds[,i], n, p) * (model$b_draws[i,2] - model$b_draws[i,1]) * t(matrix(rep(model$beta_draws[,i - num_burnin], n), p, n))
        # TODO: check betadraws t_mod matches t
    }
    
    obj$beta_draws = model$beta_draws
    return(obj)
}
