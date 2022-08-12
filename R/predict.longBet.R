# #' Get post-burnin draws from trained model
# #'
# #' @param model A trained LongBet model.
# #' @param x An input matrix of size n by p. Column order matters: continuos features should all bgo before of categorical.
# #' @param t A vector of time to be predicted
# #' @param sig_knl variance parameter for squared exponential kernel (default is based on trained model).
# #' @param lambda_knl lengthscale parameter for squared exponential kernel (default is based on trained model).
# #'
# #'
# #' @return A 3d array matrix of predicted mu with dimensions correspond to datapoints, time and iterations, respectively.
# #' @return A 3d array matrix of predicted tau with dimensions correspond to datapoints, time and iterations, respectively.
# #' @export
# predict.longBet <- function(model, x, t, sig_knl = NULL, lambda_knl = NULL) {
#     if(!("matrix" %in% class(x))) {
#         cat("Msg: input x is not a matrix, try to convert type.\n")
#         x = as.matrix(x)
#     }
#     if(!("matrix" %in% class(t))) {
#         cat("Msg: input t is not a matrix, try to convert type.\n")
#         t = as.matrix(t)
#     }
#     # get post treatment time
#     t_mod <- as.matrix(sapply(t, function(x) max(x - model$t0, 0)))

#     if(ncol(x) != model$input_var_count$x_mod) {
#         stop(paste0('Check dimensions of input matrices. The model was trained on
#         x with ', model$input_var_count$x_mod,
#         ' columns; trying to predict on x_con with ', ncol(x),' columns.'))
#     }

#     if (is.null(sig_knl)){
#       sig_knl = model$model_params$sig_knl
#     }
#     if (is.null(lambda_knl)){
#       lambda_knl = model$model_params$lambda_knl
#     }

#     print("model time residuals = ")
#     print(model$time_info$time_residuals)


#     print("x = ")
#     print(x)

#     obj1 = .Call(`_longBet_predict`, x, t_mod, model$model_list$tree_pnt_pr)
#     obj2 = .Call(`_longBet_predict`, x, t_mod, model$model_list$tree_pnt_trt)
#     obj_beta = .Call(`_longBet_predict_beta`, as.matrix(t_mod), 
#       as.matrix(model$time_info$train_t),
#       as.matrix(model$time_info$time_beta),
#       as.matrix(model$time_info$time_residuals),
#       as.matrix(model$time_info$time_diag_A),
#       as.matrix(model$time_info$time_diag_Sig),
#       sig_knl, lambda_knl)
    
#     sweeps <- ncol(model$tauhats)
#     burnin <- model$model_params$burnin

#     tauhats.adjusted <- array(NA, dim = c(nrow(y), ncol(y), sweeps - burnin))
#     muhats.adjusted <- array(NA, dim = c(nrow(y), ncol(y), sweeps - burnin))
#     seq <- (burnin+1):sweeps
#     for (i in seq) {
#         tauhats.adjusted[,, i - burnin] = model$sdy * matrix(obj2$predicted_values[,i], nrow(y), ncol(y)) * (model$b_draws[i,2] - model$b_draws[i,1])
#         tauhats.adjusted[,,i - burnin] = tauhats.adjusted[,,i - burnin] * t(matrix(rep(obj_beta$beta_draws[,i - burnin], nrow(y)), ncol(y), nrow(y)))
#         muhats.adjusted[,, i - burnin] = model$sdy * matrix(obj1$predicted_values[,i], nrow(y), ncol(y)) * (model$a_draws[i]) + model$meany
#     }

#     obj <- list(mudraws=muhats.adjusted, taudraws=tauhats.adjusted)

#     return(obj)
# }