ppi_plusplus_logistic <- function (X_l, Y_l, f_l, X_u, f_u, lhat = NULL, coord = NULL, 
                                   opts = NULL, w_l = NULL, w_u = NULL) 
{
  n <- nrow(f_l)
  N <- nrow(f_u)
  p <- ncol(X_u)
  w_l <- if (is.null(w_l)) 
    rep(1, n)
  else w_l/sum(w_l) * n
  w_u <- if (is.null(w_u)) 
    rep(1, N)
  else w_u/sum(w_u) * N
  use_u <- is.null(lhat) || lhat != 0
  ## CHANGED: data.frame(X_l, Y_l) -> data.frame(X_l)
  theta0 <- coef(glm(Y_l ~ . - 1, data = data.frame(X_l), 
                     family = binomial))
  est <- ppi_plusplus_logistic_est(X_l, Y_l, f_l, X_u, f_u, 
                                   opts = opts, lhat = lhat, coord = coord, w_l = w_l, w_u = w_u)
  stats <- logistic_get_stats(est, X_l, Y_l, f_l, X_u, f_u, 
                              w_l, w_u, use_u = use_u)
  if (is.null(lhat)) {
    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat, stats$grads_hat_unlabeled, 
                          stats$inv_hessian, clip = TRUE)
    return(ppi_plusplus_logistic(X_l, Y_l, f_l, X_u, f_u, 
                                 lhat = lhat, coord = coord, opts = opts, w_l = w_l, 
                                 w_u = w_u))
  }
  var_u <- cov(lhat * stats$grads_hat_unlabeled)
  var_l <- cov(stats$grads - lhat * stats$grads_hat)
  Sigma_hat <- stats$inv_hessian %*% (n/N * var_u + var_l) %*% 
    stats$inv_hessian
  return(list(est = est, se = sqrt(diag(Sigma_hat)/n), lambda = lhat, 
              rectifier_est = theta0 - est, var_u = var_u, var_l = var_l, 
              grads = stats$grads, grads_hat_unlabeled = stats$grads_hat_unlabeled, 
              grads_hat = stats$grads_hat, inv_hessian = stats$inv_hessian))
}

ppi_plusplus_logistic_est <- function (X_l, Y_l, f_l, X_u, f_u, lhat = NULL, coord = NULL, 
                                       opts = NULL, w_l = NULL, w_u = NULL) 
{
  n <- nrow(f_l)
  N <- nrow(f_u)
  p <- ncol(X_u)
  if (is.null(w_l)) 
    w_l <- rep(1, n)
  else w_l <- w_l/sum(w_l) * n
  if (is.null(w_u)) 
    w_u <- rep(1, N)
  else w_u <- w_u/sum(w_u) * N
  if (is.null(opts) || !("factr" %in% names(opts))) {
    opts <- list(factr = 1e-15)
  }
  ## CHANGED: data.frame(X_l, Y_l) -> data.frame(X_l)
  theta <- coef(glm(Y_l ~ . - 1, data = data.frame(X_l), 
                    family = binomial))
  theta <- matrix(theta, ncol = 1)
  lhat_curr <- ifelse(is.null(lhat), 1, lhat)
  rectified_logistic_loss <- function(theta) {
    sum(w_u * (-f_u * (X_u %*% theta) + log1pexp(X_u %*% 
                                                   theta))) * lhat_curr/N - sum(w_l * (-f_l * (X_l %*% 
                                                                                                 theta) + log1pexp(X_l %*% theta))) * lhat_curr/n + 
      sum(w_l * (-Y_l * (X_l %*% theta) + log1pexp(X_l %*% 
                                                     theta)))/n
  }
  rectified_logistic_grad <- function(theta) {
    lhat_curr/N * t(X_u) %*% (w_u * (plogis(X_u %*% theta) - 
                                       f_u)) - lhat_curr/n * t(X_l) %*% (w_l * (plogis(X_l %*% 
                                                                                         theta) - f_l)) + 1/n * t(X_l) %*% (w_l * (plogis(X_l %*% 
                                                                                                                                            theta) - Y_l))
  }
  est <- optim(par = theta, fn = rectified_logistic_loss, gr = rectified_logistic_grad, 
               method = "L-BFGS-B", control = list(factr = opts$factr))$par
  if (is.null(lhat)) {
    stats <- logistic_get_stats(est, X_l, Y_l, f_l, X_u, 
                                f_u, w_l, w_u)
    lhat <- calc_lhat_glm(stats$grads, stats$grads_hat, stats$grads_hat_unlabeled, 
                          stats$inv_hessian, clip = TRUE)
    return(ppi_plusplus_logistic_est(X_l, Y_l, f_l, X_u, 
                                     f_u, opts = opts, lhat = lhat, coord = coord, w_l = w_l, 
                                     w_u = w_u))
  }
  else {
    return(est)
  }
}