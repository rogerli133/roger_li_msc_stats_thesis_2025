# ppi_plusplus_logistic requires matrices as inputs
# data needs additional intercept column
X_l <- cbind(as.matrix(wine_exp[, -12]))
Y_l <- as.matrix(wine_exp[, 12])
f_l <- as.matrix(wine_compare[, 12])
X_u <- cbind(as.matrix(wine_nov[, -12]))
f_u <- as.matrix(wine_nov[, 12])

# stores all ppi++ outputs and intermediate values
wine_ppi <- ppi_plusplus_logistic(
  X_l = cbind(1, X_l),   # X
  Y_l = Y_l,             # Y
  f_l = f_l,             # f(X)
  X_u = cbind(1, X_u),   # ~X~
  f_u = f_u,             # f(~X~)
)

# coef and se estimates for ppi++
ppi_coef <- t(wine_ppi$est)[-1]
ppi_se <- wine_ppi$se[-1]