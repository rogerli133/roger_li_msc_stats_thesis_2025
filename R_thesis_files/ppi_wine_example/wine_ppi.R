# recreate dfs to undo the factor from earlier
wine_exp <- wine %>% 
  filter(label == "expert") %>% 
  select(-label)

wine_nov <- wine %>% 
  filter(label == "novice") %>% 
  select(-label)

# ppi_plusplus_logistic requires matrices as inputs
# data needs additional intercept column
X_l <- cbind(as.matrix(wine_exp[, -12]), 1)
Y_l <- as.matrix(wine_exp[, 12])
f_l <- as.matrix(as.numeric(as.character(wine_exp_preds)))
X_u <- cbind(as.matrix(wine_nov[, -12]), 1)
f_u <- as.matrix(as.numeric(as.character(wine_nov_preds)))

# stores all ppi++ outputs and intermediate values
wine_ppi <- ppi_plusplus_logistic(
  X_l = 1, X_l,
  Y_l = Y_l,
  f_l = f_l,
  X_u = 1, X_u, 
  f_u = f_u
)

# coef and se estimates for ppi++
ppi_coef <- t(wine_ppi$est)[-1]
ppi_se <- wine_ppi$se[-1]