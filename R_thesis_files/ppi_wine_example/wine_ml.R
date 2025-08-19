# create the ml model on novice data
wine_rf <- randomForest(
  high_quality ~ .,
  data = wine_nov, 
  importance = TRUE,
  mtry = 3,
  ntree = 1000,
  replace = TRUE,
  sampsize = 100,
  strata = high_quality)

# create predictions
wine_nov_preds <- wine_rf$predicted           # f(X~)
wine_exp_preds <- predict(wine_rf, wine_exp)  # f(X)