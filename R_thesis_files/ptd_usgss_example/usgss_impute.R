# create the ml imputation on old census data
train_x <- usgss_98 %>% select(-kids) %>% data.matrix()
train_y <- usgss_98$kids
dtrain <- lgb.Dataset(data = train_x, label = train_y)

params <- list(
  objective = "regression",  
  metric = "regression", 
  learning_rate = 0.05,
  num_leaves = 30,
  max_depth = -1, 
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq = 5
)

lgb_model <- lgb.train(
  params,
  dtrain,
  nrounds = 1000,
  valids = list(train = dtrain),
  early_stopping_rounds = 50,
  verbose = 0
)

# predict number of kids with new census data
test_x <- usgss_02 %>% select(-kids) %>% data.matrix()
usgss02_preds <- predict(lgb_model, test_x)

usgss_imp <- usgss_02
usgss_imp$kids <- usgss02_preds