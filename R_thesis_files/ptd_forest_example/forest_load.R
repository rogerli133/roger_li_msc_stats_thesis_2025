set.seed(1)

# ground truth data with matching ~D~_hard values
cover_df_all <- read_csv("data/forest_cover_all.csv")[, -1]

# map-product data
cover_df_mp <- read_csv("data/forest_cover_mp.csv")[, -1] %>% 
  rename(
    cover = cover_pred,
    pop = pop_pred,
    elev = elev_truth
  )

# ground truth data: D_easy, D_hard
cover_df_gt <- cover_df_all %>% 
  select(ends_with("truth")) %>% 
  rename(
    cover = cover_truth,
    pop = pop_truth,
    elev = elev_truth
  )

# ground truth data: D_easy, ~D~_hard
cover_df_compare <- cover_df_all %>% 
  select(cover_pred, pop_pred, elev_truth) %>% 
  rename(
    cover = cover_pred,
    pop = pop_pred,
    elev = elev_truth
  )

# ALL ground truth points (including those excluded from example); proxy for true param
cover_df_true <- read_csv("data/forest_cover_true.csv")[, -1] %>% 
  rename(
    cover = cover_truth,
    pop = pop_truth,
    elev = elev_truth
  )