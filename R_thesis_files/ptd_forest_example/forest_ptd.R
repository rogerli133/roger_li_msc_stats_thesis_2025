# ptd bootstrapping
cover_ptd <- PTD_bootstrap.glm(
  true_data_completeSamp = cover_df_gt,             # ground-truth data
  predicted_data_completeSamp = cover_df_compare,   # ground_truth easy with mp hard
  predicted_data_incompleteSamp = cover_df_mp,      # full mp data
  regFormula.glm = "cover ~ elev + pop",
  GLM_type = "linear",
  alpha = 0.05,
  B = 2000,                                         # number of bootstraps
  TuningScheme = "Optimal", 
  speedup = TRUE
)

ptd_coef <- cover_ptd$PTD_estimate[-1]
ptd_lower <- cover_ptd$PTD_Boot_CIs[-1, 1]
ptd_upper <- cover_ptd$PTD_Boot_CIs[-1, 2]