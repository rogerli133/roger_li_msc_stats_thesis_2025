# ptd bootstrapping
usgss_ptd <- PTD_bootstrap.glm(
  true_data_completeSamp = usgss_02,
  predicted_data_completeSamp = usgss_imp,
  predicted_data_incompleteSamp = usgss_98,
  regFormula.glm = "kids ~ .",
  GLM_type = "linear",
  alpha = 0.05,
  B = 2000,
  TuningScheme = "Optimal", 
  speedup = TRUE
)

# store coefficient estimates and confidence intervals
ptd_ci <- usgss_ptd$PTD_Boot_CIs[-1, ]
ptd_coef <- t(usgss_ptd$PTD_estimate)[-1]