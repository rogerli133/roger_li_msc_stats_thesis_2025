set.seed(1)

wine <- read_csv("data/wine_quality_labelled.csv")

# separate data into expert and novice dfs

# 100 x 12 df of wines tasted by both experts and novices
# responses here are from the expert
wine_exp <- wine %>% 
  filter(label == "expert") %>% 
  select(-label)

# 100 x 12 df of wines tasted by both experts and novices
# responses here are from the novices
wine_compare <- wine %>% 
  filter(label == "compare") %>% 
  select(-label)

# 3000 x 12 df of wines tasted only by novies
wine_nov <- wine %>% 
  filter(label == "novice") %>% 
  select(-label)

