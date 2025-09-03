set.seed(1)

wine <- read_csv("data/wine_quality_labelled.csv")

# separate data into expert and novice dfs
wine_exp <- wine %>% 
  filter(label == "expert") %>% 
  select(-label)

wine_nov <- wine %>% 
  filter(label == "novice") %>% 
  select(-label)

# separate the wines tasted by both expert and novices (Y1, f(X1)), ..., (Yn, f(Xn))
wine_compare <- wine %>% 
  filter(label == "compare") %>% 
  select(-label)