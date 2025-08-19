set.seed(1)

wine <- read_csv("data/wine_quality_labelled.csv")

# randomForest needs response to be a factor for classification
winedf <- wine %>% 
  mutate(high_quality = factor(high_quality))

# separate data into expert and novice dfs
wine_exp <- winedf %>% 
  filter(label == "expert") %>% 
  select(-label)

wine_nov <- winedf %>% 
  filter(label == "novice") %>% 
  select(-label)