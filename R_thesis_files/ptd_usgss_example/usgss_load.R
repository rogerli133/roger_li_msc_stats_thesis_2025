set.seed(1)

usgss_df <- read_csv("data/usgss02_labelled.csv")

usgss_98 <- usgss_df %>% 
  filter(year != 2002) %>% 
  select(-year)
usgss_02 <- usgss_df %>% 
  filter(year == 2002) %>% 
  select(-year)