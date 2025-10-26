library(haven)

data <- read_sav("../data/SESR 1 - II Group Discussion from Fall 2022 Working.sav")

filtered_data <- data %>% 
  filter(Include == 1)

rm(data)

codebook <- look_for(filtered_data)

print("'filtered_data' object created.")