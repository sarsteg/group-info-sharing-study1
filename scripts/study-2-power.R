library(simr)
library(lme4)

# 1. Create a dummy data frame representing your factorial design
n_subj <- 240 
subj_ids <- 1:n_subj
structure <- c("S1", "S2", "S3")
motivation <- c("M1", "M2")
info_type <- c("T1", "T2")
importance <- c("High", "Low")
items <- 1:5

df <- expand.grid(ResponseId=subj_ids, structure=structure, motivation=motivation, 
                  info_type=info_type, importance=importance, item=items)

# 2. Define the Fixed Effects (The Intercept and the 4-way interaction)
# other effects to 0 for a conservative test of the 4-way term
# 0.18 is intercept, 1.34 is the last interaction term
fixed <- c(0.18, rep(0, 22), 1.34) 

# 3. Random Effects (Subject variance)
rand <- 0.1 

# 4. Create the model object
model <- makeGlmer(y ~ structure*motivation*info_type*importance + (1|ResponseId), 
                   family="poisson", fixef=fixed, VarCorr=rand, data=df)

# 5. Run the power analysis
power_sim <- powerSim(model, fixed("structure:motivation:info_type:importance", "lr"), nsim=100)
print(power_sim)