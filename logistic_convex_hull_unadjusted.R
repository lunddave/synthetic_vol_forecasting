# https://stackoverflow.com/questions/51937380/fast-post-hoc-computation-using-r

# https://stackoverflow.com/questions/52124359/anova-on-r-with-different-dependent-variables

library(dplyr)

options(digits = 6)

load("/home/david/Desktop/synthetic_vol_forecasting/simulation_results/simcount_2_savetime_TueMar1420:54:122023_runtime_0.389301863577631_grid_size2592_recovery_0.667631172839506_permute_0.Rdata")

library(lmboot)

length_unique <- function(x) {return(length(unique(x)))}

#Check that things vary correctly cross vol models
check <- output %>% group_by(vol_model) %>% 
  summarise(across(everything(), length_unique),
            .groups = 'drop')  %>%
  as.data.frame()

check

# We look for the columns that have no variation, i.e. those with only 1 value.  
unique_count_df <- apply(output, 2, function(x) length(unique(x)))

# We drop these columns.
columns_we_want <- names(unique_count_df)[unique_count_df != 1]
reduced_df <- output[names(output) %in% columns_we_want]

#First specify which vol model one want to focus on
#vol_model_chosen <- reduced_df[reduced_df$vol_model == 1,]
vol_model_chosen <- reduced_df

#Look for missing values
library(Amelia)
missmap(reduced_df, main = "Missing values vs observed")

# Make the dataframe numeric
vol_model_chosen <- as.data.frame(lapply(vol_model_chosen, as.numeric))

# Now remove rows with output$M21_M22_vol_mu_delta absurdly high
vol_model_chosen <- vol_model_chosen[vol_model_chosen$M21_M22_vol_mu_delta != .5, ]

#Make vol model into factor
#reduced_df$vol_model <- as.factor(reduced_df$vol_model)

df_only_one_outcome <- cbind(vol_model_chosen[,c(1:11)], as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14))

non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]

model <- glm(non_NA[,12] ~. , family=binomial(link='logit'),data=non_NA[,c(1:5,10)])
summary(model)

model <- glm(non_NA[,12] ~.^2 , family=binomial(link='logit'),data=non_NA[,c(1:5,10)])
summary(model)

#Let's do an analysis by volatility model
library(dplyr)

## Hyperplanes
c1 <- non_NA$n > 5
c2 <- non_NA$p > 5
c3 <- non_NA$arch_param < .3
c4 <- non_NA$garch_param < .3
c5 <- non_NA$vol_shock_length > 1
c6 <- non_NA$vol_sig_noise_ratio > 4.9
dom_subset <- non_NA[c1 & c2 & c3 & c4 & c5 & c6,]
nrow(dom_subset)
mean(dom_subset$`as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14)`)

means <- non_NA %>% group_by(vol_model) %>% summarise(across(everything(), list(means)))
means <- as.data.frame(means)
means
