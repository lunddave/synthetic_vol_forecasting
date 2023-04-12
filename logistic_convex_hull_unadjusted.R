# https://stackoverflow.com/questions/51937380/fast-post-hoc-computation-using-r

# https://stackoverflow.com/questions/52124359/anova-on-r-with-different-dependent-variables

library(dplyr)

options(digits = 6)

load("/home/david/Desktop/synthetic_vol_forecasting/simulation_results/simcount_2_savetime_ThuMar3016:56:042023_runtime_4.892_hr__grid_size2880_recovery_1_permute_0.Rdata")

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

#Make vol model into factor
#reduced_df$vol_model <- as.factor(reduced_df$vol_model)

df_only_one_outcome <- cbind(vol_model_chosen[,c(1:11)], as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14))

non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]

model <- glm(non_NA[,12] ~. , family=binomial(link='logit'),data=non_NA[,c(1:5)])
summary(model)

model <- glm(non_NA[,12] ~.^2 , family=binomial(link='logit'),data=non_NA[,c(1:5)])
summary(model)

#Let's do an analysis by volatility model
library(dplyr)

par(mfrow = c(1,2))

## DOMINATING SUBSET
c1 <- non_NA$n > 4
#c1 <- TRUE
c2 <- non_NA$p > 2
c2 <- TRUE
c3 <- non_NA$arch_param < .4
#c3 <- TRUE
c4 <- non_NA$garch_param < .4
#c4 <- TRUE
c5 <- non_NA$vol_shock_length > 0
c5 <- TRUE
c6 <- non_NA$vol_sig_noise_ratio > 25
#c7 <- TRUE
dom_subset <- non_NA[c1 & c2 & c3 & c4 & c5 & c6,]
nrow(dom_subset)
p <- mean(dom_subset$`as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14)`)
p
prop.test(x=sum(dom_subset$`as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14)`), 
          n=nrow(dom_subset), conf.level=.95, correct=FALSE)

#We fix values for all parameters EXCEPT vol_sig_noise_ratio
Synth_Vol_Dominating_Proportion <- c()
sig_noise_ratio <- seq(min(non_NA$vol_sig_noise_ratio),50, 2)
for (val in sig_noise_ratio){
  c6 <- non_NA$vol_sig_noise_ratio > val
  #c7 <- TRUE
  dom_subset <- non_NA[c1 & c2 & c3 & c4 & c5 & c6,]
  p <- mean(dom_subset$`as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14)`)
  Synth_Vol_Dominating_Proportion <- c(Synth_Vol_Dominating_Proportion, p)
}
plot(x = sig_noise_ratio, y = Synth_Vol_Dominating_Proportion, col = 'red', 
     main = 'Dominating Subset', ylim = c(0,1))

## COIN-FLIP SUBSET
c1 <- non_NA$n < 20
c2 <- non_NA$p < 12
c2 <- TRUE
c3 <- non_NA$arch_param > .1
#c3 <- TRUE
c4 <- non_NA$garch_param > .1
#c4 <- TRUE
c5 <- non_NA$vol_shock_length < 2
c5 <- TRUE
c6 <- non_NA$vol_sig_noise_ratio < 10
#c7 <- TRUE
dom_subset <- non_NA[c1 & c2 & (c3 | c4) & c5 & c6,]
nrow(dom_subset)
p <- mean(dom_subset$`as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14)`)
p
prop.test(x=sum(dom_subset$`as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14)`), 
          n=nrow(dom_subset), conf.level=.95, correct=FALSE)

#We fix values for all parameters EXCEPT vol_sig_noise_ratio
Synth_Vol_Dominating_Proportion <- c()
sig_noise_ratio <- seq(min(non_NA$vol_sig_noise_ratio),50, 2)
for (val in sig_noise_ratio){
  c6 <- non_NA$vol_sig_noise_ratio > val
  #c7 <- TRUE
  dom_subset <- non_NA[c1 & c2 & c3 & c4 & c5 & c6,]
  p <- mean(dom_subset$`as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14)`)
  Synth_Vol_Dominating_Proportion <- c(Synth_Vol_Dominating_Proportion, p)
}
plot(x = sig_noise_ratio, y = Synth_Vol_Dominating_Proportion, col = 'red', 
     main = 'Coin-flip subset', ylim = c(0,1))

means <- non_NA %>% group_by(vol_model) %>% summarise(across(everything(), list(means)))
means <- as.data.frame(means)
means
