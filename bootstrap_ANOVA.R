# https://stackoverflow.com/questions/51937380/fast-post-hoc-computation-using-r

# https://stackoverflow.com/questions/52124359/anova-on-r-with-different-dependent-variables

options(digits = 3)

load("/home/david/Desktop/synthetic_vol_forecasting/simulation_results/simcount_20_savetime_Sat Dec 31 20:02:35 2022.Rdata")

library(lmboot)

# We look for the columns that have no variation, i.e. those with only 1 value.  
unique_count_df <- apply(output, 2, function(x) length(unique(x)))

# We drop these columns.
columns_we_want <- names(unique_count_df)[unique_count_df != 1]
reduced_df <- output[names(output) %in% columns_we_want]

# Make the dataframe numeric
reduced_df <- as.data.frame(lapply(reduced_df, as.numeric))

#Make vol model into factor
reduced_df$vol_model <- as.factor(reduced_df$vol_model)

df_only_one_outcome <- cbind(reduced_df[,c(1:10)], reduced_df$QL_adj1)

non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]

bootstrap_samples <- nrow(non_NA) * 1.5

myANOVA2 <- ANOVA.boot(non_NA[,11] ~. , data=non_NA[,c(1:10)], B = bootstrap_samples)
names(myANOVA2)
myANOVA2$terms
myANOVA2$`p-values` #bootstrap p-values for 2-way interactions model

ANOVA_df <- matrix(c(myANOVA2$terms[-1], as.numeric(myANOVA2$`p-values`) ), ncol = 2)
ANOVA_df <- ANOVA_df[order(ANOVA_df[,2], decreasing = FALSE),]
ANOVA_df
