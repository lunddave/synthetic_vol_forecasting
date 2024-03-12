# https://stackoverflow.com/questions/51937380/fast-post-hoc-computation-using-r

# https://stackoverflow.com/questions/52124359/anova-on-r-with-different-dependent-variables

options(digits = 6)

library(dplyr)
library("reshape")
library(dplyr)
library(ggplot2)
library(Amelia)
library(gridExtra)

# https://stackoverflow.com/questions/69054275/loading-multiple-rdata-and-binding-into-a-single-data-frame

#dev.new(width=2, height=2)

#https://gist.github.com/bannister/8002800
path <- '~/Desktop/PhD/simulation_results'
files <- list.files(path=path, pattern = "*simcount_500_savetime_MonMar1100:26:382024_runtime_2.209*.*Rdata$")
setwd(path)
results <- sapply(files, function(x) mget(load(x)), simplify = TRUE)
output <- do.call(rbind, results)
rownames(output) <- NULL
data.frame(sapply(output,class))

length_unique <- function(x) {return(length(unique(x)))}
data.frame(sapply(output,length_unique))

output <- as.data.frame(sapply(output, as.numeric)) #<- sapply is here
data.frame(sapply(output,class))



#Check that things vary correctly cross vol models
check <- output %>% group_by(vol_model) %>%
  summarise(across(everything(), length_unique),
            .groups = 'drop')  %>% as.data.frame()

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
#missmap(reduced_df, main = "Missing values vs observed")

# Make the dataframe numeric
vol_model_chosen <- as.data.frame(lapply(vol_model_chosen, as.numeric))

#Make vol model into factor
#reduced_df$vol_model <- as.factor(reduced_df$vol_model)

vol_model_chosen$success <- as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14)
df_only_one_outcome <- cbind(vol_model_chosen, vol_model_chosen$success)
names(df_only_one_outcome) <- c(names(df_only_one_outcome)[1:(length(df_only_one_outcome)-1)], 'success')

non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]

non_NA <- non_NA[non_NA$p == 15,]

lr <- glm(success ~ (dbw_loss1)/X1_norm1 , data = non_NA, family = "binomial")
summary(lr)

