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

dev.new(width=2, height=2)

#https://gist.github.com/bannister/8002800

sysname <- Sys.info()["sysname"]

if(sysname == "Darwin") {
  setwd("~/Desktop/PhD/simulation_results") # example on mac machine
} else if(sysname == "Linux") {
  setwd('~/Desktop/synthetic_vol_forecasting/simulation_results') # example on linux machine
}

path <- getwd()

files <- list.files(path=path, pattern = ".*simcount_500_savetime_MonMar1100:26:382024_runtime_2.209*.*Rdata$")
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
df_only_one_outcome <- cbind(vol_model_chosen[,1:11], vol_model_chosen$success)
names(df_only_one_outcome) <- c(names(df_only_one_outcome)[1:11], 'success')

non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]

# https://statisticsglobe.com/heatmap-in-r

non_NA <- non_NA[non_NA$p == 45,]

means <- non_NA %>% group_by(vol_shock_sd, M21_M22_vol_mu_delta) %>% summarise(prop=mean(success))
means <- as.data.frame(sapply(means, as.numeric))
means$prop <- round(means$prop, 2)

count <- non_NA %>% group_by(vol_shock_sd, M21_M22_vol_mu_delta) %>% count()
count

means$n <- count$n

ggp1 <- ggplot(means,
              aes(x = factor(vol_shock_sd), y = factor(M21_M22_vol_mu_delta), fill = prop)) +
  scale_fill_gradient(low="white",  high="red") +
  geom_tile() +
  geom_text(aes(label = paste(prop, '\n(',n,')', sep =''))) +
  guides(fill = guide_colourbar(title = "Success Proportion")) +
  ggtitle("Synthetic Volatility Forecast Outperformance of Unadjusted Forecast
          \n Each Square: Outperformance Proportion and (Simulation Count)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Volatility Shock Standard Deviation", y = "Volatility Shock Mean")

ggp1

#
# # Now we write a function
# heatmap_maker <- function(df, var1, var2){
#
#   means <- df %>% group_by(var1, var2) %>%
#     summarise(prop=mean(success))
#   means <- as.data.frame(sapply(means, as.numeric))
#
#   ggp <- ggplot(means,
#                 aes(x = factor(var1), y = factor(var2), fill = prop)) +
#     scale_fill_gradient(low="white",  high="red") +
#     geom_tile()
#
#   ggp
# }
#
# heatmap_maker(non_NA, garch_param, vol_shock_sd)


#https://gist.github.com/bannister/8002800
path <- '/home/david/Desktop/synthetic_vol_forecasting/simulation_results'
files <- list.files(path=path, pattern = ".*_120.*Rdata$")
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
library(Amelia)
#missmap(reduced_df, main = "Missing values vs observed")

# Make the dataframe numeric
vol_model_chosen <- as.data.frame(lapply(vol_model_chosen, as.numeric))

#Make vol model into factor
#reduced_df$vol_model <- as.factor(reduced_df$vol_model)

vol_model_chosen$success <- as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14)
df_only_one_outcome <- cbind(vol_model_chosen[,1:11], vol_model_chosen$success)
names(df_only_one_outcome) <- c(names(df_only_one_outcome)[1:11], 'success')

non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]

# https://statisticsglobe.com/heatmap-in-r


means <- non_NA %>% group_by(vol_shock_sd, M21_M22_vol_mu_delta) %>% summarise(prop=mean(success))
means <- as.data.frame(sapply(means, as.numeric))
means$prop <- round(means$prop, 2)

count <- non_NA %>% group_by(vol_shock_sd, M21_M22_vol_mu_delta) %>% count()
count

means$n <- count$n

ggp2 <- ggplot(means,
              aes(x = factor(vol_shock_sd), y = factor(M21_M22_vol_mu_delta), fill = prop)) +
  scale_fill_gradient(low="white",  high="red") +
  geom_tile() +
  geom_text(aes(label = paste(prop, '\n(',n,')', sep =''))) +
  guides(fill = guide_colourbar(title = "Success Proportion")) +
  ggtitle("Synthetic Volatility Forecast Outperformance of Unadjusted Forecast
          \n Each Square: Outperformance Proportion and (Simulation Count)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Volatility Shock Standard Deviation", y = "Volatility Shock Mean")


#
# # Now we write a function
# heatmap_maker <- function(df, var1, var2){
#
#   means <- df %>% group_by(var1, var2) %>%
#     summarise(prop=mean(success))
#   means <- as.data.frame(sapply(means, as.numeric))
#
#   ggp <- ggplot(means,
#                 aes(x = factor(var1), y = factor(var2), fill = prop)) +
#     scale_fill_gradient(low="white",  high="red") +
#     geom_tile()
#
#   ggp
# }
#
# heatmap_maker(non_NA, garch_param, vol_shock_sd)


#https://gist.github.com/bannister/8002800
path <- '/home/david/Desktop/synthetic_vol_forecasting/simulation_results'
files <- list.files(path=path, pattern = ".*20_80.*Rdata$")
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
df_only_one_outcome <- cbind(vol_model_chosen[,1:11], vol_model_chosen$success)
names(df_only_one_outcome) <- c(names(df_only_one_outcome)[1:11], 'success')

non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]

# https://statisticsglobe.com/heatmap-in-r


means <- non_NA %>% group_by(vol_shock_sd, M21_M22_vol_mu_delta) %>% summarise(prop=mean(success))
means <- as.data.frame(sapply(means, as.numeric))
means$prop <- round(means$prop, 2)

count <- non_NA %>% group_by(vol_shock_sd, M21_M22_vol_mu_delta) %>% count()
count

means$n <- count$n

ggp3 <- ggplot(means,
              aes(x = factor(vol_shock_sd), y = factor(M21_M22_vol_mu_delta), fill = prop)) +
  scale_fill_gradient(low="white",  high="red") +
  geom_tile() +
  geom_text(aes(label = paste(prop, '\n(',n,')', sep =''))) +
  guides(fill = guide_colourbar(title = "Success Proportion")) +
  ggtitle("Synthetic Volatility Forecast Outperformance of Unadjusted Forecast
          \n Each Square: Outperformance Proportion and (Simulation Count)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Volatility Shock Standard Deviation", y = "Volatility Shock Mean")

grid.arrange(ggp1, ggp2, ggp3, ncol=3)


#
# # Now we write a function
# heatmap_maker <- function(df, var1, var2){
#
#   means <- df %>% group_by(var1, var2) %>%
#     summarise(prop=mean(success))
#   means <- as.data.frame(sapply(means, as.numeric))
#
#   ggp <- ggplot(means,
#                 aes(x = factor(var1), y = factor(var2), fill = prop)) +
#     scale_fill_gradient(low="white",  high="red") +
#     geom_tile()
#
#   ggp
# }
#
# heatmap_maker(non_NA, garch_param, vol_shock_sd)


means$prop <- means$prop - shortest_series$prop
