# https://stackoverflow.com/questions/51937380/fast-post-hoc-computation-using-r

# https://stackoverflow.com/questions/52124359/anova-on-r-with-different-dependent-variables

options(digits = 6)

library(dplyr)
library("reshape")
library(dplyr)
library(ggplot2)
library(Amelia)
library(gridExtra)
library(ggpubr)


# https://stackoverflow.com/questions/69054275/loading-multiple-rdata-and-binding-into-a-single-data-frame

# par(mfrow = c(1,3))
# dev.new(width=2, height=2)

#https://gist.github.com/bannister/8002800
path <- '/home/david/Desktop/simulation_results'
#files <- list.files(path=path, pattern = ".*Apr16.*Rdata$")
files <- list.files(path=path, pattern = ".*SunDec1017:47:062023*.*Rdata$")

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

ggp_list = list()

for (extra_day in unique(vol_model_chosen$extra_measurement_days)){
  
  extra_day_set_chosen <- vol_model_chosen[vol_model_chosen$extra_measurement_days == extra_day,]
  
  extra_day_set_chosen$success <- as.integer(extra_day_set_chosen$QL_adj1 <= extra_day_set_chosen$QL_adj14)
  df_only_one_outcome <- cbind(extra_day_set_chosen[,1:11], extra_day_set_chosen$success)
  names(df_only_one_outcome) <- c(names(df_only_one_outcome)[1:11], 'success')
  
  non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]
  
  # https://statisticsglobe.com/heatmap-in-r
  
  means <- non_NA %>% 
    group_by(mu_eps_star, M21_M22_vol_mu_delta) %>% summarise(prop=mean(success),.groups = 'drop')
  means <- as.data.frame(sapply(means, as.numeric))
  means$prop <- round(means$prop, 2)
  
  means
  
  count <- non_NA %>% group_by(mu_eps_star, M21_M22_vol_mu_delta) %>% count()
  count
  
  means$n <- count$n
  
  title <- paste("Synthetic Volatility Forecast Outperformance of Unadjusted GARCH Forecast
          \n Each Square: Outperformance Proportion and (Simulation Count)
          \n Number of Extra Measurement Days = ",
                 extra_day, sep = '')
  
  ggp1 <- ggplot(means,
                 aes(x = factor(mu_eps_star), y = factor(M21_M22_vol_mu_delta), fill = prop)) +
    scale_fill_gradient(low="white",  high="red") +
    geom_tile() +
    geom_text(aes(label = paste(prop, '\n(',n,')', sep =''))) +
    guides(fill = guide_colourbar(title = "Success Proportion")) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "M1 Level Shock Mean", y = "M21 Volatility Shock Mean")
  
  ggp_list[[extra_day+1]] <- ggp1
  
}

#http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
# ggarrange(ggp_list[[1]], ggp_list[[2]], ggp_list[[3]], ggp_list[[4]], 
#           labels = c("A", "B", "C", "D"),
#           ncol = 2, nrow = 2)

#https://stackoverflow.com/questions/26034177/save-multiple-ggplots-using-a-for-loop
# Another option: create pdf where each page is a separate plot.

time_date <- gsub(" ", "", format(Sys.time(), "%a%b%d%X%Y"), fixed = TRUE)
pdf_file_name <- paste(
                       "/home/david/Desktop/synthetic_vol_forecasting/simulation_plots/"
                       ,'dual_level_vol_shock_'
                       ,time_date
                       ,""
                       ,".pdf"
                       ,sep="")



pdf(pdf_file_name)
for (i in 1:length(ggp_list)) {
  print(ggp_list[[i]])
}
dev.off()

#http://127.0.0.1:20841/graphics/ed675dcc-4353-4d42-bb93-c844f62500cf.png
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
path <- '/home/david/Desktop/simulation_results'
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
