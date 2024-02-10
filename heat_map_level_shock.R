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
files <- list.files(path=path, pattern = ".*ThuJan0421:45:032024*.*Rdata$")

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

  title <- paste("SVF Outperformance of Unadjusted GARCH Forecast
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



# Now we do the same as above, but without varying the number of
# extra measurement days printed in the plots

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

  title <- paste("SVF Outperformance of
          \n Unadjusted GARCH Forecast
          \n Each Square: Outperformance Proportion and (Simulation Count)", sep = '')

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

#Or if we want to save the first plot as a png, we use
time_date <- gsub(" ", "", format(Sys.time(), "%a%b%d%X%Y"), fixed = TRUE)
png_file_name <- paste(
  "/home/david/Desktop/synthetic_vol_forecasting/simulation_plots/"
  ,'dual_level_vol_shock_'
  ,time_date
  ,""
  ,".png"
  ,sep="")
png(png_file_name)
print(ggp_list[[1]])
dev.off()
