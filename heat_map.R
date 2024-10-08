# https://stackoverflow.com/questions/51937380/fast-post-hoc-computation-using-r

# https://stackoverflow.com/questions/52124359/anova-on-r-with-different-dependent-variables

options(digits = 6)

library(dplyr)
library("reshape")
library(ggplot2)
library(Amelia)
library(gridExtra)
library("scatterplot3d")

# https://stackoverflow.com/questions/69054275/loading-multiple-rdata-and-binding-into-a-single-data-frame

#https://gist.github.com/bannister/8002800

sysname <- Sys.info()["sysname"]

if(sysname == "Darwin") {
  setwd("~/Desktop/PhD/simulation_results") # example on mac machine
} else if(sysname == "Linux") {
  setwd('~/Desktop/synthetic_vol_forecasting/simulation_results') # example on linux machine
}

path <- getwd()

files <- list.files(path=path, pattern = "*_SunJun0216*.*Rdata$")
#results <- sapply(files, function(x) mget(load(x)), simplify = TRUE)
#output <- do.call(rbind, results)
load(files)
# rownames(output) <- NULL
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

print(columns_we_want)

#First specify which vol model one want to focus on
#vol_model_chosen <- reduced_df[reduced_df$vol_model == 1,]
vol_model_chosen <- reduced_df

#Look for missing values
#missmap(reduced_df, main = "Missing values vs observed")

# Make the dataframe numeric
vol_model_chosen <- as.data.frame(lapply(vol_model_chosen, as.numeric))

#Make vol model into factor
#reduced_df$vol_model <- as.factor(reduced_df$vol_model)

#Define outcomes variable(s)
vol_model_chosen$success <- as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj14)
vol_model_chosen$success_versus_mean <- as.integer(vol_model_chosen$QL_adj1 <= vol_model_chosen$QL_adj12)

df_only_one_outcome <- cbind(vol_model_chosen[,1:11],
                             vol_model_chosen$success,
                             vol_model_chosen$success_versus_mean )
names(df_only_one_outcome) <- c(names(df_only_one_outcome)[1:11], 'success', 'against_mean')

non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]

scatterplot3d(non_NA[,c(6,5,2)], pch = 16, color=non_NA$success)

success_glm <- glm(success ~.^2, data = non_NA[,c("mu_x"
                                          ,"sigma_x"
                                          , 'M21_M22_vol_mu_delta'
                                          , "mu_omega_star"
                                          ,"vol_shock_sd"
                                          ,"success")])

summary(success_glm)

against_mean_glm <- glm(against_mean ~.^2, data = non_NA[,c("mu_x"
                                                       ,"sigma_x"
                                                       , 'M21_M22_vol_mu_delta'
                                                       , "mu_omega_star"
                                                       ,"vol_shock_sd"
                                                       ,"against_mean")])

summary(against_mean_glm)

fitting_loss <- lm(log(dbw_loss1) ~.^2, data = reduced_df[,c("mu_x"
                                                            ,"sigma_x"
                                                            , 'M21_M22_vol_mu_delta'
                                                            , "mu_omega_star"
                                                            ,"vol_shock_sd"
                                                            ,"dbw_loss1")])

summary(fitting_loss)

# https://statisticsglobe.com/heatmap-in-r

apply(non_NA, 2, function(x) unique(x))[c(
                                           "mu_x"
                                           ,"sigma_x"
                                          , 'M21_M22_vol_mu_delta'
                                           , "mu_omega_star"
                                           ,"vol_shock_sd" )]

hm_generator <- function(y_input
                         , x_input
                         , outcome
                         , ylab
                         , xlab)
{

  means <- non_NA %>% group_by({{x_input}}, {{y_input}}) %>% summarise(prop=mean({{outcome}}))
  means <- as.data.frame(sapply(means, as.numeric))
  means$prop <- round(means$prop, 2)

  count <- non_NA %>% group_by({{x_input}}, {{y_input}}) %>% count()

  means$n <- count$n

  print(means)

  ggp1 <- ggplot(means,
                aes(x = factor({{x_input}}), y = factor({{y_input}}), fill = prop)) +
    scale_fill_gradient(low="white",  high="red") +
    geom_tile() +
    geom_text(aes(label = paste(prop, '\n(',n,')', sep =''))) +
    guides(fill = guide_colourbar(title = "Success Proportion")) +
    ggtitle("Forecast Performance of Adjusted Forecast
            \n Each Square: Outperformance Proportion and (Simulation Count)") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 30)
          ) +
    labs(x = xlab, y = ylab)

  time_date <- gsub(" ", "", gsub(':', '', format(Sys.time(), "%b%d_%X_%Y")), fixed = TRUE)

  file_name <- paste('~/Desktop/PhD/synthetic_vol_forecasting/simulation_plots/'
                     ,time_date
                     ,'_'
                     ,ylab
                     ,'_'
                     ,xlab
                     ,'.png'
                     ,sep = '')

  ggsave(filename=file_name, plot = ggp1, width = 8, height = 6, dpi = 300)

  dev.off()

}

## Subsets to include in paper:

## Subset 1:
non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]
non_NA <- non_NA[non_NA$mu_x == 0.125,]
non_NA <- non_NA[non_NA$sigma_x == .125,]
#non_NA <- non_NA[non_NA$M21_M22_vol_mu_delta == .2,]
non_NA <- non_NA[non_NA$mu_omega_star == .125,]
#non_NA <- non_NA[non_NA$vol_shock_sd == .125,]

hm_generator(y_input = M21_M22_vol_mu_delta
             ,x_input = vol_shock_sd
             , success
             ,expression(mu[delta])
             ,expression(sigma[u]))

#Comment: it is hard to explain why increasing vol_shock_sd increases prop

## Subset 2:
non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]
non_NA <- non_NA[non_NA$mu_x == .5,]
non_NA <- non_NA[non_NA$sigma_x == .125,]
#non_NA <- non_NA[non_NA$M21_M22_vol_mu_delta == .2,]
non_NA <- non_NA[non_NA$mu_omega_star == .125,]
#non_NA <- non_NA[non_NA$vol_shock_sd == .125,]

hm_generator(y_input = M21_M22_vol_mu_delta
             ,x_input = vol_shock_sd
             , success
             ,expression(mu[delta])
             ,expression(sigma[u]))

#Comment: by increasing mu_x, we not only increase effect of M21_M22_vol_mu_delta,
#we also see negative effect of vol_shock_sd

## Subset 3:
non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]
non_NA <- non_NA[non_NA$mu_x == 1,]
non_NA <- non_NA[non_NA$sigma_x == .125,]
#non_NA <- non_NA[non_NA$M21_M22_vol_mu_delta == .2,]
non_NA <- non_NA[non_NA$mu_omega_star == .125,]
#non_NA <- non_NA[non_NA$vol_shock_sd == .125,]

hm_generator(y_input = M21_M22_vol_mu_delta
             ,x_input = vol_shock_sd
             , success
             ,expression(mu[delta])
             ,expression(sigma[u]))

#Comment: with mu_x at .5, the phenomena that we're tracking are somewhere in
#the middle

## Subset 4:
non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]
#non_NA <- non_NA[non_NA$mu_x == .5,]
non_NA <- non_NA[non_NA$sigma_x == .125,]
#non_NA <- non_NA[non_NA$M21_M22_vol_mu_delta == .2,]
non_NA <- non_NA[non_NA$mu_omega_star == .125,]
non_NA <- non_NA[non_NA$vol_shock_sd == .125,]

hm_generator(y_input = M21_M22_vol_mu_delta
             ,x_input = mu_x
             , success
             ,expression(mu[delta])
             ,expression(mu[v]))

#Comment: mu_x is on display

## Subset 5:
non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]
#non_NA <- non_NA[non_NA$mu_x == .5,]
non_NA <- non_NA[non_NA$sigma_x == .125,]
#non_NA <- non_NA[non_NA$M21_M22_vol_mu_delta == .2,]
non_NA <- non_NA[non_NA$mu_omega_star == .125,]
non_NA <- non_NA[non_NA$vol_shock_sd == 1,]

hm_generator(y_input = M21_M22_vol_mu_delta
             ,x_input = mu_x
             , success
             ,expression(mu[delta])
             ,expression(mu[v]))

#Comment: with vol_shock_sd = 1, we get a much more muted effect.
#Also, notice the number of non-convergent simulations.

## Subset 6:
non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]
#non_NA <- non_NA[non_NA$mu_x == .5,]
non_NA <- non_NA[non_NA$sigma_x == .125,]
non_NA <- non_NA[non_NA$M21_M22_vol_mu_delta == .125,]
non_NA <- non_NA[non_NA$mu_omega_star == .125,]
#non_NA <- non_NA[non_NA$vol_shock_sd == 1,]

hm_generator(y_input = mu_x
             ,x_input = vol_shock_sd
             , success
             ,expression(mu[v])
             ,expression(sigma[u]))

#Comment: we see a strange increase in proportion with rising vol_shock_sd

# I WILL SKIP THIS ONE.

## Subset 7:
non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]
non_NA <- non_NA[non_NA$mu_x == .125,]
non_NA <- non_NA[non_NA$sigma_x == .125,]
non_NA <- non_NA[non_NA$M21_M22_vol_mu_delta == .125,]
#non_NA <- non_NA[non_NA$mu_omega_star == .125,]
#non_NA <- non_NA[non_NA$vol_shock_sd == 1,]

hm_generator(y_input = mu_omega_star
             ,x_input = vol_shock_sd
             , success
             ,expression(mu[omega^"*"])
             ,expression(sigma[u]))

## Subset 8:
non_NA <- df_only_one_outcome[complete.cases(df_only_one_outcome),]
non_NA <- non_NA[non_NA$mu_x == .125,]
non_NA <- non_NA[non_NA$sigma_x == .125,]
non_NA <- non_NA[non_NA$M21_M22_vol_mu_delta == 2,]
#non_NA <- non_NA[non_NA$mu_omega_star == .125,]
#non_NA <- non_NA[non_NA$vol_shock_sd == 1,]

hm_generator(y_input = mu_omega_star
             ,x_input = vol_shock_sd
             , success
             ,expression(mu[omega^"*"])
             ,expression(sigma[u]))

#Let's generate 3d scatterplots

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
  ggtitle("Forecast Performance of Adjusted Forecast
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
  ggtitle("Forecast Performance of Adjusted Forecast
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
