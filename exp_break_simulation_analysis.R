options(digits = 6)

library(dplyr)
library("reshape")
library(ggplot2)
library(Amelia)
library(gridExtra)
library("scatterplot3d")

sysname <- Sys.info()["sysname"]

if(sysname == "Darwin") {
  setwd("~/Desktop/PhD/synthetic_vol_forecasting/exponential_sims") # example on mac machine
} else if(sysname == "Linux") {
  setwd('~/Desktop/synthetic_vol_forecasting/simulation_results') # example on linux machine
}

load("simcount_10_savetime_TueNov2613_20_022024_runtime_0.176_hr__grid_size27_recovery_0.722_seed_1986.Rdata")

df <- as.matrix(output)
df <- data.frame(df)
names(df)

count_unique <- function(input) length(unique(input))
apply(df,2,count_unique)

df$simplex_dominates <- as.integer(as.logical(df$simplex_dominates))

# Group by mean of multiple columns
# df2 <- df %>% group_by(n
#                        ,p
#                        ,covariate_sigma
#                       ,shock_sd
#                         ,mu_delta) %>% 
#   summarise(mean_simplex=mean(simplex_dominates),
#             .groups = 'drop') %>%
#   as.data.frame()


hm_generator <- function(input_df
                         ,y_input
                         , x_input
                         , outcome
                         , ylab
                         , xlab)
{
  
  means <- input_df %>% group_by({{x_input}}, {{y_input}}) %>% summarise(prop=mean({{outcome}}))
  means <- as.data.frame(sapply(means, as.numeric))
  means$prop <- round(means$prop, 2)
  
  count <- input_df %>% group_by({{x_input}}, {{y_input}}) %>% count()
  
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
temp <- df[df$shock_sd == 1,]

hm_generator(temp
              ,y_input = mu_delta
             ,x_input = covariate_sigma
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[x]))

#Comment: 

## Subset 2:
temp <- df[df$mu_delta == .05,]

hm_generator(temp
             ,y_input = covariate_sigma
             ,x_input = shock_sd
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[epsilon]))

#Comment: 

## Subset 3:
temp <- df[df$shock_sd == 1,]

hm_generator(temp
             ,y_input = mu_delta
             ,x_input = covariate_sigma
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[x]))

#Comment: 

## Subset 4:
temp <- df[df$shock_sd == 1,]

hm_generator(temp
             ,y_input = mu_delta
             ,x_input = covariate_sigma
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[x]))

#Comment: mu_x is on display

## Subset 5:
temp <- df[df$shock_sd == 1,]

hm_generator(temp
             ,y_input = mu_delta
             ,x_input = covariate_sigma
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[x]))

#Comment: 

## Subset 6:
temp <- df[df$shock_sd == 1,]

hm_generator(temp
             ,y_input = mu_delta
             ,x_input = covariate_sigma
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[x]))

#Comment: 
