options(digits = 6)

library(dplyr)
library(ggplot2)

sysname <- Sys.info()["sysname"]

if(sysname == "Darwin") {
  setwd("~/Desktop/PhD/synthetic_vol_forecasting/exponential_sims") # example on mac machine
} else if(sysname == "Linux") {
  setwd('~/Desktop/synthetic_vol_forecasting/simulation_results') # example on linux machine
}

load("simcount_40_savetime_TueNov2614_14_562024_runtime_0.381_hr__grid_size36_recovery_0.614_seed_1986.Rdata")

df <- as.matrix(output)
df <- data.frame(df)
df$simplex_dominates <- as.integer(as.logical(df$simplex_dominates))

head(df)

apply(df,2,class)

df <- df[,-c(10)]

df <- data.frame(apply(df, 2, function(x) as.numeric(as.character(x))))

count_unique <- function(input) length(unique(input))



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
  
  file_name <- gsub('\\[|\\]','',file_name)
  
  ggsave(filename=file_name, plot = ggp1, width = 8, height = 6, dpi = 300)
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

#Comment: we can see that as the variation in the covariates increases, 
#we get better performance from the simplex.

temp <- df[df$shock_sd == 2,]
hm_generator(temp
             ,y_input = mu_delta
             ,x_input = covariate_sigma
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[x]))

#Comment: performance drops as noise increases

temp <- df[df$shock_sd == 3,]
hm_generator(temp
             ,y_input = mu_delta
             ,x_input = covariate_sigma
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[x]))

#Comment: performance drops as noise increases


## Subset 2:
temp <- df[df$mu_delta == .001,]

hm_generator(temp
             ,y_input = covariate_sigma
             ,x_input = shock_sd
             , simplex_dominates
             ,expression(sigma[x])
             ,expression(sigma[epsilon]))

temp <- df[df$mu_delta == 0.01,]

hm_generator(temp
             ,y_input = covariate_sigma
             ,x_input = shock_sd
             , simplex_dominates
             ,expression(sigma[x])
             ,expression(sigma[epsilon]))

temp <- df[df$mu_delta == .05,]

hm_generator(temp
             ,y_input = covariate_sigma
             ,x_input = shock_sd
             , simplex_dominates
             ,expression(sigma[x])
             ,expression(sigma[epsilon]))

#Comment: small mu_delta is bad because that means psi could be negative

## Subset 3:
temp <- df[df$covariate_sigma == .01,]

hm_generator(temp
             ,y_input = mu_delta
             ,x_input = shock_sd
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[epsilon]))

temp <- df[df$covariate_sigma == 1,]

hm_generator(temp
             ,y_input = mu_delta
             ,x_input = shock_sd
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[epsilon]))

temp <- df[df$covariate_sigma == 2,]

hm_generator(temp
             ,y_input = mu_delta
             ,x_input = shock_sd
             , simplex_dominates
             ,expression(mu[delta])
             ,expression(sigma[epsilon]))

