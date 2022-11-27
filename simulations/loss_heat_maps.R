# https://stackoverflow.com/questions/51937380/fast-post-hoc-computation-using-r

# https://stackoverflow.com/questions/52124359/anova-on-r-with-different-dependent-variables

library(gplots)
library(RColorBrewer)
library(dplyr)

options(scipen = 9)

load("/home/david/Desktop/synthetic_vol_forecasting/simulation_results/simcount_52022-11-21 15:10:27.Rdata")

#We inspect the data
head(output_n_sim_1)

# We look for the columns that have no variation, i.e. those with only 1 value.  
unique_count_df <- apply(output_n_sim_1, 2, function(x) length(unique(x)))

# We drop these columns.
columns_we_want <- names(unique_count_df)[unique_count_df != 1]
reduced_df <- output_n_sim_1[names(output_n_sim_1) %in% columns_we_want]

# Make boolean column intger
reduced_df$beat_unadjusted <- as.integer(as.factor(reduced_df$beat_unadjusted))



reduced_df <- reduced_df[,-c(2)]

#Make vol model into factor
reduced_df$vol_model <- as.factor(reduced_df$vol_model)

# Make the dataframe numeric
reduced_df[,-c(1)] <- as.data.frame(lapply(reduced_df[,-c(1)], as.numeric))

# Make the dataframe numeric
reduced_df$linear_comb_names <- as.factor(reduced_df$linear_comb_names)

# Take only rows with complete data (ie non NA)
reduced_df <- reduced_df[complete.cases(reduced_df),]

# We want only columns with result information in them; i.e. no columns with parameters
outcomes <- reduced_df[,c(1:4)]
#outcomes <- as.data.frame(lapply(outcomes, as.numeric))
dim(outcomes)

X_df <- reduced_df[,c(1,6:ncol(reduced_df))] 
dim(X_df)


means <- outcomes %>% group_by(linear_comb_names) %>% summarise(across(everything(), list(mean)))
means <- as.data.frame(means)
row.names(means) <- means$linear_comb_names
means$linear_comb_names <- NULL

df <- scale(means)
# Default plot
heatmap(df, scale = "none", cexCol = 1)

#What about a violation plot?
library(ggplot2)

# Most basic violin chart
p <- ggplot(outcomes, aes(x=linear_comb_names, y=QL_adj, fill=linear_comb_names)) + 
  # fill=name allow to automatically dedicate a color for each group
  scale_y_continuous(limits = quantile(outcomes$QL_adj, c(0.0, 0.7))) +
  geom_violin()

p
## end of violin

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green

# creates a 5 x 5 inch image
png("../images/heatmaps_in_r.png",    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(output_n_sim_1,
          cellnote = output_n_sim_1,  # same data set for cell labels
          main = "Quasi-Loss", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device