# https://stackoverflow.com/questions/51937380/fast-post-hoc-computation-using-r

# https://stackoverflow.com/questions/52124359/anova-on-r-with-different-dependent-variables

options(digits = 3)

load("/home/david/Desktop/synthetic_vol_forecasting/simulation_results/simulation_results/simcount_1_savetime_Wed Dec 28 01:00:21 2022.Rdata")

library(lmboot)

#We inspect the data
View(output)

# We look for the columns that have no variation, i.e. those with only 1 value.  
unique_count_df <- apply(output, 2, function(x) length(unique(x)))

# We drop these columns.
columns_we_want <- names(unique_count_df)[unique_count_df != 1]
reduced_df <- output[names(output) %in% columns_we_want]

# Make the dataframe numeric
reduced_df <- as.data.frame(lapply(reduced_df, as.numeric))

#Make vol model into factor
reduced_df$vol_model <- as.factor(reduced_df$vol_model)

# Take only rows with complete data (ie non NA)
reduced_df <- reduced_df[complete.cases(reduced_df),]

# We want only columns with result information in them; i.e. no columns with parameters
outcomes <- reduced_df[,c(2:5)]
outcomes <- as.data.frame(lapply(outcomes, as.numeric))
dim(outcomes)

X_df <- reduced_df[,c(1,6:ncol(reduced_df))] 
dim(X_df)

X_df[,duplicated(cor(X_df))]


# https://cran.r-project.org/web/packages/lmboot/lmboot.pdf

reduced_df <- reduced_df[complete.cases(cbind(reduced_df$QL_adj1, reduced_df[,c(1:10)])),]

myANOVA2 <- ANOVA.boot(reduced_df$QL_adj1 ~. , data=reduced_df[,c(1:10)], B = 40000)
names(myANOVA2)
myANOVA2$terms
myANOVA2$`p-values` #bootstrap p-values for 2-way interactions model

ANOVA_df <- matrix(c(myANOVA2$terms[-1], as.numeric(myANOVA2$`p-values`) ), ncol = 2)
ANOVA_df <- ANOVA_df[order(ANOVA_df[,2], decreasing = FALSE),]
ANOVA_df

res <- lm(reduced_df$QL_adj1 ~. , data = reduced_df[,c(1:10)] )
summary(res)

par(mfrow=c(2,2))
plot(res)
par(mfrow=c(1,1))

res <- aov(cbind(outcomes[,1], outcomes[,12]) ~. , data = X_df )
summary(res)



#Now let's do all outcome columns
res <- aov( as.matrix(outcomes) ~.^2 , data = X_df )
summary(res)

#Kruskal Wallis test
#kruskal.test( as.matrix(outcomes[,1]) ~. , data = X_df ) 
#https://stackoverflow.com/questions/35091164/r-kruskal-wallis-with-multiple-factors

#Just look at means
applied <- apply(outcomes, 2, mean, na.rm=TRUE)
df_applied <- as.data.frame(applied, ncol = 1)
df_applied

# Look at shapiro wilkes
ws <- apply(outcomes, 2, shapiro.test)
ws

