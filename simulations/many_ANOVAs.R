# https://stackoverflow.com/questions/51937380/fast-post-hoc-computation-using-r

# https://stackoverflow.com/questions/52124359/anova-on-r-with-different-dependent-variables


unique_count_df <- apply(output_n_sim_1, 2, function(x) length(unique(x)))

columns_we_want <- names(unique_count_df)[unique_count_df != 1]

reduced_df <- output_n_sim_1[names(output_n_sim_1) %in% columns_we_want]

reduced_df <- as.data.frame(lapply(reduced_df, as.numeric))

outcomes <- reduced_df[,c(6:ncol(reduced_df))]
outcomes <- as.data.frame(lapply(outcomes, as.numeric))
dim(outcomes)

X_df <- reduced_df[,c(1:5)] 
dim(X_df)

res <- aov(cbind(outcomes[,1], outcomes[,12]) ~. , data = X_df )
summary(res)

#Now let's do all outcome columns
res <- aov( as.matrix(outcomes) ~.^2 , data = X_df )
summary(res)

#Just look at means
apply(outcomes, 2, mean, na.rm=TRUE)
