


df <- as.data.frame(matrix(c('red', 'red','blue',1,1,1,'tiger','bear','elephant'),nrow = 3, byrow = F))
df

unique_count_df <- apply(df, 2, function(x) length(unique(x)))

columns_we_want <- names(unique_count_df)[unique_count_df != 1]

df[names(df) %in% columns_we_want]
