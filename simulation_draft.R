


\item Size of donor pool, n: n = c(10,30,90)
\item Size of $(\alpha, \beta)$, perhaps 6 pairs?
  \item Three kinds of shock pairs: vol only, level only, dual-shock.
\item $\mc{M}_{1}$ and $\mc{M}_{21}$ and $\mc{M}_{22}$ 
  \item 3 magnitudes of vol shock, 3 magnitudes of level shock
\item 5-6 signal to noise ratios
\item 2 lengths for level shock, 2 lengths for vol shock
\item number of seeds used in simulations

donor_pool_size <- c(10,30,50)
alpha <- c(.15,.3,.65, .75)
beta <- c(.15,.3,.65, .75)
vol_model <- seq(1,4,1)
level_model <- seq(1,4,1)
vol_shock_length <- c(1,2,3)
level_shock_length <- c(1,2)


gridd <- expand.grid(list(donor_pool_size
                          , alpha
                          , beta
                          , vol_model
                          , level_model
                          , vol_shock_length
                          , level_shock_length))

# Now we get rid of arrangements where alpha + beta >= 1
gridd_subset <- gridd[gridd$Var2 + gridd$Var3 < 1,]

# Get rid of arrangement with no vol shock and no level shock
gridd_subset <- gridd_subset[gridd_subset$Var4 < 4 || gridd_subset$Var5 < 4,]


nrow(gridd_subset)
head(gridd_subset)

