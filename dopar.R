
# https://privefl.github.io/blog/a-guide-to-parallelism-in-r/

# https://cran.r-project.org/web/packages/foreach/vignettes/foreach.html

# https://www.r-bloggers.com/2013/06/grid-search-for-free-parameters-with-parallel-computing/

library(doParallel)

cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)

foreach(i = 1:3, .combine = 'c') %dopar% {
  sqrt(i)
}
