library(sn)
# https://cran.r-project.org/web/packages/sn/sn.pdf

xi <- c(-4, 1)
Omega <- matrix(c(1,0,0,.83), byrow = T, ncol = 2)
Omega[2,1] <- Omega[1,2] <- -0.8
alpha <- c(-3.5,4)
rnd <- rmst(1000,  xi, Omega, alpha, 6)

plot(rnd[,1])
plot(density(rnd[,1]))
mean(rnd[,1])

plot(rnd[,2])
plot(density(rnd[,2]))
mean(rnd[,2])
min(rnd[,2])

cor.test(rnd[,1], rnd[,2])
