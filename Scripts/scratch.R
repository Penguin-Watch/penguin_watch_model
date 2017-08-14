#            mean     2.5%    50%   97.5%    Rhat
# beta_p     0.0065  0.0052 0.0065 0.0077 1.0001
# beta_phi   0.0017 -0.0041 0.0016 0.0076 1.0001
# mean_p     0.6180  0.5704 0.6184 0.6634 1.0003
# mean_phi   0.9922  0.9870 0.9924 0.9963 1.0020
# mu_p       0.4823  0.2837 0.4827 0.6786 1.0003
# mu_phi     4.8978  4.3322 4.8755 5.5910 1.0019
# sigma_p2   0.2091  0.1066 0.1969 0.3794 1.0000
# sigma_phi2 0.3890  0.0002 0.2790 1.4682 1.0089


require(boot)

cov <- 1:199
#mu_phi + beta_phi + eps_phi
inv.logit(4.8755 + 0.0016*cov - 0.3)


#mu_p + beta_p + eps_p
inv.logit(0.4827 + 0.0065*cov)



#priors for params

#beta_phi and beta_p
a <- rnorm(100000, 0, 10)
to.rm <- which(a > 1 | a < -1)
b <- a[-to.rm]
c <- b[1:5000] #first 5k
plot(density(c))

#mean_phi and mean_p
a2 <- runif(5000, 0, 1)
plot(density(a2))

#eps_phi and eps_p
a3 <- rnorm(10000, 0, 0.34)
hist(a3)
to.rm3 <- which(a3 > 20 | a3 < -20)
b3 <- a3[-to.rm3]
c3 <- b3[1:5000]
plot(density(b3))

