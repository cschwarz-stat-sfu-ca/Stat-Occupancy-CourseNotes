library(rstan)

setwd("C:/")

##Model 44

source('Formatting Data_GBBB44.R')

# takes number of hours. Repeat until get non-zero initials or edit inits.
model44 <- stan('Stan_model48.stan', data = data, pars = params, chains = 1, init = inits,
                iter = 2000, warmup = 1000, thin = 1)

# checking convergence
max(summary(model48)[[1]][, 'Rhat'])

post.lik <- extract(model48, 'll')

# WAIC
-2 * sum(log(apply(exp(post.lik[[1]]), 2, mean))) +
  2 * sum(apply(post.lik[[1]], 2, var))

# summarizing slope coefficients
post.beta <- extract(model48, 'beta')

rbind(apply(post.beta[[1]], 2, quantile, probs = 0.025),
      apply(post.beta[[1]], 2, mean),
      apply(post.beta[[1]], 2, quantile, probs = 0.975))
