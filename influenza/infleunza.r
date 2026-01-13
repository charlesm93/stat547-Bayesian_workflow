
rm(list = ls())
gc()
set.seed(1954)

# adjust to your setting
.libPaths("~/Rlib/")
setwd("~/Desktop/Teaching/UBC/STAT 547/course_code/influenza")

library(outbreaks)
library(ggplot2)
theme_set(theme_bw())
library(rjson)
library(cmdstanr)

library(bayesplot)
library(posterior)
library(loo)
library(tidyverse)

source("tools_is.r")

#####################################################################
## Prep data

# plot the data, obtained from outbreaks
ggplot(data = influenza_england_1978_school) + 
  geom_point(mapping = aes(x = date, y = in_bed)) + 
  labs(y = "Number of students in bed")

# create a data list to be passed to Stan
cases <- influenza_england_1978_school$in_bed
N <- 763;
n_days <- length(cases) 
t <- seq(0, n_days, by = 1)
t0 = 0 
t <- t[-1]

#initial conditions
i0 <- 1
s0 <- N - i0
r0 <- 0
y0 <- c(s0, i0, r0)

data_sir <- list(n_days = n_days,
                 y0 = y0,
                 t0 = t0,
                 ts = t,
                 N = N,
                 cases = cases)

# optional: saved process data set
write(toJSON(data_sir), "data/data_sir.json")

#####################################################################
## Run and fit Stan

# define starting distribution based on prior
init <- function() {
  list(beta = abs(rnorm(1, mean = 2, sd = 1)),
       gamma = abs(rnorm(1, mean = 0.4, sd = 0.5)),
       phi_inv = rexp(1, rate = 5))
}

# transpile (translate Stan to C++ and then compile)
mod <- cmdstan_model("model/sir_demo.stan")

n_chains <- 4
fit <- mod$sample(data = data_sir,
                  chains = n_chains,
                  parallel_chains = n_chains,
                  init = init,
                  save_warmup = TRUE)

#####################################################################
## Check the inference
pars <- c("gamma", "beta", "T", "R0")
fit$summary(variables = pars)

bayesplot::mcmc_trace(fit$draws(inc_warmup = TRUE),
                      n_warmup = 1000, pars = pars)
bayesplot::mcmc_dens_overlay(fit$draws(), pars = pars)



# Extract posterior predictive checks
pred_cases <- as.matrix(
  as_draws_df(fit$draws(variables = c("pred_cases"))))[, -(15:17)]

bayesplot::ppc_ribbon(y = data_sir$cases, yrep = pred_cases, 
                      x = data_sir$ts, y_draw = "point") + 
  theme_bw() +
  ylab("cases") + xlab("days")

#####################################################################
## Run same model with a negative Binomial likelihood

mod <- cmdstan_model("model/sir_nb.stan")

fit_nb <- mod$sample(data = data_sir,
                     chains = n_chains,
                     parallel_chains = n_chains,
                     init = init,
                     save_warmup = TRUE)

pars <- c("gamma", "beta", "phi_inv", "T", "R0")
fit_nb$summary(variables = pars)

pred_cases_poisson <- as.matrix(
  as_draws_df(fit_poisson$draws(variables = c("pred_cases"))))[, -(15:17)]

bayesplot::ppc_ribbon(y = data_sir$cases, yrep = pred_cases_poisson, 
                      x = data_sir$ts, y_draw = "point") + 
  theme_bw() +
  ylab("cases") + xlab("days")

#####################################################################
# compute PSIS-loo estimate

log_lik_draws <- fit$draws("log_lik")
loo_estimate <- loo(log_lik_draws, r_eff = relative_eff(log_lik_draws))

print(loo_estimate)
print(loo_estimate$diagnostics)


log_lik_draws_nb <- fit_nb$draws("log_lik")
loo_estimate_nb <-
  loo(log_lik_draws_nb, r_eff = relative_eff(log_lik_draws_nb))

print(loo_estimate_nb)
