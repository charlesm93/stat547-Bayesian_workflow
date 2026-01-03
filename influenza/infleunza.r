
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
mod <- cmdstan_model("model/sir_poisson.stan")

n_chains <- 4
fit <- mod$sample(data = data_sir,
                  chains = n_chains,
                  parallel_chains = n_chains,
                  init = init,
                  save_warmup = TRUE)

#####################################################################
## Check the inference
pars <- c("gamma", "beta", "R0", "T")
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
## Run same model with a Poisson likelihood

mod <- cmdstan_model("model/sir_poisson.stan")

fit_poisson <- mod$sample(data = data_sir,
                          chains = n_chains,
                          parallel_chains = n_chains,
                          init = init,
                          save_warmup = TRUE)

fit_poisson$summary(variables = pars)

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


log_lik_draws_poisson <- fit_poisson$draws("log_lik")
loo_estimate_poisson <-
  loo(log_lik_draws_poisson, r_eff = relative_eff(log_lik_draws_poisson))

print(loo_estimate_poisson)
print(loo_estimate)


#####################################################################
## run SIR model with custom tolerance for the ODE solver

mod <- cmdstan_model("model/sir_tol.stan")

tol <- 1e-2
data_sir$tol <- tol
fit_tol <- mod$sample(data = data_sir,
                      chains = n_chains,
                      parallel_chains = n_chains,
                      init = init,
                      save_warmup = TRUE)

fit_tol$time()

log_lik <- fit_tol$draws("log_lik")

fit_tol$summary(variables = pars)
loo_estimate <- loo(log_lik, r_eff = relative_eff(log_lik))
print(loo_estimate)

log_ratios <- fit_tol$draws("log_ratios")

psis_fit <- psis(log_ratios, r_eff =  relative_eff(log_ratios))
psis_fit

# Correct Monte Carlo samplers, using importance weights.
# Only works if the log ratios don't go to 0!
is_summary(fit_tol, pars, psis_fit, log_ratios)



###############################################################################
## Additional run for MCMC handbook chapter
library(latex2exp)

# tolerance is validated via IS, see above.
tol <- 1e-2
data_sir$tol <- tol

seed_list <- c(1954, 1974, 1990, 1998, 2014)
for (seed_x in seed_list) {
  fit_tol <- mod$sample(data = data_sir, chains = n_chains,
                        parallel_chains = n_chains,
                        init = init, save_warmup = TRUE,
                        seed = seed_x)
  
  fit_tol$save_object(paste0("deliv/sir_seed", seed_x, ".fit.RDS"))
} 

pars <- c("recovery_time")
saved_results <- matrix(NA, nrow = length(seed_list), ncol = 4)
saved_results_sub <- matrix(NA, nrow = length(seed_list), ncol = 4)
saved_diagnostics <- matrix(NA, nrow = length(seed_list), ncol = 2)
saved_diagnostics_sub <- matrix(NA, nrow = length(seed_list), ncol = 2)

n_sub_sample <- 10

for (i in 1:length(seed_list)) {
  fit_saved <- readRDS(paste0("deliv/sir_seed", seed_list[i], ".fit.RDS"))
  fit_summary <- fit_saved$summary(variables = pars)
  saved_results[i, ] <- c(fit_summary$median,
                          fit_summary$q5,
                          fit_summary$q95,
                          fit_summary$ess_bulk)
  
  saved_diagnostics[i, ] <- c(fit_summary$rhat,
                              fit_summary$ess_bulk)
  
  draws_sub <- fit_saved$draws(variables = pars[1])[1:n_sub_sample, , ]
  fit_summary <- summarise_draws(draws_sub)
  saved_results_sub[i, ] <- c(fit_summary$median,
                              fit_summary$q5,
                              fit_summary$q95,
                              fit_summary$ess_bulk)
  
  saved_diagnostics_sub[i, ] <- c(fit_summary$rhat,
                                  fit_summary$ess_bulk)
}

saved_results
saved_results_sub

saved_diagnostics
saved_diagnostics_sub

###############################################################################
## plot diagnostics

# Let's retweak the posterior function a little bit.

mcmc_rhat_cust <- function(data, ..., size = NULL) {
  graph <- ggplot(data = data,
                  mapping = aes(x = .data$value, 
                                y = .data$parameter,
                                color = .data$rating)) + 
    geom_segment(mapping = aes(yend = .data$parameter,
                               xend = ifelse(min(.data$value) < 1, 1, -Inf)),
                 na.rm = TRUE, linewidth = 1) +
    geom_point(size = 5) +
    scale_color_manual(values = c("orange", "lightblue"),
                       labels = unname(TeX(c("$\\hat{R} > 1.01", "$\\hat{R} < 1.01")))) +
    bayesplot_theme_get()
  if (min(data$value) < 1) {
    graph <- graph + vline_at(1, color = "gray", linewidth = 1)
  }

  brks <- c(1.01)
  graph <- graph + vline_at(brks, color = "gray", linetype = 2, linewidth = 0.5) +
    # scale_color_manual(values = c("orange", "light blue")) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ylab(" ") +
    labs(x = expression(hat(R)))

  graph

  # graph + diagnostic_points(size)
  # brks <- set_rhat_breaks(data$value)
  # graph + diagnostic_points(size) + vline_at(brks[-1], color = "gray", 
  #                                            linetype = 2, linewidth = 0.25) + labs(y = NULL, x = expression(hat(R))) + 
  #   scale_fill_diagnostic("rhat") + scale_color_diagnostic("rhat") + 
  #   scale_x_continuous(breaks = brks, expand = c(0, 0.01)) + 
  #   scale_y_discrete(expand = c(0.025, 0)) + yaxis_title(FALSE) + 
  #   yaxis_text(FALSE) + yaxis_ticks(FALSE)
}

mcmc_rhat_cust(data_rhat)

saved_rhat_all <- c(saved_diagnostics_sub[, 1], saved_diagnostics[, 1])
saved_rhat_exp <- c(rep("N = 40", nrow(saved_diagnostics)),
                    rep("N = 4000", nrow(saved_diagnostics)))
saved_rhat_exp <- saved_rhat_exp[order(saved_rhat_all)]

data_rhat <- mcmc_rhat_data(c(saved_diagnostics_sub[, 1], saved_diagnostics[, 1]))
for (i in 1:nrow(data_rhat)) {
  if (data_rhat$value[i] <= 1.01) {
    data_rhat$rating[i] <- factor("low")
    data_rhat$description[i] <- "hat(R) <= 1.01"
  } else if (data_rhat$value[i] > 1.01) {
    data_rhat$rating[i] <- factor("high")
    data_rhat$description[i] <- "hat(R) > 1.01"    
  } else {
    data_rhat$rating[i] <- factor("high")
    data_rhat$description[i] <- "hat(R) > 1.05" 
  }
}
data_rhat$rating <- factor(data_rhat$rating, levels = c("high", "low")) 
data_rhat$experiment <- factor(saved_rhat_exp)

# a bit of a hack but saves me time
data_rhat$parameter <- c(1:2, 1:5, 3:5)

mcmc_rhat_cust(data_rhat) + facet_wrap(~experiment) +
  theme(legend.position = c(0.9, 0.25)) +
  theme(text = element_text(size = 20)) +
  labs(color="")

