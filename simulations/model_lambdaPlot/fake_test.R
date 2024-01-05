library(tidyverse)
library(rstan)
library(posterior)
library(ggdist)
  
# fake data
n = 1000

tibble(
  temp = runif(n, -3, 20)
) |>
bind_cols(
  beta = 0.1,
  Int = 0.2,
  sigma_int = -2
  beta_sigma = 0.15
) |>
mutate(
  err = rnorm(n(), 0, sigma),
  db_origin = ifelse(temp < 4, 1, 2),
  offset = ifelse(origin == 'qc', 0.2, 0),
  sigma = exp(sigma_int + beta_sigma * temp),
  y = rnorm(n(), temp * beta + Int + offset, sigma),
) ->
dat


data_stan <- with(
  dat,
  list(
    N = nrow(dat),
    lambda_log = y,
    MAT = temp,
    db_origin = db_origin
  )
)

out <- stan(
  data = data_stan,
  file = "simulations/model_lambdaPlot/model.stan",
  chains = 2,
  cores = 2
)


out |>
  as_draws_df() |>
  select(!c('lp__', contains('.'), contains('y_pred'))) |>
  mutate(draw = row_number()) |>
  pivot_longer(cols = !draw, names_to = 'par') ->
post_pars

post_pars |>
  pivot_wider(names_from = par, values_from = value) |>
  group_by(draw) |>
  expand_grid(temp = seq(-3, 20, 0.5)) |>
  mutate(
    lambda_fia = rnorm(n(), beta_mean * temp + `Int[1]`, exp(`sigma_int[1]` + beta_sigma * temp)),
    lambda_qc = rnorm(n(), beta_mean * temp + `Int[2]`, exp(`sigma_int[2]` + beta_sigma * temp))
  ) |>
  select(temp, contains('lambda')) |>
  pivot_longer(cols = contains('lambda')) |>
  ggplot() +
  aes(temp, value) +
  facet_wrap(~name) +
  geom_point(
    data = dat |>
      mutate(
        name = case_match(db_origin, 1 ~ 'lambda_fia', 2 ~ 'lambda_qc'),
      ),
    aes(temp, log(lambda)),
    alpha = 0.3
  ) +
  stat_lineribbon(.width = 0.9, alpha = 0.6)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate the effect of unbalanced explanatory variable on the variance of y

library(tidyverse)
library(rstan)
library(posterior)
library(ggdist)

tibble(
  temp = rnorm(5e3, 2, 10)
) |>
mutate(
  temp = abs(temp)
) |>
bind_cols(
  int_mean = 0.8,
  beta_mean = -0.016,
  int_sigma = -0.0002,
  beta_sigma = -0.02
) |>
mutate(
  lambda = rnorm(n(), int_mean + beta_mean * temp, exp(int_sigma + beta_sigma * temp))
) ->
fake_data

png('growth_abibal.png', width = 1200, height = 700, units = "px", res = 200) 

fake_data |>
  ggplot() + 
  aes(temp, lambda) + 
  geom_point() +
  theme_classic() + 
  labs(x = 'Mean annual temperature (°C)', y = expression(lambda)) + geom_hline(yintercept = 0, linetype = 2, color = 'red')  
dev.off()

data_stan <- with(
  fake_data,
  list(
    N = nrow(fake_data),
    lambda_log = lambda,
    MAT = temp,
    db_origin = rep(1, nrow(fake_data))
  )
)

out <- stan(
  data = data_stan,
  file = "simulations/model_lambdaPlot/model.stan",
  chains = 2,
  cores = 2
)


out |>
  as_draws_df() |>
  select(!c('lp__', contains('.'), contains('y_pred'))) |>
  mutate(draw = row_number()) |>
  pivot_longer(cols = !draw, names_to = 'par') ->
post_pars

post_pars |>
  pivot_wider(names_from = par, values_from = value) |>
  group_by(draw) |>
  expand_grid(temp = seq(min(fake_data$temp), max(fake_data), 0.5)) |>
  mutate(
    lambda = rnorm(n(), beta_mean * temp + `inter[1]`, exp(`sigma_inter[1]` + beta_sigma * temp))
  ) |>
  ggplot() +
  aes(temp, lambda) +
  geom_point(
    data = fake_data,
    aes(temp, lambda),
    alpha = 0.3
  ) +
  stat_lineribbon(.width = 0.9, alpha = 0.6)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use skewed normal distribution

library(tidyverse)
library(rstan)
library(posterior)
library(ggdist)

tibble(
  temp = rnorm(5e3, 2, 10)
) |>
mutate(
  temp = abs(temp)
) |>
bind_cols(
  int_mean = 0.014,
  beta_mean = 0.00007,
  int_sigma = -5.5,
  beta_sigma = 0,
  alpha = -7
) |>
mutate(
  lambda = sn::rsn(
    n(),
    int_mean + beta_mean * temp,
    exp(int_sigma + beta_sigma * temp),
    alpha
  )
) ->
fake_data

fake_data |>
  ggplot() + 
  aes(temp, lambda) + 
  geom_point() +
  theme_classic() + 
  labs(x = 'Mean annual temperature (°C)', y = expression(lambda)) + geom_hline(yintercept = 0, linetype = 2, color = 'red')  


data_stan <- with(
  fake_data,
  list(
    N = nrow(fake_data),
    lambda_log = lambda,
    MAT = temp,
    db_origin = rep(1, nrow(fake_data))
  )
)

out <- stan(
  data = data_stan,
  file = "simulations/model_lambdaPlot/model.stan",
  chains = 2,
  cores = 2
)


out |>
  as_draws_df() |>
  select(!c('lp__', contains('.'), contains('y_pred'))) |>
  mutate(draw = row_number()) |>
  pivot_longer(cols = !draw, names_to = 'par') ->
post_pars

post_pars |>
  pivot_wider(names_from = par, values_from = value) |>
  group_by(draw) |>
  expand_grid(temp = seq(min(fake_data$temp), max(fake_data), 0.5)) |>
  mutate(
    lambda = rnorm(n(), beta_mean * temp + `inter[1]`, exp(`sigma_inter[1]` + beta_sigma * temp))
  ) |>
  ggplot() +
  aes(temp, lambda) +
  geom_point(
    data = fake_data,
    aes(temp, lambda),
    alpha = 0.3
  ) +
  stat_lineribbon(.width = 0.9, alpha = 0.6)
