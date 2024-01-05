

md_out = read_cmdstan_csv(dir('test', full.names = T)) 

md_out$post_warmup_draws |>
  as_draws_df() |>
  select(!lp__) |>
  pivot_longer(
    cols = c(contains('_'), inter, alpha),
    names_to = 'par'
  ) ->
pars

pars |>
  ggplot() +
  aes(.iteration, value) +
  aes(color = factor(.chain)) +
  facet_wrap(~par, scales = 'free') +
  geom_line()

pars |>
  # filter(.chain == 2) |>
  pivot_wider(
    names_from = 'par'
  ) |>
  group_by(.draw) |>
  expand_grid(
    temp = seq(min(data_i$bio_01_mean), max(data_i$bio_01_mean), length.out=200)
  ) |>
  mutate(
    Mean = temp * beta_mean + inter,
    sigma = exp(temp * beta_sigma + sigma_inter),
    lambda = sn::rsn(n = n(), xi = Mean, omega = sigma, alpha = alpha)
  ) |>
  ggplot() +
  aes(temp, lambda) +
  geom_point(
    data = data_stan |> as.data.frame(),
    aes(MAT, lambda_log),
    alpha = 0.5, size = 0.5
  ) +
  ggdist::stat_lineribbon(.width = 0.9, alpha = 0.7) +
  theme_classic() +
  geom_hline(yintercept = 0)

