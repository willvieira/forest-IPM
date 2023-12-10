# Compute individual size variability due to parameter uncertainty

library(tidyverse)

sim_path <- file.path('simulations', 'lambda_plot')

# First step is to select two species with high and low paramater variability
# I will use the 500 parameters draw from the sensitivity_analysis v3 simulation
# load parameters and output


readRDS(file.path('simulations', 'sensAnalysis_V3', 'output_complete.RDS')) |>
  mutate(sim = 1:n()) |>
  left_join(
    read_csv(file.path('simulations', 'sensAnalysis_V3', 'simulation_pars.csv')) |>
      mutate(sim = 1:n()) |>
      select(species_id, comp, clim, replication, sim)
  ) |>
  select(!sim) ->
lambdas

# plot Coefficient of variation to grep two species
lambdas |>
  select(c(species_id, contains('.'))) |>
  pivot_longer(
    cols = contains('.'),
    names_to = 'par'
  ) |>
  group_by(species_id, par) |>
  reframe(
    cv = abs(sd(value)/mean(value))
  ) |>
  group_by(species_id) |>
  reframe(
    cv_total = sum(cv),
    cv_median = median(cv),
    cv_sd = sd(cv)
  ) |>
  ggplot() +
  aes(cv_total, cv_median, label = species_id) +
  geom_point() +
  ggrepel::geom_text_repel(size = rel(1.8))

# The extreme species are "18032ABIBAL" and "27821NYSSYL"
# I will use the posterior distribution of 100 draws to compute how large
# and small individuals vary in size with the paramaters

spIds <- c('18032ABIBAL', '27821NYSSYL')

# load IPM functions
source('R/vital_rates.R')
source('R/kernel.R')
source('R/matrix_image.R')
source('R/BasalArea_competition.R')


# Prepare parameters for both species
lambdas |>
  filter(species_id %in% spIds) |>
  filter(comp == 'low' & clim == 'center') |>
  select(species_id, contains('.')) |>
  group_by(species_id) |>
  slice_sample(n = 100) |> 
  group_split() ->
pars

# function to tranform rowwise tibble of parameters in list format
pars_to_list <- function(dat) {
    dat |>
      select(contains('.')) |>
      pivot_longer(cols = everything()) |>
      mutate(
        vr = str_replace(name, '\\..*', ''),
        par = str_replace(name, paste0(vr, '.'), '')
      ) |>
      select(!name) |>
      group_split(vr) %>%
      # fuck tidyverse for not wanting to implement an argument to keep a named list
      set_names(map_chr(., ~.x$vr[1])) |>
      map(
        ~.x |>
          select(!vr) |>
          pivot_wider(names_from = par) |>
          as_vector()
      )
}

# Get parameters to calculate initial state
pars |>
  map(~.x |> filter(growth.Lmax == max(growth.Lmax))) |>
  map(pars_to_list) ->
pars_init

N_init = list()  
N_init[[1]] = init_pop(pars_init[[1]], L = 127, h = 1, N = 0)
N_init[[2]] = init_pop(pars_init[[2]], L = 127, h = 1, N = 0)

size = list()
size[[1]] = sample(127:pars_init[[1]]$growth['Lmax'], 10)
size[[2]] = sample(127:pars_init[[2]]$growth['Lmax'], 20)

out_dt = tibble(); count = 1
for(sp in 2) {

  for(N in c(5, 10, 15))
  {
    
    dbh_size <- sample(127:pars_init[[1]]$growth['Lmax'], N)

    for(Sigma in seq(5, 12, 1))
    {
      for(i in 1:100)
      {
        N_sp = rnorm(N, dbh_size, Sigma)

        N_pop = dbh_to_sizeDist(dbh = N_sp, N_init[[sp]])

        # prepare parameters
        pars_sp_i = pars_to_list(pars[[sp]][i, ])
        
        
        K = mkKernel(
          Nvec_intra = N_pop,
          Nvec_inter = N_init[[sp]],
          delta_time = 1,
          plotSize = 300,
          Temp = pars_sp_i|>map(~.x['optimal_temp']) |>unlist()|>mean(na.rm = TRUE),
          Prec = pars_sp_i|>map(~.x['optimal_prec']) |>unlist()|>mean(na.rm = TRUE),
          pars = pars_sp_i,
          plot_random = rep(0, 3)
        )
        out_dt <- rbind(
          out_dt,
          tibble(
            Sp = sp,
            n = N,
            sigma = Sigma,
            rep = i,
            BA = size_to_BAplot(N_pop, plot_size = 300),
            lambda = max(Re(eigen(K$K)$values))
          )
        )
        cat('   Progress ', round(count/(2*3*8*100) * 100, 0), '%\r')
        count = count + 1
      }
    }
  }
}

out_dt |>
  ggplot() +
  aes(factor(sigma), log(lambda)) +
  facet_wrap(Sp~n, scales = 'free') +
  ggdist::stat_pointinterval()

out_dt |>
  ggplot() +
  aes(factor(sigma), log(BA)) +
  facet_wrap(Sp~n, scales = 'free') +
  ggdist::stat_pointinterval()

out_dt |>
  group_by(Sp, n, sigma) |>
  reframe(
    cv_lambda = sd(lambda)/mean(lambda),
    cv_BA = sd(BA)/mean(BA)
  ) |>
  ggplot() +
  aes(sigma, cv_lambda, color = factor(n)) +
  facet_wrap(~Sp, scales = 'free') +
  geom_point() +
  geom_smooth()

out_dt |>
  group_by(Sp, n, sigma) |>
  reframe(
    cv_lambda = sd(lambda)/mean(lambda),
    cv_BA = sd(BA)/mean(BA)
  ) |>
  ggplot() +
  aes(sigma, cv_BA, color = factor(n)) +
  facet_wrap(~Sp, scales = 'free') +
  geom_point() +
  geom_smooth()

out_dt |>
  ggplot() +
  aes(BA, factor(sigma), fill = factor(Sp)) +
  facet_wrap(~n, scales = 'free') +
  ggridges::geom_density_ridges2()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run IPM to equilibrium and quantify SD of BAeq from the 100 post draws

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


time_max = 200

# Get parameters to calculate initial state
pars |>
  map(~.x |> filter(growth.Lmax == max(growth.Lmax))) |>
  map(pars_to_list) ->
pars_init

# Create initial condition (small pop) and het competition (null)
N_init = list()  
N_init[[1]] = init_pop(pars_init[[1]], L = 127, h = 1, N = 1)
N_init[[2]] = init_pop(pars_init[[2]], L = 127, h = 1, N = 1)
N_comp <- N_init
N_comp[[1]]$Nvec = rep(0, length(N_comp[[1]]$Nvec))
N_comp[[2]]$Nvec = rep(0, length(N_comp[[2]]$Nvec))

mat = list()
for(Sp in 1:2) {

  sp_mat = list()
  for(rep in 1:100) {
    
    N_sp <- N_init[[Sp]]

    pars_sp <- pars_to_list(pars[[Sp]][rep, ])

    # generate initial Kernel
    K0 <- mkKernel(
      Nvec_intra = N_init[[Sp]],
      Nvec_inter = N_comp[[Sp]],
      delta_time = 1,
      plotSize = 300,
      Temp = pars_sp_i|>map(~.x['optimal_temp']) |>unlist()|>mean(na.rm = TRUE),
      Prec = pars_sp_i|>map(~.x['optimal_prec']) |>unlist()|>mean(na.rm = TRUE),
      pars = pars_sp,
      plot_random = rep(0, 3)
    )

    # matrix to save size distribution
    ntmat = matrix(0, nrow = length(N_sp$Nvec), ncol = time_max)
    ntmat[, 1] <- N_sp$Nvec
    nochange = 0
    for(Time in 2:time_max)
    {
      # update the state
      ntmat[, Time] <- K0$K %*% ntmat[, Time - 1]
      N_sp$Nvec <- ntmat[, Time]

      # update the kernel
      K0 <- mkKernel(
        Nvec_intra = N_sp,
        Nvec_inter = N_comp[[Sp]],
        delta_time = 1,
        plotSize = 300,
        Temp = pars_sp_i|>map(~.x['optimal_temp']) |>unlist()|>mean(na.rm = TRUE),
        Prec = pars_sp_i|>map(~.x['optimal_prec']) |>unlist()|>mean(na.rm = TRUE),
        pars = pars_sp,
        plot_random = rep(0, 3)
      )

      cat(' Time step', Time, 'of', time_max, '(', round(Time/time_max * 100, 1), '%)\r')
      
    }

    sp_mat[[rep]] <- ntmat
  }

  mat[[Sp]] <- sp_mat
}

get_ba = function(Npop, value, plot_size) {
  Npop$Nvec = value
  size_to_BAplot(Npop, plot_size = 300)
}






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run IPM to equilibrium and quantify SD of BAeq from the 100 post draws

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get parameters to calculate initial state
pars |>
  map(~.x |> filter(growth.Lmax == max(growth.Lmax))) |>
  map(pars_to_list) ->
pars_init

N_init = list()  
N_init[[1]] = init_pop(pars_init[[1]], L = 127, h = 1, N = 0)
N_init[[2]] = init_pop(pars_init[[2]], L = 127, h = 1, N = 0)

size = list()
size[[1]] = sample(127:pars_init[[1]]$growth['Lmax'], 10)
size[[2]] = sample(127:pars_init[[2]]$growth['Lmax'], 20)

out_dt = tibble(); count = 1
for(sp in 1:2) {

  for(N in c(5, 10, 15))
  {
    
    dbh_size <- sample(127:pars_init[[1]]$growth['Lmax'], N)

    for(Sigma in seq(5, 35, 5))
    {
      for(i in 1:100)
      {
        N_sp = rnorm(N, dbh_size, Sigma)

        N_pop = dbh_to_sizeDist(dbh = N_sp, N_init[[sp]])

        # prepare parameters
        # pars_sp_i = pars_to_list(pars[[sp]][i, ])
        
        
        # K = mkKernel(
        #   Nvec_intra = N_pop,
        #   Nvec_inter = N_init[[sp]],
        #   delta_time = 1,
        #   plotSize = 300,
        #   Temp = pars_sp_i|>map(~.x['optimal_temp']) |>unlist()|>mean(na.rm = TRUE),
        #   Prec = pars_sp_i|>map(~.x['optimal_prec']) |>unlist()|>mean(na.rm = TRUE),
        #   pars = pars_sp_i,
        #   plot_random = rep(0, 3)
        # )
        out_dt <- rbind(
          out_dt,
          tibble(
            Sp = sp,
            n = N,
            sigma = Sigma,
            rep = i,
            BA = size_to_BAplot(N_pop, plot_size = 300)#,
            # lambda = max(Re(eigen(K$K)$values))
          )
        )
        cat('   Progress ', round(count/(2*3*8*100) * 100, 0), '%\r')
        count = count + 1
      }
    }
  }
}

out_dt |>
  ggplot() +
  aes(factor(sigma), BA) +
  facet_wrap(Sp~n, scales = 'free') +
  ggdist::stat_pointinterval()

out_dt |>
  group_by(Sp, n, sigma) |>
  reframe(
    cv_BA = sd(BA)/mean(BA)
  ) |>
  ggplot() +
  aes(sigma, cv_BA, color = factor(n)) +
  facet_wrap(~Sp, scales = 'free') +
  geom_point() +
  geom_smooth()

out_dt |>
  ggplot() +
  aes(BA, factor(sigma), fill = factor(Sp)) +
  facet_wrap(~n, scales = 'free') +
  ggridges::geom_density_ridges2()
    