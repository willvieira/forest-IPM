#' Get Species-Specific Posterior Parameters
#'
#' Load posterior parameter samples for a given species from the local RDS files.
#'
#' @param sp Character. Species ID (e.g. "ABIBAL").
#' @param method Character. Either \code{"mean"} (posterior mean) or
#'   \code{"random"} (a single random posterior draw).
#' @param model Character of length 1 or 3. Model complexity for
#'   growth, mortality, and recruitment respectively. Accepted values:
#'   \code{"intcpt"}, \code{"intcpt_plot"}, \code{"intcpt_plot_comp"},
#'   \code{"intcpt_plot_comp_clim"}. Default is the full model.
#' @param path Character. Path to the directory containing the processed
#'   posterior RDS files. Default is \code{file.path("data", "output_sim_processed")}.
#'
#' @return A named list with elements \code{growth}, \code{mort}, \code{rec},
#'   and \code{sizeIngrowth}, each a named numeric vector of posterior parameter
#'   values.
#' @export
getPars_sp <- function(
  sp,
  method,
  model = 'intcpt_plot_comp_clim',
  path = file.path('data', 'output_sim_processed')
){
  # deal for `model` parameter
  if(length(model) == 1) {
    model = rep(model, 3)
  }else if(length(model) != 3) {
    stop('`method` parameter must be a character of length 1 or 3.')
  }

  for(MD in model)
    if(!MD %in% c('intcpt', 'intcpt_plot', 'intcpt_plot_comp', 'intcpt_plot_comp_clim')) {
      stop(
        paste0('model parameter should be one of the following character:\n', paste0(paste0('-', c('intcpt', 'intcpt_plot', 'intcpt_plot_comp', 'intcpt_plot_comp_clim')), collapse='\n'))
      )
    }
  
  # define full path to load parameters for each vital rate
  pars_dir <- setNames(
    paste0(
      path, '/', c('growth', 'mort', 'recruit', 'sizeIngrowth'), '/', c(model, 'time_truc'), '/posterior_pop_', sp, '.RDS'
    ),
    c('growth', 'mort', 'rec', 'sizeIngrowth')
  )

  # export pars
  if(method == 'mean') {
    pars_dir |>
      map(
        ~ readRDS(.x) |>
          filter(.data$par != 'lp__') |>
          group_by(.data$par) |>
          reframe(value = mean(.data$value)) |>
          pivot_wider(names_from = 'par') |>
          unlist(use.names = TRUE)
      )
  }else if(method == 'random') {
    pars_dir |>
      map(
        ~ readRDS(.x) |>
          filter(.data$par != 'lp__') |>
          filter(.data$iter == sample(1:4000, 1)) |>
          select(!.data$iter) |>
          pivot_wider(names_from = 'par') |>
          unlist(use.names = TRUE)
      )
  }else{
    stop('`param_method` parameter must be either `mean` or `random`')
  }
}

# Function to tranform rowwise tibble of parameters in list format
# Which is the current way of passing parameters to the IPM
pars_to_list <- function(pars)
{
  pars |>
    select(contains('.')) |>
    pivot_longer(cols = everything()) |>
    mutate(
      vr = str_replace(.data$name, '\\..*', ''),
      par = str_replace(.data$name, paste0(.data$vr, '.'), ''),
      # some code uses 'recruit' instead of 'rec'
      vr = case_match(
        .data$vr,
        'recruit' ~ 'rec',
        .default = .data$vr
      )
    ) |>
    select(!.data$name) |>
    group_split(.data$vr) %>%
    # fuck tidyverse for not wanting to implement an argument to keep a named list
    set_names(map_chr(., ~.x$vr[1])) |>
    map(
      ~.x |>
        select(!.data$vr) |>
        pivot_wider(names_from = 'par') |>
        unlist(use.names = TRUE)
    )
}
