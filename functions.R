
#Used in curve fitting
logistic <- function(x) 1/(1+exp(-x))


#Measures how well a simulated curve matches the historical data
distance_f <- function(pandemic, f) {
  rel_df <- fluScenarioR::relative_incidence_df %>%
    dplyr::filter(!(Pandemic == "2009" & Date > lubridate::as_date("2009-11-15"))) %>%
    dplyr::filter(Pandemic == pandemic) |>
    dplyr::rename(y = value)# |> dplyr::mutate(Date = floor(as.numeric(Date)))
  
  function(par) {   #a function to produce a curve fitting score
    f(par) %>% dplyr::group_by(Date) |>
      dplyr::summarise(value = sum(value)) |>
      dplyr::ungroup() |>
      dplyr::mutate(value = value/max(value)) |>
      dplyr::inner_join(rel_df, by = "Date") |>
      dplyr::mutate(dist = abs(value - y)) -> df
    
    
    sum(df$dist)
  }
}

# Function that allows you to change mixing rates
contact_mask_f <- function(peak, startDate, age_limits, contact_group = "default") {
  yr <- lubridate::year(startDate)
  
  # peak is the number of peaks
  # Setup default cm. Note to use dates holiday dates etc with the same year as the given startDate
  # TODO: fix it so we can give any startDate and use holidays relative from then etc.
  if (peak == 2) {
    cw <- fluScenarioR::create_contact_mask_df(dates = c(startDate, lubridate::as_date(paste0(yr, c("-07-20", "-08-31")))),
                                               scale = c(1, -1, 1), i = c(2), j = c(2), age_limits = age_limits)
    cw_id <- which(cw$value == -1)
  } else if (peak == 3) {
    if (contact_group == "paediatric")
      df <- data.frame(i = 2, j = 2)
    else if (contact_group == "young" || contact_group == "two")
      df <- data.frame(i = c(2,3), j = c(2,3))
    else
      df <- tidyr::expand(data.frame(), i = 1:6, j = 1:6)
    
    cw <- fluScenarioR::create_contact_mask_df(dates = c(startDate, lubridate::as_date(c(paste0(yr, c("-06-29", "-09-21", "-11-02")),
                                                                                         #cw <- fluScenarioR::create_contact_mask_df(dates = c(startDate, lubridate::as_date(c(paste0(yr, c("-07-06", "-09-14", "-11-09")),
                                                                                         paste0(yr+1, "-01-10")))),
                                               scale = c(1, -1, 1, -2, 1), i = df$i, j = df$j,
                                               age_limits = age_limits)
    if (contact_group == "two")
      cw <- cw %>% dplyr::mutate(value = ifelse(lubridate::month(Date) == 11, -2, value))
    cw_id1 <- which(cw$value == -1)
    cw_id2 <- which(cw$value == -2)
  } else {
    stop("Unsupported number of peaks")
  }
  
  function(par) {
    # Change the masks as needed..
    # If it doesn't work for 3 peaks then try to reduce mixing in the whole population between the second and
    # third wave.
    if (peak == 2)
      cw$value[cw_id] <- logistic(par[1])
    else {
      cw$value[cw_id1] <- logistic(par[1])
      cw$value[cw_id2] <- logistic(par[2])
    }
    cw
  }
}

#assembles the model and it's parameters
incidence_f <- function(pandemic, peak, severity, panId, n = 1, startDate, data_lst, contact_group = "default", ...) {
  if (missing(peak)) {
    if (pandemic == "1918")
      peak <- 3
    else
      peak <- 2
  }
  if (missing(severity)) {
    if (peak == 3)
      severity <- "High"
    else
      severity <- "Low"
  }
  load("scenario_model_parameter_df.rda")
  load("pandemic_df.rda")
  if (missing(panId)) {
    panId <- 1
    if (peak == 2 && severity == "Low")
      panId <- 2
    mp_df <- pandemic_df %>% dplyr::left_join(scenario_model_parameter_df) %>%
      dplyr::filter(PanID == panId)
  } else {
    mp_df <- pandemic_df %>% dplyr::left_join(scenario_model_parameter_df) %>%
      dplyr::filter(PanID == panId)
  }
  
  epi_par <- list(
    R0 = mp_df %>% dplyr::filter(Parameter == "R0") %>% dplyr::select(value),
    latent_period = mp_df %>% dplyr::filter(Parameter == "Latent period") %>% dplyr::select(value),
    infectious_period = mp_df %>% dplyr::filter(Parameter == "Infectious period") %>% dplyr::select(value),
    susc = mp_df %>% dplyr::filter(Parameter == "Susceptibility") %>% dplyr::select(value)
  ) %>% purrr::imap(function(df, k) df$value)
  
  if (missing(data_lst))
    data_lst <- fluScenarioR::loadData("2009")
  if (!missing(pandemic)) {
    rel_df <- fluScenarioR::relative_incidence_df %>%
      dplyr::filter(!(Pandemic == "2009" & Date > lubridate::as_date("2009-11-15"))) %>%
      dplyr::filter(Pandemic == pandemic)
    data_lst$dates <- rel_df$Date
  }
  
  if (!missing(startDate)) {
    data_lst$dates <- seq(startDate, by = 7, length.out = 52)
  }
  
  data_lst$cm_df <- data_lst$cm_df %>% dplyr::mutate(value = ifelse(as.numeric(AgeGroup) >= 3 | as.numeric(AgeGroup2) >= 3, value*0.469, value))
  
  cm_f <- contact_mask_f(peak, data_lst$dates[1], data_lst$age_limits, contact_group = contact_group)
  function(par) {
    cw <- cm_f(par)
    if (peak == 2) {
      data_lst$initial_infected <- rep(1000*logistic(par[2]), length(data_lst$initial_infected))
      #epi_par$R0 <- 1.33
    } else if (peak == 3) {
      data_lst$initial_infected <- rep(1000*logistic(par[3]), length(data_lst$initial_infected))
    }
    delay <- 0
    if (n > 1) {
      delay <- logistic(par[length(par)])*lubridate::duration(4, units = "week")/n
      #qs::qsave(delay, file = "~/Downloads/delay.qs")
    }
    arg <- list(n = 1, R0 = epi_par$R0, susceptibility = epi_par$susc, data_lst = data_lst,
                infection_delays = c(epi_par$latent_period, epi_par$infectious_period),
                contact_mask = cw, ..., startPandemic = 0, include_novacc = F)
    if (!is.null(arg$duration) && arg$duration == 0)
      arg$dosesPerWeek <- NULL
    suppressMessages(suppressWarnings(capture.output(
      ode <- do.call(fluScenarioR::rModel, arg), type = c("output", "message"))))
    inc_df <- ode$incidence %>%
      dplyr::mutate(value = value/n) %>%
      dplyr::select(-Week)
    if (n > 1) {
      for (i in 2:n) {
        #lubridate::as_date(lubridate::as_date("1919-01-02") - lubridate::duration(delay, units = "week"))
        sDate <- lubridate::as_date(arg$data_lst$dates[1] + delay)
        arg$data_lst$dates <- c(sDate, arg$data_lst$dates[arg$data_lst$dates > sDate])
        if ("vacc_delay" %in% names(arg))
          arg$vacc_delay <- arg$vacc_delay - (as.numeric(delay, units = "week"))
        suppressMessages(suppressWarnings(capture.output(
          ode <- do.call(fluScenarioR::rModel, arg), type = c("output", "message"))))
        #saveRDS(ode, file = "~/Downloads/ode.rds")
        inc_df <- inc_df %>%
          dplyr::full_join(ode$incidence %>% dplyr::rename(nv = value) %>% dplyr::select(-Week),
                           by = c("np", "Date", "RiskGroup", "AgeGroup")) %>%
          dplyr::mutate(value = replace(value, is.na(value), 0), nv = replace(nv, is.na(nv), 0)) %>%
          dplyr::mutate(value = value + nv/n) %>% dplyr::select(-nv)
      }
    }
    inc_df
  }
}

incidence_region_f <- function(pandemic, peak, severity, panId, n = 1, startDate, data_lst, contact_group = "default", ...) {
  if (missing(peak)) {
    if (pandemic == "1918")
      peak <- 3
    else
      peak <- 2
  }
  if (missing(severity)) {
    if (peak == 3)
      severity <- "High"
    else
      severity <- "Low"
  }
  load("assets/scenario_model_parameter_df.rda")
  load("assets/pandemic_df.rda")
  if (missing(panId)) {
    panId <- 1
    if (peak == 2 && severity == "Low")
      panId <- 2
    mp_df <- pandemic_df %>% dplyr::left_join(scenario_model_parameter_df) %>%
      dplyr::filter(PanID == panId)
  } else {
    mp_df <- pandemic_df %>% dplyr::left_join(scenario_model_parameter_df) %>%
      dplyr::filter(PanID == panId)
  }
  
  epi_par <- list(
    R0 = mp_df %>% dplyr::filter(Parameter == "R0") %>% dplyr::select(value),
    latent_period = mp_df %>% dplyr::filter(Parameter == "Latent period") %>% dplyr::select(value),
    infectious_period = mp_df %>% dplyr::filter(Parameter == "Infectious period") %>% dplyr::select(value),
    susc = mp_df %>% dplyr::filter(Parameter == "Susceptibility") %>% dplyr::select(value)
  ) %>% purrr::imap(function(df, k) df$value)
  
  if (missing(data_lst))
    data_lst <- fluScenarioR::loadData("2009")
  if (!missing(pandemic)) {
    rel_df <- fluScenarioR::relative_incidence_df %>%
      dplyr::filter(!(Pandemic == "2009" & Date > lubridate::as_date("2009-11-15"))) %>%
      dplyr::filter(Pandemic == pandemic)
    data_lst$dates <- rel_df$Date
  }
  
  if (!missing(startDate)) {
    data_lst$dates <- seq(startDate, by = 7, length.out = 52)
  }
  
  data_lst$cm_df <- data_lst$cm_df %>% dplyr::mutate(value = ifelse(as.numeric(AgeGroup) >= 3 | as.numeric(AgeGroup2) >= 3, value*0.469, value))
  
  cm_f <- contact_mask_f(peak, data_lst$dates[1], data_lst$age_limits, contact_group = contact_group)
  function(par) {
    cw <- cm_f(par)
    if (peak == 2) {
      data_lst$initial_infected <- rep(1000*logistic(par[2]), length(data_lst$initial_infected))
      #epi_par$R0 <- 1.33
    } else if (peak == 3) {
      data_lst$initial_infected <- rep(1000*logistic(par[3]), length(data_lst$initial_infected))
    }
    delay <- 0
    if (n > 1) {
      delay <- logistic(par[length(par)])*lubridate::duration(4, units = "week")/n
    }
    arg <- list(n = 1, R0 = epi_par$R0, susceptibility = epi_par$susc, data_lst = data_lst,
                infection_delays = c(epi_par$latent_period, epi_par$infectious_period),
                contact_mask = cw, ..., startPandemic = 0,
                maxUptake = 0.9, include_novacc = F)
    if (!is.null(arg$duration) && arg$duration == 0)
      arg$dosesPerWeek <- NULL
    suppressMessages(suppressWarnings(capture.output(
      ode <- do.call(fluScenarioR::rModel, arg), type = c("output", "message"))))
    inc_df <- ode$incidence %>%
      dplyr::mutate(value = value/n, Region = 1) %>%
      dplyr::select(-Week)
    if (n > 1) {
      for (i in 2:n) {
        #lubridate::as_date(lubridate::as_date("1919-01-02") - lubridate::duration(delay, units = "week"))
        sDate <- lubridate::as_date(arg$data_lst$dates[1] + delay)
        arg$data_lst$dates <- c(sDate, arg$data_lst$dates[arg$data_lst$dates > sDate])
        if ("vacc_delay" %in% names(arg))
          arg$vacc_delay <- arg$vacc_delay - (as.numeric(delay, units = "week"))
        suppressMessages(suppressWarnings(capture.output(
          ode <- do.call(fluScenarioR::rModel, arg), type = c("output", "message"))))
        inc_df <- inc_df %>%
          dplyr::bind_rows(ode$incidence %>% dplyr::select(-Week) %>%
                             dplyr::mutate(value = value/n, Region = i))
      }
    }
    inc_df
  }
}

as_cpp_array <- function(v) {
  matrix(v, ncol = 1)
}

collect_variables_f <- function(pandemic, peak, severity, behaviour, panId, n = 1, startDate, data_lst, contact_group = "default", hospital_rates = c(), mortality_rates = c(), long = F, max_scenario = 1000, daily = F,
                                minimal_irate = 2000/56000000,...) {
  if (missing(peak)) {
    if (pandemic == "1918")
      peak <- 3
    else
      peak <- 2
  }
  if (missing(severity)) {
    if (peak == 3)
      severity <- "High"
    else
      severity <- "Low"
  }
  if (missing(behaviour)) {
    behaviour <- 1
  }
  load("scenario_model_parameter_df.rda")
  load("pandemic_df.rda")
  if (missing(panId)) {
    panId <- 1
    if (peak == 2 && severity == "Low")
      panId <- 2
    mp_df <- pandemic_df %>% dplyr::left_join(scenario_model_parameter_df) %>%
      dplyr::filter(PanID == panId)
  } else {
    mp_df <- pandemic_df %>% dplyr::left_join(scenario_model_parameter_df) %>%
      dplyr::filter(PanID == panId)
  }
  
  epi_par <- list(
    R0 = mp_df %>% dplyr::filter(Parameter == "R0") %>% dplyr::select(value),
    latent_period = mp_df %>% dplyr::filter(Parameter == "Latent period") %>% dplyr::select(value),
    infectious_period = mp_df %>% dplyr::filter(Parameter == "Infectious period") %>% dplyr::select(value),
    susc = mp_df %>% dplyr::filter(Parameter == "Susceptibility") %>% dplyr::select(value)
  ) %>% purrr::imap(function(df, k) df$value)
  
  if (missing(data_lst))
    data_lst <- fluScenarioR::loadData("2009")
  if (!missing(pandemic)) {
    rel_df <- fluScenarioR::relative_incidence_df %>%
      dplyr::filter(!(Pandemic == "2009" & Date > lubridate::as_date("2009-11-15"))) %>%
      dplyr::filter(Pandemic == pandemic)
    data_lst$dates <- rel_df$Date
  }
  
  if (!missing(startDate)) {
    data_lst$dates <- seq(startDate, by = 7, length.out = 52)
  }
  
  #stopifnot(data_lst$no_age_groups == 6, "Currently the cpp model is hardcoded for 12 groups in total (6 age groups, 2 risk groups)")
  
  data_lst$cm_df <- data_lst$cm_df %>% dplyr::mutate(value = ifelse(as.numeric(AgeGroup) >= 3 | as.numeric(AgeGroup2) >= 3, value*0.469, value))
  
  cm_f <- contact_mask_f(peak, data_lst$dates[1], data_lst$age_limits, contact_group = contact_group)
  
  pars <- list()
  pars$max_scenario <- max_scenario
  pars$times <- as.numeric(data_lst$dates)
  pars$N <- sum(data_lst$age_groups$value)
  pars$minimal_infected <- pars$N*minimal_irate
  seq(1, long)
  
  if (long) {
    stopifnot(pars$times[2] - pars$times[1] == 7)
    pars$times <- c(pars$times, seq(max(pars$times) + 7, by = 7, length.out = 52))
  }
  
  if (daily)
    pars$times <- seq(min(pars$times), max(pars$times))
  
  function(par) {
    stopifnot(n == 5)
    cw <- cm_f(par)
    delay <- 0
    if (peak == 2) {
      data_lst$initial_infected <- rep(1000*logistic(par[2]), length(data_lst$initial_infected))
      delay <- logistic(par[3])*lubridate::duration(4, units = "week")/n
      #epi_par$R0 <- 1.33
    } else if (peak == 3) {
      data_lst$initial_infected <- rep(1000*logistic(par[3]), length(data_lst$initial_infected))
      delay <- logistic(par[4])*lubridate::duration(4, units = "week")/n
    }
    # Make sure we are done with par
    par <- NULL
    
    arg <- list(n = 1, R0 = epi_par$R0, susceptibility = epi_par$susc, data_lst = data_lst,
                infection_delays = c(epi_par$latent_period, epi_par$infectious_period),
                contact_mask = cw, ..., startPandemic = 0, include_novacc = F)
    if (!is.null(arg$duration) && arg$duration == 0)
      arg$dosesPerWeek <- NULL
    
    
    #suppressMessages(suppressWarnings(capture.output(
    #   ode <- do.call(fluScenarioR::rModel, arg), type = c("output", "message"))))
    #inc_df <- ode$incidence %>%
    #  dplyr::mutate(value = value/n) %>%
    #k  dplyr::select(-Week)
    #qs::qsave(arg, file = paste0(fn_id, "incidence_f.qs"))
    
    nregions <- 5
    #fn_id <<- fn_id + 1
    args <- arg
    ddelay <- as.numeric(delay, units = "day")
    pars$region_times_base <- c(0, cumsum(rep(ddelay, 4)))
    # transmissibility
    cm <- args$data_lst$cm_df %>%
      unnest_contact_matrix() %>%
      spread_contact_matrix(args$data_lst$dates[1] - 1)
    cm <- cm[1, 1:args$data_lst$no_age_groups, 1:args$data_lst$no_age_groups]
    pars$transmissibility <- fluEvidenceSynthesis::as_transmission_rate(args$R0, cm,
                                                                        args$data_lst$age_groups$value,
                                                                        duration = args$infection_delays[2])
    
    pars$gamma1 <- 1/args$infection_delays[1]
    pars$gamma2 <- 1/args$infection_delays[2]
    if (is.null(args$efficacy)) {
      pars$ve <- rep(0, 12*nregions)
    } else {
      pars$ve <- rep(rep(args$efficacy, 2), 5)
    }
    pars$susc <- rep(args$susceptibility, 2)
    pars$region_times <- pars$times[1] + pars$region_times_base
    pars$I0 <- args$data_lst$initial_infected
    pars$pop <- args$data_lst$pop_df$value
    
    # Contact list and dates
    cm_arr <- args$contact_mask %>% dplyr::rename(scale = value) %>%
      dplyr::full_join(args$data_lst$cm_df) %>% dplyr::mutate(value = value*scale) %>%
      dplyr::select(-scale) %>%
      unnest_contact_matrix() %>%
      spread_contact_matrix()
    for (i in 1:dim(cm_arr)[1]) {
      pars$cm_lst[[i]] <- cm_arr[i,,]
    }
    pars$contact_dates <- args$contact_mask$Date |> unique() |> sort() |> as.numeric()
    
    # Vacc list
    if (is.null(args$dosesPerWeek)) {
      vc <- list(object = empty_vaccination_calendar(6, 2))
    } else {
      vc <- rVaccCalendar(min(args$data_lst$dates), max(args$data_lst$dates),
                          data_lst$pop_df, dosesPerWeek = args$dosesPerWeek,
                          startDate = min(args$data_lst$dates) + lubridate::duration(args$vacc_delay, units = "weeks"),
                          maxUptake = args$maxUptake,
                          data_lst = args$data_lst, efficacy = pars$ve, order = args$order, duration = args$duration)
    }
    vacc_lst <- list()
    nregions <- 5
    vacc_lst[[1]] <- as_cpp_array(rep(0, 12*nregions))
    
    for (i in 1:nrow(vc$object$calendar)) {
      vacc_lst[[i+1]] <- as_cpp_array(rep(vc$object$calendar[i,1:12], nregions))
    }
    
    pars$vacc_lst <- vacc_lst
    pars$vacc_times <- c(pars$times[1], vc$object$dates |> as.numeric())
    
    sev <- 0
    if (severity == "High")
      sev <- 1
    
    pars$gamma_m <- 1/3.0 # Upper estimate pandemic 2009 @camp_risk_2010
    if (length(hospital_rates) == 0) {
      #message(hospital_rates)
      if (behaviour != 1)
        stop("Need hospitalisation rates to run behavioural models")
      pars$hospital <- as_cpp_array(rep(0, 12))
    } else {
      #pars$mortality <- as_cpp_array(rep(0, 12))
      pars$hospital <- hospital_rates
    }
    
    pars$gamma_mort <- 1/11.0 # Needs reference
    if (length(mortality_rates) == 0) {
      #message(hospital_rates)
      pars$mortality <- as_cpp_array(rep(0, 12))
    } else {
      #pars$mortality <- as_cpp_array(rep(0, 12))
      pars$mortality <- mortality_rates
    }
    
    #message(pars$mortality)
    
    #qs::qsave(pars, file = paste0("~/Downloads/", fn_id, "pars.qs"))
    pars$severity <- sev
    pars$behaviour <- behaviour - 1
    pars
  }
}

incidence_cpp_f <- function(collect_f) {
  function(params) {
    pars <- collect_f(params)
    minimal_irate <- dplyr::last(params)/56000000
    pars$minimal_infected <- pars$N*minimal_irate
    
    foo(pars$transmissibility, pars$gamma1, pars$gamma2, pars$gamma_m, pars$gamma_mort, pars$ve,
        pars$hospital, pars$mortality, as_cpp_array(pars$susc), pars$cm_lst, pars$contact_dates,
        pars$vacc_lst, pars$vacc_times, as_cpp_array(pars$I0),
        as_cpp_array(pars$pop), pars$region_times, pars$times,
        pars$severity, pars$behaviour, pars$max_scenario, pars$minimal_infected)
  }
}

distance_cpp_f <- function(collect_f, data_v, base_scale) {
  function(params) {
    pars <- collect_f(params)
    data_scale <- logistic(params[length(params) - 1])/base_scale
    minimal_irate <- dplyr::last(params)/56000000
    pars$minimal_infected <- pars$N*minimal_irate
    
    foo_abs(data_v, data_scale, pars$transmissibility, pars$gamma1,
            pars$gamma2, pars$gamma_m, pars$gamma_mort, pars$ve, pars$hospital,
            pars$mortality, as_cpp_array(pars$susc), pars$cm_lst, pars$contact_dates,
            pars$vacc_lst, pars$vacc_times, as_cpp_array(pars$I0), as_cpp_array(pars$pop),
            pars$region_times, pars$times, pars$severity, pars$behaviour, pars$max_scenario, pars$minimal_infected) -> ll
    ll + 0.5*abs(2000-dplyr::last(params))
  }
}
