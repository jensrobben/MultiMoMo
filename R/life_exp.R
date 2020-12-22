#' Calculate period and cohort life expectancy
#'
#' @description This function calculates the estimated period and/or cohort life
#' expectancies based on the closed, fitted and simulated mortality rates, as
#' retreived from the function \code{\link{close_mortality_rates}}.
#'
#' @param le_yv The vector of years on which estimations are made.
#' @param le_type Character vector, specifying the type of life expectancy:
#' \code{c("per")}, \code{c("coh")} or \code{c("per", "coh")}.
#' @param le_ages The vector of ages for which life expectancy estimations are made.
#' @param sim_qxt_cl The closed mortality rates to close. An array with dimensions
#' (closed ages, years, simulations).
#' @param parallel Logical. Use parallel processing to speed up calculations?
#'
#' @details
#' The argument \code{sim_qxt_cl} should be the output of the function
#' \code{\link{close_mortality_rates}}
#'
#' @return A list containing
#' \itemize{
#' \item $per A list of length(le_ages) containing matrices of dimension
#' length(le_yv) x (number of simulations) that contain the estimated period
#' life expectancies.
#' \item $coh A list of length(le_ages) containing matrices of dimension
#' length(le_yv) x (number of simulations) that contain the estimated cohort
#' life expectancies.
#' }
#'
#' @examples
#' lst   <- MultiMoMo::european_mortality_data
#' dat_M <- lst$M
#' dat_F <- lst$F
#' xv    <- 0:90
#' yv = yvSPEC <- 1970:2018
#' Countries   <- names(dat_M$UNI)
#' CountrySPEC <- "BE"
#' fit_M <- fit_li_lee(xv, yv, yvSPEC, CountrySPEC, dat_M, "NR", TRUE, FALSE)
#' fit_F <- fit_li_lee(xv, yv, yvSPEC, CountrySPEC, dat_F, "NR", TRUE, FALSE)
#'
#' arima_spec  <- list(K.t_M = "RWD", k.t_M = "AR3.1", K.t_F = "RWD", k.t_F = "AR5.0")
#' n_ahead     <- 52 + 120
#' n_sim       <- 1000
#' est_method  <- "PORT"
#' proj_par    <- project_parameters(fit_M, fit_F, n_ahead, n_sim, arima_spec, est_method)
#' proj_rates  <- project_mortality_rates(fit_M, fit_F, proj_par)
#'
#' kannisto_nages <- 30
#' kannisto_nobs  <- 11
#' close_rates_M  <- close_mortality_rates(yvSPEC, proj_rates$Male, kannisto_nages, kannisto_nobs)
#' close_rates_F  <- close_mortality_rates(yvSPEC, proj_rates$Female, kannisto_nages, kannisto_nobs)
#'
#' le_yv   <- 1970:2070
#' le_ages <- c(0,65)
#' le_type <- c("per", "coh")
#' le_M    <- life_exp(le_yv, le_type, le_ages, close_rates_M)
#'
#'
#' @importFrom parallel detectCores
#' @importFrom snow makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar%
#'
#' @export







life_exp <- function(le_yv, le_type, le_ages, sim_qxt_cl, parallel = TRUE){
  # General information retreived from input objects
  dimension <- dim(sim_qxt_cl)
  simul1    <- sim_qxt_cl[,,1]
  xv        <- as.numeric(rownames(simul1))
  yvv       <- as.numeric(colnames(simul1))
  n_sim     <- dimension[3]
  n_xv      <- length(xv)

  # Check if enough simulated mortality rates in the future
  if("coh" %in% le_type)
    if(! all(seq(le_yv[1], tail(le_yv,1)+tail(xv,1)) %in% yvv))
      stop("Too few years in 'sim_qxt_cl' to estimate life expectancies for all years in 'le_yv'.")

  if(! all(le_yv %in% yvv))
    stop("Too few years in 'sim_qxt_cl.")

  # Create list to store cohort and/or period life expectancy per method
  mat           <- matrix(, nrow = length(le_yv), ncol = n_sim)
  dimnames(mat) <- list(le_yv, 1:n_sim)
  outp          <- rep(list(mat), length(le_ages))
  names(outp)   <- paste0("Age_", le_ages)
  output        <- rep(list(outp), length = length(le_type))
  names(output) <- le_type

  # Death_rates
  sim_mxt_cl <- -log(1-sim_qxt_cl)

  # Make grid for years and ages to consider
  grid_year_age <- expand.grid(le_yv, le_ages)

  # Parallel execution of computations
  nc <- ifelse(parallel == FALSE, 2, parallel::detectCores())
  cl <- snow::makeCluster(nc-1)
  doSNOW::registerDoSNOW(cl)
  le <- foreach::foreach(t=1:nrow(grid_year_age), .verbose = FALSE,
                                .export = c('sim_mxt_cl', 'diags','grid_year_age', 'le_type','le_calc','diags')) %dopar% {
    le_calc(sim_mxt_cl, age = grid_year_age[t,2], year = grid_year_age[t,1], le_type)}
  snow::stopCluster(cl)

  names(le) <- sapply(1:nrow(grid_year_age),
                      function(x) paste0(grid_year_age[x,], collapse = "_"))

  # Store results
  for(a in le_ages)
    for (t in le_yv)
      for(type in le_type)
        output[[type]][[paste0("Age_",a)]] <- t(sapply(le_yv, function(t)
          le[[paste0(t,"_",a)]]$per, simplify = 'array'))

  output
}


#' @keywords internal
le_calc <- function(sim_mxt_cl, age, year, type){
  # Indicators
  age_ind  <- age + 1
  year_ind <- which(as.numeric(colnames(sim_mxt_cl[,,1])) %in% year)
  n_sim    <- dim(sim_mxt_cl)[3]
  n_xv     <- dim(sim_mxt_cl)[1]

  # Constant term in cohort/period LE
  constant <- lapply(1:n_sim, function(sim) (1-exp(-sim_mxt_cl[age_ind,year_ind,sim]))/
                       sim_mxt_cl[age_ind,year_ind,sim])

  # Cohort life expectancy
  if(c("coh") %in% type){
    kvalue     <- age_ind - year_ind
    bandvector <- lapply(1:n_sim, function(sim) diags(sim_mxt_cl[,,sim], which = kvalue))
    n.band     <- length(bandvector[[1]])
    start.band <- ifelse(kvalue>0, year_ind, age_ind)
    eproj      <- sapply(1:n_sim, function(sim) sum(cumprod(exp(-bandvector[[sim]][start.band:(n.band-1)]))*
                                                      (1-exp(-bandvector[[sim]][(start.band + 1):n.band]))/
                                                      bandvector[[sim]][(start.band + 1):n.band]) + constant[[sim]])} else
                                                        eproj <- NULL
  # Period life expectancy
  if(c("per") %in% type){
    eproj2 <- sapply(1:n_sim, function(sim) sum(cumprod(exp(-sim_mxt_cl[age_ind:(n_xv-1), year_ind, sim]))*
                                                  (1-exp(-sim_mxt_cl[(age_ind+1):n_xv, year_ind, sim]))/
                                                  sim_mxt_cl[(age_ind+1):n_xv, year_ind, sim]) + constant[[sim]])} else
                                                    eproj2 <- NULL

  list("coh" = eproj, "per" = eproj2)

}



#' @keywords internal
diags <- function(x, which = 0){
  if(which < -(ncol(x)-1) | which > nrow(x) - 1)
    stop("Impossible which-value")
  s <- sign(which)
  k <- abs(which)
  ind <- (k+1):min(ncol(x),k+nrow(x))
  t <- if(s < 0) as.numeric(t(x))[ind + ncol(x)*(0:(length(ind)-1))] else
      as.numeric(t(x))[1:min(nrow(x)-k,ncol(x)) + ncol(x)*(k:min(nrow(x)-1,ncol(x)+k-1))]
  t
}

