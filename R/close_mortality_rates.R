#' Closing the projected mortality rates
#'
#' @description This function closes the fitted and simulated mortality
#' rates using the method of Kannisto.
#'
#' @param yv The vector of years.
#' @param sim_qxt The mortality rates to close. An array with dimensions
#' (age, years, simulations).
#' @param kannisto_nages The extrapolation ages used for closing the mortality
#' rates.
#' @param kannisto_nobs The number of observations used to extrapolate.
#' @param parallel Logical. Use parallel processing to speed up the calculations?
#'
#' @details
#' The argument \code{sim_qxt} should be the output of the function
#' \code{\link{project_mortality_rates}}
#'
#' @return An array of dimension (age, years, simulations) reporting the closed
#' mortality rates.
#'
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
#' n_ahead     <- 50
#' n_sim       <- 10000
#' est_method  <- "PORT"
#' proj_par    <- project_parameters(fit_M, fit_F, n_ahead, n_sim, arima_spec, est_method)
#' proj_rates  <- project_mortality_rates(fit_M, fit_F, proj_par)
#'
#' kannisto_nages <- 30
#' kannisto_nobs  <- 11
#' close_rates_M  <- close_mortality_rates(yvSPEC, proj_rates$Male, kannisto_nages, kannisto_nobs)
#' close_rates_F  <- close_mortality_rates(yvSPEC, proj_rates$Female, kannisto_nages, kannisto_nobs)
#'
#' @importFrom stats lm.fit
#' @importFrom parallel detectCores
#' @importFrom snow makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach %dopar%
#'
#' @export

close_mortality_rates <- function(yv, sim_qxt, kannisto_nages, kannisto_nobs, parallel = FALSE){

  # Make sure to deal with an array
  if(length(dim(sim_qxt)) != 3)
    sim_qxt <- array(sim_qxt, c(dim(sim_qxt),1), dimnames = dimnames(sim_qxt))

  # From mortality rates to forces of mortality
  mxt <- -log(1-sim_qxt)

  # Some information retreived from objects
  xv      <- as.numeric(rownames(mxt[,,1]))
  yvv     <- as.numeric(colnames(mxt[,,1]))
  n_ahead <- length(yvv) - length(yv)
  n_sim   <- dim(mxt)[3]
  n_age   <- length(xv)

  # Number of cores
  nc <- ifelse(parallel == FALSE, 2, parallel::detectCores())

  # Some upperbound on mtx
  mxt[which(mxt >= 1, arr.ind = T)] <- 0.99

  # Split fitted and simulated mxt
  mxt_hat <- array(mxt[,as.character(yv),1], dim = c(n_age, length(yv), 1))
  mxt_sim <- mxt[,as.character(yvv[-c(1:length(yv))]),]


  ## Closing the simulations
  if(n_ahead > 0){
    # Enlarge mxt array for kannisto extrapolation
    mxt_sim_all             <- array(dim = c(n_age + kannisto_nages, n_ahead, n_sim))
    mxt_sim_all[1:n_age,,]  <- mxt_sim

    # Extend using kannisto for each year and simulation
    kannisto_obsages         <- tail(xv, kannisto_nobs)
    kannisto_extrapolateages <- tail(xv, 1) + 1:kannisto_nages

    # Parallel code for closing simulated forces of mortality
    cl <- snow::makeCluster(nc-1)
    doSNOW::registerDoSNOW(cl)
    close_sim <- foreach::foreach(t=1:n_ahead, .export = c('kannisto')) %dopar% {
      kannisto(t, mxt_sim, n_sim, kannisto_nobs, kannisto_extrapolateages, kannisto_obsages)
    }
    snow::stopCluster(cl)

    # Store into mxt_all array
    arr <- sapply(close_sim, function(x) x, simplify = 'array')
    mxt_sim_all[n_age + 1:kannisto_nages,,] <- aperm(arr, c(1,3,2))} else
      mxt_sim_all <- NULL



  ## Closing fitted mxt

  # Enlarge mhat array for kannisto extrapolation
  mxt_hat_all            <- array(dim = c(n_age + kannisto_nages, length(yv), 1))
  mxt_hat_all[1:n_age,,] <- mxt_hat

  # Parallel code for closing fitted forces of mortality
  cl <- snow::makeCluster(nc-1)
  doSNOW::registerDoSNOW(cl)
  close_fit <- foreach::foreach(t=1:length(yv), .export = c('kannisto')) %dopar% {
    kannisto(t, mxt_hat, n_sim = 1, kannisto_nobs, kannisto_extrapolateages,kannisto_obsages)
                                      }
  snow::stopCluster(cl)

  # Store into mxt_all array
  arr <- sapply(close_fit, function(x) x, simplify = 'array')
  mxt_hat_all[n_age + 1:kannisto_nages,,] <- aperm(arr, c(1,3,2))

  ## Concatenate all results
  mxt_new <- array(0,dim=c(n_age + kannisto_nages, length(yvv), n_sim))
  dimnames(mxt_new) <- list(xv[1]:(tail(xv,1)+kannisto_nages), yvv, 1:n_sim)

  mxt_new[,as.character(yv),] <- mxt_hat_all[,,1]
  mxt_new[,as.character(yvv[! yvv %in% yv]),] <- mxt_sim_all

  # Transform to mortality rates
  1 - exp(-mxt_new)}



#' @keywords internal

kannisto <- function(t, mxt, n_sim, kannisto_nobs, kannisto_extrapolateages, kannisto_obsages){
  #obtain mx used for extrapolation
  kannisto_mx_obs <- lapply(1:n_sim, function(x) tail(mxt[,t,x], kannisto_nobs))

  #apply logit transformation
  kannisto_mx_obs <- lapply(1:n_sim, function(x) log(kannisto_mx_obs[[x]]/(1-kannisto_mx_obs[[x]])))

  #apply the regression and obtain phi1 and phi2
  mat <- matrix(cbind("int" = 1,"Age" = kannisto_obsages),ncol = 2)
  kannisto_lm   <- lapply(1:n_sim, function(x) lm.fit(mat, kannisto_mx_obs[[x]]))
  kannisto_phi1 <- lapply(1:n_sim, function(x) exp(kannisto_lm[[x]]$coefficient[1]))
  kannisto_phi2 <- lapply(1:n_sim, function(x) kannisto_lm[[x]]$coefficient[2])

  # Calculate mx and put into array
  sapply(1:n_sim, function(x) kannisto_phi1[[x]]*exp(kannisto_phi2[[x]]*kannisto_extrapolateages)/
           (1+kannisto_phi1[[x]]*exp(kannisto_phi2[[x]]*kannisto_extrapolateages)))}




