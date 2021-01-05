#' Closing the observed death rations
#'
#' @description This functions closes the observed death rates of the country of interest.
#'
#' @param deaths_obs The death counts.
#' @param exp_obs The exposures numbers.
#' @param kannisto_nages The extrapolation ages used for closing the observed mortality
#' rates.
#' @param kannisto_nobs The number of observations used to extrapolate.
#'
#' @details The argument \code{deaths_obs} should be a matrix (year x ages) of
#' the deaths of your country of interest, e.g. as provided by the function
#' \code{\link{get_mortality_data}}.
#'
#' The same holds for the argument \cdoe{exp_obs}.
#'
#' @examples
#' lst   <- MultiMoMo::european_mortality_data
#' dat_M <- lst$M
#' dat_F <- lst$F
#'
#' deaths_obs <- dat_M$UNI$BE$dtx
#' exp_obs    <- dat_M$UNI$BE$etx
#' kannisto_nages <- 30
#' kannisto_nobs  <- 11
#' close_obs_m    <- close_obs_death_rates(deaths_obs, exp_obs, kannisto_nages, kannisto_nobs)
#'
#' @export

close_obs_death_rates <- function(deaths_obs, exp_obs, kannisto_nages, kannisto_nobs){

  # Age/Years
  dimn <- dimnames(deaths_obs)
  xv   <- as.numeric(dimn[[2]])
  yv   <- as.numeric(dimn[[1]])

  n_age   <- length(xv)
  n_years <- length(yv)

  # Observed central death rates
  m_obs <- t(deaths_obs/exp_obs)

  # Upperbound on observed central death rates
  m_obs[which(m_obs >= 1, arr.ind = T)] <- 0.99


  # Enlarge array for kannisto extrapolation
  mhat_all <- array(dim = c(n_age + kannisto_nages, n_years))
  mhat_all[1:n_age, 1:n_years] <- m_obs[1:n_age, 1:n_years]

  # Extend using kannisto for each year and simulation
  kannisto_obsages <- tail(xv, kannisto_nobs)
  kannisto_extrapolateages <- tail(xv, 1) + 1:kannisto_nages

  for(t in 1:n_years){
    mhat_all[n_age + 1:kannisto_nages,t] <- kannisto(t, array(m_obs,c(dim(m_obs),1)),
                                                     1, kannisto_nobs,
                                                     kannisto_extrapolateages,
                                                     kannisto_obsages)}
  dimnames(mhat_all) <- list(0:120, yv)

  mhat_all
}


