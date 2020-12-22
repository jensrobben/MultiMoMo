#' Project mortality rates
#'
#' \loadmathjax
#' @description This function project mortality rates using the simulations that
#' were set up for the joint time series \mjeqn{(K_t^M, \kappa_t^M, K_t^F,
#' \kappa_t^F)}{ASCII representation}.
#'
#' @param fit_M The calibrated Li-Lee model for males.
#' @param fit_F The calibrated Li-Lee model for females.
#' @param proj_par The (simulations of) the projected time dependent parameters
#'   \mjeqn{(K_t^M, \kappa_t^M, K_t^F, \kappa_t^F)}{ASCII representation}.
#'
#' @details
#' The arguments \code{fit_M} and \code{fit_F} should be the outputs of the function
#' \code{\link{fit_li_lee}}.
#'
#' The argument \code{proj_par} should be the output of the function
#' \code{\link{project_parameters}}.
#'
#' @return A list containing the following objects
#' \itemize{
#'   \item The projected/simulated male mortality rates ($Male).
#'   \item The projected/simulated female mortality rates ($Female).
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
#' arima_spec <- list(K.t_M = "RWD", k.t_M = "AR3.1", K.t_F = "RWD", k.t_F = "AR5.0")
#' n_ahead    <- 50
#' n_sim      <- 10000
#' est_method <- "PORT"
#' proj_par   <- project_parameters(fit_M, fit_F, n_ahead, n_sim, arima_spec, est_method)
#' proj_rates <- project_mortality_rates(fit_M, fit_F, proj_par)
#'
#'
#' @export

project_mortality_rates <- function(fit_M, fit_F, proj_par){
  # Retreive information from input objects
  xv       <- as.numeric(colnames(fit_M$m.tx))
  yv       <- as.numeric(rownames(fit_M$m.tx))
  n_year 	 <- length(yv)
  n_age 	 <- length(xv)
  n_sim    <- ncol(proj_par$K.t_M)
  n_ahead  <- nrow(proj_par$K.t_M) - n_year

  # Define objects to store country specific mortality
  m.xt_M             <- array(NA, dim = c(n_age, (n_year + n_ahead), n_sim))
  m.xt_M[,1:n_year,] <- t(fit_M$m.tx)
  m.xt_F             <- array(NA, dim = c(n_age, (n_year + n_ahead), n_sim))
  m.xt_F[,1:n_year,] <- t(fit_F$m.tx)

  # Project forces of mortality
  m.xt <- lapply((n_year + 1):(n_year + n_ahead), FUN = function(t){
    x1 <- m.xt_M[,n_year,]*exp(fit_M$B.x %o% (proj_par$K.t_M[t,] - fit_M$K.t[n_year]) +
                                 fit_M$b.x %o% (proj_par$k.t_M[t,] - fit_M$k.t[n_year]))
    x2 <- m.xt_F[,n_year,]*exp(fit_F$B.x %o% (proj_par$K.t_F[t,] - fit_F$K.t[n_year]) +
                                 fit_F$b.x %o% (proj_par$k.t_F[t,] - fit_F$k.t[n_year]))
    list("Male" = x1, "Female" = x2)})

  m.xt_M[,(n_year + 1):(n_year + n_ahead),] <- aperm(sapply(1:n_ahead, function(t)
    m.xt[[t]]$Male, simplify = 'array'), c(1,3,2))
  m.xt_F[,(n_year + 1):(n_year + n_ahead),] <- aperm(sapply(1:n_ahead, function(t)
    m.xt[[t]]$Female, simplify = 'array'), c(1,3,2))

  dimnames(m.xt_M) <- list(xv, yv[1]:(tail(yv,1) + n_ahead), 1:n_sim)
  dimnames(m.xt_F) <- list(xv, yv[1]:(tail(yv,1) + n_ahead), 1:n_sim)

  # Switch to mortality rates
  q.xt_M <- 1-exp(-m.xt_M)
  q.xt_F <- 1-exp(-m.xt_F)

  list("Male" = q.xt_M, "Female" = q.xt_F)}
