#' Estimate the time series of the time dependent parameters in the calibrated Li-Lee model and
#' make projections
#'
#' \loadmathjax
#' @description
#' This function estimates the parameters in the  time series specifications
#' for the time dependent parameters in the calibrated Li-Lee model, namely
#' \mjeqn{K_t}{ASCII representation} and \mjeqn{\kappa_t}{ASCII representation},
#' jointly for males and females.
#'
#' @param fit_M The calibrated Li-Lee model for males.
#' @param fit_F The calibrated Li-Lee model for females.
#' @param n_ahead The number of years to project in the future.
#' @param n_sim The number of simulations/projections.
#' @param arima_spec The ARIMA time series specifications for \mjeqn{(K_t^M, \kappa_t^M, K_t^F, \kappa_t^F)}{ASCII representation}.
#' @param est_method The estimation method of the parameters in the time series specifications.
#'
#' @details
#' The arguments \code{fit_M} and \code{fit_F} are the outputs of the function
#' \code{\link{fit_li_lee}}.
#'
#' The argument \code{arima_spec} is a list specifying the type of time series used
#' for each time dependent parameter in the Li-Lee model and is of the form:
#' \deqn{list("K.t_M" = , "k.t_M" =  , "K.t_F" = , "k.t_F" = ).}
#' The capital \code{K.t} refers to the time dependent parameter in the Lee-Carter model
#' for the common trend, i.e. \mjeqn{K_t}{ASCII representation}. The small \code{k.t} refers
#' to the time dependent parameter \mjeqn{\kappa_t}{ASCII representation}
#' in the Lee-Carter model for the country-specific deviation from the common trend.
#' We jointly model the time series dynamics for men and women by assuming a
#' multivariate Gaussian distribution with mean \mjeqn{(0,0,0,0)}{ASCII represenation} and covariance matrix
#' \mjeqn{C}{ASCII represenation} for the error terms
#' \mjeqn{(\epsilon_t^M, \delta_t^M, \epsilon_t^F, \delta_t^F)}{ASCII representation}
#' of the multivariate time series \mjeqn{(K_t^M, \kappa_t^M, K_t^F, \kappa_t^F)}{ASCII representation}.
#' Possible choices for the time series dynamics are \code{RWD} (e.g. for \mjeqn{K_t}{ASCII representation})
#' or \code{ARk.j}, where \code{k} refers to an AR(k) process, an auto-regressive process
#' with lag \code{k}. Further \mjeqn{j \in \{0,1\}}{ASCII representation} and refers
#' to an AR process without or with an intercept respectively.
#'
#' The argument \code{est_method} specifies the estimation method of the parameters in
#' the time series specifications, given in \code{arima_spec}. A first option is \code{SUR},
#' referring to seemingly unrelated regression, that uses the function \code{\link[systemfit]{systemfit}}.
#' This method only works for a limited amount of cases. In particular, all time series
#' must jointly start at the same time t. For example, this means that in the example below,
#' we start the joint time series estimation at year 1975 (since an \code{AR5.0} is used).
#' A second option is \code{PORT}. This is a self-written objective function that
#' maximizes the log-likelihood of the multivariate Gaussian distribution for
#' \mjeqn{(\epsilon_t^M, \delta_t^M, \epsilon_t^F, \delta_t^F)}{ASCII represenation} using the function
#' \code{\link[stats]{nlminb}} in the \code{stats}-package. In the example below, the joint estimation
#' consists of three parts. For t = 1971-1972 only the \code{RWD} processes are considered leading to
#' a bi-variate normal log-likelihood. Then for t = 1973-1974 the two \code{RWD} processes together with
#' \code{k.t_M} are used, leading to a 3-variate normal log-likelihood. From 1975 on, we deal with the
#' four processes (4-variate normal log-likelihood). The sum of these three likelihoods is maximized
#' at once using PORT routines.
#'
#' @return A list containing the following objects
#' \itemize{
#'   \item The estimated parameters in the time series for \mjeqn{K_t^M}{ASCII representation}:
#'         $coef_KtM
#'   \item The estimated parameters in the time series for \mjeqn{\kappa_t^M}{ASCII representation}:
#'         $coef_ktM
#'   \item The estimated parameters in the time series for \mjeqn{K_t^F}{ASCII representation}:
#'         $coef_KtF
#'   \item The estimatedparameters in the time series for \mjeqn{\kappa_t^F}{ASCII representation}:
#'         $coef_ktF
#'   \item The estimated covariance matrix \mjeqn{C}{ASCII representation}:
#'         $cov_mat
#'   \item The \mjeqn{K_t^M}{ASCII representation} simulations: $K.t_M
#'   \item The \mjeqn{\kappa_t^M}{ASCII representation} simulations: $k.t_M
#'   \item The \mjeqn{K_t^F}{ASCII representation} simulations: $K.t_F
#'   \item The \mjeqn{\kappa_t^F}{ASCII representation} simulations: $k.t_F
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
#' proj       <- project_parameters(fit_M, fit_F, n_ahead, n_sim, arima_spec, est_method)
#'
#' @importFrom stats nlminb as.formula
#' @importFrom systemfit systemfit
#' @importFrom stats rnorm
#'
#' @export


project_parameters <- function(fit_M, fit_F, n_ahead, n_sim, arima_spec, est_method){

  # Some checks
  if(! est_method %in% c("SUR","PORT"))
    stop("est_method should be either 'SUR' or 'PORT'.")

  if(! is.list(arima_spec))
    stop("arima_spec should be a list (see ?project_parameters).")

  if(length(arima_spec) != 4)
    stop("arima_spec should be a list of length 4 specifying the time series dynamics of
         (K.t_M, k.t_M, K.t_F, k.t_F).")

  # Time series dynamics and specifications of (K.t_M, k.t_M, K.t_F, k.t_F)
  ts_list <- list(K.t_M = fit_M$K.t, k.t_M = fit_M$k.t, K.t_F = fit_F$K.t, k.t_F = fit_F$k.t)
  m_list  <- list(K.t_M = arima_spec$K.t_M, k.t_M = arima_spec$k.t_M,
                  K.t_F = arima_spec$K.t_F, k.t_F = arima_spec$k.t_F)

  # Time series information
  TS_n   <- length(ts_list)		# number of time series
  yv     <- as.numeric(rownames(fit_M$m.tx)) # calibration period yvSPEC
  nY     <- length(yv)	# number of years in calibration period
  lags   <- sapply(m_list, function(x) ifelse(x == "RWD", 1, as.numeric(substr(x,start = 3,stop = 3))))
  max.lag <- max(lags)
  min.lag <- min(lags)

  # Number of parameters in each time series
  nPar.ts <- sapply(m_list, function(x) ifelse(x=="RWD",1,as.numeric(substr(x,start=3,stop=3)) +
                                                 ifelse(substr(x,start=5,stop=5) == "1",1,0)))
  nPar <- Reduce('+', nPar.ts) # Total
  nPar.vec <- unlist(sapply(1:TS_n, function(t) rep(t,nPar.ts[t]))) # indices

  # More checks
  if(! all(lags == lags[[1]]) & est_method == "SUR")
    stop("The 'SUR' estimation method does not work for varying lags. Take 'PORT'.")

  # Estimate time series models jointly
  if(est_method == "PORT")
    kappa_est <- port_estimation(fit_M, fit_F, ts_list, m_list) else
      kappa_est <- sur_estimation(fit_M, fit_F, ts_list, m_list)

  # Project time dependent parameters / simulations
  output <- list()

  for(i in 1:TS_n){
    if(n_sim != 0)
      output[[i]] <- array(dim = c(nY + n_ahead, n_sim))
    if(n_sim == 0)
      output[[i]] <- array(dim = c(nY + n_ahead, 1))
    output[[i]][1:nY,] <-  ts_list[[i]]
  }

  if(n_ahead != 0){
    if(n_sim == 0){
      errors <- matrix(0, nrow = TS_n, ncol = n_ahead)
      for(j in 1:TS_n){
        namej <- paste0("coef_", gsub("[^0-9A-Za-z///' ]", "" , names(m_list)[j], ignore.case = TRUE))
        output[[j]][(nY+1):(nY + n_ahead),1] <- simulate_ts(ts = ts_list[[j]],
                                                            name = m_list[[j]],
                                                            coeff = kappa_est[[namej]],
                                                            errors  = errors[j,])}}
    if(n_sim > 0){
      eval <- function(){
        obj <- matrix(NA, nrow = TS_n, ncol = n_ahead)
        errors <- t(mvrnorm(n = n_ahead, mu = rep(0, TS_n), Sigma = kappa_est$cov_mat))
        f <-  for(j in 1:TS_n){
          namej <- paste0("coef_", gsub("[^0-9A-Za-z///' ]", "" , names(m_list)[j], ignore.case = TRUE))
          obj[j,] <- simulate_ts(ts = ts_list[[j]], name = m_list[[j]], coeff = kappa_est[[namej]], errors  = errors[j,])}
        obj}
      sim <- replicate(n_sim, eval())

      for (s in 1:4)
        output[[s]][(nY + 1):(nY + n_ahead),] <- sim[s,,]}}


  kappa_sim <- kappa_est
  kappa_sim[["K.t_M"]] <- output[[1]]
  kappa_sim[["k.t_M"]] <- output[[2]]
  kappa_sim[["K.t_F"]] <- output[[3]]
  kappa_sim[["k.t_F"]] <- output[[4]]
  kappa_sim
}

#' @keywords internal
port_estimation <- function(fit_M, fit_F, ts_list, m_list){
  # Time series information
  TS_n   <- length(ts_list)		# number of time series
  yv     <- as.numeric(rownames(fit_M$m.tx)) # calibration period yvSPEC
  nY     <- length(yv)	# number of years in calibration period
  lags   <- sapply(m_list, function(x) ifelse(x == "RWD", 1, as.numeric(substr(x,start = 3,stop = 3))))
  max.lag <- max(lags)
  min.lag <- min(lags)

  # Number of parameters in each time series
  nPar.ts <- sapply(m_list, function(x) ifelse(x=="RWD",1,as.numeric(substr(x,start=3,stop=3)) +
                                                 ifelse(substr(x,start=5,stop=5) == "1",1,0)))
  nPar <- Reduce('+', nPar.ts) # Total
  nPar.vec <- unlist(sapply(1:TS_n, function(t) rep(t,nPar.ts[t]))) # indices


  ##  Yt-vector, Xt-matrices: Yt = Xt Phi + Zt

  # Construct Yt-vectors
  make_Yt <- function(t){
    ts    <- names(m_list)[t] # Kind of time series: K.t_M, K.t_F, k.t_M of k.t_F
    model <- m_list[[t]] # Model name: RWD, AR1.0, AR1.1, ...
    fit <- if(grepl("M", ts, fixed = T)) fit_M else fit_F
    if(model == "RWD")
      c(rep(NA, lags[t] - min.lag), (fit[[substr(ts,1,3)]])[(lags[t]+1):nY] -
          (fit[[substr(ts,1,3)]])[(lags[t]):(nY-1)]) else
        c(rep(NA, lags[t] - min.lag), (fit[[substr(ts,1,3)]])[(lags[t]+1):nY])}

  Y <- t(sapply(1:TS_n, function(t) make_Yt(t)))
  Y <- sapply(1:ncol(Y), function(x) as.numeric(na.omit(Y[,x])), simplify = F)

  # Construct Xt-matrices
  make_Xt <- function(t){
    mat <- matrix(NA, ncol = nPar, nrow = 0)
    for (s in 1:TS_n){
      ts <- names(m_list)[s] # Kind of time series: K.t_M, K.t_F, k.t_M of k.t_F
      model <- m_list[[s]] # Model name: RWD, AR1.1
      fit <- if(grepl("M", ts, fixed = T)) fit_M else fit_F
      L <- if(s == 1) 0 else sum(nPar.ts[1:(s-1)]) # Number of zeros that precede
      intercept <- if(substr(model,5,5) == "1") 1 else NULL
      if(t >= lags[s])
        V <- if(model == "RWD") 1 else c(intercept,fit[[substr(ts,1,3)]][t:(t-lags[s]+1)]) else
          next
      row.s <- c(rep(0,L), V, rep(0, nPar - sum(nPar.ts[1:s])))
      mat <- rbind(mat, row.s)
    }
    mat
  }

  X <- sapply(1:length(Y), function(t) make_Xt(t), simplify = F)

  # Construct initial starting vector for PORT optimalization
  bin1 <- unlist(sapply(1:TS_n, function(t) if(m_list[[t]] == "RWD") 0 else{
    int <- if(substr(m_list[[t]],5,5) == "1") 0 else NULL
    c(int,1,rep(0,nPar.ts[t]-1-length(int)))
  }))
  CC <- diag(rep(1,TS_n))
  bin2 <- chol(CC)[lower.tri(chol(CC),diag=T)]
  bin <- c(bin1, bin2)

  # Optimization of objective function
  nll <- nlminb(start = bin, objective = obj_f_port, gradient = NULL, hessian = NULL, X = X, Y = Y, N = nPar.vec,
                control = list(trace = F, eval.max = 10000, iter.max = 10000, rel.tol = 1e-10))

  # Save parameter estimates time series
  regression_results <- list()
  regression_results$coef_KtM <- nll$par[1:nPar][nPar.vec==1]
  regression_results$coef_ktM <- nll$par[1:nPar][nPar.vec==2] # Parameters time series
  regression_results$coef_KtF <- nll$par[1:nPar][nPar.vec==3] # Parameters time series
  regression_results$coef_ktF <- nll$par[1:nPar][nPar.vec==4] # Parameters time series

  # Save estimate covariance matrix C
  mat <- matrix(0, nrow = TS_n, ncol = TS_n)
  mat[lower.tri(mat,diag=T)] <- nll$par[-c(1:nPar)]
  regression_results$cov_mat <- mat%*%t(mat) # Covariance matrix of residuals Z_t = (epsilon_M, delta_M, epsilon_F, delta_F)
  dimnames(regression_results$cov_mat) <- list(names(m_list), names(m_list))

  regression_results}

#' @keywords internal
obj_f_port  <- function(b, X, Y, N){
  nr.ts <- sapply(1:length(X), function(t) nrow(X[[t]]))
  TS_n  <- max(nr.ts)
  nPar  <- length(b) - sum(TS_n:1) # Parameters in time series, not in C
  Phi   <- b[1:nPar] # Params (theta_M, theta_F, c_M, alpha_M, c_F alpha_F)
  C     <- matrix(0, ncol = TS_n, nrow = TS_n, byrow = F) # create covariance matrix structure
  C[lower.tri(C,diag = T)] <- b[-c(1:nPar)] #
  C <- C%*%t(C) # we get positive definite matrix now (should be to take inverse and to have positive determinant)

  obj <- list()
  for (j in 1:length(unique(nr.ts))){
    ind <- which(nr.ts == unique(nr.ts)[j])
    ts.ind <- unique(N[! apply(X[[ind[1]]],2,function(x) all(x==0))])
    Csub <- C[ts.ind,ts.ind]
    obj[[j]] <- 0.5*tr(solve(Csub) %*% Reduce("+",lapply(ind, function(t)
      (Y[[t]] - X[[t]]%*%Phi)%*%t(Y[[t]] - X[[t]]%*%Phi)))) + (length(ind)/2)*log(det(Csub)) +
      0.5*(length(ind))*length(ts.ind)*log(2*pi)
  }

  Reduce('+', obj)
}

#' @keywords internal
sur_estimation <- function(fit_M, fit_F, ts_list, m_list){
  # Time series information
  TS_n   <- length(ts_list)		# number of time series
  yv     <- as.numeric(rownames(fit_M$m.tx)) # calibration period yvSPEC
  nY     <- length(yv)	# number of years in calibration period

  # Get input (string) systemfit
  formula.list <- list()
  formula.list <- create_formula_string(ts_list,m_list) #added -> generalization

  # Perform SUR regression
  fit_sur <- systemfit(formula = formula.list, method="SUR", maxiter = 10000, methodResidCov = "noDfCor")

  # Save results
  regression_results <- list()
  regression_results$coef_KtM <- unname(fit_sur$eq[[1]]$coefficients)
  regression_results$coef_ktM <- unname(fit_sur$eq[[2]]$coefficients)
  regression_results$coef_KtF <- unname(fit_sur$eq[[3]]$coefficients)
  regression_results$coef_ktF <- unname(fit_sur$eq[[4]]$coefficients)
  regression_results$cov_mat  <- fit_sur$residCov
  dimnames(regression_results$cov_mat) <- list(names(m_list),names(m_list))

  regression_results}

#' @keywords internal
create_formula_string <- function(ts_list, m_list){
  # time series information
  lags   <- sapply(m_list, function(x) ifelse(x == "RWD", 1, as.numeric(substr(x,start = 3,stop = 3))))
  max.lag <- max(lags)
  min.lag <- min(lags)

  intercept <- sapply(m_list, function(x) ifelse(substr(x,start = 5,stop = 5) == "1", "1", "-1"))

  formula.list=list()

  for (k in 1:length(m_list)){
    index <- k
    ts_n <- length(ts_list[[index]])
    if(m_list[[index]] == "RWD"){
      formula.list[[index]] = as.formula(paste("ts_list","[[",index,"]][",
                                               max.lag+1,":",ts_n,"] - ts_list[[",index,"]][",
                                               max.lag,":",ts_n-1,"] ~ 1",sep=""))}
    if(m_list[[index]] != "RWD"){
      form <- paste("ts_list[[",index,"]][",max.lag+1,":",ts_n,"] ~ ",intercept[index],sep="")
      for (l in 1:lags[k]){
        form <- paste(form," + ts_list[[",index,"]][",max.lag+1-l,":",ts_n-l,"]",sep="")
      }
      formula.list[[index]] <- as.formula(form)}}

  return(formula.list)

}

#' @keywords internal
simulate_ts <- function(ts, name, coeff, errors){

  if(name == "RWD")
    ts_sim <- stats::filter(coeff + errors, filter = c(1), init = tail(ts,1), method = "recursive")

  if(name != "RWD"){
    intercept <- ifelse(substr(name,5,5) == "1",T, F) #in reverse time order k_(t+1) = c + a_1 * k_t + a_2 *k_(t-1) + delta_(t+1)
    x      <- if(intercept) coeff[1] + errors else errors
    filter <- if(intercept) coeff[-1] else coeff
    lag    <- as.numeric(substr(name, start = 3, stop = 3)) #lag of AR process
    init   <- rev(tail(ts, lag)) #initial kappa_t, in reverse time order
    ts_sim <- stats::filter(x = x, filter = filter, init = init, method = "recursive")
  }

  ts_sim}

#' @keywords internal
mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06)
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1)
    drop(X)
  else t(X)
}

#' @keywords internal
tr <- function (m)
{
  if (!is.matrix(m) | (dim(m)[1] != dim(m)[2]))
    stop("m must be a square matrix")
  return(sum(diag(m), na.rm = TRUE))
}

