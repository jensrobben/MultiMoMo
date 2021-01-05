#' Plot simulated life expectancies
#'
#' @description This function plots the simulated period and/or cohort life
#' expectancies for a particular age, as retreived from the function
#' \code{\link{life_exp}}, as a function of time.
#'
#' @param le_yv The vector of years.
#' @param age The age.
#' @param le_sim The simulated period and/or cohort life expectancies.
#' @param sex Character vector. "Male" or "Female".
#' @param type Character vector. Options are "per" and/or "coh".
#' @param quantile The quantile of the confidence bound, e.g. "99%".
#' @param m_obs Optional. The observed and closed death rates of the country of
#'   interest.
#'
#'
#' @details
#' The argument \code{le_sim} should be the output of the function
#' \code{\link{life_exp}}
#'
#' The argument \code{m_obs} should be the output of the function
#' \code{\link{close_obs_death_rates}}
#'
#' @return A ggplot representing the evolution and the uncertainty of the period
#' and/or cohort life expectancy for a particular age and gender.
#'
#' @examples
#' lst   <- MultiMoMo::european_mortality_data
#' xv    <- 0:90
#' yv = yvSPEC <- 1988:2018
#' Countries   <- names(lst$M$UNI)
#' CountrySPEC <- "BE"
#'
#' lst   <- get_data_specific_range(xv, yv, yvSPEC, lst, Countries, CountrySPEC)
#' dat_M <- lst$M
#' dat_F <- lst$F
#'
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
#' le_yv   <- 1988:2070
#' le_ages <- c(0,65)
#' le_type <- c("per", "coh")
#' le_M    <- life_exp(le_yv, le_type, le_ages, close_rates_M)
#'
#' m_obs_cl <- close_obs_death_rates(dat_M$UNI$BE$dtx, dat_M$UNI$BE$etx, kannisto_nages, kannisto_nobs)
#' plot_life_exp(le_yv, 0, le_M, "Male", le_type, "99%", m_obs = m_obs_cl)
#'
#'
#' @import dplyr
#'
#' @export



plot_life_exp <- function(le_yv, age, le_sim, sex, type, quantile, m_obs = NULL){

  # Age indicator
  age_ind <- paste0("Age_",age)

  # Life expectancy simulation
  le_per <- if("per" %in% type) as.data.frame(le_sim[["per"]][[age_ind]]) else NULL
  le_coh <- if("coh" %in% type) as.data.frame(le_sim[["coh"]][[age_ind]]) else NULL

  # Add column for plotting
  le_per[,"Year"] <- as.numeric(rownames(le_per))
  le_coh[,"Year"] <- as.numeric(rownames(le_coh))

  # calculate observed period life expectancy
  le_yv_per <- as.numeric(dimnames(m_obs)[[2]])
  le_obs <- if("per" %in% type) life_exp(le_yv_per, c("per"), age, 1-exp(-m_obs), FALSE)

  # quantile
  nq <- as.numeric(gsub("[^0-9A-Za-z///' ]","" , quantile ,ignore.case = TRUE))
  qq <- c((1 - nq/100)/2, 0.50, (1 + nq/100)/2)

  # Make object ready for plotting
  le_per <- if("per" %in% type)
    as.data.frame(le_per %>% gather("key","Value",-Year))
  le_coh <- if("coh" %in% type)
    as.data.frame(le_coh %>% gather("key","Value",-Year))

  # Find quantiles period life expectancies
  if("per" %in% type){
    le_per <- le_per %>% group_by(Year) %>% summarise(qv = quantile(Value, probs = qq), qt = qq)
    le_per[,"min"] <- rep(pull(le_per[which(le_per[,"qt"] == min(qq)),],"qv"), each = 3)
    le_per[,"max"] <- rep(pull(le_per[which(le_per[,"qt"] == max(qq)),],"qv"), each = 3)}

  # Find quantiles cohort life expectancies
  if("coh" %in% type){
    le_coh <- le_coh %>% group_by(Year) %>% summarise(qv = quantile(Value, probs = qq), qt = qq)
    le_coh[,"min"] <- rep(pull(le_coh[which(le_coh[,"qt"] == min(qq)),],"qv"), each = 3)
    le_coh[,"max"] <- rep(pull(le_coh[which(le_coh[,"qt"] == max(qq)),],"qv"), each = 3)}

  # Join data frames
  le <- rbind(if("per" %in% type) data.frame(le_per,"type" = "per"),
              if("coh" %in% type) data.frame(le_coh, "type" = "coh"))


  p <- ggplot(le, aes(x = Year, y = qv)) + geom_line(aes(group = interaction(qt,type))) +
    geom_ribbon(aes(ymin = min, ymax = max, group = interaction(qt,type), fill = type), alpha = 0.2) +
    theme_bw(base_size = 20) + xlab("Year") + ylab(bquote(e[.(age)])) + ggtitle(paste0(sex,'s')) +
    scale_fill_discrete(name = "")

  if("per" %in% type & ! is.null(m_obs)){
    df_per_le <- le_obs[["per"]][[paste0("Age_",age)]]
    df_per_le <- data.frame("Value" = unname(df_per_le[,1]),
                            "Year" = as.numeric(rownames(df_per_le)))
    p <- p + geom_point(mapping = aes(x = Year, y = Value, col = "dots"),
                   data = df_per_le) +
      scale_color_manual(name = "", values = c("dots" = "black"), labels = c("Obs. le")) +
      theme(legend.position = 'bottom')
  }

  p

}


