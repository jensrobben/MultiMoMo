#' Plot the Poisson Pearson residuals of the multipopulation Li-Lee model
#'
#' @description This function returns a level plot of the Poisson Pearson residuals
#'   of the fitted multipopulation Li-Lee model. These Pearson residuals are defined as
#'   \deqn{(d_{x,t} - \hat{mu}_{x,t}E_{x,t})/\sqrt{\hat{mu}_{x,t}E_{x,t}},}
#'   and are visualized as a function of age (\code{xv}) and time (\code{yv}).
#'   The greener the Pearson residuals are coloured, the closer they are to zero
#'   and hence the better the estimated deaths from the Li-Lee model correspond
#'   to the observed deaths.
#'
#' @param xv The vector of ages.
#' @param yv The vector of years (for your country of interest).
#' @param fit The fitted Li-Lee model.
#' @param country_spec The user code of your country of interest.
#' @param sex The gender of interest (\code{"Male"} or \code{"Female"})
#'
#' @return A level plot showing the Poisson Pearson residuals.
#'
#' @details The object \code{fit} must be the output of the function \code{\link{fit_li_lee}}.
#'
#' @examples
#' lst   <- MultiMoMo::european_mortality_data
#' dat_M <- lst$M
#' xv    <- 0:90
#' yv = yvSPEC <- 1970:2018
#' Countries   <- names(dat_M$UNI)
#' country_spec <- "BE"
#' fit_M <- fit_li_lee(xv, yv, yvSPEC, country_spec, dat_M, "NR", TRUE, FALSE)
#' plot_residuals_li_lee(xv, yvSPEC, fit_M, country_spec, "Male")
#'
#' @importFrom tidyr gather
#'
#' @export



plot_residuals_li_lee  <- function(xv, yv, fit, country_spec, sex){

  if (length(fit$k.t) != length(yv))
    stop("The vector yv should have the same length as the length of the time series k.t.")

  if (length(fit$A.x) != length(xv))
    stop("The vector xv should have the same length as the length of the time series A.x, B.x, a.x and b.x.")

  if (! sex %in% c("Male", "Female"))
    stop("The argument sex must be either 'Male' or 'Female'.")

  if(! country_spec %in% MultiMoMo::country_codes[,"User_code"])
    stop("Wrong user code. Check MultiMoMo::country_codes to look up the correct user code for your country")

  # Make data ready for plotting
  res      <- as.data.frame(res)
  res$Year <- as.numeric(rownames(res))
  res      <- res %>% gather("Age","Value",-Year)
  res$Age  <- as.numeric(as.character(res$Age))

  ggplot(res, aes(x = Age, y = Year, fill = Value)) + geom_tile() + theme_bw(base_size = 20) +
    scale_fill_gradient2(midpoint = 0, high = "red2", low = "blue2") +
    scale_x_continuous(expand = c(0, 0),breaks = seq(0,90,10)) + scale_y_continuous(expand = c(0, 0))
}

