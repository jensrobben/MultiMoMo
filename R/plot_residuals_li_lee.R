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
#' @param CountrySPEC The user code of your country of interest.
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
#' CountrySPEC <- "BE"
#' fit_M <- fit_li_lee(xv, yv, yvSPEC, CountrySPEC, dat_M, "NR", TRUE, FALSE)
#' plot_residuals_li_lee(xv, yvSPEC, fit_M, CountrySPEC, "Male")
#'
#' @importFrom lattice levelplot
#' @importFrom grDevices rainbow
#'
#' @export



plot_residuals_li_lee  <- function(xv, yv, fit, CountrySPEC, sex){

  if (length(fit$k.t) != length(yv))
    stop("The vector yv should have the same length as the length of the time series k.t.")

  if (length(fit$A.x) != length(xv))
    stop("The vector xv should have the same length as the length of the time series A.x, B.x, a.x and b.x.")

  if (! sex %in% c("Male", "Female"))
    stop("The argument sex must be either 'Male' or 'Female'.")

  if(! CountrySPEC %in% MultiMoMo::country_codes[,"User_code"])
    stop("Wrong user code. Check MultiMoMo::country_codes to look up the correct user code for your country")

  # Setting up some scales to nicely visualize the Pearson residuals
  nx     <- floor((length(xv)-1)/10)
  ny     <- floor((length(yv)-1)/10)
  scales <- list(y = list(at = seq(0.5, ny*10 + 0.5, by = 10), labels = seq(yv[1], yv[length(yv)], by = 10)),
                 x = list(at = seq(0.5, nx*10 + 0.5, by = 10), labels = seq(xv[1], xv[length(xv)], by = 10)))

  res  <- fit$residuals
  rg   <- max(c(-floor(min(res)),ceiling(max(res))))

  par(cex.main = 20, cex.lab = 1.5, cex.axis = 1.5)
  plot(levelplot(t(res), border = NA, border.lty = 1, border.lwd = 1, at = c(seq(-rg,rg, length.out = 16)),
                 col.regions = rainbow(16, start = 0, end = 0.6),
                 colorkey = list(height = 0.8, width = 0.75, axis.line = list(col="black")),
                 scales = scales, aspect = 'fill', ylab = list(cex = 1.5, label = "Year"),
                 xlab = list(cex = 1.5, label = "Age"),
                 main = list(cex = 2, label = paste("Pearson residuals:", sex))))
}

