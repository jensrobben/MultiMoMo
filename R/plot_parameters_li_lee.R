#' Plot the parameters of the multipopulation (Li-Lee) model
#'
#' @description This function plots the estimated parameters in the multipopulation Li-Lee model.
#'
#' @param xv The vector of ages.
#' @param yv The vector of years (for your country of interest).
#' @param fit The estimated Li-Lee model.
#' @param sex The gender of interest (\code{"Male"} or \code{"Female"})
#' @param country_spec The user code of your country of interest.
#' @param method The method used to fit the Li-Lee mortality model (currently only \code{"NR"} working).
#' @param type Which type of parameter estimates should be visualized? Either one of the common
#'   parameter estimates (\code{"A.x"}, \code{"B.x"} or \code{"K.t"}), one of the country-specific
#'   estimates (\code{"a.x"}, \code{"b.x"} or \code{"k.t"}) or all of them combined in one plot (\code{"ALL"})
#'
#' @details The object \code{fit} must be the output of the function \code{\link{fit_li_lee}}.
#'
#' @return A plot showing the parameter estimates of the fitted Li-Lee model.
#'
#' @examples
#' lst   <- MultiMoMo::european_mortality_data
#' dat_M <- lst$M
#' xv    <- 0:90
#' yv = yvSPEC <- 1970:2018
#' Countries   <- names(dat_M$UNI)
#' country_spec <- "BE"
#' fit_M <- fit_li_lee(xv, yv, yvSPEC, country_spec, dat_M, "NR", TRUE, FALSE)
#' plot_parameters_li_lee(xv, yvSPEC, fit_M, "Male", country_spec, "NR", "ALL")
#'
#' @import ggplot2
#' @importFrom plyr mapvalues
#' @importFrom graphics par
#' @importFrom stats ts
#'
#' @export



plot_parameters_li_lee <- function(xv, yv, fit, sex, country_spec, method, type){

  if (method != "NR")
    stop("Only the Newton-Raphson method ('NR') works for now.")

  if (length(fit$k.t) != length(yv))
    stop("The vector yv should have the same length as the length of the time series k.t.")

  if (length(fit$A.x) != length(xv))
    stop("The vector xv should have the same length as the length of the time series A.x, B.x, a.x and b.x.")

  if (! sex %in% c("Male", "Female"))
    stop("The argument sex must be either 'Male' or 'Female'.")

  if(! country_spec %in% MultiMoMo::country_codes[,"User_code"])
    stop("Wrong user code. Check MultiMoMo::country_codes to look up the correct user code for your country")

  param_names <- c("A.x", "B.x", "K.t", "a.x", "b.x", "k.t")

  if (! type %in% c(param_names, "ALL"))
    stop("Choose one of the valid parameter estimates: 'A.x', 'B.x', 'K.t', 'a.x', 'b.x', 'k.t' or 'ALL'")

  # k.t and K.t are possibly of different length
  yv.kt <- yv
  yv.Kt <- (tail(yv,1) - length(fit$K.t) + 1):tail(yv,1)


  xrange        <- list(xv, xv, yv.Kt, xv, xv, yv.kt)
  names(xrange) <- param_names
  if (type != "ALL")
    xrange <- xrange[[type]]

  smain <- list(bquote("Common"~.(sex):~A[x]), bquote("Common"~.(sex):~B[x]), bquote("Common"~.(sex):~K[t]),
                bquote(.(country_spec)~ .(sex): ~alpha[x]^(.(country_spec))),
                bquote(.(country_spec)~ .(sex): ~beta[x]^(.(country_spec))),
                bquote(.(country_spec)~ .(sex): ~kappa[t]^(.(country_spec))))
  names(smain) <- param_names

  if (type != "ALL")
    smain <- smain[[type]]

  if (type == "ALL")
    sxlab <- list("Age", "Age", "Year", "Age", "Age", "Year") else
      sxlab <- if(grepl("t",type)) "Year" else "Age"

  if (type != "ALL"){
    par(cex.main = 2, cex.lab = 1.5, cex.axis = 1.5)
    plot(ts(fit[[type]],start = xrange[1]), main = smain, xlab = sxlab, ylab = "", type = "l",
         ylim = c(min(fit[[type]]), max(fit[[type]])), lwd = 2)}

  if (type == "ALL"){
    par(mfrow = c(2,3), cex.main = 2, cex.lab = 1.5, cex.axis = 1.5)
    for (t in 1:length(param_names))
      plot(ts(fit[[param_names[t]]],start = (xrange[[t]])[1]), main = smain[[t]], xlab = sxlab[[t]], ylab = "", type = "l",
           ylim = c(min(fit[[param_names[t]]]), max(fit[[param_names[t]]])), lwd = 2)
    par(mfrow=c(1,1))}}
