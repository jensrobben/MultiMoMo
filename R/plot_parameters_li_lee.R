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
#' plot_parameters_li_lee(xv, yvSPEC, fit_M, "Male", country_spec, "NR", "k.t")
#'
#' @import ggplot2
#' @importFrom plyr mapvalues
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
    df <- data.frame("x" = xrange, "y" = fit[[type]])
    p <- ggplot(df, aes(x = x, y = y)) + geom_line(col = "black", size = 1.1) + xlab(sxlab) +
      ylab("") + theme_bw(base_size = 20) + ggtitle(smain)
    print(p)}

  if (type == "ALL"){
    g_legend <- function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)}

    plots <- list()

    for (t in 1:length(param_names)){
      df <- data.frame("x" = xrange[[t]], "y" = fit[[param_names[[t]]]])
      plots[[t]] <- ggplot(df, aes(x = x, y = y)) + geom_line(col = "black", size = 1.05) +
        xlab(sxlab[t]) + ylab("") + theme_bw(base_size = 15) + ggtitle(smain[[t]])}

    ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]],
              ncol=3, nrow=2)}
}
