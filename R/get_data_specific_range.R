#' Adapt the year and/or age range of the multi-population mortality dataset
#'
#' @description This function adapts the year and/or age range of the mortality dataset,
#' as obtained from the function \code{\link[MultiMoMo]{get_mortality_data}},
#' for each country in the multi-population model.
#'
#' @param xv The new vector of ages.
#' @param yv The new vector of years.
#' @param yv_spec The new vector of years for the country of interest.
#' @param data The multi-population mortality dataset.
#' @param country The vector of countries.
#' @param country_spec The country of interest.
#'
#' @details The argument \code{data} must be (in the same form as) the output
#' from the function \code{\link[MultiMoMo]{get_mortality_data}}. See
#' \code{\link[MultiMoMo]{european_mortality_data}} for such an example mortality
#' dataset.
#'
#' @return A list containing the following updated objects:
#'   \itemize{
#'   \item $M (male data), $F (female data)
#'   \item $UNI (individual mortality data: $BE (Belgium), $NL (Netherlands), $AT (Austria), ...), $ALL (combined data)
#'   \item $dtx (deaths), $etx (exposures), $wa (weights).
#'   }
#'
#' @examples
#' xv       <- 0:90
#' yv       <- 1988:2018
#' yv_spec  <- 1988:2018
#' data     <- MultiMoMo::european_mortality_data
#' country  <- names(data$M$UNI)
#' country_spec <- "BE"
#'
#' data_new <- get_data_specific_range(xv, yv, yv_spec, data, country, country_spec)
#'
#' @export



get_data_specific_range  <- function(xv, yv, yv_spec, data, country, country_spec){
  # The number of countries
  nc <- length(country)

  # Object to store results
  data_sep <- rep(list(NA),2)
  names(data_sep) <- names(data)

  # Check vector country
  if(! all(country %in% names(data[[1]][[1]])))
    stop("There is at least one country not contained in the mortality dataset.")

  # The actual code
  for(s in names(data)){
    for(type in names(data[[s]])){
      if(type == "UNI")
        for(c in country){
          if(c == country_spec) yvv <- yv_spec else yvv <- yv
            data[[s]][[type]][[c]][["dtx"]] <- data[[s]][[type]][[c]][["dtx"]][
              as.character(yvv), as.character(xv)]
            data[[s]][[type]][[c]][["etx"]] <- data[[s]][[type]][[c]][["etx"]][
              as.character(yvv), as.character(xv)]
            data[[s]][[type]][[c]][["wa"]] <- data[[s]][[type]][[c]][["wa"]][
              as.character(yvv), as.character(xv)]}
      if(type == "ALL"){
        data[[s]][[type]][["dtx"]] <- data[[s]][[type]][["dtx"]][
          as.character(yv), as.character(xv)]
        data[[s]][[type]][["etx"]] <- data[[s]][[type]][["etx"]][
          as.character(yv), as.character(xv)]
        data[[s]][[type]][["wa"]] <- data[[s]][[type]][["wa"]][
          as.character(yv), as.character(xv)]
      }
    }
  }
  data
}


