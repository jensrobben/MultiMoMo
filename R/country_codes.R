#' User, HMD and Eurostat codes
#'
#' @description Dataset including the user, HMD and Eurostat codes for the available countries. The country
#'   names are available as well.
#'
#' @docType data
#'
#' @usage data(country_codes)
#'
#' @keywords datasets
#'
#' @format A data frame with 59 rows and 4 variables:
#' \describe{
#' \item{Country}{The name of the country}
#' \item{User_code}{The user code name of the country}
#' \item{HMD_code}{The HMD code name of the country}
#' \item{Eurostat_code}{The Eurostat code name of the country}
#' }
#'
#' @source
#'   \url{https://ec.europa.eu/eurostat/} \cr
#'
#'   \url{https://www.mortality.org/}
#'
#' @examples
#'   data(country_codes)
#'   user_code <- attr(country_codes, "User_code")
#'   hmd_code  <- attr(country_codes, "HMD_code")
#'   euro_code <- attr(country_codes, "Eurostat_code")
#'
#'
"country_codes"
