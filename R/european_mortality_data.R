#' Example of a mortality dataset
#'
#' @description Dataset (matrix) including male and female deaths, exposures and weights for 14
#'   different European countries (Austria, Belgium, Denmark, Finland, France, Germany, Iceland,
#'   Ireland, Luxembourg, Netherlands, Norway, Sweden, Switzerland, United Kindom). Data is
#'   recorded for an age range \code{xv = 0:90} and year range \code{yv = yvSPEC = 1970:2018}.
#'
#' @docType data
#'
#' @usage data(european_mortality_data)
#'
#' @keywords datasets
#'
#' @format A list object containing:
#' \itemize{
#' \item $M (male data) and $F (female data).
#' \item $UNI (individual data per country: $AT, $BE, $DK, ...) and $ALL (combined European data)
#' \item $dtx (deaths), $etx (exposures) and $wa (weights). Matrices of size 49 x 91.
#' }
#'
#' @examples
#'   data(european_mortality_data)
#'
#'   # Male deaths of Belgium
#'   deaths_M_BE <- european_mortality_data$M$UNI$BE$dtx
#'
#'   # Female exposures (total combined European countries)
#'   exp_F       <- european_mortality_data$F$ALL$etx
#'
#'   # Male weights of Iceland
#'   wa_M_IS     <- european_mortality_data$M$UNI$IS$wa
#'
#'
"european_mortality_data"
