#' Get the available countries within the package
#'
#' @description This function shows the available countries that can be used in the
#'   multipopulation mortality model. Further, their user codes are given as well as the
#'   available year range for each country. Information is substracted from the data sources
#'   \href{https://www.mortality.org/}{HMD} and \href{https://ec.europa.eu/eurostat/}{Eurostat}.
#'
#' @param data_magec The dataset of the period death counts on Eurostat. If \code{NULL},
#' this dataset is downloaded from \href{https://ec.europa.eu/eurostat/}{Eurostat}.
#' @param data_mager The dataset of the cohort death counts on Eurostat. If \code{NULL},
#' this dataset is downloaded from \href{https://ec.europa.eu/eurostat/}{Eurostat}.
#' @param data_pjan The dataset of the population numbers at 1 January of each year on Eurostat. If \code{NULL},
#' this dataset is downloaded from \href{https://ec.europa.eu/eurostat/}{Eurostat}.
#'
#' @return A dataframe containing the available country names, the corresponding HMD, Eurostat
#'   and user country labels as well as the avaible year range for each country.
#'
#' @details This function can be used to find your desired country set and year range to construct
#'   the multipopulation mortality model.
#'
#'   A special case: Germany. Here before 1990 (German unification), only data from West-Germany
#'   is considered.
#'
#' @examples df_country <- get_country_codes()
#'
#' @importFrom utils read.csv tail
#' @importFrom dplyr %>% left_join
#' @importFrom stats na.omit setNames
#'
#' @export


get_country_codes <- function(data_magec = NULL, data_mager = NULL, data_pjan = NULL){
  # Get the possible countries at the HMD
  file_hmd <- MultiMoMo::countries

  # Special case: Germany
  wGer <- file_hmd[which(file_hmd[,"Country"] == "Germany: West Germany"),]
  eGer <- file_hmd[which(file_hmd[,"Country"] == "Germany: East Germany"),]

  # Only work with Total population of country (except Germany)
  ind_tp   <- duplicated(sapply(file_hmd[, "HMD"], function(x) substr(x, 1,3)))
  file_hmd <- file_hmd[!ind_tp,]

  # Alpha 3 code
  file_hmd[, "Alpha_3"] <- sapply(file_hmd[,"HMD"], function(x) substr(x, 1, 3))

  # Get Alpha 2 code
  ISO_df <- ISOcodes::ISO_3166_1
  ISO_df <- ISO_df[,c("Alpha_2", "Alpha_3")]
  file   <- file_hmd %>% left_join(ISO_df, by = c("Alpha_3" = "Alpha_3"))

  # Change definition of Germany
  old_ger <- "Germany: Total population"
  new_ger <- "Germany (until 1990 former FRG)"
  ind_ger <- which(file[,"Country"] == old_ger)
  file[ind_ger, c("Country", "FY")] <- c(new_ger, wGer[,"FY"])

  # Change codes a bit (difference with european commision)
  old <- c("GR", "GB")
  new <- c("EL", "UK")
  file[,"Alpha_2"] <- plyr::mapvalues(file[,"Alpha_2"], old, new)

  # Years HMD
  years_hmd <- lapply(1:nrow(file), function(x) {
    gap <- file[x,"Gap_FY"]:file[x,"Gap_EY"]
    rg  <- file[x,"FY"]:file[x,"EY"]
    if(all(gap != 0))
      rg <- rg[-which(rg %in% gap)]
    rg})
  names(years_hmd) <- file[,"Alpha_2"]

  # Eurostat data
  dat1 <- if(is.null(data_magec)) get_data_eurostat("demo_magec") else data_magec
  dat2 <- if(is.null(data_mager)) get_data_eurostat("demo_mager") else data_mager
  dat3 <- if(is.null(data_pjan))  get_data_eurostat("demo_pjan")  else data_pjan

  # Remove rows with NaN values
  dat1 <- subset(dat1, values >= 0)
  dat2 <- subset(dat2, values >= 0)
  dat3 <- subset(dat3, values >= 0)

  # Country codes
  codes <- intersect(intersect(unique(dat1$geo), unique(dat2$geo)), unique(dat3$geo))
  codes <- sapply(codes, function(x) ifelse(nchar(x) != 2 & x !="DE_TOT", NA, x), USE.NAMES = F)
  codes <- sort(na.omit(codes))

  # Available years eurostat
  years_euro <- lapply(codes, function(x) {
    c_magec <- unique(subset(dat1, geo %in% x)$time)
    c_mager <- unique(subset(dat2, geo %in% x)$time)
    c_pjan  <- unique(subset(dat3, geo %in% x)$time)
    c_int   <- Reduce(intersect, list(c_magec, c_mager))

    check   <- sapply(c_int, function(x) x %in% c_pjan & (x + 1) %in% c_pjan)
    c_int   <- c_int[check]
    c_int
  })
  names(years_euro) <- codes

  # Combine years hmd and eurostat
  keys    <- unique(c(names(years_hmd), names(years_euro)))
  comb    <- function(x, y) sort(unique(c(x,y)))
  years_c <- setNames(mapply(comb, years_hmd[keys], years_euro[keys]), keys)

  # Take longest consecutive sequence
  years_f <- lapply(years_c, function(x) {
    s <- split(x, cumsum(c(TRUE, diff(x) != 1)))
    l <- sapply(s, length)
    s <- s[[tail(which(l == max(l)),1)]]})

  # Make final summary
  ISO_df    <- ISOcodes::ISO_3166_1
  names_c   <- mapvalues(names(years_f), new, old)
  ind_m     <- match(names_c, ISO_df$Alpha_2)
  countries <- ISO_df$Name[ind_m]

  # Some adaptations
  cc <- c("DE", "DE_TOT", "FX", "XK")
  countries[which(names_c %in% cc)] <- c("Germany (until 1990 former FRG)", "Germany including former GDR",
                                         "France (Metropolitan)", "Kosovo")
  names_c <- mapvalues(names_c, old, new)

  in_hmd  <- names_c %in% names(years_hmd)
  in_euro <- names_c %in% names(years_euro)
  range   <- lapply(names_c, function(x) range(years_f[[x]]))
  FY      <- sapply(range, function(x) x[1])
  EY      <- sapply(range, function(x) x[2])

  ma   <- match(names_c, file[,"Alpha_2"])
  hmd  <- file[ma, "HMD"]
  euro <- sapply(1:length(names_c), function(x) ifelse(in_euro[x], names_c[x], NA))


  # Table
  Table <- data.frame("Country" = countries, "in_HMD" = in_hmd, "in_Eurostat" = in_euro,
                      "HMD_code" = hmd, "Eurostat_code" = euro, "FY" = FY, "EY" = EY,
                      "User_code" = names_c)

  rank  <- order(Table[,"Country"])
  Table <- Table[rank,]
  rownames(Table) <- NULL

  # Adapt definition of France -> HMD reports metropoliton france from 1946 on
  ind.fra <- which(Table$Country == "France")
  ind.metrofra <- which(Table$Country == "France (Metropolitan)")
  Table[ind.fra,] <- data.frame("France", FALSE, TRUE, NA, "FR", range(years_euro$FR)[1],
                                range(years_euro$FR)[2], "FR_TOT")
  Table[ind.metrofra,] <- data.frame("France (Metropolitan)", TRUE, TRUE, "FRATNP","FX",
                                     range(c(years_hmd$FR, years_euro$FX))[1],
                                     range(c(years_hmd$FR, years_euro$FX))[2], "FR")

  Table

}

