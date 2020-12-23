#' Get the multipopulation mortality dataset
#'
#' @description This function returns the country-specific and combined `death counts` and `exposures` for a selection
#'   of countries from the Human Mortality Database (\href{https://www.mortality.org/}{HMD}) and/or
#'   \href{https://ec.europa.eu/eurostat/}{Eurostat} and for both males and females.
#'
#' @param xv Numeric. The vector of ages.
#' @param yv The vector of years.
#' @param yv_spec The vector of years for the country of interest.
#' @param countries The vector of countries.
#' @param country_spec The country of interest.
#' @param username The username of your HMD account.
#' @param password The password of your HMD acccount.
#'
#' @details This function downloads and/or calculates the death counts and exposures for an
#'  age range \code{xv}, a time range \code{yv} and for a specified set of countries \code{countries}.
#'  The main data source is HMD. The second one is Eurostat for the years for which there is no
#'  data available at the HMD. From the HMD database this functions downloads the tables
#'  `Deaths` and `Exposure to risk` in 1x1 format. In the case of Germany (\code{DE}), we consider
#'  West-Germany instead of whole Germany for the period before the German unification
#'  (...-1989). At the Eurostat database we download the databases `demo_mager` (period deaths),
#'  `demo_magec` (cohort deaths) and `demo_pjan` (population size). Period exposures are
#'  calculated from these datasets according to the HMD protocol.
#'
#'  The argument \code{countries} should contain the user code labels of your countries of
#'  interest. The available countries and their corresponding `user codes` are given in the dataset
#'  \code{\link{country_codes}}, provided in this package.
#'
#'  Further, you have to register on \href{https://www.mortality.org/}{HMD} and provide
#'  the function with your \code{username} and \code{password} such that the function can
#'  automatically download from HMD.
#'
#'  Remark that the year ranges \code{yv} and \code{yv_spec} should not be the same. Typically
#'  a longer year range for yv is chosen and yv_spec consists of some more recent years.
#'  When we download data for your country of interest (\code{country_spec}), we consider
#'  the year range \code{min(yv,yv_spec):tail(yv_spec,1)} since the older data is also used for estimating
#'  the common mortality trend.
#'
#' @return A list containing:
#'   \itemize{
#'   \item $M (male data), $F (female data)
#'   \item $UNI (individual mortality data: $BE (Belgium), $NL (Netherlands), $AT (Austria), ...), $ALL (combined data)
#'   \item $dtx (deaths), $etx (exposures), $wa (weights).
#'   }
#'
#' @examples
#' \dontrun{
#' xv <- 0:90
#' yv <- 1970:2015
#' yv_spec <- 1985:2018
#' countries <- c("BE", "IT", "UK", "DE")
#' country_spec <- "BE"
#' username <- ""
#' password <- ""
#' data <- get_mortality_data(xv, yv, yv_spec, countries, country_spec, username, password)
#' }
#'
#' @importFrom dplyr filter
#' @importFrom curl has_internet
#' @importFrom RCurl getURL
#' @importFrom utils read.table tail
#'
#' @export

get_mortality_data <-  function(xv, yv, yv_spec, countries, country_spec, username, password){

  # Only continue with an internet connection
  if(! has_internet())
    stop(paste0("No internet connection. You need an internet connection to download mortality ",
                "data from the internet."))

  # Bulk download of Eurostat mortality data
  message("Downloading mortality data from Eurostat: deaths, cohorts and population")
  deaths_eurostat  <- get_data_eurostat("demo_magec")
  cohorts_eurostat <- get_data_eurostat("demo_mager")
  pop_eurostat     <- get_data_eurostat("demo_pjan")

  # Only continue with valid country labels
  df_country <- get_country_codes(data_magec = deaths_eurostat, data_mager = cohorts_eurostat,
                                  data_pjan = pop_eurostat)
  if(! all(is.element(countries, df_country[, "User_code"])))
    stop(paste0("Invalid countries or country labels: ", paste(countries[which(! countries %in% df_country[,"User_code"])], collapse = ', '),
                ". ", "Please consult the function 'get_country_codes()' to take valid countries ",
                "or the correct user-labels for your chosen countries."))

  # Only continue with numeric yv and yv_spec
  yv     <- as.numeric(yv)
  yv_spec <- as.numeric(yv_spec)

  # Only continue when country_spec is included in countries
  if(! country_spec %in% countries)
    stop(paste0("You need to include 'country_spec' into the 'countries' argument."))

  # Only continue with valid yv and yv_spec
  if(! (all(diff(yv) == 1) & all(diff(yv_spec) == 1) & all(diff(xv)) == 1))
    stop(paste0("There may not be any missing years in the time ranges 'yv' and 'yv_spec'",
                "or missing ages in the age range 'xv'."))

  # Check valid time range
  sub_df   <- subset(df_country, User_code %in% countries)
  ind_spec <- which(sub_df[,"User_code"] == country_spec)
  first_y  <- max(sub_df[-ind_spec,"FY"])
  last_y   <- min(sub_df[-ind_spec,"EY"])
  begin_y  <- min(yv)
  end_y    <- max(yv)
  b_yspec  <- min(yv_spec)
  e_yspec  <- max(yv_spec)
  if(begin_y < first_y | end_y > last_y)
    stop(paste0("The age range yv should be within the range: ", first_y,"-",last_y,". ",
                "Have a look at the function 'get_country_codes()'."))

  if(b_yspec < sub_df[ind_spec,"FY"] | e_yspec > sub_df[ind_spec,"EY"])
    stop(paste0("The age range yv_spec should be within the range: ", sub_df[ind_spec,"FY"],"-",sub_df[ind_spec,"EY"],". ",
                "Have a look at the function 'get_country_codes()'."))

  if(xv[1] != 0)
    stop("The age range should start with the age of 0.")


  # Subset about information about the predefined countries
  df_country <- subset(df_country, User_code %in% countries)

  # Mortality data must be filtered out
  deaths_eurostat  <- dplyr::filter(deaths_eurostat, age %in% 0:99 &
                                      geo %in% df_country[,"User_code"] & sex %in% c("M","F"))
  cohorts_eurostat <- dplyr::filter(cohorts_eurostat, age %in% 0:99 &
                                      geo %in% df_country[,"User_code"] & sex %in% c("M","F"))
  pop_eurostat     <- dplyr::filter(pop_eurostat, age %in% 0:99 &
                                      geo %in% df_country[,"User_code"] & sex %in% c("M","F"))

  # Get the datasets ordered by age
  deaths_eurostat  <- deaths_eurostat[order(deaths_eurostat$age), ]
  cohorts_eurostat <- cohorts_eurostat[order(cohorts_eurostat$age), ]
  pop_eurostat     <- pop_eurostat[order(pop_eurostat$age), ]

  # Create empty lists to store results for males and females later on
  lijst_M <- list()
  lijst_F <- list()
  yv0     <- yv
  yv_spec0 <- yv_spec

  # Link username and password to download from HMD
  userpwd <- paste(username, ":", password, sep = "")

  # Iterating over the different countries
  message("Iterating over the different countries:")
  for (c in countries){
    # Period
    if (c == country_spec) yv <- yv_spec0 else yv <- yv0

    ### Get the HMD and Eurostat country label
    c_hmd <- subset(df_country, User_code %in% c)[,"HMD_code"]
    c_eur <- subset(df_country, User_code %in% c)[,"Eurostat_code"]

    end_year <- min(yv) - 1

    if(!is.na(c_hmd)){
      ### Deaths from HMD
      path_deaths <- paste0("https://www.mortality.org/hmd/", c_hmd, "/STATS/", "Deaths_1x1.txt")
      txt         <- getURL(path_deaths, userpwd = userpwd)
      con         <- textConnection(txt)
      deaths_hmd  <- data.frame(try(read.table(con, skip = 2, header = TRUE, na.strings = "."), TRUE))
      close(con)
      .check_file(deaths_hmd)
      start_year  <- deaths_hmd$Year[1]; end_year <- tail(deaths_hmd$Year,1)

      # Put deaths data in right format for males/females
      matrix_deaths_M <- matrix(deaths_hmd$Male, byrow = T, nrow = length(start_year:end_year))
      matrix_deaths_F <- matrix(deaths_hmd$Female, byrow = T, nrow = length(start_year:end_year))
      dimnames(matrix_deaths_M) = dimnames(matrix_deaths_F) <- list(start_year:end_year,
                                                                    0:(ncol(matrix_deaths_M)-1))
      matrix_deaths_M <- matrix_deaths_M[, as.character(xv)]
      matrix_deaths_F <- matrix_deaths_F[, as.character(xv)]

      ### Exposures from HMD
      path_exp <- paste("https://www.mortality.org/hmd/", c_hmd, "/STATS/", "Exposures_1x1.txt", sep = "")
      txt <- getURL(path_exp, userpwd = userpwd)
      con <- textConnection(txt)
      exp_hmd <- data.frame(try(read.table(con, skip = 2, header = TRUE, na.strings = "."), TRUE))
      close(con)
      .check_file(exp_hmd)

      # Put exposure data in right format for males/females
      matrix_exp_M <- matrix(exp_hmd$Male, byrow = T, nrow = length(start_year:end_year))
      matrix_exp_F <- matrix(exp_hmd$Female, byrow = T, nrow = length(start_year:end_year))
      dimnames(matrix_exp_M) = dimnames(matrix_exp_F) <- list(start_year:end_year,
                                                              0:(ncol(matrix_exp_M)-1))
      matrix_exp_M <- matrix_exp_M[, as.character(xv)]
      matrix_exp_F <- matrix_exp_F[, as.character(xv)]

      # In the case of Germany -> Include data from West-Germany before the German unification (1990)
      if (c == "DE"){
        # Get deaths on the HMD from West-Germany
        path_deaths <- paste("https://www.mortality.org/hmd/", "DEUTW", "/STATS/", "Deaths_1x1.txt", sep = "")
        txt         <- getURL(path_deaths, userpwd = userpwd)
        con         <- textConnection(txt)
        deaths_hmd  <- data.frame(try(read.table(con, skip = 2, header = TRUE, na.strings = "."), TRUE))
        close(con)
        .check_file(deaths_hmd)
        start_year  <- deaths_hmd$Year[1]
        end_yr      <- tail(deaths_hmd$Year,1)

        # Put deaths data in right format for males and females
        deaths_hmd       <- dplyr::filter(deaths_hmd, Year %in% c(start_year:end_yr))
        matrix_wdeaths_M <- matrix(deaths_hmd$Male, byrow = T, nrow = length(start_year:end_yr))
        matrix_wdeaths_F <- matrix(deaths_hmd$Female, byrow = T, nrow = length(start_year:end_yr))
        dimnames(matrix_wdeaths_M) = dimnames(matrix_wdeaths_F) <- list(start_year:end_yr,
                                                                        0:(ncol(matrix_wdeaths_M)-1))
        matrix_deaths_M <- rbind(matrix_wdeaths_M[as.character(deaths_hmd$Year[1]:1989),
                                                  as.character(xv)], matrix_deaths_M)
        matrix_deaths_F <- rbind(matrix_wdeaths_F[as.character(deaths_hmd$Year[1]:1989),
                                                  as.character(xv)], matrix_deaths_F)

        # Get exposures on HMD from West-Germany
        path_exp <- paste("https://www.mortality.org/hmd/", "DEUTW", "/STATS/", "Exposures_1x1.txt", sep = "")
        txt <- getURL(path_exp, userpwd = userpwd)
        con <- textConnection(txt)
        exp_hmd <- data.frame(try(read.table(con, skip = 2, header = TRUE, na.strings = "."), TRUE))
        close(con)
        .check_file(exp_hmd)

        #Put exposures data in right format for males/females
        matrix_wexp_M <- matrix(exp_hmd$Male, byrow = T, nrow = length(start_year:end_yr))
        matrix_wexp_F <- matrix(exp_hmd$Female, byrow = T, nrow = length(start_year:end_yr))
        dimnames(matrix_wexp_M) = dimnames(matrix_wexp_F) <- list(start_year:end_yr,
                                                                  0:(ncol(matrix_wexp_M)-1))
        matrix_exp_M <- rbind(matrix_wexp_M[as.character(exp_hmd$Year[1]:1989), as.character(xv)],
                              matrix_exp_M)
        matrix_exp_F <- rbind(matrix_wexp_F[as.character(exp_hmd$Year[1]:1989), as.character(xv)],
                              matrix_exp_F)
      }}

    if (end_year < tail(yv,1) & !is.na(c_eur)){

      ### Filtering Eurostat data (deaths, cohort deaths, population sizes)
      deaths_eurostat_M  <- dplyr::filter(deaths_eurostat, time %in% c((end_year+1):tail(yv,1)) &
                                            geo == c_eur, sex == "M")
      deaths_eurostat_F  <- dplyr::filter(deaths_eurostat, time %in% c((end_year+1):tail(yv,1)) &
                                            geo == c_eur, sex == "F")
      cohorts_eurostat_M <- dplyr::filter(cohorts_eurostat, time %in% c((end_year+1):tail(yv,1)) &
                                            geo == c_eur, sex == "M")
      cohorts_eurostat_F <- dplyr::filter(cohorts_eurostat, time %in% c((end_year+1):tail(yv,1)) &
                                            geo == c_eur, sex == "F")
      pop_eurostat_M     <- dplyr::filter(pop_eurostat, time %in% c((end_year+1):tail(yv+1,1)) &
                                            geo == c_eur, sex == "M")
      pop_eurostat_F     <- dplyr::filter(pop_eurostat, time %in% c((end_year+1):tail(yv+1,1)) &
                                            geo == c_eur, sex == "F")

      ### Put Eurostat data in right format -> (t,x)-matrix
      age_seq <- sort(unname(unlist(unique(deaths_eurostat_M[,"age"]))))

      deaths_eurostat_M  <- .transform_euro(deaths_eurostat_M, age_seq)
      deaths_eurostat_F  <- .transform_euro(deaths_eurostat_F, age_seq)
      cohorts_eurostat_M <- .transform_euro(cohorts_eurostat_M, age_seq)
      cohorts_eurostat_F <- .transform_euro(cohorts_eurostat_F, age_seq)
      pop_eurostat_M     <- .transform_euro(pop_eurostat_M, age_seq)
      pop_eurostat_F     <- .transform_euro(pop_eurostat_F, age_seq)

      ### Calculating the exposures
      exp_eurostat_M <- exp_eurostat_F <- matrix(NA, nrow = length((end_year+1):tail(yv,1)),
                                                 ncol = length(xv))
      dimnames(exp_eurostat_F) = dimnames(exp_eurostat_M) <- list((end_year+1):tail(yv,1), xv)
      n <- nrow(exp_eurostat_M)

      # Age 0
      exp_eurostat_F[,"0"] <- 0.5 * (pop_eurostat_F[1:n, "0"] + pop_eurostat_F[2:(n+1), "0"]) +
        1/6*(cohorts_eurostat_F[1:n, "0"] - 0.5*cohorts_eurostat_F[1:n, "1"])
      exp_eurostat_M[,"0"] <- 0.5 * (pop_eurostat_M[1:n,"0"] + pop_eurostat_M[2:(n+1),"0"]) +
        1/6*(cohorts_eurostat_M[1:n,"0"] - 0.5*cohorts_eurostat_M[1:n,"1"])

      # Remaining ages of xv
      exp_eurostat_F[,2:length(xv)] <-
        0.5 * (pop_eurostat_F[1:n, 2:length(xv)] + pop_eurostat_F[2:(n+1), 2:length(xv)]) +
        1/6*(0.5*cohorts_eurostat_F[1:n, 2:length(xv)] - 0.5*cohorts_eurostat_F[1:n, 3:(length(xv)+1)])
      exp_eurostat_M[,2:length(xv)] <-
        0.5 * (pop_eurostat_M[1:n, 2:length(xv)] + pop_eurostat_M[2:(n+1), 2:length(xv)]) +
        1/6*(0.5*cohorts_eurostat_M[1:n, 2:length(xv)] - 0.5*cohorts_eurostat_M[1:n, 3:(length(xv)+1)])

      ### Merge the Eurostat and HMD data
      matrix_deaths_M <- if(exists("matrix_deaths_M"))
        rbind(matrix_deaths_M, deaths_eurostat_M[,(xv+1)]) else
          deaths_eurostat_M[,(xv+1)]
      matrix_deaths_F <- if(exists("matrix_deaths_F"))
        rbind(matrix_deaths_F, deaths_eurostat_F[,(xv+1)]) else
          deaths_eurostat_F[,(xv+1)]
      matrix_exp_M <- if(exists("matrix_exp_M"))
        rbind(matrix_exp_M, exp_eurostat_M) else
          exp_eurostat_M
      matrix_exp_F <- if(exists("matrix_exp_F"))
        rbind(matrix_exp_F, exp_eurostat_F) else
          exp_eurostat_F

      dimnames(matrix_deaths_M) = dimnames(matrix_deaths_F) = dimnames(matrix_exp_M) =
        dimnames(matrix_exp_F) <- list(as.numeric(rownames(matrix_deaths_M)[1]):tail(yv,1), xv)}

    ## Adapt range
    matrix_deaths_M <- matrix_deaths_M[as.character(min(c(yv0,yv)):tail(yv,1)),]
    matrix_deaths_F <- matrix_deaths_F[as.character(min(c(yv0,yv)):tail(yv,1)),]
    matrix_exp_M    <- matrix_exp_M[as.character(min(c(yv0,yv)):tail(yv,1)),]
    matrix_exp_F    <- matrix_exp_F[as.character(min(c(yv0,yv)):tail(yv,1)),]


    ### Define weight matrix
    wa_M = wa_F <- matrix(data = 1, ncol = dim(matrix_deaths_M)[2], nrow = dim(matrix_deaths_M)[1])
    ind_M <- which(is.na(matrix_deaths_M), arr.ind = TRUE)
    ind_F <- which(is.na(matrix_deaths_F), arr.ind = TRUE)
    wa_M[ind_M] <- 0
    wa_F[ind_F] <- 0
    dimnames(wa_M) = dimnames(wa_F) <- dimnames(matrix_deaths_M)

    ### Avoid problem of taking logarithm of zero
    matrix_deaths_M[matrix_deaths_M == 0] <- 0.01
    matrix_deaths_F[matrix_deaths_F == 0] <- 0.01

    ### Put results in list object
    lijst_M[[c]] <- list('dtx' = matrix_deaths_M, 'etx' = matrix_exp_M, 'wa' = wa_M)
    lijst_F[[c]] <- list('dtx' = matrix_deaths_F, 'etx' = matrix_exp_F, 'wa' = wa_F)
    message(c)}

  ### prepare data all populations together:
  # dtxALL = sum of all dtx(i)
  # etxALL = sum of all etx(i)
  # waALL = 1 if etx(i) and dtx(i) are available for all populations (i), = 0 otherwise

  dtxALL_M = dtxALL_F <- matrix(0, nrow = length(yv0), ncol = length(xv))	# initialize matrix
  etxALL_M = etxALL_F <- matrix(0, nrow = length(yv0), ncol = length(xv))	# initialize matrix
  waALL_M  = waALL_F  <- matrix(1, nrow = length(yv0), ncol = length(xv))	# initialize matrix

  for(i in 1:length(countries)){
    dtxALL_M <- dtxALL_M + lijst_M[[i]]$dtx[as.character(yv0), ]
    dtxALL_F <- dtxALL_F + lijst_F[[i]]$dtx[as.character(yv0), ]
    etxALL_M <- etxALL_M + lijst_M[[i]]$etx[as.character(yv0), ]
    etxALL_F <- etxALL_F + lijst_F[[i]]$etx[as.character(yv0), ]
    waALL_M  <- (waALL_M == TRUE) & (lijst_M[[i]]$wa[as.character(yv0), ] == TRUE)
    waALL_F  <- (waALL_F == TRUE) & (lijst_F[[i]]$wa[as.character(yv0), ] == TRUE)}

  dataALL_M = list(dtx = dtxALL_M, etx = etxALL_M, wa = waALL_M*1)
  dataALL_F = list(dtx = dtxALL_F, etx = etxALL_F, wa = waALL_F*1)

  lst_M <- list(UNI = lijst_M, ALL = dataALL_M)
  lst_F <- list(UNI = lijst_F, ALL = dataALL_F)
  list("M" = lst_M, "F" = lst_F)
}


#' @keywords internal
.fac_to_num <- function(x) as.numeric(as.character(x))

#' @keywords internal
.check_file <- function(x) {
  if (class(x) == "try-error")
    stop("File not found at www.mortality.org")}

#' @keywords internal
.transform_euro <- function(x,xv){
  x <- x[order(x$time),]
  X <- matrix(x$values, byrow = T, ncol = length(xv))
  colnames(X) <- xv
  rownames(X) <- unique(x$time)
  X
}

