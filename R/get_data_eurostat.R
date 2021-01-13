#' Bulk download from Eurostat of mortality data
#'
#' @description This function downloads mortality data from
#'   \href{https://ec.europa.eu/eurostat/}{Eurostat}. Little preprocessing is done: only ages
#'   up to 99 are considered and data is put in an other format.
#'
#' @param code The code name of the dataset to download.
#'
#' @return A dataframe containing the available country names, the corresponding HMD, Eurostat
#'   and user country labels as well as the avaible year range for each country.
#'
#' @details The avaible code names are \code{demo_magec} (period deaths), \code{demo_mager}
#'   (cohort deaths) and \code{demo_pjan} (population sizes).
#'
#' @examples data <- get_data_eurostat("demo_magec")
#'
#' @importFrom dplyr %>% left_join pull filter arrange
#' @importFrom RCurl url.exists
#' @importFrom plyr mapvalues
#' @importFrom tidyr separate gather
#' @importFrom readr read_tsv cols col_character
#' @importFrom data.table as.data.table
#'
#' @export


get_data_eurostat <- function(code){

  # Only continue with valid code name
  if(! code %in% c("demo_magec","demo_mager","demo_pjan"))
    stop("The code name must be either 'demo_magec', 'demo_mager' or 'demo_pjan'.")

  # URL of file
  base <- "https://ec.europa.eu/eurostat/"
  ext  <- paste0("estat-navtree-portlet-prod/BulkDownloadListing?sort=1&file=data%2F",
                 code,".tsv.gz")
  url  <- paste0(base,ext)

  # Check if url exists
  if(! url.exists(url))
    stop(paste0("The URL ", url, " does not exist anymore. Please contact the maintainer",
                " of this package to make the necessary changes."))

  # Download URL
  dest_file <- tempfile()
  on.exit(unlink(dest_file))
  utils::download.file(url, dest_file, quiet = TRUE)
  data <- read_tsv(gzfile(dest_file), na = ":", col_types =
                            cols(.default = col_character()))

  # Make the dataset ready
  coln  <- unlist(strsplit(colnames(data)[1],"[\\,]"))
  coln  <- coln[-length(coln)]
  data  <- separate(data, col = colnames(data)[1], into = coln, sep = ",", convert = FALSE)
  data  <- gather(data, "time", "values", -(1:length(coln)))
  #data  <- dplyr::filter(data, !is.na(values))

  # Numeric/factor columns
  data$values <- gsub("[^0-9.-]+", "", data$values)
  data$values <- as.numeric(data$values)
  data$time   <- as.numeric(data$time)
  data$sex    <- factor(data$sex, levels = c("M", "F", "T"))

  # Age
  age_select <- c("Y_LT1", paste0("Y",1:99))
  data       <- subset(data, age %in% age_select)
  data$age   <- as.numeric(mapvalues(data$age, age_select, 0:99))
  data       <- as.data.table(data)
  data       <- data.frame(data[order(data$time, data$age, data$sex, data$geo),])
  data
}

