#' The evolution of the life expectancy of a newborn
#'
#' @description Graph of the evolution of the life expectancy of a male and female newborn.
#'  Data is downloaded from the Human Mortality Database (\href{https://www.mortality.org/}{HMD}).
#'  Therefore, only countries available in the HMD are visualized in this graph.
#'
#' @param Country The set of countries.
#' @param username The username of your HMD account.
#' @param password The password of your HMD account.
#'
#' @source \\href{https://www.mortality.org/}
#'
#' @keywords ggplots
#'
#' @import ggplot2
#' @importFrom grDevices colors
#' @importFrom RCurl getURL
#' @importFrom curl has_internet
#' @importFrom utils read.table
#' @importFrom cowplot background_grid
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom ggpubr ggarrange
#'
#' @examples
#' \dontrun{
#'   Country   <- c("BE", "IT", "DE", "FR", "NL")
#'   username  <- ""
#'   password  <- ""
#'   plot_e0_hmd(Country, username, password)}
#'
#' @export


plot_e0_hmd <- function(Country, username, password){

  # Only Countries available in the HMD can be used
  df      <- MultiMoMo::country_codes
  ind     <- which(df[,"User_code"] %in% Country)
  Country <- df[ind,"HMD_code"]

  # A good colour range is essential
  nc           <- length(Country)
  colour_field <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "gold1",
    "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
    "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
    "yellow3", "darkorange4", "brown")
  if(nc > length(colour_field))
    stop("Take a maximum of 25 countries.")
  colours      <- colour_field[seq_along(Country)]

  # Password is needed to download data from HMD
  userpwd <- paste(username, ":", password, sep = "")

  # Download life expectancy data from HMD to visualize
  data_e0 <- list()
  message("Downloading data from HMD for each country")
  for(i in 1:nc){
    url_c   <- paste0("https://former.mortality.org/hmd/", Country[i], "/STATS/", "E0per.txt")
    txt     <- getURL(url_c, userpwd = userpwd)
    con     <- textConnection(txt)
    e0_hmd  <- data.frame(try(read.table(con, skip = 2, header = TRUE, na.strings = "."), TRUE))
    close(con)
    .check_file(e0_hmd)

    if(Country[i] == "DEUTNP"){
      url_c   <- paste0("https://www.mortality.org/hmd/", "DEUTW", "/STATS/", "E0per.txt")
      txt     <- getURL(url_c, userpwd = userpwd)
      con     <- textConnection(txt)
      e0_hmd0 <- data.frame(try(read.table(con, skip = 2, header = TRUE, na.strings = "."), TRUE))
      e0_hmd0 <- subset(e0_hmd0, Year < 1990)
      e0_hmd  <- rbind(e0_hmd0, e0_hmd)
      close(con)
    }

    matrix  <- subset(e0_hmd, Year >= 1920)
    data_e0 <- append(data_e0, list(matrix))
  }

  names(data_e0)  <-  Country

  # Make the figure
  x         <- unlist(unname(sapply(data_e0, function(x) x[,"Year"])))
  Countries <- unlist(unname(sapply(1:length(Country),
                                    function(x) rep(Country[x], nrow(data_e0[[x]])))))
  value_F   <- unlist(unname(sapply(data_e0, function(x) x[,"Female"])))
  value_M   <- unlist(unname(sapply(data_e0, function(x) x[,"Male"])))
  df_F      <- data.frame(x, Countries, value_F)
  df_M      <- data.frame(x, Countries, value_M)

  # Same y range
  ymin <- min(value_M, value_F)
  ymax <- max(value_M, value_F)

  d_F <- ggplot(df_F, aes(x = x, y = value_F, group = Countries, colour = Countries)) +
    ylim(ymin, ymax) + geom_line() + theme_bw(base_size = 20) + ggtitle("Females") +
    scale_color_manual(name = "", breaks = Country, values = colours, labels = df[ind, "User_code"]) +
    background_grid(major = "xy", minor = "none") + xlab("Age") + ylab(expression(e["0,t"])) +
    guides(col = guide_legend(nrow = 1, override.aes = list(size = 1.5))) +
    theme(legend.position = 'bottom')


  d_M <- ggplot(df_M, aes(x = x, y = value_M, group = Countries, colour = Countries) ) +
    ylim(ymin, ymax) + geom_line() + theme_bw(base_size = 20) + ggtitle("Males") +
    scale_color_manual(name = "", breaks = Country, values = colours,  labels = df[ind,"User_code"]) +
    background_grid(major = "xy", minor = "none") + xlab("Age") + ylab(expression(e["0,t"])) +
    theme(legend.position = "bottom") +
    guides(col = guide_legend(nrow = 1, override.aes = list(size = 1.5)))

  g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

  mylegend <- g_legend(d_M)

  ggarrange(d_M, d_F, ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom",
                        legend.grob = mylegend)

}
