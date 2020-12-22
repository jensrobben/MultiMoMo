#' Graph of the combined male and female deaths and exposures
#'
#' @description This function constructs a graph (ggplot) of the combined male and female deaths
#'   over the different countries (\code{Country}). The same is done for the combined exposures.
#'   On top of that a dotted line is added to show the deaths and exposures for the country
#'   of interest (\code{CountrySPEC}).
#'
#' @param Country The vector of countries.
#' @param CountrySPEC The country of interest.
#' @param data_M The male mortality data, containing the deaths and exposures for each country.
#' @param data_F The female mortality data, containing the deaths and exposures for each country.
#'
#' @details The input parameters \code{data_M} and \code{data_F} must be (in the same format as)
#'   the output of the function \code{get_mortality_data()}.
#'
#' @keywords ggplots
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom ggpubr ggarrange
#' @importFrom tidyr gather
#'
#' @examples
#' xv          <- 0:90
#' yv          <- 1970:2018
#' yvSPEC      <- 1970:2018
#' username    <- "jens.robben@live.be"
#' password    <- "AtletieK1996"
#' Country     <- c("FR", "BE", "NL", "LU")
#' CountrySPEC <- "BE"
#' data        <- get_mortality_data(xv, yv, yvSPEC, Country, CountrySPEC, username, password)
#' data_M      <- data$M
#' data_F      <- data$F
#' plot_combined_deaths_exposures(Country, CountrySPEC, data_M, data_F)
#'
#'
#' @export

plot_combined_deaths_exposures <- function(Country, CountrySPEC, data_M, data_F){

  # Age range
  xv <- as.numeric(colnames(data_M[[1]][[1]][[1]]))

  # A good colour range is essential
  nc           <- length(Country)
  colour_field <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "gold1",
                    "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
                    "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise", "green1", "yellow4",
                    "yellow3", "darkorange4", "brown")
  if(nc > length(colour_field))
    stop("Take a maximum of 25 countries.")

  col_country       <- colour_field[seq_along(Country)]
  col_legend        <- col_country
  names(col_legend) <- Country

  # Data of deaths must be in specific format to plot
  comb_deaths_male   <- lapply(Country, function(x) colSums(data_M$UNI[[x]]$dtx, na.rm = FALSE, dims = 1))
  comb_deaths_female <- lapply(Country, function(x) colSums(data_F$UNI[[x]]$dtx, na.rm = FALSE, dims = 1))
  names(comb_deaths_female) = names(comb_deaths_male) <- Country
  comb_deaths_total  <- lapply(Country, function(x)
    matrix(comb_deaths_male[[x]] + comb_deaths_female[[x]], nrow = 1, ncol = length(xv)))

  dat           <- data.frame(do.call(rbind, comb_deaths_total))
  dimnames(dat) <- list(Country, seq(from = xv[1], to = max(xv)))
  dat[,"row"]   <- rownames(dat)
  dat2          <- dat %>% gather(key = "Age", value = "Deaths", -row)
  dat2[,"Age"]  <- as.numeric(dat2[,"Age"])

  lab           <- names(sort(rowSums(dat[,-which(colnames(dat) == "row")]), decreasing = F))
  val           <- col_legend[lab]
  dat2[,"row"]  <- factor(dat2[,"row"], levels = lab)

  # Plot combined deaths
  death_plot <- ggplot(dat2, aes(x = .data$Age, y = .data$Deaths, fill = factor(.data$row))) +
    geom_bar(stat = "identity", width = 1) + xlab("\nAge") +  ylab("Deaths\n") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Deaths") +
    scale_fill_manual(name = "", values = val, labels = lab) + theme_bw(base_size = 20) +
    scale_x_discrete(breaks = unique(c(seq(xv[1], max(xv) - 1, by = 10), max(xv) - 1))) +
    geom_point(data = filter(dat2, .data$row == CountrySPEC),
               aes(x = .data$Age, y = .data$Deaths), show.legend = F)


  # Data of exposures must be in specific format to plot
  comb_exp_male   <- lapply(Country, function(x) colSums(data_M$UNI[[x]]$etx, na.rm = FALSE, dims = 1))
  comb_exp_female <- lapply(Country, function(x) colSums(data_F$UNI[[x]]$etx, na.rm = FALSE, dims = 1))
  names(comb_exp_female) = names(comb_exp_male) <- Country
  comb_exp_total  <- lapply(Country, function(x)
    matrix(comb_exp_male[[x]] + comb_exp_female[[x]], nrow = 1, ncol = length(xv)))

  dat           <- data.frame(do.call(rbind, comb_exp_total))
  dimnames(dat) <- list(Country, seq(from = xv[1], to = max(xv)))
  dat[,"row"]   <- rownames(dat)
  dat2          <- dat %>% gather(key = "Age", value = "Exp", -row)
  dat2[,"Age"]  <- as.numeric(dat2[,"Age"])
  dat2[,"row"]  <- factor(dat2[,"row"], levels = lab)

  # Plot the combined exposures
  exposure_plot <- ggplot(dat2, aes(x = .data$Age, y = .data$Exp, fill = factor(.data$row))) +
    geom_bar(stat = "identity", width = 1) + xlab("\nAge") + ylab("Exposures\n") +
    scale_fill_manual(name = "", labels = lab, values = val) + ggtitle("Exposures") +
    scale_x_discrete(breaks = unique(c(seq(xv[1], max(xv) - 1, by = 10), max(xv) - 1))) +
    theme(legend.position = "bottom") + theme_bw(base_size = 20) +
    theme(legend.direction = "horizontal") +
    geom_point(data = filter(dat2, .data$row == CountrySPEC),
               aes(x = .data$Age, y = .data$Exp), show.legend = F) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 0)))

  g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}


  ggarrange(death_plot + theme_bw(base_size = 20) + theme(legend.position = "none"),
                    exposure_plot + theme_bw(base_size = 20) + theme(legend.position = "none"),
                    ncol=2, nrow=1 ,common.legend = TRUE, legend="bottom",
                                               legend.grob = g_legend(exposure_plot))

}
