#' @title Function to calculate the amplification efficiency
#'
#' @description This is a custom version of the "EasyqpcR::slope function.
#'  It has been adapted to run on recent versions of R, and to
#'  calculate the standard error of the efficiency (not only the efficiency).
#'
#' @param data data.frame containing row datas (genes in columns, samples in
#' rows, Cq values).
#'
#' @param q numeric, cDNA dilution values.
#'
#' @param r numeric, number of qPCR replicates.
#'
#' @param na.rm logical, indicating whether NA values should be stripped before
#' the computation proceeds.
#'
#' @import dplyr
#' @import tidyr
#' @import broom
#'
#' @author
#' Sylvain Le pape <sylvain.le.pape@univ-poitiers.fr> (original function author)
#'
#' Daniel S. Rodríguez-Leon (author of easyqpcr2)
#'
#' @usage
#' amp_efficiency(data, q, r, na.rm = FALSE)
#'
#' @returns
#'  A tibble data frame containing the amplification efficiency calculation
#'  results for each gene.
#'
#'  intercept Intercept of the dilution curve
#'
#'  slope Slope of the dilution curve.
#'
#'  se.slope Standard error of the slope of the dilution curve.
#'
#'  efficiency Primer amplification efficiency.
#'
#'  se.efficiency Standard error of the primer amplification efficiency.
#'
#' @references Jan Hellemans, Geert Mortier, Anne De Paepe, Frank Speleman and
#' Jo Vandesompele. qBase relative quantification framework and software for
#' management and automated analysis of real-time quantitative PCR data. Genome
#' Biology 2007, 8:R19 (doi:10.1186/gb-2007-8-2-r19).
#'
#'
#' @export
amp_efficiency  <-  function(data, q, r, na.rm = FALSE) {

  # Get the number of columns
  n <- ncol(data)
  # Get column index for the genes
  i <- 2:n

  if (length(i) < 2) {
    i <- 2
  }
  if (length(i) >= 2) {
    data[, i] <- as.data.frame(data[, i])
    # Calculate the mean Cq from the replicates
    am <- (aggregate(data[, i]
                     , by = list(rep(1:(nrow(data[, i]) / r)
                                     , each = r))
                     , mean
                     , na.rm = TRUE))
    # Remove Group column
    am <- data.frame(am[, 2:ncol(am)])
    am1 <- as.matrix(am)
  } else {
    data <- as.data.frame(data[,2])
    # Calculate the mean Cq from the replicates
    am <- (aggregate(data
                     , by = list(rep(1:(nrow(data) / r)
                                     , each = r))
                     , mean
                     , na.rm = TRUE))
    # Remove Group column
    am <- data.frame(am[, 2:ncol(am)])
    am1 <- as.matrix(am)

  }

  # Fit a linear model predicting the mean Cq values as a function of the
  # know amount of cDNA
  fit <- lm(am1 ~ log10(q))

  # Extract the intercept, slope, and std. error of the slope
  intercept <- fit |>
    tidy() |>
    filter(term == "(Intercept)") |>
    pull(estimate)

  names(intercept) <- fit |>
    tidy() |>
    filter(term == "(Intercept)") |>
    pull(response)

  E <- bind_cols(gene = fit |>
                   tidy() |>
                   filter(term == "(Intercept)") |>
                   pull(response)
                 , intercept = fit |>
                   tidy() |>
                   filter(term == "(Intercept)") |>
                   pull(estimate)
                 , fit |>
                   tidy() |>
                   filter(term == "log10(q)") |>
                   select(estimate, std.error)) |>
    mutate(slope = estimate
           , se.slope = std.error
           , .keep = "unused") |>
    mutate(efficiency = 10^(-1 / slope)
           , se.efficiency = ((efficiency) * log(10) * se.slope / slope^2))

  return(E)
}

