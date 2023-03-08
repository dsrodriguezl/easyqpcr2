#' @title amp_efficiency
#'
#' @description This is a custom version of the "slope" function of the EasyqpcR package.
#' It has been adapted to run on recent versions of R, and to calculate the
#' standard error of the efficiency (not only the efficiency)
#'
#' @param data
#' @param q
#' @param r
#' @param na.rm
#'
#' @import dplyr
#' @import tidyr
#' @import broom
#'
#' @return list
#'
#' @export

amp_efficiency  <-  function(data, q, r, na.rm = FALSE) {

  # Get the number of columns
  n <- ncol(data)
  # Get column index for the genes
  i <- 2:n

  # if (length(i) >= 2) {
    i <- 2:n
  # } else {i <- 2}

  # if (length(i) >= 2) {
    data[,i] <- as.data.frame(data[,i])
    # Calculate the mean Cq from the replicates
    am <- (aggregate(data[,i]
                     , by = list(rep(1:(nrow(data[,i])/r) , each=r))
                     , mean
                     , na.rm=TRUE))
    # Remove Group column
    am <- data.frame(am[,2:ncol(am)])
    am1 <- as.matrix(am)
  # } else {
  #   data <- as.data.frame(data[,2])
  #   # Calculate the mean Cq from the replicates
  #   am <- (aggregate(data
  #                    , by=list(rep(1:(nrow(data)/r), each=r))
  #                    , mean
  #                    , na.rm=TRUE))
  #   # Remove Group column
  #   am <- data.frame(am[, 2:ncol(am)])
  #   am1 <- as.matrix(am)
  #
  # }

  # if (length(i) >= 2) {
    # Fit a linear model predicting the mean Cq values as a function of the
    # know amount of cDNA
    fit <- lm(am1 ~ log10(q))

    # Extract the intercept, slope, and std. error of the slope
    # intercept <- coef(fit)[1,]
    intercept <- fit |>
      tidy() |>
      filter(term == "(Intercept)") |>
      pull(estimate)

    names(intercept) <- fit |>
      tidy() |>
      filter(term == "(Intercept)") |>
      pull(response)

    # slope <- coef(fit)[2,]
    # slope <- fit |>
    #   tidy() |>
    #   filter(term == "log10(q)") |>
    #   select(estimate, std.error) |>
    #   mutate() |>
    #   as.data.frame()

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

