#' @title Calculate the NRQs values for each run
#'
#' @description Function to automatically calculate the NRQs values for each run.
#'
#' @param run.data.df
#' data frame or tibble containing qPCR Cq values data for each gene per sample
#' per run per well. It requires a column per gene, with, as well as columns
#' to indicate the Well, Run, and Sample in which a given Cq value was measured.
#'
#' @param n_replicates Number of replicates per sample within runs
#'
#' @param reference_genes
#' Character vector with the name(s) of the reference gene(s)
#'
#' @param cals_identifier
#' Character string indicating a text pattern to identify the calibrators by
#' their names in the Sample column of the data frame, via stringr::str_detect.
#' Default is set to "Cal".
#'
#' @param amp_efficiencies
#' Object returned by the amp_efficiency function.
#'
#' @param cal.factors.list
#' Object returned by calibration_factors function.
#'
#' @param n_control
#' Number of control samples.
#' The function assumes that this does not apply to the calculations it has to
#'  make, thus its default value is NA.
#'
#' @param trace
#' logical, print additional information.
#'
#' @param geo
#' logical, to scale to your control group, the function will use the
#'  geometrical mean if TRUE or the arithmetic mean if FALSE.
#'
#' @param na.rm
#' a logical value indicating whether NA values should be stripped before the
#'  computation proceeds.
#'
#' @import dplyr
#' @import tidyr
#' @import purrr
#'
#' @export
nrqs <- function(run.data.df
                 , n_replicates
                 , reference_genes
                 , cals_identifier = "Cal"
                 , amp_efficiencies
                 , cal.factors.list
                 , n_control = NA
                 , trace = F
                 , geo = T
                 , na.rm = T) {
  # Get the label of the runs
  runs <- run.data.df$Run |>
    unique()

  # Sets the data frame structure to the one required by normalize_data
  tmp_df <- run.data.df |>
    select(-(Well:Run)) |>
    as.data.frame()

  # Defines n_control as 5, if it was set to NA in the function call (default)
  # Note: This does not affect the calibration factors calculation,
  # but the original EasyqpcR function requires it to work.
  # I did not changed this to allow using the custom function to be used
  # in the case of a control vs treatment analysis, although that was not my
  # case
  if (is.na(n_control)) {
    n_control <- 5
  }

  nrqs.list <- list()
  for (run in runs) {
    run_calibration_factor <- cal.factors.list |>
      pluck(run)

    run_cal <- run.data.df |>
      filter(Run == run
             , Sample |> stringr::str_detect(cals_identifier)) |>
      pull(Sample) |>
      unique()

    if (length(run_cal) > 1) {
      first_cal <- 1 + (run.data.df |>
                          filter(!stringr::str_detect(Sample
                                                      , cals_identifier)) |>
                          nrow()/n_replicates)
      last_cal <- run.data.df |>
        filter(!Sample %in% run_cal) |>
        nrow()/n_replicates

      run_cal <- c(first_cal:last_cal)
    }

    if (length(run_cal) == 1) {
      run_cal <- run.data.df |>
        filter(!Sample %in% run_cal) |>
        nrow()/n_replicates
    }

    nrqs.list[[run]] <- tmp_df |>
      normalize_data(
        # Set the number of replicates
        r = n_replicates
        # Set amplification efficiencies of the genes
        , E = amp_efficiencies$efficiency
        # Set the standard error of the amplification efficiencies of the genes
        , Eerror = amp_efficiencies$se.efficiency
        # Calculate the number of samples within the data frame
        , nSpl = nrow(tmp_df)/n_replicates
        # Set the number of reference genes
        , nbRef = length(reference_genes)
        # Set the column position(s) of the reference gene(s),
        # assuming they are the first genes in the data frame
        , Refposcol = 1:length(reference_genes)
        # Set the number of control samples
        , nCTL = n_control
        # Sets the calibration factor for the corresponding run
        , CF = run_calibration_factor
        # Sets the row positions of the calibrator of the corresponding run
        , CalPos = run_cal
        , trace = trace
        , geo = geo
        , na.rm = na.rm) |>
      pluck(2)
  }
  nrqs.list
}
