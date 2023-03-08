#' @title nrqs
#'
#' @description Function to automatically calculate the NRQs values for each run
#'
#' @param run.data.df
#' @param n_replicates
#' @param reference_genes
#' @param cals_identifier
#' @param amp_efficiencies
#' @param cal.factors.list
#' @param n_control
#' @param trace
#' @param geo
#' @param na.rm
#'
#' @import dplyr
#' @import tidyr
#' @import purrr
#'
#' @export

nrqs <- function(run.data.df
                 # Number of replicates per sample within runs
                 , n_replicates
                 # Vector with the name(s) of the reference gene(s)
                 , reference_genes
                 # String with a pattern to identify
                 # the calibrators by their names in the Sample
                 # column of the data frame
                 , cals_identifier = "Cal"
                 # Object returned by the amp_efficiency function
                 , amp_efficiencies
                 # Calibration factors
                 # It assumes this variable does not apply,
                 # as the function purpose is to calculate them
                 , cal.factors.list
                 # Number of control samples
                 # As my data does not have control/treatment
                 # groups, I made the function assuming that
                 # this does not apply to the calculations it
                 # has to make
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
