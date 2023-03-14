#' @title Calculate the calibration factors of each run.
#'
#' @description  Function to automatically get the calibration
#' factors of each run.
#'
#' It relies on the original EasyqpcR::calData function.
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
#' Object returned by the amp_efficiency function
#'
#' @param CF
#' Calibration factors. The function assumes that this variable does not apply,
#' as the function purpose is to calculate them.
#'
#' This is a parameter of EasyqpcR::nrmData, where it would be numeric (or
#' object if you have used the EasyqpcR::calData function), values of
#' the calibration factors for each gene (follow the same order of the genes).
#'
#' calibration_factors can take defined values for CF, in case you have
#' calculated them externally (e.g. using the original EasyqpcR::calData).
#' In that case it will just pass the specified CF to normalize_data.
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
#' @return list
#'
#' @export

calibration_factors <- function(run.data.df
                                , n_replicates
                                , reference_genes

                                , cals_identifier = "Cal"
                                , amp_efficiencies
                                , CF = NA
                                , n_control = NA
                                , trace = F
                                , geo = T
                                , na.rm = T) {
  # Get the label of the runs
  runs <- run.data.df$Run |>
    unique()

  # Sets the data frame structure to the one required by normalize_data
  tmp_df <- run.data.df |>
    select(-(contains("Well"):contains("Run"))) |>
    as.data.frame()

  # If the CF is set to NA in the function call (default)
  # CF is then re-defined as a numeric vector full of 1s, with a length equal
  # to the number of genes
  if (is.na(CF)) {
    CF <- rep(1
              , ncol(tmp_df |>
                       select(-contains("Sample"))))
  }

  # Defines n_control as 5, if it was set to NA in the function call (default)
  # Note: This does not affect the calibration factors calculation,
  # but the original EasyqpcR::nrmData function requires it to work.
  # I did not changed this to allow using the custom function to be used
  # in the case of a control vs treatment analysis, although that was not the
  # analysis that drived my interest on developing the package.
  if (is.na(n_control)) {
    n_control <- 5
  }

  # Get the position of the first calibrator
  # Note: Assumes the calibrators are all at the end of the data frame
  first_cal <- 1 + (run.data.df |>
                      filter(!stringr::str_detect(get("Sample")
                                                  , cals_identifier)) |>
                      nrow() / n_replicates)

  # Get the position of the last calibrator
  # Note: Assumes the calibrators are all at the end of the data frame
  last_cal <- run.data.df |>
    nrow() / n_replicates

  cal.factors.list <- list()
  for (run in runs) {
    cal.factors.list[[run]] <- tmp_df |>
      normalize_data(
        # Set the number of replicates
        r = n_replicates
        # Set amplification efficiencies of the genes
        , E = amp_efficiencies$efficiency
        # Set the standard error of the amplification efficiencies of the genes
        , Eerror = amp_efficiencies$se.efficiency
        # Calculate the number of samples within the data frame
        , nSpl = nrow(tmp_df) / n_replicates
        # Set the number of reference genes
        , nbRef = length(reference_genes)
        # Set the column position(s) of the reference gene(s),
        # assuming they are the first genes in the data frame
        , Refposcol = 1:length(reference_genes)
        # Set the number of control samples
        , nCTL = n_control
        # Sets the calibration factors
        , CF = CF
        # Sets the row positions of the calibrators
        , CalPos = c(first_cal:last_cal)
        , trace = trace
        , geo = geo
        , na.rm = na.rm) |>
      pluck(3) |>
      slice(as.numeric(run)) |>
      calData()
  }
  cal.factors.list
}
