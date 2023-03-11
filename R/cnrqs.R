#' @title Calculate CNRQs values for the whole qPCR experiment
#' (batch / run series)
#'
#' @description Function to obtain the CNRQs values for the whole qpcr
#' experiment (batch / run series)
#'
#' @param nrqs.list
#' Object returned by nrqs function.
#'
#' @param run.data.df
#' data frame or tibble containing qPCR Cq values data for each gene per sample
#' per run per well. It requires a column per gene, with, as well as columns
#' to indicate the Well, Run, and Sample in which a given Cq value was measured.
#'
#' @import dplyr
#' @import tidyr
#' @import purrr
#'
#' @return list
#'
#' @export

cnrqs <- function(nrqs.list, run.data.df) {

  cnrqs.list <- list()
  for (run in names(nrqs.list)) {
    run_nrqs <- nrqs.list |>
      pluck(run)

    cnrqs.list[[run]] <- run_nrqs |>
      filter(rownames(run_nrqs) %in%
               (run.data.df |>
                  filter(Run == run) |>
                  pull(Sample)))
  }
  cnrqs <- cnrqs.list |>
    list_rbind()
  cnrqs
}
