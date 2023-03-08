#' @title cnrqs
#'
#' @description Function to obtain the CNRQs values for the whole qpcr
#' experiment (batch / run series)
#'
#' @param nrqs.list
#' @param run.data.df
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
}
