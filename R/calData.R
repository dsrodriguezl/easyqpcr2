#' @title Calculation of calibration factors
#'
#' @description Original EasyqpcR function.
#' This function determines the calibration factors (CF) using the method
#' described in Hellemans et al. (2007).
#'
#' @param data data.frame containing the NRQs of your calibrator(s) for each
#' gene obtained by the nrmData function of this package.
#'
#' @details This function is necessary for comparing different quantitative
#' real-time PCR runs to reduce the inter-run variability (Hellemans et al.
#' (2007)). Then, the results obtained have to be included in an R object and
#' then be inputed in the nrmData function (see the vignette for more
#' informations).
#'
#' @usage calData(data)
#'
#' @returns This function returns the calibration factor associated to each gene
#'  for the whole runs.
#'
#' @author
#' Sylvain Le pape <sylvain.le.pape@univ-poitiers.fr> (original function author)
#'
#' Daniel S. Rodríguez-Leon (author of easyqpcr2)
#'
#' @references Jan Hellemans, Geert Mortier, Anne De Paepe, Frank Speleman and
#' Jo Vandesompele. qBase relative quantification framework and software for
#' management and automated analysis of real-time quantitative PCR data. Genome
#' Biology 2007, 8:R19 (doi:10.1186/gb-2007-8-2-r19).
#'
#' @export

calData <-function(data)
{

    n<-ncol(data)
    i<-2:n
    datas<-data[,i]
    a<-sapply(data,prod,margin=1)
    b<-a^(1/nrow(data))

    return(b)

}
