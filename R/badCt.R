#' @title Evaluation of the qPCR technical replicates
#'
#' @description Original EasyqpcR function.
#' This function allows you to evaluate your qPCR technical
#' replicates, you only need to define the threshold (according to Hellemans et
#'  al. (2007), 0.5 is a good threshold value), the dataset and the number of
#'  technical replicates you have done. I recommand you to use the gWidgets
#'  package to easily exclude the failed replicates.
#'
#' @param data data.frame containing row datas (genes in columns, samples in
#' rows, Cq values).
#' @param r numeric, number of qPCR replicates.
#' @param threshold numeric, the maximal variation between your qPCR replicates.
#' @param na.rm a logical value indicating whether NA values should be stripped
#' before the computation proceeds.
#'
#' @details To facilitate the use of the function, I suggest you to use the
#' gWidgets package as described in the vignette.
#'
#' @usage badCt(data, r, threshold, na.rm = FALSE)
#'
#' @returns This function returns the position (sample position and column
#'  position) where the variation between qPCR replicates is superior to the
#'   threshold value.
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
#' @examples
#'
#' data(qPCR_run1)
#'
#' badCt(data=qPCR_run1, r=3, threshold=0.3, na.rm=TRUE)
#'
#' @export
#'
badCt <- function(data, r, threshold, na.rm=FALSE)
{

    if (!is.data.frame(data) & !is.matrix(data))
    stop("'relData' has to of class matrix or data.frame")

    n <- ncol(data)

    i <- 2:n

    if (is.numeric(i))

        ctmean<-aggregate(data[,i], by=list(rep(1:(nrow(data[,i])/r), each=r)),
        mean, na.rm=na.rm)
        rownames(ctmean) <- data[seq(1, nrow(data), by=r), 1]
        ctmean1 <- ctmean[,2:n]

        SD <- aggregate(data[,i], by=list(rep(1:(nrow(data[,i])/r), each=r)),
        std.error, na.rm=na.rm)
        rownames(SD) <- data[seq(1, nrow(data), by=r), 1]
        SD1 <- SD[,2:n]

        ctmax <- aggregate(data[,i],
        by=list(rows=rep(1:(nrow(data[,i])/r), each=r)), max, na.rm=na.rm)
        rownames(ctmax) <- data[seq(1, nrow(data), by=r), 1]
        ctmax1 <- ctmax[,2:n]

        ctmin <- aggregate(data[,i],
        by=list(rows=rep(1:(nrow(data[,i])/r), each=r)), min, na.rm=na.rm)
        rownames(ctmin) <- data[seq(1, nrow(data), by=r), 1]
        ctmin1 <- ctmin[,2:n]

        ctdrep <- ctmax-ctmin
        rownames(ctdrep) <- data[seq(1, nrow(data), by=r), 1]
        ctdrep1 <- ctmax1-ctmin1

        ctbadrep <- which(ctdrep>threshold, arr.ind=TRUE)


    return(list('Bad replicates localization'=ctbadrep,
    'Mean of the Cq'=ctmean1, 'Standard error of the Cq'=SD1))


}
