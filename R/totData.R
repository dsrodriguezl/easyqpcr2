#' @title Aggregation of qPCR biological replicates and data transformation
#'
#' @description Original EasyqpcR function.
#' This function aggregates qPCR biological replicates and calculates the main
#' parameters : mean (arithmetic or geometric), the standard deviation and the
#' standard error from your biological replicates of your experience.
#' This function has an algorithm published by Willems et al. (2008) which
#' performs a standardization procedure that can be applied to data sets that
#' display high variation between biological replicates. This enables proper
#' statistical analysis and drawing relevant conclusions. The procedure is not
#' new, it has been used in microarray data analysis and is based on a series
#' of sequential corrections, including log transformation, mean centering,
#' and autoscaling.
#'
#' @param data data.frame containing row datas (genes in columns, samples
#' in rows, Cq values).
#'
#' @param r numeric, number of qPCR replicates.
#'
#' @param geo logical, the function will use the geometrical mean of your
#' biological replicates if TRUE or the arithmetic mean if FALSE.
#'
#' @param logarithm logical, the NRQs will be log-transformed.
#'
#' @param base numeric, the logarithmic base (2 or 10).
#'
#' @param transformation logical, if TRUE, the transformation procedure for
#' highly variable biological replicates (but with the same tendency) will be
#' done.
#'
#' @param nSpl numeric, the number of samples.
#'
#' @param linear logical, after the transformation procedure done, your raw
#' data will be normalized (anti-log-transformed).
#'
#' @param na.rm logical, indicating whether NA values should be stripped before
#' the computation proceeds.
#'
#' @details The standardization procedure used in this function (if TRUE for
#' the transformation argument) is based on the article of Willems et al.
#' (2008). This function perform successively thEerror operations :
#' log-transformation of your raw data, mean of your log-transformed data for
#' each biological replicate, mean centering for each biological replicate,
#' standard deviation of each mean-centered biological replicate, autoscaling
#' of your data, i.e., your mean-centered data for each biological replicate
#' will be divided by the standard deviation of the mean-centered biological
#' replicate and then multiplicated by the mean of the standard deviation of
#' all the biological replicates.
#' For more information for the way to use this
#' function, please see the vignette
#'
#' @usage
#' totData(data, r, geo = TRUE, logarithm = TRUE, base, transformation = TRUE,
#'  nSpl, linear = TRUE, na.rm = na.rm)
#'
#' @returns
#' Mean of your qPCR runs.
#' The geometric (if TRUE for geo) or arithmetic mean of your biological
#' replicates.
#'
#' Standard deviations of your qPCR runs.
#' The standard deviation of your biological replicates.
#'
#' Standard errors of your qPCR runs.
#' The standard error of your biological replicates.
#'
#' Transformed data.
#' If TRUE for transformation, your raw data will be transformed by the algorithm
#' of Willems et al. (2008).
#'
#' Reordered transformed data.
#' The transformed data reordered by rowname.
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
#'
totData  <- function(data, r, geo=TRUE, logarithm=TRUE, base,
                     transformation=TRUE, nSpl, linear=TRUE, na.rm=na.rm)
{

    colnames_removing_prefix <- function(df, prefix)
    {
        names <- colnames(df)
        indices <- (substr(names, 1, nchar(prefix))==prefix)
        names[indices] <- substr(names[indices], nchar(prefix)+1,
                                 nchar(names[indices]))
    }

    if (nchar(colnames(data[1]))<=28) {

        colnames(data) <- colnames(data)

    } else {

        colnames(data) <- colnames_removing_prefix(data,
        "NRQs.normalized.to.control.")
    }

    x <- data
    x2 <- x[order(rownames(x)),]

    if (transformation) {

        if(base==2) {x <- log2(x)}
	    else  {x <- log10(x)}

        meancentered <- colMeans(x, na.rm=na.rm)

        mc <- rbind(rep(meancentered, each=(nrow(x))))
        mc1 <- as.data.frame(matrix(mc, ncol=ncol(x)))
        mc2 <- x-mc1

        expsd <- aggregate(mc2, by=list(rep(1:(nrow(x2)/nSpl), each=nSpl)), sd,
      	              na.rm=na.rm)
        expsd <- expsd[, 2:ncol(expsd)]
        expsd1 <- expsd[rep(1:nrow(expsd), each=nrow(expsd)),]

        expsd2 <- colMeans(expsd1)
        expsd3 <- rbind(rep(expsd2, each=(nrow(x))))
        expsd4 <- as.data.frame(matrix(expsd3, ncol=ncol(x)))

        autoscaling <- mc2/expsd4
        autoscalingstd <- autoscaling*expsd4
        colnames_removing_prefix(autoscalingstd, "NRQs.normalized.to.control.")
        autoscalingstd2 <- autoscalingstd[order(rownames(autoscalingstd)),]

        totsd <- aggregate(autoscalingstd2,
    	                by=list(rep(1:(nrow(autoscalingstd2)/r), each=r)), sd,
		                na.rm=na.rm)
        totsd <- totsd[, 2:ncol(totsd)]
        rownames(totsd) <- rownames(autoscalingstd2[seq(1, nrow(autoscalingstd2),
		                            by=r),])

        totse <- aggregate(x2, by=list(rep(1:(nrow(autoscalingstd2)/r), each=r)),
                            sd, na.rm=na.rm)
        totse <- totse[, 2:ncol(totse)]
        totse <- totse/sqrt(r)
        rownames(totse) <- rownames(autoscalingstd2[seq(1, nrow(autoscalingstd2),
	                                by=r),])

        totmean <- aggregate(autoscalingstd2,
                        by=list(rep(1:(nrow(autoscalingstd2)/r), each=r)), mean,
                        na.rm=na.rm)
    totmean <- totmean[, 2:ncol(totmean)]
    rownames(totmean) <- rownames(autoscalingstd2[seq(1, nrow(autoscalingstd2),
	                                by=r),])

    if (linear)

    {
        if (base==2) {

            autoscalingstd2 <- 2^autoscalingstd2
            autoscalingstd <- 2^autoscalingstd

        } else {autoscalingstd2 <- 10^autoscalingstd2
                autoscalingstd <- 10^autoscalingstd}

    } else {autoscalingstd2 <- autoscalingstd2
            autoscalingstd <- autoscalingstd}

    }

    if (geo)
    {
        totmean <- aggregate(x2, by=list(rep(1:(nrow(x2)/r), each=r)), prod,
                           na.rm=na.rm)
        totmean <- totmean[, 2:ncol(totmean)]
        rownames(totmean) <- rownames(x2[seq(1, nrow(x2), by=r),])

        colnames_removing_prefix(totmean, "NRQs.normalized.to.control.")
        totmean <- totmean^(1/r)

    } else {

        totmean <- aggregate(x2, by=list(rep(1:(nrow(x2)/r), each=r)), mean,
                           na.rm=na.rm)
        totmean <- totmean[, 2:ncol(totmean)]
        rownames(totmean) <- rownames(x2[seq(1, nrow(x2), by=r),])}

    if (logarithm)

    {
        if (base==2) {totmean <- log2(totmean)}

        else {totmean <- log10(totmean)}

    }

    else {x2 <- x2}

    totsd <- aggregate(x2, by=list(rep(1:(nrow(x2)/r), each=r)), sd,
                        na.rm=na.rm)
    totsd <- totsd[, 2:ncol(totsd)]
    rownames(totsd) <- rownames(x2[seq(1, nrow(x2), by=r),])

    totse <- aggregate(x2, by=list(rep(1:(nrow(x2)/r), each=r)), sd,
                        na.rm=na.rm)
    totse <- totse[, 2:ncol(totse)]
    totse <- totse/sqrt(r)
    rownames(totse) <- rownames(x2[seq(1, nrow(x2), by=r),])

    if (transformation) {

        return(list('Mean of your qPCR runs'=totmean, 'Standard deviations of
        your qPCR runs'=totsd, 'Standard errors of your qPCR runs'=totse,
        'Transformed data'=autoscalingstd,
	    'Reordered transformed data'=autoscalingstd2))

    } else {

        return(list('Mean of your qPCR runs'=totmean, 'Standard deviations of
        your qPCR runs'=totsd, 'Standard errors of your qPCR runs'=totse))

    }
}
