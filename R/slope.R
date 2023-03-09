#' @title Function to calculate the amplification efficiency
#'
#' @description Original EasyqpcR function.
#' This function calculates the amplification efficiency from
#' classical qPCR dilution experiment using the Cq values. The Cq values are
#' plotted against the logarithmized concentration (or dilution) values, a
#' linear regression line is fit and the efficiency calculated by
#' E = 10^-1/slope.
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
#' @seealso You can also see the qpcR package with the calib and calib2
#' functions.
#'
#' @usage
#'
#' slope(data, q, r, na.rm = FALSE)
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
#' data(Efficiency_calculation)
#'
#' slope(data = Efficiency_calculation, q=c(1000, 100 ,10, 1, 0.1), r=3,
#'  na.rm=TRUE)
#'
#' @export

slope  <-  function(data, q, r, na.rm=FALSE)
{
    n <- ncol(data)
    i <- 2:n

    if (length(i) >= 2) {

        i <- 2:n

    } else {i <- 2}

    if (length(i) >= 2) {

        data[,i] <- as.data.frame(data[,i])
        am <- (aggregate(data[,i], by=list(rep(1:(nrow(data[,i])/r), each=r)),
                        mean, na.rm=TRUE))
        am <- data.frame(am[,2:ncol(am)])
        am1 <- as.matrix(am)

    } else {

        data <- as.data.frame(data[,2])
        am <- (aggregate(data, by=list(rep(1:(nrow(data)/r), each=r)), mean,
                na.rm=TRUE))
        am <- data.frame(am[, 2:ncol(am)])
        am1 <- as.matrix(am)

    }

    if (length(i) >= 2) {

        intercept <- coef(lm(am1~log10(q)))[1,]
        slope <- coef(lm(am1~log10(q)))[2,]

    } else {

        intercept <- coef(lm(am1~log10(q)))[1]
        slope <- coef(lm(am1~log10(q)))[2]

    }

    q1 <- as.matrix(q)

    q2 <- matrix(c(rep(q1, length(i))), byrow=TRUE)
    q3 <- matrix(log10(q2), ncol=length(i))


    slope1 <- as.data.frame(slope, byrow=TRUE)
    colnames(slope1) <- "Slope"


    if (length(i) >= 2) {

        z <- rbind(t(slope1), q3)
        colnames(z) <- colnames(am)
        combss  <-  combn(seq_len(nrow(z)),  2)
        matt <- data.matrix(z)
        z1 <- (matt[combss[1,],] * matt[combss[2,],])
        z1 <- z1[1:nrow(am),]

    } else {

        z <- rbind(t(slope1), q3)
        colnames(z) <- colnames(am)
        combss  <-  combn(seq_len(nrow(z)),  2)
        matt <- data.matrix(z)
        z1 <- (matt[combss[1,],] * matt[combss[2,],])
        z1 <- as.data.frame(z1)
        z1 <- z1[1:nrow(am),]

    }

    intercepta <- as.matrix(intercept)
    intercept1 <- intercepta[rep(1:(nrow(intercepta)), each=length(q))]
    intercept2 <- matrix(intercept1, byrow=TRUE)
    intercept3 <- matrix(intercept2, ncol=length(i))

    z2 <- z1+intercept3

    E <- as.data.frame(10^(-1/slope))
    colnames(E) <- "E"
    rownames(E) <- c(colnames(am1))

    return(list('Efficiency'=E))

}
