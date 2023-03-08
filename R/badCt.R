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
