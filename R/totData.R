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
