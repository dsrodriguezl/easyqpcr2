nrmData  <- function(data, r, E, Eerror, nSpl, nbRef, Refposcol, nCTL, CF, 
                     CalPos, trace=FALSE, geo=FALSE, na.rm=na.rm) 
        
{
        
    n <- ncol(data)
    i <- 2:n
        
    data1 <- data[,i]
    ctmean <- aggregate(data[,i], by=list(rep(1:(nrow(data[,i])/r), each=r)), 
                        mean, na.rm=na.rm)
    rownames(ctmean) <- data[seq(1, nrow(data), by=r), 1]
    ctmean1 <- ctmean[,2:n]
        
    ctmeanall <- colMeans(ctmean[,i], na.rm=na.rm)
    ctmeanall <- as.data.frame(ctmeanall)
    ctmeanall1 <- ctmeanall[rep(1:nrow(ctmeanall), each=nrow(ctmean)),]
    E1 <- matrix(rep(E, each=nrow(ctmean)), nrow=nrow(ctmean), byrow=FALSE)
    RQ <- E1^(ctmeanall1-ctmean1)     
        
    rootgeoM <- function(data, E){data^(1/nbRef)}
        
    geoM <- mdply(RQ[,Refposcol], 'prod', MARGIN=1)
    geoM <- geoM[, ncol(geoM)]
    geoM <- as.data.frame(geoM)
    rq1 <- rootgeoM(geoM, E)
    rq1 <- as.data.frame(rq1)
    rownames(rq1) <- rownames(ctmean)
    colnames(rq1) <- "NF"
    data1 <- data[,i]
           
    ctmean <- aggregate(data[,i], by=list(rep(1:(nrow(data[,i])/r), each=r)), mean, 
                      na.rm=na.rm)
    ctmean <- ctmean[,i]
    ctmean <- as.data.frame(ctmean)
    ctmean2 <- ctmean[rep(1:nrow(ctmean), each=r),]
        
        
    u <- data1-ctmean2
    u <- as.data.frame(u)
    u1 <- u^2
    u2 <- aggregate(u1, by=list(rep(1:(nrow(u1)/r), each=r)), sum, na.rm=na.rm)
    u2 <- u2[,i]
        
    u3 <- sqrt((1/(r*(r-1)))*u2)
    u3 <- as.data.frame(u3)    
        
    u4 <- aggregate(data1, by=list(rep(1:(nrow(data1)/r), each=r)), sd, na.rm=na.rm)
    u4 <- u4[,i]
    u4 <- as.data.frame(u4)
             
    ub <- aggregate(u, by=list(rep(1:(nrow(u)/r), each=r)), mean, na.rm=na.rm)
    ub <- ub[,i]
        
    r2 <- rq1[rep(1:(nSpl), each=ncol(data1)),]
    r3 <- matrix(r2, ncol=length(i), byrow=TRUE)
    r3 <- r3[,Refposcol]
        
    Eerror1 <- matrix(rep(Eerror, each=nrow(ctmean)), nrow=nrow(ctmean), byrow=FALSE)
        
        
    w1 <- r3*sqrt( (Eerror1[,Refposcol]*sqrt(r)*((ctmeanall1-ctmean1)[,Refposcol])/E1[,Refposcol])^2+((log(E1[,Refposcol])*u4[,Refposcol])^2))

    if (nbRef == 1) {
        y <- (r3/sqrt(r)) * sqrt(rowSums(as.data.frame(w1/(nbRef * r3))^2))
        y1 <- as.data.frame(y)
    }
    else {
        y <- (r3/sqrt(r)) * sqrt(rowSums((w1/(nbRef * r3))^2))
        y1 <- as.data.frame(y[, 2])
    }

    rownames(y1) <- rownames(rq1)
    colnames(y1) <- "SE(NF)"
    yy2 <- cbind(rq1, y1)
                
    r4 <- rq1[rep(1:(nrow(rq1)), each=ncol(data1)),]
    r5 <- matrix(r4, ncol=length(i), byrow=TRUE)
        
    y2 <- sqrt(r)*(y1[rep(1:(nrow(y1)), each=ncol(data1)),])
    if (nbRef == 1) {
        y3 <- matrix(as.matrix(y2), ncol = length(i), byrow = TRUE)
    }
    else {
        y3 <- matrix(y2, ncol = length(i), byrow = TRUE)
    }
        
    w2 <- (RQ)*sqrt( (Eerror1*sqrt(r)*(ctmeanall1-ctmean1)/E1)^2+((log(E1)*u4)^2)) 
        
    CF1 <- rbind(rep(CF, each=(nrow(RQ))))
    CF1 <- as.data.frame(matrix(CF1, ncol=ncol(RQ)))
        
    o <- as.data.frame((RQ/r5)/CF1)
    rownames(o) <- rownames(rq1)
    colnames(o) <- colnames(ctmean)
              
    if (geo) { 
            
        if (nCTL >= 2) {

            o1 <- (colProds(o[1:nCTL,], na.rm=na.rm))^(1/length(1:nCTL))
            o1 <- as.data.frame(o1) 

        } else {
                            
    		o1 <- (o[1:nCTL,])
    		o1 <- as.data.frame(o1)
                
        }

    } else {
            
        if (nCTL >= 2){
                
            o1 <- colMeans(o[1:nCTL,], na.rm=na.rm)
            o1 <- as.data.frame(o1)
                
        } else {
                            
            o1 <- (o[1:nCTL,])
            o1 <- as.data.frame(o1)
                
        }
    }
        
    if (nCTL >=2) {

        o1 <- o1

    } else {

        o1 <- as.numeric(o1)

    }

    o1 <- as.data.frame(o1)
    o5 <- (o1[rep(1:(ncol(data1)), each=nrow(ctmean1)),])
    o6 <- as.data.frame(matrix(o5, ncol=length(i), byrow=FALSE))
    rownames(o6) <- rownames(o)
    colnames(o6) <- colnames(o)
    o4 <- o/o6
        
    t <- (o4/sqrt(r))*sqrt(((y3/r5)^2)+((w2/RQ)^2))
        
    o2906 <- o[CalPos,]
        
    if (trace) { 
          
        message(paste("The Standard Error of Normalized Relative Quantity for", 
            "each gene and each sample scaled to control is:",  "",  sep="\n"))
        print(t)
            
        message(paste("The Standard Error of Relative Quantity for",  
            "each gene and each sample is:",  "", sep="\n"))
        print(w2)
            
        message(paste("The Standard Error of Normalization Factor for each", 
            "sample is:",  "", sep="\n"))
    	print(y1)
            
    }
        
    return(list('NRQs normalized to control'=o4, 'NRQs'=o, 
    'NRQs of your calibrator for this run'=o2906))
}
