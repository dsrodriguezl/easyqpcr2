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
