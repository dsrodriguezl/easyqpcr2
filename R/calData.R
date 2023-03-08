calData <-function(data) 
{
        
    n<-ncol(data)
    i<-2:n
    datas<-data[,i]
    a<-sapply(data,prod,margin=1)
    b<-a^(1/nrow(data))
        
    return(b)
        
}
