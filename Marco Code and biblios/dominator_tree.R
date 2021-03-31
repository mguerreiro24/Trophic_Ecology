library("igraph")

##tets<-rbind(c(0, 0, 1, 0),c(1, 0, 0, 1),c(1, 0, 0, 0), c(1, 0, 0, 0))
##
##
##

RootFW <- function(m){
		n<- nrow(m)
		n <- n + 1
		##
		b.2 <- m
		b.1 <- rbind(rep(0,ncol(m)),b.2)
		b <- cbind(rep(0,nrow(b.1)),b.1)
		##
		a <- b
		a[1,] <- apply(a,2,sum)
		a[1,1]<-0
		##
		for (i in 2:n){
			if (a[1,i]>0){
				a[1,i]<-0
				}
			##
			else{
				a[1,i]<-1
				}
			}
		##
		colnames(a) <- as.character(0:(n-1))
		rownames(a) <- as.character(0:(n-1))
		##
		return(a)
}

##
##
##
##

Dom <- function(MAT){
	N <- nrow(MAT)
	DM <- matrix(rep(1,N^2), nrow = N)
	DM[1,2:N] <- 0
	##
	Diff <- 1 
	while(Diff > 0){
		Diff <- 0
		for (i in 2:N){
			tmp <- rep(1,N)	
			for (j in 1:N){
				if(MAT[j,i]==1){ 
					tmp <- tmp * DM[j,]
				}
			}
			tmp[i] <- 1
			if (sum(abs(tmp - DM[i,]))!= 0) { 
				Diff <- 1
			} 
			DM[i,] <- tmp
		}
	}
	##
	colnames(DM) <- c(0:(N-1))
	rownames(DM) <- c(0:(N-1))
	##
	return(DM) 
}

##
##
##
##

DomTree <- function(DOM){
	N <- nrow(DOM)
	DT <- matrix(rep(0,N^2), nrow = N)
	##
	for (j in 1:N){ 
		tmp <- rep(0,N)
		##
		for (k in 1:N){
			if (j!=k){
				tmp <- DOM[j,]
				tmp[k] <- 1
				if (sum(abs(tmp-DOM[k,]))== 0){ 
					DT[j,k] <- 1
				}
			}
		}
	}
	##
	colnames(DT) <- c(0:(N-1))
	rownames(DT) <- c(0:(N-1))
	##
	return(DT)
}
