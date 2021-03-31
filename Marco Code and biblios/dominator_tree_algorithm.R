setwd("D:/ECOLOGY/2020-06-08---dominator_tree/dominator_tree_algorithm")

##
##

## (1) data preparation -- Prince William Sound food web

PWS <- read.table(file = "PWS.txt", header = FALSE, sep = "\t", colClasses = rep("numeric",3))
n_max <- max(c(PWS[,1], PWS[,2]))
tt <- nrow(PWS)

PWS_AM <- matrix(rep(0,n_max^2), nrow = n_max)
PWS_B <- matrix(rep(0,n_max^2), nrow = n_max)

for(i in 1:tt){
	prey <- PWS[i,1]
	pred <- PWS[i,2]
	weight <- PWS[i,3]
	PWS_B[prey,pred] <- 1
	PWS_AM[prey,pred] <- weight
}

PWS_IN <- apply(PWS_AM,2,sum)
PWS_G <- matrix(rep(0,n_max^2), nrow = n_max)

for(i in 1:n_max)
	for(j in 1:n_max){
		if(PWS_AM[i,j]!=0)PWS_G[i,j] <- PWS_AM[i,j]/PWS_IN[j]
}

PWS_B_5 <- PWS_B
PWS_B_10 <- PWS_B
PWS_B_20 <- PWS_B

for(i in 1:n_max)
	for(j in 1:n_max){
		##
		## removing links with relative importance <= 0.05
		if(PWS_G[i,j] > 0 & PWS_G[i,j] <= 0.05)PWS_B_5[i,j] <- 0
		##
		## removing links with relative importance <= 0.10
		if(PWS_G[i,j] > 0 & PWS_G[i,j] <= 0.10)PWS_B_10[i,j] <- 0
		##
		## removing links with relative importance <= 0.20
		if(PWS_G[i,j] > 0 & PWS_G[i,j] <= 0.20)PWS_B_20[i,j] <- 0
}

write.table(PWS_B, file = "PWS_B.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(PWS_B_5, file = "PWS_B_5.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(PWS_B_10, file = "PWS_B_10.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(PWS_B_20, file = "PWS_B_20.txt", sep = " ", col.names = FALSE, row.names = FALSE, quote = FALSE)

##
##

## (2) the rooted food web (RFW)

## the input for the function rooted_FW is a binary matrix (i.e. adjacency matrix)
## such matrix depicts predator-prey interactions (i.e. who eats whom, with prey in
## the rows and predators in the columns)
##
## the output of this function is a binary matrix equivalent to the input matrix,
## but with one row and one column more -- the main objective is adding the root node
## to the adjaceny matrix used as an input; the root node is connected to the
## primary producers and has no energy/matter flows pointing to it (i.e. it represents
## the external source of energy on which the food web depends)

rooted_FW <- function(M){
		n <- nrow(M)
		n_r <- n + 1
		##
		M_s1 <- rbind(rep(0,ncol(M)),M)
		M_r <- cbind(rep(0,n_r),M_s1)
		##
		IN <- apply(M,2,sum)
		##
		M_r[1,which(IN==0)+1] <- 1
		##
		colnames(M_r) <- c(0:n)
		rownames(M_r) <- c(0:n)
		##
		return(M_r)
}

RFW_PWS_B <- rooted_FW(PWS_B)
RFW_PWS_B_5 <- rooted_FW(PWS_B_5)
RFW_PWS_B_10 <- rooted_FW(PWS_B_10)
RFW_PWS_B_20 <- rooted_FW(PWS_B_20)


## (3) the dominator matrix (DM)

## this function returns a dominator matrix in which all of the dominators of each
## node are summarized (i.e. both immediate and non immediate dominators)

RFW_Findley <- matrix(c(0,1,1,0,0,
0,0,0,1,0,
0,1,0,0,1,
0,1,0,0,0,
0,1,0,0,0), nrow = 5, byrow = TRUE)

dominator_matrix <- function(M){
		##
		n <- nrow(M)
		DM <- matrix(rep(1,n^2), nrow = n)
		DM[1,2:n] <- 0
		##
		diff <- 1
		while(diff > 0){
			diff <- 0
			for(i in 2:n){
				vett <- rep(1,n)
				for(j in 1:n){
					if(M[j,i]==1){
						vett <- vett * DM[j,]
					}
				}
				##
				vett[i] <- 1
				##
				if(sum(abs(vett - DM[i,])) != 0){
					diff <- 1
				}
				##
				DM[i,] <- vett
			}
		}
		##
		colnames(DM) <- c(0:(n-1))
		rownames(DM) <- c(0:(n-1))
		##
		return(DM)
}

DM_Findley <- dominator_matrix(RFW_Findley)


## (4) dominator tree (DT)

## this function returns a binary matrix with immediate dominators only
## the input is the dominator matrix computed with the previous function
## the output is a binary matrix that corresponds to the dominator tree
## the dominator tree is a graph with N nodes (N-1 species or trophospecies plus
## the root node) and with N-1 links

dominator_tree <- function(M){
		##
		n <- nrow(M)
		DT <- matrix(rep(0,n^2), nrow = n)
		##
		for(j in 1:n){
			vett <- rep(0,n)
			for(k in 1:n){
				if(j != k){
					vett <- M[j,]
					vett[k] <- 1
					if(sum(abs(vett - M[k,])) == 0){
						DT[j,k] <- 1
					}
				}
			}
		}
		##
		colnames(DT) <- c(0:(n-1))
		rownames(DT) <- c(0:(n-1))
		##
		return(DT)
} 

DT_Findley <- dominator_tree(DM_Findley)

##
##

####################################
##								  ##
##   error & attack sensitivity   ##
##								  ##
####################################

## nodes in the rooted food web (example for Lake Findley) ---> DM_Findley
N <- nrow(DM_Findley)

## [-1] indicates that all nodes are considered except the root
## DO_n returns the fraction of nodes dominated by each node
DO_n <- (apply(DM_Findley,2,sum)[-1]-1)/(N-1)

## the average of DO_n corresponds to the error sensitivity
## below I reported three equivalent formulations
mean(DO_n)
sum(DO_n)/(N-1)
sum(apply(DM_Findley,2,sum)[-1]-1)/(N-1)^2
##
ES <- sum(DO_n)/(N-1)
ES

## ---> standard deviation is different from zero when computed for the number of secondary
## extinctions triggered by all possible primary extinctions (i.e. the only case of sd = 0
## is when primary extinctions of all nodes trigger the same number of secondary extinctions)
## ---> the standard deviation cannot be calculated for error or attack sensitivity as
## these are single numbers and will always have sd = 0
sd((apply(DM_Findley,2,sum)[-1]-1)/(N-1))
sd(DO_n)

## attack sensitivity is the largest fraction of secondary extinctions triggered by a primary extinction
AS <- max(DO_n)
AS
