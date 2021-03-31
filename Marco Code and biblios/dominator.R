## (1) the rooted food web (RFW)

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


## (2) the dominator matrix (DM)

## this function returns a dominator matrix in which all of the dominators of each
## node are summarized (i.e. both immediate and non immediate dominators)

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


## (3) dominator tree (DT)

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
