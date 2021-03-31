rm(list=ls(all=TRUE))

setwd("D:/People/Miguel_Fernandes_Guerreiro/R_Allesina")

##
##

## load the three functions for constructing the dominator tree
source("dominator.R")

## load the "igraph" library
## library("igraph")

## example with a hypothetical food web composed of seven taxa and the following directed links

MAT <-  matrix(c(0,0,1,1,0,0,0,
0,0,0,1,1,0,0,
0,0,0,0,0,1,1,
0,0,0,0,0,0,1,
0,0,0,0,0,0,0,
0,0,0,0,0,0,0,
0,0,0,0,0,0,0), byrow = TRUE, nrow = 7)
colnames(MAT) <- rownames(MAT) <- c("A","B","C","D","E","F","G")

## construction of the rooted food web
MAT_RFW <- rooted_FW(MAT)
colnames(MAT_RFW) <- rownames(MAT_RFW) <- c("root","A","B","C","D","E","F","G")
MAT_RFW

## assembly of the dominator matrix
MAT_DM <- dominator_matrix(MAT_RFW)
colnames(MAT_DM) <- rownames(MAT_DM) <- c("root","A","B","C","D","E","F","G")
MAT_DM

## dominator matrix without the "root" node
MAT_DMNR <- MAT_DM[-1,-1]

## attack sensitivity ---> largest fraction of food web taxa disconnected from the root
AS <- max(apply(MAT_DMNR,2,sum))/nrow(MAT_DMNR)

## error sensitivity ---> average number of taxa disconnected from the root because of a random removal
ES <- mean(apply(MAT_DMNR,2,sum)/nrow(MAT_DMNR))
## it is analogous to ---> sum(apply(MAT_DMNR,2,sum)/nrow(MAT_DMNR)^2)


##
##
##
##


## (1) Gulf of California food web ---> analysis based on the full topology

GC_db <- read.table(file = "ESM2---Olmo_Gilabert_2019.txt", header = FALSE,
colClasses = c("character","character","numeric"))

## name of all nodes
all_nodes <- unique(c(GC_db[,1],GC_db[,2]))

## total number of nodes
l_nodes <- length(all_nodes)

## preparing an empty matrix to be filled with info from the edgelist
GC_M <- matrix(rep(0,l_nodes^2), nrow = l_nodes)
rownames(GC_M) <- colnames(GC_M) <- all_nodes

## filling in the matrix
for(i in 1:nrow(GC_db)){
	GC_M[GC_db[i,1],GC_db[i,2]] <- GC_db[i,3]
}

## checking that the matrix includes feeding preferences normalized (column sum = 1) ---> condition verified
## apply(GC_M,2,sum)

## conversion of the feeding preferences matrix into a binary matrix (i.e. only with presence/absence of links)
all_elements <- rep(0,l_nodes^2)
all_elements[which(GC_M!=0)] <- 1
GC_B <- matrix(all_elements, byrow = FALSE, nrow = l_nodes)

## construction of the rooted food web
GC_RFW <- rooted_FW(GC_B)
colnames(GC_RFW) <- rownames(GC_RFW) <- c("root",colnames(GC_M))

## assembly of the dominator matrix
GC_DM <- dominator_matrix(GC_RFW)
colnames(GC_DM) <- rownames(GC_DM) <- c("root",colnames(GC_M))

## dominator matrix without the "root" node
GC_DMNR <- GC_DM[-1,-1]

## error sensitivity ---> average number of taxa disconnected from the root because of a random removal
ES_GC <- mean(apply(GC_DMNR,2,sum)/nrow(GC_DMNR))
## it is analogous to ---> sum(apply(GC_DMNR,2,sum)/nrow(GC_DMNR)^2)

## attack sensitivity ---> largest fraction of food web taxa disconnected from the root
AS_GC <- max(apply(GC_DMNR,2,sum))/nrow(GC_DMNR)


##
##
##
##


## Gulf of California food web ---> progressive removal of all links below a given threshold
## the procedure continues irrespective of the original food web fragmentation into subgraphs ---> based on bottom-up approach

## conversion of the feeding preferences matrix into a binary matrix (i.e. only with presence/absence of links)
all_elements <- rep(0,l_nodes^2)
all_elements[which(GC_M > 0)] <- 1
GC_B <- matrix(all_elements, byrow = FALSE, nrow = l_nodes)
diag(GC_B) <- rep(0,nrow(GC_B))

## total number of primary producers in the food web
prim_prod_el <- which(apply(GC_B,2,sum)==0)
prim_prod <- length(prim_prod_el)
prim_prod_max <- max(prim_prod_el)

## total number of nodes in the original food web
tot_nodes <- nrow(GC_B)

## k <- 1
tot_steps <- 100

## condition to verify that ES and AS are calculated on the entire food web
## (i.e. no taxa disconnected from the root)
for(k in 1:tot_steps){
	##
	print(k/tot_steps)
	##
	## construction of the rooted food web
	GC_RFW <- rooted_FW(GC_B)
	## colnames(GC_RFW) <- rownames(GC_RFW) <- c("root",colnames(GC_M))
	##
	## assembly of the dominator matrix
	GC_DM <- dominator_matrix(GC_RFW)
	## colnames(GC_DM) <- rownames(GC_DM) <- c("root",colnames(GC_M))
	##
	## dominator matrix without the "root" node
	## GC_DMNR <- GC_DM[-1,-1]
	GC_DMNR <- GC_DM
	##
	## total number of nodes in the giant components
	nodi_giant <- nrow(GC_DMNR) - 1
	##
	## error sensitivity ---> average number of taxa disconnected from the root because of a random removal
	## ES_GC <- sum(apply(GC_DMNR[-1,-1],2,sum)+(tot_nodes - nrow(GC_B)))/tot_nodes^2
	ES_NN <- sum(apply(GC_DMNR[-1,-1],2,sum))
	ES_GC <- mean(apply(GC_DMNR[-1,-1],2,sum))/nodi_giant
	##
	## attack sensitivity ---> largest fraction of food web taxa disconnected from the root
	AS_NN <- max(apply(GC_DMNR[-1,-1],2,sum))
	AS_GC <- (max(apply(GC_DMNR[-1,-1],2,sum))+(tot_nodes - nrow(GC_B)))/tot_nodes
	## AS_GC <- max(apply(GC_DMNR[-1,-1],2,sum))/nodi_giant
	##
	##
	{
	if(k == 1){
		NN_v <- nodi_giant
		##
		ES_v <- ES_GC
		ES_n <- ES_NN
		##
		AS_v <- AS_GC
		AS_n <- AS_NN
		}
	else{
		NN_v <- c(NN_v,nodi_giant)
		##
		ES_v <- c(ES_v,ES_GC)
		ES_n <- c(ES_n,ES_NN)
		##
		AS_v <- c(AS_v,AS_GC)
		AS_n <- c(AS_n,AS_NN)
		}
	}
	##
	##
	thresh <- k/tot_steps
	##
	all_elements <- rep(0,l_nodes^2)
	all_elements[which(GC_M > thresh)] <- 1
	GC_B1 <- matrix(all_elements, byrow = FALSE, nrow = l_nodes)
	##
	diffe <- which(apply(GC_B1,2,sum) == 0)[which(which(apply(GC_B1,2,sum) == 0) > prim_prod_max)]
	GC_B2 <- GC_B1
	##
	while(length(diffe) > 0){
		GC_B2 <- GC_B2[-diffe,-diffe]
		diffe <- which(apply(GC_B2,2,sum) == 0)[which(which(apply(GC_B2,2,sum) == 0) > prim_prod_max)]
	}
	##
	GC_B <- GC_B2
}

plot(AS_v,col = "red", ylim = c(0,1), xlab = "diet fraction", ylab = "fraction nodes removed", type = "b", pch = 19)
lines(ES_v,col = "blue", type = "b", pch = 19)
legend(0.1, 1, legend = c("AS", "ES"), col = c("red", "blue"), lty = 1:2, cex = 1, box.lty = 0)

whole_db <- cbind(NN_v, ES_v, ES_n, AS_v, AS_n)


##
##
##
##


## function that performs the removal scenario to test error and attack sensitivity according to Allesina and Bodini (2004)
## Allesina, S. and Bodini, A., 2004. Who dominates whom in the ecosystem? Energy flow
## bottlenecks and cascading extinctions. Journal of Theoretical Biology 230, 351-358.

removal_allesina <- function(nome, matrice, passi){
	##
	##
	## identification of all feeding links
	l_nodes <- nrow(matrice)
	all_elements <- rep(0,l_nodes^2)
	all_elements[which(matrice > 0)] <- 1
	##
	## construction of the binary food web
	GC_B <- matrix(all_elements, byrow = FALSE, nrow = l_nodes)
	diag(GC_B) <- rep(0,nrow(GC_B))
	##
	## total number of primary producers in the food web
	prim_prod_el <- which(apply(GC_B,2,sum)==0)
	prim_prod <- length(prim_prod_el)
	prim_prod_max <- max(prim_prod_el)
	##
	## total number of nodes in the original food web
	tot_nodes <- nrow(GC_B)
	##
	##
	for(k in 1:passi){
		##
		## print(k/passi)
		##
		## construction of the rooted food web
		GC_RFW <- rooted_FW(GC_B)
		##
		## assembly of the dominator matrix
		GC_DM <- dominator_matrix(GC_RFW)
		GC_DMNR <- GC_DM
		##
		## total number of nodes in the giant components, excluding the "root"
		nodi_giant <- nrow(GC_DMNR) - 1
		##
		## error sensitivity ---> average number of taxa disconnected from the giant component "root" due to random removal
		ES_NN <- sum(apply(GC_DMNR[-1,-1],2,sum))
		ES_GC <- mean(apply(GC_DMNR[-1,-1],2,sum))/nodi_giant
		##
		## attack sensitivity ---> largest fraction of taxa disconnected from the "root" (value based on initial food web size)
		AS_NN <- max(apply(GC_DMNR[-1,-1],2,sum))
		AS_GC <- (max(apply(GC_DMNR[-1,-1],2,sum))+(tot_nodes - nrow(GC_B)))/tot_nodes
		##
		##
		{
		if(k == 1){
			NN_v <- nodi_giant
			##
			ES_v <- ES_GC
			ES_n <- ES_NN
			##
			AS_v <- AS_GC
			AS_n <- AS_NN
			}
		else{
			NN_v <- c(NN_v,nodi_giant)
			##
			ES_v <- c(ES_v,ES_GC)
			ES_n <- c(ES_n,ES_NN)
			##
			AS_v <- c(AS_v,AS_GC)
			AS_n <- c(AS_n,AS_NN)
			}
		}
		##
		##
		thresh <- k/passi
		##
		all_elements <- rep(0,l_nodes^2)
		all_elements[which(matrice > thresh)] <- 1
		GC_B1 <- matrix(all_elements, byrow = FALSE, nrow = l_nodes)
		##
		diffe <- which(apply(GC_B1,2,sum) == 0)[which(which(apply(GC_B1,2,sum) == 0) > prim_prod_max)]
		GC_B2 <- GC_B1
		##
		while(length(diffe) > 0){
			GC_B2 <- GC_B2[-diffe,-diffe]
			diffe <- which(apply(GC_B2,2,sum) == 0)[which(which(apply(GC_B2,2,sum) == 0) > prim_prod_max)]
		}
		##
		GC_B <- GC_B2
	}
	##
	##
	x11()
	plot(AS_v,col = "red", ylim = c(0,1), main = nome, xlab = "diet fraction",
	ylab = "fraction nodes removed", type = "b", pch = 19)
	lines(ES_v,col = "blue", type = "b", pch = 19)
	legend(0.1, 1, legend = c("AS", "ES"), col = c("red", "blue"), lty = 1:2, cex = 1, box.lty = 0)
	##
	whole_db <- cbind(NN_v, ES_v, ES_n, AS_v, AS_n)
	return(whole_db)
}


#
#Guerreiro removal
#taking combinatorial approach

removal_guerreiro <- function(nome, matrice, passi){
	##
	##
	## identification of all feeding links
	l_nodes <- nrow(matrice)
	all_elements <- rep(0,l_nodes^2)
	all_elements[which(matrice > 0)] <- 1
	##
	## construction of the binary food web
	GC_B <- matrix(all_elements, byrow = FALSE, nrow = l_nodes)
	diag(GC_B) <- rep(0,nrow(GC_B))
	##
	## total number of primary producers in the food web
	prim_prod_el <- which(apply(GC_B,2,sum)==0)
	prim_prod <- length(prim_prod_el)
	prim_prod_max <- max(prim_prod_el)
	##
	## total number of nodes in the original food web
	tot_nodes <- nrow(GC_B)
	##
	##
	for(k in 1:passi){
		##
		## print(k/passi)
		##
		## construction of the rooted food web
		GC_RFW <- rooted_FW(GC_B)
		##
		## assembly of the dominator matrix
		GC_DM <- dominator_matrix(GC_RFW)
		GC_DMNR <- GC_DM
		##
		## total number of nodes in the giant components, excluding the "root"
		nodi_giant <- nrow(GC_DMNR) - 1
		##
		## error sensitivity ---> average number of taxa disconnected from the giant component "root" due to random removal
		ES_NN <- sum(apply(GC_DMNR[-1,-1],2,sum))
		ES_GC <- mean(apply(GC_DMNR[-1,-1],2,sum))/nodi_giant
		##
		## attack sensitivity ---> largest fraction of taxa disconnected from the "root" (value based on initial food web size)
		AS_NN <- max(apply(GC_DMNR[-1,-1],2,sum))
		AS_GC <- (max(apply(GC_DMNR[-1,-1],2,sum))+(tot_nodes - nrow(GC_B)))/tot_nodes
		##
		##
		{
		if(k == 1){
			NN_v <- nodi_giant
			##
			ES_v <- ES_GC
			ES_n <- ES_NN
			##
			AS_v <- AS_GC
			AS_n <- AS_NN
			}
		else{
			NN_v <- c(NN_v,nodi_giant)
			##
			ES_v <- c(ES_v,ES_GC)
			ES_n <- c(ES_n,ES_NN)
			##
			AS_v <- c(AS_v,AS_GC)
			AS_n <- c(AS_n,AS_NN)
			}
		}
		##
		#########################
		thresh <- k/passi
		##
		all_elements <- rep(0,l_nodes^2)
		all_elements[which(matrice > thresh)] <- 1
		GC_B1 <- matrix(all_elements, byrow = FALSE, nrow = l_nodes)
		##
		diffe <- which(apply(GC_B1,2,sum) == 0)[which(which(apply(GC_B1,2,sum) == 0) > prim_prod_max)]
		GC_B2 <- GC_B1
		##
		while(length(diffe) > 0){
			GC_B2 <- GC_B2[-diffe,-diffe]
			diffe <- which(apply(GC_B2,2,sum) == 0)[which(which(apply(GC_B2,2,sum) == 0) > prim_prod_max)]
		#########################
		}
		##
		GC_B <- GC_B2
	}
	##
	##
	x11()
	plot(AS_v,col = "red", ylim = c(0,1), main = nome, xlab = "diet fraction",
	ylab = "fraction nodes removed", type = "b", pch = 19)
	lines(ES_v,col = "blue", type = "b", pch = 19)
	legend(0.1, 1, legend = c("AS", "ES"), col = c("red", "blue"), lty = 1:2, cex = 1, box.lty = 0)
	##
	whole_db <- cbind(NN_v, ES_v, ES_n, AS_v, AS_n)
	return(whole_db)
}


gulf_california <- removal_allesina(nome = "Gulf of California", matrice = GC_M, passi = 100)


##
##
##
##


## import data on all weighted networks and calculation of feeding preferences

nomi_folders <- c("balticsea", "caspiansea", "chesa", "cryscon", "cypressdry",
"evergladesdry", "floridadry", "goose", "huizache", "mangrovedry", "marks",
"mondego", "narra", "northsea", "prince")

l_folders <- length(nomi_folders)

tutti_w <- l_folders + 1
tutti_net <- as.list(rep(NA,tutti_w))

for(i in 1:l_folders){
	##
	## set the working directory
	percorso <- paste("D:/People/Miguel_Fernandes_Guerreiro/R_Allesina/", nomi_folders[i], sep = "")
	setwd(percorso)
	##
	## import the food web and its size
	FW <- read.table(file = "data.txt", header = FALSE, sep = "\t")
	##
	FW_l <- nrow(FW)
	##
	## total amount of energy/matter consumed by each compartment
	FW_cs <- apply(FW,2,sum)
	##
	## creation of the matrix with feeding preferences
	FW_G <- matrix(rep(0,FW_l^2), nrow = FW_l)
	for(k in 1:FW_l){
		for(j in 1:FW_l){
			if(FW[k,j]!=0)FW_G[k,j] <- FW[k,j]/FW_cs[j]
		}
	}
	##
	## storing the matrix with feeding preferences in the list
	tutti_net[[i]] <- FW_G
}

GC_M_nonames <- GC_M
colnames(GC_M_nonames) <- rownames(GC_M_nonames) <- NULL
tutti_net[[length(tutti_net)]] <- GC_M_nonames
names(tutti_net) <- c(nomi_folders,"california")

## for(i in 1:length(tutti_net))print(apply(tutti_net[[i]],2,sum))

## cleaning of the dataset to remove consumers that do not feed from other compartments in the system
## (i.e. those that from a pure assessment of the initial topology would be assumed to be primary producers)
## ---> they are present and problematic in all food webs but "chesa", "cryscon", "floridadry"


##
##
##
##


## complete analysis for attack and error sensitivity in the 11 weighted networks ---> Allesina's algorithm

removal_allesina_OUT <- as.list(rep(NA,tutti_w))

nomi_full <- c("western Baltic Sea", "Caspian Sea", "Chesapeake Bay", "Crystal River Creek (control)",
"Cypress (dry season)", "Everglades (dry season)", "Floriday Bay (dry season)", "Goose Creek Bay",
"Huizache–Caimanero Lagoon", "Mangroves (dry season)", "St. Marks River", "Mondego Estuary",
"Narragansett Bay", "southern North Sea", "Prince William Sound", "Gulf of California")

sim_steps <- 100

for(k in 1:tutti_w){
	removal_allesina_OUT[[k]] <- removal_allesina(nome = nomi_full[k],
	matrice = tutti_net[[k]], passi = sim_steps)
}

save.image("D:/People/Miguel_Fernandes_Guerreiro/R_Allesina/2021-02-10---output_v0.RData")
