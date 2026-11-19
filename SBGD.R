#################################################
# Sequential Build of Genetic Diallel Matrices  #
# last runned in        27 of June, 2007        #
# With Dominance                                #
# defines the size of the diallel               #
# m: number of parents                          #
# c: number of possible crossings               #
# r: number of repeats of each crossing         #
#################################################


#####################
# Reading functions #
#####################

source("design.R")
source("relatd.R")

m     <- 16
r     <- 1
c     <- m*(m-1)/2

######################################################################
# Building genetic Relationship matrices to the Whole set of crosses #
######################################################################

genot  <- diallel(m,1)
cross  <- genot$cross

correl <- related(cross,m)

####################################################
# Assumption that only data on crosses are available 
# (no data on parents)
A      <- correl$A[(m+1):(m+c),(m+1):(m+c)] # Numerator Relationship
D      <- correl$D[(m+1):(m+c),(m+1):(m+c)] # Dominance Relationship

######################
# MIGUEL
# Adição de uma testemunha para ter t=121

At <- matrix(0, nrow=c+1, ncol=c+1)
At[1:c, 1:c] <- A
At[c+1, c+1] <- 1 # Diagonal da testemunha
A <- At

Dt<- matrix(0, nrow=c+1, ncol=c+1)
Dt[1:c, 1:c] <- D
Dt[c+1, c+1] <- 1
D <- Dt

############################

write(A,file="A10_Miguel.txt",ncol=c,sep=",")
write(D,file="D10_Miguel.txt",ncol=c,sep=",")


