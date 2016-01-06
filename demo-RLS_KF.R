# demo for RLS-KF algorithm of Drug-Target-Interaction (DTI) prediction

#setwd("To Your Directory Including All Required Files")

rm(list = ls())

# libraries for compiling C++ codes used in this work
library(Rcpp) 
library(RcppArmadillo)

# library for checking positive semi-definite
library(matrixcalc)
 
# libraries for calAUPR
library(MESS)
library(pracma)
library(ROCR) 
library(Bolstad2)

# sourceCpp
sourceCpp("fastKgipMat.cpp")
sourceCpp("fastKF.cpp")
sourceCpp("fastSolve.cpp")

# source
source("RLS_KF.R")
source("calAUPR.R")

###################################################################################################
# You just modify the partfn to different data sets
# file name to be used: nr, gpcr, ic, e
partfn = "nr"
# Take long time for Enzyme dataset
if (partfn == "e") cat("Need several hours to finish the big [Enzymes] data set, please be patient!\n")
###################################################################################################

switch(partfn,
	nr = {
	  # y
	  yFn <- paste0(partfn, "_admat_dgc.txt")
	  y <- read.table(yFn)
	  # simmatCompd
	  simCompdFn <- paste0(partfn, "_simmat_dc.txt")
	  simmatCompd <- read.table(simCompdFn)
	  # simmatTarget
	  simTargetFn <- paste0(partfn, "_simmat_dg.txt")
	  simmatTarget <- read.table(simTargetFn)
	  # convert into matrix
	  y <- as.matrix(y)
	  simmatCompd <- as.matrix(simmatCompd)
	  simmatTarget <- as.matrix(simmatTarget)
	  # check matrix symmetric
	  if (!isSymmetric(simmatCompd)) simmatCompd <- (simmatCompd + t(simmatCompd))/2
	  # check matrix positive semi-definite 
	  epsilon <- 0.1
	  while (!is.positive.semi.definite(simmatCompd)) simmatCompd <- simmatCompd + epsilon * diag(nrow(simmatCompd))
	  # check matrix symmetric
	  if (!isSymmetric(simmatTarget)) simmatTarget <- (simmatTarget + t(simmatTarget))/2
	  # check matrix positive semi-definite 
	  epsilon <- 0.1
	  while (!is.positive.semi.definite(simmatTarget)) simmatTarget <- simmatTarget + epsilon * diag(nrow(simmatTarget))
	},
	gpcr = {
	  # y
	  yFn <- paste0(partfn, "_admat_dgc.txt")
	  y <- read.table(yFn)
	  # simmatCompd
	  simCompdFn <- paste0(partfn, "_simmat_dc.txt")
	  simmatCompd <- read.table(simCompdFn)
	  # simmatTarget
	  simTargetFn <- paste0(partfn, "_simmat_dg.txt")
	  simmatTarget <- read.table(simTargetFn)
		# convert into matrix
	  y <- as.matrix(y)
	  simmatCompd <- as.matrix(simmatCompd)
	  simmatTarget <- as.matrix(simmatTarget)
		# check matrix symmetric
	  if (!isSymmetric(simmatCompd)) simmatCompd <- (simmatCompd + t(simmatCompd))/2
	  # check matrix positive semi-definite
	  epsilon <- 0.1
	  while (!is.positive.semi.definite(simmatCompd)) simmatCompd <- simmatCompd + epsilon * diag(nrow(simmatCompd))
	  # check matrix symmetric
	  if (!isSymmetric(simmatTarget)) simmatTarget <- (simmatTarget + t(simmatTarget))/2
	  # check matrix positive semi-definite 
	  epsilon <- 0.1
	  while (!is.positive.semi.definite(simmatTarget)) simmatTarget <- simmatTarget + epsilon * diag(nrow(simmatTarget))
	},
  ic = {
		# y
		yFn <- paste0(partfn, "_admat_dgc.txt")
		y <- read.table(yFn)
		# simmatCompd
		simCompdFn <- paste0(partfn, "_simmat_dc.txt")
		simmatCompd <- read.table(simCompdFn)
		# simmatTarget
		simTargetFn <- paste0(partfn, "_simmat_dg.txt")
		simmatTarget <- read.table(simTargetFn)
		# convert into matrix
	  y <- as.matrix(y)
	  simmatCompd <- as.matrix(simmatCompd)
	  simmatTarget <- as.matrix(simmatTarget)
		# check matrix symmetric
	  if (!isSymmetric(simmatCompd)) simmatCompd <- (simmatCompd + t(simmatCompd))/2
	  # check matrix positive semi-definite 
	  epsilon <- 0.1
	  while (!is.positive.semi.definite(simmatCompd)) simmatCompd <- simmatCompd + epsilon * diag(nrow(simmatCompd))
	  # check matrix symmetric
	  if (!isSymmetric(simmatTarget)) simmatTarget <- (simmatTarget + t(simmatTarget))/2
	  # check matrix positive semi-definite 
	  epsilon <- 0.1
	  while (!is.positive.semi.definite(simmatTarget)) simmatTarget <- simmatTarget + epsilon * diag(nrow(simmatTarget))
  },
  e = {
		# y
		yFn <- paste0(partfn, "_admat_dgc.txt")
		y <- read.table(yFn)
		# simmatCompd
		simCompdFn <- paste0(partfn, "_simmat_dc.txt")
		simmatCompd <- read.table(simCompdFn)
		# simmatTarget
		simTargetFn <- paste0(partfn, "_simmat_dg.txt")
		simmatTarget <- read.table(simTargetFn)
		# convert into matrix
		y <- as.matrix(y)
		simmatCompd <- as.matrix(simmatCompd)
		simmatTarget <- as.matrix(simmatTarget)
		# check matrix symmetric
		if (!isSymmetric(simmatCompd)) simmatCompd <- (simmatCompd + t(simmatCompd))/2
		# check matrix positive semi-definite 
		epsilon <- 0.1
		while (!is.positive.semi.definite(simmatCompd)) simmatCompd <- simmatCompd + epsilon * diag(nrow(simmatCompd))
		# check matrix symmetric
		if (!isSymmetric(simmatTarget)) simmatTarget <- (simmatTarget + t(simmatTarget))/2
		# check matrix positive semi-definite 
		epsilon <- 0.1
		while (!is.positive.semi.definite(simmatTarget)) simmatTarget <- simmatTarget + epsilon * diag(nrow(simmatTarget))
  },
  stop("partfn should be one of: {nr, gpcr, ic, enzyme}\n")
)

# parameters
gamma0   = 1    
lambda   = 1    
nfold    = 10
# k for KF
numNeig  = 4    
# t for KF
numIter  = 2

# number of replicated runs
nreps <- 10 
AUC_ave <- vector(length = nreps)
AUPR_ave <- vector(length = nreps)
AUC_max <- vector(length = nreps)
AUPR_max <- vector(length = nreps)
prob_ave <- NULL

for (i_rep in 1:nreps) {
  cat(i_rep, "/", nreps, "\n")
	cat("data set:", partfn, "\n")
	cat("k =", numNeig, "\n")
	cat("t =", numIter, "\n")
	flush.console()
	# (1) Prediction based on the target similarity
	yTarget <- y
	numRows <- nrow(yTarget)
	numCols <- ncol(yTarget)
	predBasedTarget <- matrix(0, nr = numRows, nc = numCols)
	myColPrediction <- matrix(0, nr = numRows, nc = 1)

	# Kernel matrix from similarity matrix
	k4simmat <- simmatTarget

	# Cross-validation folds
	lenSeg <- ceiling(numRows / nfold)
	incomplete <- nfold * lenSeg - numRows     
	complete <- nfold - incomplete                  
	inds <- matrix(c(sample(1:numRows), rep(NA, incomplete)), nrow = lenSeg, byrow = TRUE)
	folds <- lapply(as.data.frame(inds), function(x) c(na.omit(x)))

	# Main fold prediction function
	for (i in 1:numCols) {
		currY <- yTarget[, i]
		# Used when currY are all zeros: inferred currY = ytrMat %*% currSim2
		currSim2 <- simmatCompd[, i]
		for (j in 1:nfold) {
			idxTe <- folds[[j]]
			idxTr <- setdiff(1:numRows, idxTe)
			ytrMat <- yTarget
			# Put current 'test set' to zeros
			ytrMat[idxTe, i] <- 0
			currYsum <- sum(ytrMat[, i])
			# Prediction for each fold
			myColPrediction[idxTe, 1] <- RLS_KF( 
				currY        = currY,
				ytrMat       = ytrMat,
				currYsum     = currYsum,
				currSim2     = currSim2,
				idxTr        = idxTr,
				idxTe        = idxTe,
				simmat       = simmatTarget,
				simmat2      = simmatCompd, 
				gamma0       = gamma0, 
				numNeig      = numNeig, 
				numIter      = numIter,
				lambda       = lambda, 
				k4simmat     = k4simmat)
		} 
		predBasedTarget[, i] <- myColPrediction
	}

	# (2) Prediction based on the compound similarity
	yCompd <- t(y)
	numRows <- nrow(yCompd)
	numCols <- ncol(yCompd)
	predBasedCompd <- matrix(0, nr = numRows, nc = numCols)
	myColPrediction <- matrix(0, nr = numRows, nc = 1)

	# Kernel matrix from similarity matrix
	k4simmat <- simmatCompd

	# Cross-validation folds
	lenSeg <- ceiling(numRows / nfold)
	incomplete <- nfold * lenSeg - numRows     
	complete <- nfold - incomplete                  
	inds <- matrix(c(sample(1:numRows), rep(NA, incomplete)), nrow = lenSeg, byrow = TRUE)
	folds <- lapply(as.data.frame(inds), function(x) c(na.omit(x)))

	# Main fold prediction function
	for (i in 1:numCols) {
		currY <- yCompd[, i]
		# Used when currY are all zeros: inferred currY = ytrMat %*% currSim2
		currSim2 <- simmatTarget[, i]
		for (j in 1:nfold) {
			idxTe <- folds[[j]]
			idxTr <- setdiff(1:numRows, idxTe)
			ytrMat <- yCompd
			# Put current 'test set' to zeros
			ytrMat[idxTe, i] <- 0
			currYsum <- sum(ytrMat[, i])
			# Prediction for each fold
			myColPrediction[idxTe, 1] <- RLS_KF( 
				currY      = currY,
				ytrMat     = ytrMat,
				currYsum   = currYsum,
				currSim2   = currSim2,
				idxTr      = idxTr,
				idxTe      = idxTe,
				simmat     = simmatCompd,
				simmat2    = simmatTarget, 
				gamma0     = gamma0, 
				numNeig    = numNeig, 
				numIter    = numIter,
				lambda     = lambda, 
				k4simmat   = k4simmat)
		}  
		predBasedCompd[, i] <- myColPrediction
	}

	# (3) statistics
	yLabel <- as.vector(y)

	# (3-1) based on average
	predProb_ave <- (predBasedTarget + t(predBasedCompd))/2
	finalPred <- as.vector(predProb_ave)
	# used for nrEXt
	prob_ave <- cbind(prob_ave, finalPred)
	statRes_ave <- calAUPR(yLabel, finalPred)
	AUC_ave[i_rep] <- statRes_ave[, "auc"]
	AUPR_ave[i_rep] <- statRes_ave[, "aupr"]
	# (3-2) based on maximum
	predBasedCompd <- t(predBasedCompd)
	predProb_max <- ifelse(predBasedTarget > predBasedCompd, predBasedTarget, predBasedCompd)
	finalPred <- as.vector(predProb_max)
	statRes_max <- calAUPR(yLabel, finalPred)
  AUC_max[i_rep] <- statRes_max[, "auc"]
	AUPR_max[i_rep] <- statRes_max[, "aupr"]
}
# (4) show results
# (4-1) average result
meanAUC_ave <- mean(AUC_ave)
sdAUC_ave <- sd(AUC_ave)
meanAUPR_ave <- mean(AUPR_ave)
sdAUPR_ave <- sd(AUPR_ave)
# (4-2) maximum result
meanAUC_max <- mean(AUC_max)
sdAUC_max <- sd(AUC_max)
meanAUPR_max <- mean(AUPR_max)
sdAUPR_max <- sd(AUPR_max)

# output the main results
cat("Current data set is: ", partfn, "\n")
cat("mean AUC and sd based on average: ", meanAUC_ave, "+/-", sdAUC_ave, "\n")
cat("mean AUPR and sd based on average: ", meanAUPR_ave, "+/-", sdAUPR_ave, "\n")

cat("mean AUC and sd based on maximum: ", meanAUC_max, "+/-", sdAUC_max, "\n")
cat("mean AUPR and sd based on maximum: ", meanAUPR_max, "+/-", sdAUPR_max, "\n")
