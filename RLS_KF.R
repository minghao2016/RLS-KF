# The algorithm (RLS-KF) [1] combines the NII [2] and KF [3] technique to perform drug-target-interaction (DTI) prediction.
# [1] M Hao, et al., Analytica Chimica Acta (submit)
# [2] JP Mei, et. al., Bioinformatics 29 (2013) 238-245.
# [3] B Wang, et. al., Nat. Methods 11 (2014) 333-337.

# Main function
RLS_KF <- function(
  currY        = currY,
  ytrMat       = ytrMat,        ## Y matrix after put zero to test set
  currYsum     = currYsum,      ## sum of current used Y column
  currSim2     = currSim2,      ## used to infer current Y values when they are all zeros
  idxTr        = idxTr,
  idxTe        = idxTe,
  simmat       = simmatTarget, 
  simmat2      = simmatCompd,
  gamma0       = 1,      ## gamma0, for Kgip (default 1)
  numNeig      = 3,      ## number of neighbours (default 3)
  numIter      = 2,      ## number of iterators, default 2
  lambda       = 1,      ## lambda, for RLS
  k4simmat     = k4simmat) 
{
	# Gaussian kernel matrix
  Kgip <- fastKgipMat(ytrMat, gamma0)
  # Kernel Fusion (KF)
	K <- fastKF(Kgip, k4simmat, numNeig, numIter)
  
  
	# NII
  if (currYsum == 0) {
    # Infer current Y column
    currY_inferred <- ytrMat %*% currSim2
    # Put 'test set' to zero
		currY_inferred[idxTe] <- 0
    # Normalize
    currY <- (currY_inferred - min(currY_inferred)) / (max(currY_inferred) - min(currY_inferred))
    # Training set
		K1 <- K[idxTr, idxTr, drop = FALSE]
    yTr <- currY[idxTr]
    # Test set
    K2 <- K[idxTe, idxTr, drop = FALSE]

    numTrainSet <- length(idxTr)
    
    # RLS model
    model <- fastSolve(K1, yTr, numTrainSet, lambda)
    
    # Predicted values for test set in each fold
    myFoldPrediction <- matrix(0, nrow = length(idxTe), ncol = 1)
                
    # Check if current test set is new.
    # e.g., for new coming target, which is no interaction with known drugs.
    # OR, for new coming drug, which is no interaction with known targets.
    hasZeros <- which(apply(ytrMat[idxTe, ], 1, sum) == 0)
        
    if (length(hasZeros) == 0) {
      # None of test set is new
      myFoldPrediction[, 1] <- K2 %*% model
    } else if (length(hasZeros) == length(idxTe)) {
      # All test sets are new: prediction based on the similarity itself rather than K2
      myFoldPrediction[, 1] <- simmat[idxTe, idxTr, drop = FALSE] %*% model
    } else {
      # Some of test set are new, and others are not
      # (1) Test set are new
      idxZeros <- idxTe[hasZeros]
      myFoldPrediction[hasZeros, 1] <- simmat[idxZeros, idxTr, drop = FALSE] %*% model
      # (2) Test set are not new
      idxNotZeros <- setdiff(1:length(idxTe), hasZeros)
      myFoldPrediction[idxNotZeros, 1] <- K2[idxNotZeros, , drop = FALSE] %*% model
    }
  } else { # no NII
    K1 <- K[idxTr, idxTr, drop = FALSE]
    yTr <- currY[idxTr]
    K2 <- K[idxTe, idxTr, drop = FALSE]
    numTrainSet <- length(idxTr)
    
    model <- fastSolve(K1, yTr, numTrainSet, lambda)
    myFoldPrediction <- matrix(0, nrow = length(idxTe), ncol = 1)

    hasZeros <- which(apply(ytrMat[idxTe, ], 1, sum) == 0)
        
    if (length(hasZeros) == 0) {
      myFoldPrediction[, 1] <- K2 %*% model
    } else if (length(hasZeros) == length(idxTe)) {
      myFoldPrediction[, 1] <- simmat[idxTe, idxTr, drop = FALSE] %*% model
    } else {
      idxZeros <- idxTe[hasZeros]
      myFoldPrediction[hasZeros, 1] <- simmat[idxZeros, idxTr, drop = FALSE] %*% model
      idxNotZeros <- setdiff(1:length(idxTe), hasZeros)
      myFoldPrediction[idxNotZeros, 1] <- K2[idxNotZeros, , drop = FALSE] %*% model
    }
  }
  return(myFoldPrediction)
}
