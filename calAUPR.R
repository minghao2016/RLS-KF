# function to calculate the AUCR and AUC
calAUPR <- function(obsLabel, predProb) {
  unique_labels <- unique(obsLabel)
	if (length(unique_labels) != 2) stop("The first argument 'obsLabel' should be two classes!\n")
	# calculate AUC using ROCR
	# library('ROCR')
  pred <- ROCR::prediction(predProb, obsLabel)
  perf <- ROCR::performance(pred, "auc")
  auc <- as.numeric(perf@y.values)
	
  # Calculate AUPR using ROCR
  perf <- ROCR::performance(pred, 'rec', 'prec')
  Precision <- unlist(perf@x.values)
  Recall <- unlist(perf@y.values)

  # library('MESS')
  aupr_spline <- try(MESS::auc(Recall, Precision, type = 'spline'), silent = TRUE)

	# Save the result
	statRes <- matrix(0, nrow = 1, ncol = 2)
	colnames(statRes) <- c("auc", "aupr")

  if (class(aupr_spline) == 'try-error') {
    # library('Bolstad2')
    # uses Simpson's rule: numerical integrating, solve the area
    aupr_simpson <- Bolstad2::sintegral(Recall, Precision)$int
		statRes[, "auc"] <- auc
		statRes[, "aupr"] <- aupr_simpson
  } else {
    statRes[, "auc"] <- auc
		statRes[, "aupr"] <- aupr_spline
  }
	return(statRes)
}






