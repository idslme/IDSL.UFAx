UFAx_score_coefficient_corrector <- function(annotated_molf_address, maxNEME, Score_coeff, number_processing_threads = 1) {
  ##
  file_names <- gsub(".Rdata$", "", annotated_molf_address)
  annotated_molf <- IDSL.IPA::loadRdata(annotated_molf_address)
  colnames_annotated_molf <- colnames(annotated_molf)
  ##
  PCS <- as.numeric(annotated_molf[, 9])
  RCS <- as.numeric(annotated_molf[, 13])
  NEME <- as.numeric(annotated_molf[, 8])
  R13C_PL <- as.numeric(annotated_molf[, 10])
  R13C_IP <- as.numeric(annotated_molf[, 11])
  size_IP <- as.numeric(annotated_molf[, 2])
  ##
  IdentificationScore <- identification_score(Score_coeff, size_IP, PCS, RCS, NEME, maxNEME, R13C_PL, R13C_IP)
  ##
  IPApeaks <- as.numeric(annotated_molf[, 1])
  xDiff <- c(0, which(abs(diff(IPApeaks)) > 0), length(size_IP))
  ##
  annotated_molf_updated_call <- function(j) {
    x_p <- (xDiff[j] + 1):xDiff[j + 1]
    order_A <- order(IdentificationScore[x_p])
    A <- annotated_molf[x_p[order_A], ]
    A[, 14] <- seq(1:length(x_p))
    return(A)
  }
  ##
  if (number_processing_threads == 1) {
    ##
    annotated_molf_updated <- do.call(rbind, lapply(1:(length(xDiff) - 1), function(i) {
      annotated_molf_updated_call(i)
    }))
    ##
  } else {
    osType <- Sys.info()[['sysname']]
    ##
    if (osType == "Windows") {
      clust <- makeCluster(number_processing_threads)
      registerDoParallel(clust)
      ##
      annotated_molf_updated <- foreach(i = 1:(length(xDiff) - 1), .combine = 'rbind', .verbose = FALSE) %dopar% {
        annotated_molf_updated_call(i)
      }
      ##
      stopCluster(clust)
      ##
    } else if (osType == "Linux") {
      ##
      annotated_molf_updated <- do.call(rbind, mclapply(1:(length(xDiff) - 1), function(i) {
        annotated_molf_updated_call(i)
      }, mc.cores = number_processing_threads))
      ##
      closeAllConnections()
    }
  }
  ##
  annotated_molf_updated <- data.frame(annotated_molf_updated)
  colnames(annotated_molf_updated) <- colnames_annotated_molf
  rownames(annotated_molf_updated) <- NULL
  ##
  save(annotated_molf_updated, file = paste0(file_names, "_updated.Rdata"))
  write.csv(annotated_molf_updated, file = paste0(file_names, "_updated.csv"), row.names = FALSE)
}