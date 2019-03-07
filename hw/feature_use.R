feature.use <- function(urerf){
  fUse <- integer()
  for(i in 1:length(urerf$forest)){
    for(j in 1:length(urerf$forest[[i]]$matA)){
      if(!is.null(urerf$forest[[i]]$matA[[j]])){
        numFeatures <- length(urerf$forest[[i]]$matA[[j]])/2
        for(feature in 0:(numFeatures-1)){
          if(is.na(fUse[urerf$forest[[i]]$matA[[j]][feature*2+1]])){
            fUse[urerf$forest[[i]]$matA[[j]][feature*2+1]] <- 0
          }
          fUse[urerf$forest[[i]]$matA[[j]][feature*2+1]] <- fUse[urerf$forest[[i]]$matA[[j]][feature*2+1]] +1
        }
      }
    }
  }
  print("Feature Use")
  print (fUse)
}

