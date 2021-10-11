pre <- function(x.pred,x.true,average = "macro"){
  x.true <- factor(x.true)
  x.pred <- factor(x.pred)
  x.pred <- factor(x.pred, levels = levels(x.true))
  conf <- table(x.pred, x.true)
  K <- nrow(conf)
  if(K == 2){
    TP <- conf[1,1]
    FP <- conf[1,2]
    precision <- TP / (TP + FP)
  }else if(K > 2){
    if(average == "micro"){
      precision <-  sum(diag(conf)) / sum(conf)
    }
    else if(average == "macro"){
      prei = 0
      for(k in 1:K){
        TP = conf[k,k]
        FP = sum(conf[k,]) - TP
        if(TP + FP !=0){
          prei = prei + TP / (TP + FP)
        }
      }
      precision = prei / K
    }
    else{
      stop("The option for average method is either micro or macro.")
    }
  }
  else{
    stop("A confusion matrix at least has two rows.")
  }
  return(precision)
}

rec <- function(x.pred, x.true, average="macro"){
  x.true <- factor(x.true)
  x.pred <- factor(x.pred)
  x.pred <- factor(x.pred, levels = levels(x.true))
  conf <- table(x.pred, x.true)
  K <- nrow(conf)
  if(K == 2){
    TP <- conf[1,1]
    FN <- conf[2,1]
    recall <- TP / (TP + FN)
  }else if(K > 2){
    if(average == "micro"){
      recall <-  sum(diag(conf)) / sum(conf)
    }else if(average == "macro"){
      reci <- 0
      for(k in 1:K){
        TP = conf[k,k]
        FN = sum(conf[,k]) - TP
        if(TP + FN != 0){
          reci <- reci + TP / (TP + FN)
        }
      }
      recall <- reci / K
    }else{
      stop("The option for average method is either micro or macro.")
    }
  }else{
    stop("A confusion matrix at least has two rows.")
  }
  return(recall)
}

Fscore <- function(x.pred, x.true, average = "macro", beta = 1){
    precision <- pre(x.pred, x.true, average)
    recall <- rec(x.pred, x.true, average)
    Fscore <- (beta^2 + 1) * precision * recall / (beta^2 * precision + recall)
    return(Fscore)
}

errRate <- function(x.pred, x.true){
  x.true <- factor(x.true)
  x.pred <- factor(x.pred)
  x.pred <- factor(x.pred, levels = levels(x.true))
  conf <- table(x.pred, x.true)
  K <- nrow(conf)
  if(K > 2){
    erri <- 0
    for(k in 1: K){
      TP <- conf[k,k]
      FP <- sum(conf[k,]) - TP
      FN <- sum(conf[,k]) - TP
      erri <- erri + (FP + FN) / sum(conf)
    }
    errRate <- erri / K
  }else{
    stop("Error Rate is an evaluation measure for multi-class classification.")
  }
  return(errRate)
}

accuracy <- function(x.pred, x.true){
  x.true <- factor(x.true)
  x.pred <- factor(x.pred)
  x.pred <- factor(x.pred, levels = levels(x.true))
  conf <- table(x.pred, x.true)
  K <- nrow(conf)
  if(K > 2){
    accui <- 0
    for(k in 1: K){
      TP <- conf[k,k]
      FP <- sum(conf[k,]) - TP
      FN <- sum(conf[,k]) - TP
      TN <- sum(conf) - TP - FP - FN
      accui <- accui + (TP + TN) / sum(conf)
    }
    accuracy <- accui / K
  }else{
    stop("Error Rate is an evaluation measure for multi-class classification.")
  }
  return(accuracy)
}
