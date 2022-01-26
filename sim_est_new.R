rm(list=ls())
#setwd('c:/Study/Stratified-mult-GGM/Codes/source codes/')
source('jsem.R')
source('Generator.R')
source('l1LS_Main.R')
source('Objval.R')
source('JMLE.R')

library(glasso)
library(parallel)

##### Function to calculate evaluation metrics
evaluate = function(TP,TN,FP,FN){
  SEN = TP/(TP+FN)
  SPE = TN/(TN+FP)
  MCC = (TP*TN - FP*FN)/(sqrt(TP+FP)*sqrt(TP+FN)*sqrt(TN+FP)*sqrt(TN+FN))
  F1 = 2*TP/(2*TP+FP+FN)
  c(SEN,SPE,MCC,F1)
}

##### Common wrapper function
get.outputs = function(n=100, subnetSize.X=rep(10,2), subnetSize.E=rep(10,2),
                       sparsity.B=5, sparsity.Theta=5, K=5, nrep=50, filename=NULL){
  
  ## Set up some quantities
  group = rbind(
    c(1, 2),
    c(1, 4),
    c(3, 2),
    c(3, 4),
    c(5, 2)
  )  # grouping pattern
  p = sum(subnetSize.X)
  q = sum(subnetSize.E)
  
  loopfun = function(rep){
    set.seed(1e3*rep)
    
    ## Generate data *******************************************************
    # **********************************************************************
    X.layer = GenerateLayer(n, subnetSize.X, group, D=1, sparsity=sparsity.Theta/p)
    E.layer = GenerateLayer(n, subnetSize.E, group, D=1, sparsity=sparsity.Theta/q)
    
    ## generate group structure for coef array
    B0.group.array = array(0, c(p,q,K))
    g = 1
    for(i in 1:p){
      for(j in 1:q){
        B0.group.array[i,j,] = g
        g = g+1
      }
    }
    B0.array = CoefArray(B0.group.array)
    Theta0.array = array(0, c(q,q,K))
    for(k in 1:K){
      Theta0.array[,,k] = with(E.layer,
                               diag(diag(Omega[[k]])^(-0.5)) %*% Omega[[k]] %*% diag(diag(Omega[[k]])^(-0.5)))
    }
    
    ## make Y-layer
    Y.layer = E.layer
    for(k in 1:K){
      Y.layer$data[[k]] = X.layer$data[[k]] %*% B0.array[,,k] + E.layer$data[[k]]
    }
    
    ##### Given: X.list, Y.list, B.groups, Theta.groups
    Y.list = lapply(Y.layer$data, as.matrix)
    Y.indices = Y.layer$indices
    Theta.groups = Y.layer$groups
    X.list = lapply(X.layer$data, as.matrix)
    
    Theta.group.array = array(0, c(q,q,K))
    for(j in 1:q){
      Theta.group.array[j,-j,] = Y.layer$groups[[j]]
    }
    
    ## Obtain JMMLE fit ****************************************************
    # **********************************************************************
    ## tune JMMLE model
    lambda.vec = sqrt(log(p)/n) * seq(1.8, 0.4, -0.2)
    model.list = vector("list", length(lambda.vec))
    nlambda = length(lambda.vec)
    
    ## get all models
    loopfun1 = function(m){
      cat("Performing JMMLE for lambda =",lambda.vec[m],".....\n")
      jmmle.1step(Y.list, Y.indices, X.list, B.group.array=B0.group.array, Theta.groups=Theta.groups,
                  lambda = lambda.vec[m],
                  gamma = sqrt(log(q)/n) * seq(1, 0.4, -0.1),
                  init.option=1, tol=1e-3, VERBOSE=F)
      cat("done\n")
    }
    #model.list <- mclapply(1:nlambda, loopfun1, mc.cores=nlambda)
    model.list <- lapply(1:nlambda, loopfun1)
    
    ## calculate HBIC
    hbic.vec = rep(NA, nlambda)
    for(m in 1:nlambda){
      jmle.model = model.list[[m]]
      
      if(class(jmle.model)=="list"){ ## if no error in training the model
        SSE.vec = rep(0,K)
        hbic.pen.vec = rep(0,K)
        
        for(k in 1:K){
          nk = nrow(Y.list[[k]])
          Theta.k = jmle.model$Theta_refit$Theta[[k]]
          for(j in 1:q)
          {
            Theta.k[j,j] = 0
          }
          SSE.vec[k] = sum(diag(crossprod((Y.list[[k]] - X.list[[k]] %*%
                                             jmle.model$B.refit[,,k]) %*% (diag(1,q) - Theta.k))))/nk
          hbic.pen.vec[k] = log(log(nk))*log(q*(q-1)/2)/nk * sum(Theta.k != 0)/2 +
            log(log(nk))*log(p*q)/nk * sum(jmle.model$B.refit[,,k] != 0)
        }
        hbic.vec[m] = sum(SSE.vec) + sum(hbic.pen.vec) 
      }
    }
    
    ## select best model
    jmmle.model = model.list[[which.min(hbic.vec)]]
    
    ## Calculate metrics ***************************************************
    # **********************************************************************
    Theta_new.array = array(0, c(q,q,K))
    for(k in 1:K){
      Theta_new.array[,,k] = jmmle.model$Theta_refit$Theta[[k]]
    }
    
    TP.B = sum(B0.array != 0 & jmmle.model$B.refit != 0)
    TN.B = sum(B0.array == 0 & jmmle.model$B.refit == 0)
    FN.B = sum(B0.array != 0) - TP.B
    FP.B = sum(B0.array == 0) - TN.B
    TP.Theta = sum(Theta0.array != 0 & Theta_new.array != 0)
    TN.Theta = sum(Theta0.array == 0 & Theta_new.array == 0)
    FP.Theta = sum(Theta0.array != 0) - TP.Theta
    FN.Theta = sum(Theta0.array == 0) - TN.Theta
    cat("=============\nReplication",rep,"done!\n=============\n")
    rbind(c(evaluate(TP.B,TN.B,FP.B,FN.B),sqrt(sum((B0.array - jmmle.model$B.refit)^2)/sum(B0.array^2))),
          c(evaluate(TP.Theta,TN.Theta,FP.Theta,FN.Theta),
            sqrt(sum((Theta0.array - Theta_new.array)^2)/sum(Theta0.array^2)))
          )
  }
  
  # out.mat = mclapply(1:nrep, loopfun, mc.cores=8)
  out.mat = lapply(1:nrep, loopfun)
  # this mclapply sometimes gives errors for some lambdas
  # which stops corresponding cores.
  # *USE lapply HERE or wrap loopfun inside try() before using mclapply*
  if(is.null(filename)){
    filename = paste0("est_n",n,"p",p,"q",q,".rds")
  }
  saveRDS(out.mat, file=filename) # saves outputs as .rds file. read using readRDS()
}

##### a small simulation setup
get.outputs(n = 100, subnetSize.X = c(10, 10), subnetSize.E = c(10, 10))
