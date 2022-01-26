##########################################################
#    function: JSEM
## Purpose: to use joint group lasso to estimate multiple graphical models, 
##          given the known grouping structure. Thresholding corrects the group
##          misspecification 
# trainX,       # a data matrix sorted by group label, n by p
# trainY,       # a n-dim vector indicating model membership of the observed data (trainX)
# index,        # a K by p matrix of node-specific grouping pattern.
# lambda,       # penalty parameter; allow one parameter for every regression 
# delta1 = NULL,  #Thresholding parameter at the group level
# delta2 = NULL,  #thresholding parameter with group
## Output:
##    out: a list which contains the estimated adjacency matrices.
JSEM <- function(
  trainX,       
  trainY,        
  index,       
  lambda,    
  delta1 = NULL,
  delta2 = NULL,
  eps = 1e-06
){
  p = dim(trainX)[2]
  K = length(unique(trainY))
  
  lambda.grpreg = rep(lambda, p)
  
  # list observed data by model
  design = vector("list", K)
  Ahat = vector("list", K)
  for (k in 1:K){
    design[[k]] = trainX[which(trainY == k), ]
    Ahat[[k]] = matrix(0, p, p)
  }
  
  ## Start the loop for each node
  for (i in 1:p) {
    #cat("i= ", i, "\n")
    ## Create a block diagonal matrix for the variables
    list.x = vector("list", K)
    list.y = vector("list", K)
    for (k in 1:K) {
      list.x[[k]] = design[[k]][, -i]
      list.y[[k]] = design[[k]][, i]
    }
    #Duplicate columns for each k and add these extra columns to the end of the rearranged design X;
    #In the meantime, add the variable indices to the list of all variables.
    #Need to reorganize the estimate in the end.
    X = as.matrix(bdiag(list.x))
    new_X = X
    Y = unlist(list.y)
    
    ## To get the index label 
    myindex = sort.variables(index[[i]])
    X = X[, myindex$x.index]
    
    fit = grpreg(X, Y, myindex$v.index, family = "gaussian", penalty = "grLasso", lambda = lambda.grpreg[i])
    coeff = fit$beta[-1, ]    
    
    if (!is.null(delta1)){
      coeff = thr.group(coeff, myindex$v.index, delta1)
    }
    if (!is.null(delta2)){
      coeff = thr.coord(coeff, myindex$v.index, delta2)
    }
    
    tmp = matrix(coeff, nrow = K)
    for (j in 1:dim(tmp)[2]){
      tmp[, j] = tmp[myindex$b.index[, j], j]
    }
    
    for (k in 1:K){
      Ahat[[k]][i, -i] = tmp[k,]
    }     
  }
  
  # Symmetrize the thresholded coefficients; 
  # Get the symmetric adjacency matrices; 
  Theta = lapply(Ahat, FUN = symmetrize)
  Ahat = lapply(Theta, function(u) {u[abs(u)>eps] = 1; return(u)})
  
  return(list(Ahat = Ahat, Theta=Theta, lambda = lambda.grpreg))
}  

##########################################################
##  sel.lambda.jsem.R 
# Purpose: to select the tuning parameter for JSEM
sel.lambda.jsem <- function(
  trainX,
  testX, 
  trainY, 
  testY, 
  index, 
  lambda 
){
  p = dim(trainX)[2]
  K = length(unique(trainY))
  Info <- vector("list", K)
  
  # Find the sample size for each category
  n <- rep(0, K)
  for (k in 1:K){
    n[k] <- nrow(trainX[which(trainY == k), ]) 
  }
  
  N <- length(lambda)
  bic.score <- rep(0, N)
  likelihood <- rep(0, N)
  cat('Tuning JSEM model for',N,'tuning parameter values-\n')
  for (j in 1:N){
    cat("-", j)
    Ahat <- JSEM(trainX, trainY, index, lambda=lambda[j])$Ahat
    for (k in 1:K){
      Info[[k]] = zeroInd(Ahat[[k]], 1)$zeroArr
    }
    
    #In the second step of Joint group lasso, we use 0.1*log(p)/n as the default penalty for each graph.
    fit <- multi.glasso(trainX, trainY, lambda = 0.1*log(p)/n, zero = Info, BIC = T)
    bic.score[j] <- sum(fit$BIC)
    
    for (k in 1:K){
      data <- testX[which(testY == k), ]      
      empcov <- cov(data) 
      while (kappa(empcov) > 1e+2){
        empcov = empcov + 0.05 * diag(p)      
      }   
      likelihood[j] = likelihood[j] + matTr(empcov %*% fit$Omega[[k]]) - log(det(fit$Omega[[k]]))
    }
  }
  cat('\n')
  
  out <- list(BIC = bic.score, likelihood = likelihood)
  return(out)
}

##########################################################
## function: stabsel.jsem.R
## purpose: to perform selection of the estimated network based on stability criteria
## Refs: Meinshausen and Buhlmann - Stability selection - JRSSB - 2010
## Arguments:
##   X: a list of data matrices
##   cnt: the number of subsampling
##   lastar: the oracle lambda found in Glasso
##
## Outputs: A list of selected matrices and the actual number of replicates tried.
stabsel.jsem <- function(X, index, cnt, lastar) {
  K = 1
  p = ncol(X)
  if (is.null(dim(X))) {
    K = length(X)
    p = ncol(X[[1]])
  }  
  n = lapply(X, nrow)
  
  X1 = vector("list",K)
  X2 = vector("list", K)
  sel.mat = vector("list", K)
  for (k in 1:K){
    sel.mat[[k]] = matrix(0, p, p)
  }
  count = 0
  for (i in 1:cnt) {
    model.1 = NULL 
    model.2 = NULL 
    for (k in 1:K){
      ind.1 = sample(seq(1, n[[k]]), n[[k]]/2, F)
      ind.2 = seq(1, n[[k]])[match(seq(1, n[[k]]), ind.1, 0) == 0]
      X1[[k]] = X[[k]][ind.1, ]
      X2[[k]] = X[[k]][ind.2, ]
      model.1 = c(model.1, rep(k, length(ind.1)))
      model.2 = c(model.2, rep(k, length(ind.2)))
    }
    tmp.1 = try(JSEM(trainX=do.call(rbind, X1), trainY=model.1, index, lambda=lastar))
    tmp.2 = try(JSEM(trainX=do.call(rbind, X2), trainY=model.2, index, lambda=lastar))
    
    if (inherits(tmp.1, "try-error") || inherits(tmp.2, "try-error")){
      warning("There might be some error!")
      next;
    }
    
    for (k in 1:K){
      sel.mat[[k]] = sel.mat[[k]] + tmp.1$Ahat[[k]] + tmp.2$Ahat[[k]]
    }
    
    count = count + 1
  }
  
  return(list(mat = sel.mat, count = count))
}

##--------------------------------------------\
##              sort.variables
##--------------------------------------------\
# Input: the node-sepcific group index
# e.g. 
# index = cbind(do.call(cbind, rep(list(c(1,1,2,2)), 10)), do.call(cbind, rep(list(c(1,2,1,2)), 10)))
# dimension of index: K by p
# Output: 
# v.index, index indicating group membership of each variables
# g.index, index for graphs
# x.index, index for columns of the design matrix X
# b.index, index for recovering the beta coefficients to the correct order
##--------------------------------------------\
sort.variables <- function(
  index # the group index matrix of p-1 X K
){
  K = nrow(index)
  p = ncol(index) + 1
  
  len = apply(index, 2, function(x)  length(unique(x)) ) 
  
  g.index = matrix(rep(1:K, p-1), nrow = K, ncol = p-1)
  x.index = order(c(t(do.call(rbind, rep(list(1:ncol(index)), K)))))
  
  #initialize the variable index
  v.index = index
  for (j in 2:ncol(index)) {
    v.index[, j] = v.index[, j] + cumsum(len)[j-1]
  }
  v.index = c(v.index)
  
  # re-order the variable index so that they are monotone
  new.order = order(v.index)
  v.index = v.index[new.order]
  x.index = x.index[new.order]
  g.index = g.index[new.order]
  b.index = index
  for (j in 1:ncol(index)) {
    b.index[, j] = order(order(index[, j]))
  }
  
  res = list(v.index = v.index, g.index = g.index, x.index = x.index, b.index = b.index)
  return(res)
}

sf.net <- function(
  p,  # number of variables
  m = NULL, 
  rho=1 
){
  # generate a graph
  g <- barabasi.game(n = p, power=rho, m = m, directed = F, algorithm = "psumtree")
  adjm <- as.matrix(get.adjacency(g))
  d <- graph.density(g)
  return(list(A = adjm, g = g, density = d))
}

pd <- function(A, zeta=0.1){
  if (sum(A != t(A)) > 0){
    stop("This method only works for symmetric A!")
  }
  
  p <- dim(A)[1]
  diag(A) <- rep(0, p)
  diag(A) <- abs(min(eigen(A)$values)) + zeta
  Ainv <- chol2inv(chol(A))
  Ainv <- cov2cor(Ainv)
  A <- chol2inv(chol(Ainv))
  return(list(A = A, Ainv = Ainv))
}

# Symmetrize a matrix
symmetrize <- function(A, eps = 1e-06){
  A <- (A + t(A))/2
  A[abs(A)<eps] <- 0
  diag(A) <- 0
  return(A)
}

##--------------------------------------------\
## function: zeroInd
##--------------------------------------------\
## purpose: Get the indices for which we have external information
## output: 
##    zeroArr: the 2-column indices for which we will zero out
##    zeroMat: the known 0's as a matrix.
##    oneMat:  the known 1's
##--------------------------------------------\
zeroInd <- function(Amat, r, eps=1e-06){
  if (!isSymmetric(Amat)){
    stop("This method only works for symmetric matrix!")
  }
  p <- dim(Amat)[1]
  oneMat <- matrix(0, p, p)
  zeroMat <- matrix(0, p, p)
  
  one.pos <- which(abs(Amat)>=eps, arr.ind = TRUE)
  zero.pos <- which(abs(Amat)<eps, arr.ind = TRUE)
  
  zero.pos <- zero.pos[which(zero.pos[,1] > zero.pos[,2]) ,]
  sel.zero <- sample(seq(1, dim(zero.pos)[1]), r * dim(zero.pos)[1], replace = FALSE) 
  zeroMat[zero.pos[sel.zero, ]] <- 1
  zeroMat <- zeroMat + t(zeroMat)  
  zeroArr <- zero.pos[sel.zero, ]
  
  out <- list()
  out$zeroArr = zeroArr
  out$zeroMat = zeroMat
  
  if (dim(one.pos)[1] == 0){
    warning("The matrix is zero!")
    out$oneMat = matrix(0, p, p)
  } else 
  {
    one.pos <- one.pos[which(one.pos[,1] > one.pos[,2]) ,]
    if (is.null(dim(one.pos))){
      one.pos = matrix(one.pos, nrow = 1)
    }
    
    sel.one <- sample(seq(1, dim(one.pos)[1]), r * dim(one.pos)[1], replace = FALSE) 
    oneMat[one.pos[sel.one, ]] <- 1
    oneMat <- oneMat + t(oneMat)
    diag(oneMat) <- 0
    
    out$oneMat = oneMat  
  }
  
  return(out)  
}
#To compute the trace of a matrix
matTr <- function(z) sum(diag(z))

##--------------------------------------------\
#  multi.glasso
##--------------------------------------------\
#Purpose: to estimate multiple adjacency matrices using graphical lasso
# Input: 
# trainX,      a data matrix sorted by group label, n by p
# trainY,      a n-dim vector indicating model membership of trainX
# lambda,      a scalar as the penalty parameter in each glasso problem
# zero = NULL,  entries of inverse covariance matrix to be constrained to zero. a list of matrices indicating the constraints for glasso.
#              It's mainly used for our method.
# BIC = FALSE,  whether to calculate the bic.score.
# eps = 1e-06

# Output: estimated adjacency matrices; and precision matrices
##--------------------------------------------\
multi.glasso <- function(
  trainX,    
  trainY,      
  lambda,     
  zero = NULL, 
  BIC = FALSE,
  eps = 1e-06
){
  p = dim(trainX)[2]
  K = length(unique(trainY))
  n = as.numeric(table(trainY))
  
  #penalty needed for glasso
  if (length(lambda)==K) {rho = lambda} else { 
    rho = rep(lambda, K)}
  
  #Initialize the estimated precision, partial correlation and adjacency matrix
  Omega.hat = vector("list", K)
  Theta = vector("list", K)
  Ahat = vector("list", K)
  
  #Whether there are entries that need to be constrained to zero
  if (is.null(zero)){
    zero = rep(list(zero), K)
  }
  
  # if (max(sapply(zero, length)) == p*(p-1)){
  #   stop("One or more matrices are constrained to be zero")
  # }  
  
  bic.score = rep(0, K)
  
  for (k in 1:K) {
    Ahat[[k]] = matrix(0, p, p)
    data <- trainX[which(trainY == k), ]
    empcov <- cov(data) #empirical cov
    
    if(length(zero[[k]])<p*(p-1)){
      
      while (kappa(empcov) > 1e+2){
        empcov = empcov + 0.05 * diag(p)
      }
      
      fit <- glasso(empcov, rho = rho[k], zero = zero[[k]], penalize.diagonal=FALSE, maxit = 30)
      
      Omega.hat[[k]] = (fit$wi + t(fit$wi))/2
      Theta[[k]] <- diag(diag(Omega.hat[[k]])^(-0.5)) %*% Omega.hat[[k]] %*% diag(diag(Omega.hat[[k]])^(-0.5))
      Ahat[[k]][abs(Omega.hat[[k]])>eps] = 1
      diag(Ahat[[k]]) = 0
      
      if(BIC){
        bic.score[k] = matTr(empcov %*% Omega.hat[[k]]) - log(det(Omega.hat[[k]])) + log(n[k]) * sum(Ahat[[k]])/(2*n[k])
      }
    } else{ # in case all edges are zero in k-th adjacency matrix
      Omega.hat[[k]] = diag(1/(diag(empcov)-0.05))
      Theta[[k]] = diag(1, p)
    }
  }
  
  out = list(Omega = Omega.hat, Theta = Theta, Adj = Ahat, BIC = bic.score, lambda = lambda)
  return(out)
}


