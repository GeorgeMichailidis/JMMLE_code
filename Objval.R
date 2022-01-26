library(gdata)
squaredError = function(Y.list, X.list, Theta.array, B.array){
  # this function calculates the total squared error for all list entries
  
  err.list = list()
  K = dim(B.array)[3]
  for(k in 1:K){
    nk = nrow(X.list[[k]])
    qk = ncol(Y.list[[k]])
    Ek = Y.list[[k]] - X.list[[k]] %*% B.array[,,k]
    Tk = diag(1,qk) - Theta.array[,,k]
    err.list[[k]] = sum(diag(crossprod(Ek %*% Tk)))/ nk
  }
  do.call(sum, err.list)
}

Obj = function(Y.list, X.list, Theta.array, B.array,
               Theta.group.array, B.group.array, lambda, gamma){
  # this function calculates the objective function
  
  # calculate penalties
  # For each group index in the group array, collect corresponding elements in the main array ...
  # and sum their l2 norms
  unique.Theta.groups = unique(as.numeric(Theta.group.array))
  Theta.norm = 0
  for(g in unique.Theta.groups){
    Theta.norm = Theta.norm + sqrt(sum(Theta.array[which(Theta.group.array==g, arr.ind=T)]^2))
  }
  
  unique.B.groups = unique(as.numeric(B.group.array))
  B.norm = 0
  for(h in unique.B.groups){
    B.norm = B.norm + sqrt(sum(B.array[which(B.group.array==h, arr.ind=T)]^2))
  }
  
  squaredError(Y.list, X.list, Theta.array, B.array) + lambda*Theta.norm + gamma*B.norm
}

BICfunc = function(Y,X,Theta,B)
{	#BIC = -log(det(Theta)) + tr(S%*%Theta) + sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n + (log(n)/n)*(nonzeros in upperTriangle(Theta) + nonzeros in B)
	#And we prefer smaller BIC value
	n = nrow(X)
	PENs = log(n)/n*(sum(abs(upperTriangle(Theta))>1e-6) + sum(abs(B)>1e-6))
	FIT = sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n - log(det(Theta)) + sum(diag(t(Y-X%*%B) %*% (Y-X%*%B) %*% Theta))/n 
	BIC = FIT + PENs
	return(BIC)
}



		
		
				

