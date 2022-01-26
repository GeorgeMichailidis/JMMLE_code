library(glmnet)
library(huge)

l1LS_Main = function(Y,X,skeleton.hat,lambda=0.02,initializer="Lasso"){ 
## choose initializer between "Lasso" and "Ridge"
	
	n = nrow(X)
	p1 = ncol(X)
	p2 = ncol(Y)
	
	Existing.edges = diag(1:p1) %*% skeleton.hat
	
	# distribute the regression to multiple nodes
	if (initializer=="Lasso")
	{
		output.list_LS = foreach(j=1:p2)%dopar%{
		
			B_j = rep(0,p1)
	
			if (length(which(Existing.edges[,j]!=0))==0)
				B_j = rep(0,p1)
			else if (length(which(Existing.edges[,j]!=0))==1)
				B_j[which(Existing.edges[,j]!=0)] = lm(Y[,j]~X[,which(Existing.edges[,j]!=0)]+0)$coef
			else{	
				temp = glmnet(X[,which(Existing.edges[,j]!=0)],Y[,j],intercept=FALSE)
				B_j[which(Existing.edges[,j]!=0)] = predict(temp,s=lambda,type="coefficients")[-1]
			}
			
			B_j
		}
	}
	if (initializer=="Ridge")
	{
		output.list_LS = foreach(j=1:p2)%dopar%{
		
			B_j = rep(0,p1)
			if (length(which(Existing.edges[,j]!=0))==0){
				B_j = rep(0,p1);
			}
			else if (length(which(Existing.edges[,j]!=0))==1){
				B_j[which(Existing.edges[,j]!=0)] = lm(Y[,j]~X[,which(Existing.edges[,j]!=0)]+0)$coef
			}
			else {
				temp = glmnet(X[,which(Existing.edges[,j]!=0)],Y[,j],intercept=FALSE,alpha=0)
				B_j[which(Existing.edges[,j]!=0)] = predict(temp,s=lambda,type="coefficients")[-1]
			}
			B_j
		}
	}
	
	## collect the result
	B.est = array(0,c(p1,p2))
	for (j in 1:p2){
		B.est[,j] = output.list_LS[[j]]
	}
		
	ResMat = Y - X %*% B.est
	Theta0 = as.matrix(huge(ResMat,sqrt(log(p2)/n),method="glasso",verbose=FALSE)$icov[[1]])
	
	return(list(B0=B.est,Theta0=Theta0))	
}
