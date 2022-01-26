## ****************************************************#
## The Joint Multiple Multi-level Estimation (JMMLE) Method
## ****************************************************#

#****************************************************#
# One-step version: first stabilize Beta, then optimize Theta and Beta once
# USE THIS FOR PRACTICAL PURPOSES.
#****************************************************#
jmmle.1step = function(Y.list, # K-length list of Y-matrices, each n X q
                       Y.indices=NULL,
                       X.list, # K-length list of X-matrices, each n X p
                       B.group.array, # array of group indices for B, size p X q X K
                       Theta.groups, # q-length list of Y-matrices, each K X (q-1)
                       B_init.array=NULL, # initializer for B arrays
                       Theta_init.array=NULL, # initializer for theta arrays
                       init.gamma=NULL,
                       # if you supply Theta_init.array, then also
                       # supply the tuning parameter for the JSEM part
                       init.option=1, # default initializer option
                       lambda=NULL, # B estimation tuning parameter- scalar
                       gamma=NULL, # vector of Theta estimation tuning parameters
                       refit.B=TRUE, # refit B at the end of final iteration?
                       tol=1e-4, # tolerance for convergence
                       maxit=20, # max iterations
                       eps=1e-6, # threshold to set entries of Theta matrices zero
                       VERBOSE=TRUE # print output or not
){
    
    #****************************************************#
    # Define and initialize some quantities
    #****************************************************#
    arraydims = dim(B.group.array)
    p = arraydims[1]; q = arraydims[2]; K = arraydims[3]
    n = nrow(Y.list[[1]])
    Theta.group.array = array(0, c(q,q,K))
    for(j in 1:q){
        Theta.group.array[j,-j,] = Theta.groups[[j]]
    }
    
    ## default values of arguments if they are NULL
    if(is.null(Y.indices)){
        Y.indices = list()
        for(k in 1:K){
            Y.indices[[k]] = rep(k, nrow(Y.list[[k]]))
        }
    }
    if(is.null(lambda)){
        lambda = .5 * sqrt(log(p)/n)
    }
    if(is.null(gamma)){
        gamma = sqrt(log(q)/n) * seq(1, 0.1, -0.1)
    }
    
    #****************************************************#
    # Initialization of iterates
    #****************************************************#
    
    ## Option 0: initial values supplied. Nothing to do
    if(!is.null(B_init.array) & !is.null(Theta_init.array)){
        if(VERBOSE){
            cat("Initial values detected. Skipping initialization\n")
        }
        Ehat.list = list()
        for(k in 1:K){
            Ehat.list[[k]] = Y.list[[k]] - X.list[[k]] %*% B_init.array[,,k]
        }
        gamma.min = init.gamma
        
        init.option = 0
    }
    
    ## Option 1: initialize B, then run JSEM to initialize Theta
    if(init.option==1){
        ## initialize B from separate analysis
        if(VERBOSE){
            cat("Initializing B array\n")
        }
        skeleton.hat = array(1, c(p,q)); # assuming all directed edges possibly exist
        Existing.edges = diag(1:p) %*% skeleton.hat
        B_init.array = array(0, dim=c(p,q,K))
        Ehat.list = list()
        
        for(k in 1:K){
            LS = l1LS_Main(Y.list[[k]], X.list[[k]], skeleton.hat=skeleton.hat,
                           lambda=lambda,initializer="Lasso")
            B_init.array[,,k] = LS$B0
            Ehat.list[[k]] = Y.list[[k]] - X.list[[k]] %*% LS$B0
        }
        
        ## initialize Theta
        if(VERBOSE){
            cat("Initializing Theta array: ")
        }
        init.bic.jsem <- sel.lambda.jsem(do.call(rbind, Ehat.list), do.call(rbind, Ehat.list),
                                         unlist(Y.indices), unlist(Y.indices),
                                         Theta.groups, lambda=gamma)
        gamma.min = gamma[which.min(init.bic.jsem$BIC)]
        init.jsem.model = JSEM(do.call(rbind, Ehat.list), unlist(Y.indices),
                               Theta.groups, lambda=gamma.min)
        Theta_init.array = array(0, c(q,q,K))
        for(k in 1:K){
            Theta_init.array[,,k] = init.jsem.model$Theta[[k]]
        }
    }
    
    ## Option 2: initialize Theta only by running JSEM once
    if(init.option==2){
        ## initialize Theta
        if(VERBOSE){
            cat("Initializing Theta array: ")
        }
        init.bic.jsem <- sel.lambda.jsem(do.call(rbind, Y.list), do.call(rbind, Y.list),
                                         unlist(Y.indices), unlist(Y.indices),
                                         Theta.groups, lambda=gamma)
        gamma.min = gamma[which.min(init.bic.jsem$BIC)]
        init.jsem.model = JSEM(do.call(rbind, Y.list), unlist(Y.indices),
                               Theta.groups, lambda=gamma.min)
        Theta_init.array = array(0, c(q,q,K))
        for(k in 1:K){
            Theta_init.array[,,k] = init.jsem.model$Theta[[k]]
        }
        B_init.array = array(0, dim=c(p,q,K))
        Ehat.list = Y.list
    }
    
    #****************************************************#
    # Alternating algorithm
    #****************************************************#
    
    ## make long X and Y matrices
    Y = do.call(rbind, Y.list)
    X = as.matrix(do.call(bdiag, X.list))
    
    # initialize
    Normdiff = c()
    Objfunc = Obj(Y.list, X.list, Theta_init.array, B_init.array,
                  Theta.group.array, B.group.array,
                  lambda=lambda, gamma=.5 * sqrt(log(q)/n))
    iter = 0; CONVERGE=FALSE; CONVERGE1=FALSE; refit.B=TRUE;
    B_new.array = B_init.array
    Theta_new.array = Theta_init.array
    jsem.model = NULL
    
    # # store bic values from JSEM models
    # bic.mat = init.bic.jsem$BIC
    
    ## start with the alternating procedure
    if(VERBOSE){
        cat('-----\n')
    }
    while(!CONVERGE){
        iter = iter + 1;
        B_old.array = B_new.array
        Theta_old.array = Theta_new.array
        
        if(VERBOSE){
            cat("Iteration ", iter, ":\n")
        }
        # Updating B
        if(VERBOSE){
            cat("Updating B array\n")
        }
        for(j in 1:q){
            
            # make long vector or errors for j-th column for all k
            Et.j.list = list()
            for(k in 1:K){
                Et.j.list[[k]] = Ehat.list[[k]][,-j] %*% Theta_old.array[-j,j,k]
            }
            Ehat.theta.j = unlist(Et.j.list)
            rm(Et.j.list)
            
            # build model
            temp = grpreg(X, Y[,j] + Ehat.theta.j, unlist(as.numeric(B.group.array[,j,])),
                          family="gaussian", penalty="grLasso", lambda=lambda)
            B_new.array[,j,] = matrix(temp$beta[-1], ncol=K, byrow=F)
            
            ## refit if necessary
            if(refit.B){
                # make long vector or errors for j-th column for all k
                temp1 = temp$beta[-1]
                B.j.support = which(abs(temp1)>1e-6)
                if (length(B.j.support)>0){
                    # build model
                    temp2 = lm(Y[,j] + Ehat.theta.j~X[,B.j.support]+0)
                    
                    # save non-zero coefs
                    temp1[B.j.support] = temp2$coef
                    B_new.array[,j,] = matrix(temp1, ncol=K, byrow=F)
                }
            }
            
            # now update j-th column of all K matrices in B_new, and Ehat
            for(k in 1:K){
                Ehat.list[[k]][,j] = Y.list[[k]][,j] - X.list[[k]] %*% as.matrix(B_new.array[,j,k], ncol=1)
            }
        }
        
        # If Beta is stablized, updating Theta then B once then break
        if (iter >=10 | sqrt(sum(B_new.array - B_old.array)^2)<0.1){
            # Update Theta
            if(VERBOSE){
                cat("Updating Theta array: ")
            }
            bic.jsem <- sel.lambda.jsem(do.call(rbind, Ehat.list), do.call(rbind, Ehat.list),
                                        unlist(Y.indices), unlist(Y.indices),
                                        Theta.groups, lambda=gamma)
            gamma.min = gamma[which.min(bic.jsem$BIC)]
            jsem.model = JSEM(do.call(rbind, Ehat.list), unlist(Y.indices),
                              Theta.groups, lambda=gamma.min)
            Theta_new.array = array(0, c(q,q,K))
            for(k in 1:K){
                Theta_new.array[,,k] = jsem.model$Theta[[k]]
            }
            
            # Update B
            B_old.array = B_new.array
            if(VERBOSE){
                cat("Updating B array\n")
            }
            for(j in 1:q){
                
                # make long vector or errors for j-th column for all k
                Et.j.list = list()
                for(k in 1:K){
                    Et.j.list[[k]] = Ehat.list[[k]][,-j] %*% Theta_new.array[-j,j,k]
                }
                Ehat.theta.j = unlist(Et.j.list)
                rm(Et.j.list)
                
                # build model
                temp = grpreg(X, Y[,j] + Ehat.theta.j, unlist(as.numeric(B.group.array[,j,])),
                              family="gaussian", penalty="grLasso", lambda=lambda)
                B_new.array[,j,] = matrix(temp$beta[-1], ncol=K, byrow=F)
                
                ## refit if necessary
                if(refit.B){
                    # make long vector or errors for j-th column for all k
                    temp1 = temp$beta[-1]
                    B.j.support = which(abs(temp1)>1e-6)
                    if (length(B.j.support)>0){
                        # build model
                        temp2 = lm(Y[,j] + Ehat.theta.j~X[,B.j.support]+0)
                        
                        # save non-zero coefs
                        temp1[B.j.support] = temp2$coef
                        B_new.array[,j,] = matrix(temp1, ncol=K, byrow=F)
                    }
                }
                
                # now update j-th column of all K matrices in B_new, and Ehat
                for(k in 1:K){
                    Ehat.list[[k]][,j] = Y.list[[k]][,j] - X.list[[k]] %*% as.matrix(B_new.array[,j,k], ncol=1)
                }
                CONVERGE1 = TRUE # Set CONVERGE = 1 to break outer loop
            }
        } else{
            Theta_new.array = Theta_old.array
        }
        
        # check convergence
        Objfunc[iter+1] = Obj(Y.list, X.list, Theta_new.array, B_new.array,
                              Theta.group.array, B.group.array,
                              lambda=lambda, gamma=gamma.min)
        Normdiff[iter] = sqrt(sum(B_new.array - B_old.array)^2)/sqrt(sum(B_new.array^2)) +
            sqrt(sum(Theta_new.array - Theta_old.array)^2)/sqrt(sum(Theta_new.array^2))
        # if (iter == 1){
        #   Norm_diff = Normfunc[1]
        # }
        # else{
        #   Norm_diff = Normfunc[iter] - Normfunc[iter-1]
        #   Obj_diff = (Objfunc[iter] - Objfunc[iter-1])/Objfunc[iter-1]
        # }
        Obj_diff = Objfunc[iter+1]/Objfunc[iter] - 1
        
        # convergence criterion value
        nd = abs(Normdiff[iter])
        if(iter>1){
            nd = c(nd, abs(rev(diff(Normdiff))[1]))
        }
        if(iter>2){
            nd = c(nd, abs(rev(diff(Normdiff,2))[1]))
        }
        if(VERBOSE){
            cat("Norm_diff =",round(nd,4),'Obj_diff',round(abs(Obj_diff),5),'\n-----\n')
        }
        CONVERGE = (CONVERGE1 | (min(nd)<tol))
        
        if (iter == maxit){
            if(VERBOSE){
                cat("Max iterations reached.",'\n')
            }
            break;
        }
    }
    
    if(CONVERGE){
        if(VERBOSE){
            cat("Converged after",iter,"iterations.\n")
        }
    } else{
        warning("algorithm didn't converge.\nRequired epsilon for convergence is ", tol,
                " while current value is ", round(min(nd), 4))
    }
    
    ## If refitting of B matrices hasn't been done inside the loop then refit in the end
    B_refit.array = B_new.array
    if(!refit.B){
        if(VERBOSE){
            cat("Refitting B\n")
        }
        
        for(j in 1:q){
            
            # make long vector or errors for j-th column for all k
            B.j.support = which(abs(B_new.array[,j,])>1e-6)
            if (length(B.j.support)>0){
                Et.j.list = list()
                for(k in 1:K){
                    Et.j.list[[k]] = Ehat.list[[k]][,-j] %*% Theta_old.array[-j,j,k]
                }
                Ehat.theta.j = unlist(Et.j.list)
                rm(Et.j.list)
                
                # build model
                temp = lm(Y[,j] + Ehat.theta.j~X[,B.j.support]+0)
                
                # save non-zero coefs
                temp1 = B_refit.array[,j,]
                temp1[B.j.support] = temp$coef
                B_refit.array[,j,] = temp1
            }
        }
    }
    
    ## Refit to get Omega
    if(is.null(jsem.model)){
        Ahat = list()
        for(k in 1:K){
            Ahat[[k]] = matrix(0, q, q)
            Ahat[[k]][which(abs(Theta_new.array[,,k])>eps, arr.ind=T)] = 1
            diag(Ahat[[k]]) = 0
        }
    } else{
        Ahat = jsem.model$Ahat
    }
    
    Info = list()
    for (k in 1:K){
        Info[[k]] = zeroInd(Ahat[[k]], 1)$zeroArr
    }
    Theta_refit = multi.glasso(do.call(rbind, Ehat.list), unlist(Y.indices), gamma.min, Info)
    
    ## return
    return(list(B.refit=B_refit.array, # a p X q X K array of B estimates
                Theta_refit=Theta_refit
                # a list of outputs that consists of
                # 1. Omega.hat: K-length list of Omega (neightborhood coef) matrices
                # 2. Theta: K-length list of Theta (precision) matrices
                # 3. Ahat: K-length list of adjacency matrices
                # 4. bic.score: K-length vector of bic values for each precision matrix fit
                # 5. lambda: tuning parameter for glasso fits
    ))
    
}

#****************************************************#
# Full version: first stabilize Beta, then optimize Theta and Beta until convergence
# *DO NOT* USE THIS FOR PRACTICAL PURPOSES. VERY SLOW.
#****************************************************#
jmmle = function(Y.list, Y.indices=NULL, X.list,
                 B.group.array, Theta.groups,
                 B_init.array=NULL, Theta_init.array=NULL, init.gamma=NULL, init.option=1,
                 lambda=NULL, gamma=NULL,
                 refit.B=TRUE, tol=1e-4, maxit=20, eps=1e-6, verbose=TRUE){
    
    #****************************************************#
    # Define and initialize some quantities
    #****************************************************#
    arraydims = dim(B.group.array)
    p = arraydims[1]; q = arraydims[2]; K = arraydims[3]
    n = nrow(Y.list[[1]])
    Theta.group.array = array(0, c(q,q,K))
    for(j in 1:q){
        Theta.group.array[j,-j,] = Theta.groups[[j]]
    }
    
    ## default values of arguments if they are NULL
    if(is.null(Y.indices)){
        Y.indices = list()
        for(k in 1:K){
            Y.indices[[k]] = rep(k, nrow(Y.list[[k]]))
        }
    }
    if(is.null(lambda)){
        lambda = .5 * sqrt(log(p)/n)
    }
    if(is.null(gamma)){
        gamma = sqrt(log(q)/n) * seq(1, 0.1, -0.1)
    }
    
    #****************************************************#
    # Initialization of iterates
    #****************************************************#
    
    ## Option 0: initial values supplied. Nothing to do
    if(!is.null(B_init.array) & !is.null(Theta_init.array)){
        if(VERBOSE){
            cat("Initial values detected. Skipping initialization\n")
        }
        Ehat.list = list()
        for(k in 1:K){
            Ehat.list[[k]] = Y.list[[k]] - X.list[[k]] %*% B_init.array[,,k]
        }
        gamma.min = init.gamma
        init.option = 0
    }
    
    ## Option 1: initialize B, then run JSEM to initialize Theta
    if(init.option==1){
        ## initialize B from separate analysis
        cat("Initializing B array\n")
        skeleton.hat = array(1, c(p,q)); # assuming all directed edges possibly exist
        Existing.edges = diag(1:p) %*% skeleton.hat
        B_init.array = array(0, dim=c(p,q,K))
        Ehat.list = list()
        
        for(k in 1:K){
            LS = l1LS_Main(Y.list[[k]], X.list[[k]], skeleton.hat=skeleton.hat,
                           lambda=lambda,initializer="Lasso")
            B_init.array[,,k] = LS$B0
            Ehat.list[[k]] = Y.list[[k]] - X.list[[k]] %*% LS$B0
        }
        
        ## initialize Theta
        cat("Initializing Theta array: ")
        init.bic.jsem <- sel.lambda.jsem(do.call(rbind, Ehat.list), do.call(rbind, Ehat.list),
                                         unlist(Y.indices), unlist(Y.indices),
                                         Theta.groups, lambda=gamma)
        gamma.min = gamma[which.min(init.bic.jsem$BIC)]
        init.jsem.model = JSEM(do.call(rbind, Ehat.list), unlist(Y.indices),
                               Theta.groups, lambda=gamma.min)
        Theta_init.array = array(0, c(q,q,K))
        for(k in 1:K){
            Theta_init.array[,,k] = init.jsem.model$Theta[[k]]
        }
    }
    
    ## Option 2: initialize Theta only by running JSEM once
    if(init.option==2){
        ## initialize Theta
        cat("Initializing Theta array: ")
        init.bic.jsem <- sel.lambda.jsem(do.call(rbind, Y.list), do.call(rbind, Y.list),
                                         unlist(Y.indices), unlist(Y.indices),
                                         Theta.groups, lambda=gamma)
        gamma.min = gamma[which.min(init.bic.jsem$BIC)]
        init.jsem.model = JSEM(do.call(rbind, Y.list), unlist(Y.indices),
                               Theta.groups, lambda=gamma.min)
        Theta_init.array = array(0, c(q,q,K))
        for(k in 1:K){
            Theta_init.array[,,k] = init.jsem.model$Theta[[k]]
        }
        B_init.array = array(0, dim=c(p,q,K))
        Ehat.list = Y.list
    }
    
    #****************************************************#
    # Alternating algorithm
    #****************************************************#
    
    ## make long X and Y matrices
    Y = do.call(rbind, Y.list)
    X = as.matrix(do.call(bdiag, X.list))
    
    # initialize
    Normdiff = c()
    Objfunc = Obj(Y.list, X.list, Theta_init.array, B_init.array,
                  Theta.group.array, B.group.array,
                  lambda=lambda, gamma=.5 * sqrt(log(q)/n))
    iter = 0; CONVERGE=FALSE; refit.B=TRUE; update.counter=0;
    updateTheta = FALSE; # we don't update Theta until B is stabilized a bit
    B_new.array = B_init.array
    Theta_new.array = Theta_init.array
    jsem.model = NULL
    
    # # store bic values from JSEM models
    # bic.mat = init.bic.jsem$BIC
    
    ## start with the alternating procedure
    cat('-----\n')
    while(!CONVERGE){
        iter = iter + 1;
        B_old.array = B_new.array
        Theta_old.array = Theta_new.array
        
        cat("Iteration ", iter, ":\n")
        # Updating B
        cat("Updating B array\n")
        for(j in 1:q){
            
            # make long vector or errors for j-th column for all k
            Et.j.list = list()
            for(k in 1:K){
                Et.j.list[[k]] = Ehat.list[[k]][,-j] %*% Theta_old.array[-j,j,k]
            }
            Ehat.theta.j = unlist(Et.j.list)
            rm(Et.j.list)
            
            # build model
            temp = grpreg(X, Y[,j] + Ehat.theta.j, unlist(as.numeric(B.group.array[,j,])),
                          family="gaussian", penalty="grLasso", lambda=lambda)
            B_new.array[,j,] = matrix(temp$beta[-1], ncol=K, byrow=F)
            
            ## refit if necessary
            if(refit.B){
                # make long vector or errors for j-th column for all k
                temp1 = temp$beta[-1]
                B.j.support = which(abs(temp1)>1e-6)
                if (length(B.j.support)>0){
                    # build model
                    temp2 = lm(Y[,j] + Ehat.theta.j~X[,B.j.support]+0)
                    
                    # save non-zero coefs
                    temp1[B.j.support] = temp2$coef
                    B_new.array[,j,] = matrix(temp1, ncol=K, byrow=F)
                }
            }
            
            # now update j-th column of all K matrices in B_new, and Ehat
            for(k in 1:K){
                Ehat.list[[k]][,j] = Y.list[[k]][,j] - X.list[[k]] %*% as.matrix(B_new.array[,j,k], ncol=1)
            }
        }
        
        # # update array of E
        # for(k in 1:K){
        #   Ehat.list[[k]] = Y.list[[k]] - X.list[[k]] %*% B_new.array[,,k]
        # }
        
        # Updating Theta
        if (iter >=10 | sqrt(sum(B_new.array - B_old.array)^2)<0.1 | updateTheta == TRUE){
            updateTheta = TRUE; # once we start updating Theta, we just start from now on;
            cat("Updating Theta array: ")
            bic.jsem <- sel.lambda.jsem(do.call(rbind, Ehat.list), do.call(rbind, Ehat.list),
                                        unlist(Y.indices), unlist(Y.indices),
                                        Theta.groups, lambda=gamma)
            gamma.min = gamma[which.min(bic.jsem$BIC)]
            jsem.model = JSEM(do.call(rbind, Ehat.list), unlist(Y.indices),
                              Theta.groups, lambda=gamma.min)
            Theta_new.array = array(0, c(q,q,K))
            for(k in 1:K){
                Theta_new.array[,,k] = jsem.model$Theta[[k]]
            }
            update.counter = update.counter + 1
            # bic.mat = rbind(bic.mat, bic.jsem$BIC)
        } else{
            Theta_new.array = Theta_old.array
        }
        
        # check convergence
        Objfunc[iter+1] = Obj(Y.list, X.list, Theta_new.array, B_new.array,
                              Theta.group.array, B.group.array,
                              lambda=lambda, gamma=gamma.min)
        Normdiff[iter] = sqrt(sum(B_new.array - B_old.array)^2)/sqrt(sum(B_new.array^2)) +
            sqrt(sum(Theta_new.array - Theta_old.array)^2)/sqrt(sum(Theta_new.array^2))
        # if (iter == 1){
        #   Norm_diff = Normfunc[1]
        # }
        # else{
        #   Norm_diff = Normfunc[iter] - Normfunc[iter-1]
        #   Obj_diff = (Objfunc[iter] - Objfunc[iter-1])/Objfunc[iter-1]
        # }
        Obj_diff = Objfunc[iter+1]/Objfunc[iter] - 1
        
        # convergence criterion value
        nd = abs(Normdiff[iter])
        if(iter>1){
            nd = c(nd, abs(rev(diff(Normdiff))[1]))
        }
        if(iter>2){
            nd = c(nd, abs(rev(diff(Normdiff,2))[1]))
        }
        cat("Norm_diff =",round(nd,4),'Obj_diff',round(abs(Obj_diff),5),'\n-----\n')
        CONVERGE = (min(nd)<tol)
        
        if (iter == maxit){
            cat("Max iterations reached.",'\n')
            break;
        }
    }
    
    if(CONVERGE){
        cat("Converged after",iter,"iterations.\n")
    } else{
        warning("algorithm didn't converge.\nRequired epsilon for convergence is ", tol,
                " while current value is ", round(min(nd), 4))
    }
    
    ## If refitting of B matrices hasn't been done inside the loop then refit in the end
    B_refit.array = B_new.array
    if(!refit.B){
        cat("Refitting B\n")
        for(j in 1:q){
            
            # make long vector or errors for j-th column for all k
            B.j.support = which(abs(B_new.array[,j,])>1e-6)
            if (length(B.j.support)>0){
                Et.j.list = list()
                for(k in 1:K){
                    Et.j.list[[k]] = Ehat.list[[k]][,-j] %*% Theta_old.array[-j,j,k]
                }
                Ehat.theta.j = unlist(Et.j.list)
                rm(Et.j.list)
                
                # build model
                temp = lm(Y[,j] + Ehat.theta.j~X[,B.j.support]+0)
                
                # save non-zero coefs
                temp1 = B_refit.array[,j,]
                temp1[B.j.support] = temp$coef
                B_refit.array[,j,] = temp1
            }
        }
    }
    
    ## Refit to get Omega
    # cat("Getting Omega 1\n")
    if(is.null(jsem.model)){
        Ahat = list()
        for(k in 1:K){
            Ahat[[k]] = matrix(0, q, q)
            Ahat[[k]][which(abs(Theta_new.array[,,k])>eps, arr.ind=T)] = 1
            diag(Ahat[[k]]) = 0
        }
    } else{
        Ahat = jsem.model$Ahat
    }
    
    # cat("Getting Omega 2\n")
    Info = list()
    for (k in 1:K){
        Info[[k]] = zeroInd(Ahat[[k]], 1)$zeroArr
    }
    # cat("Getting Omega 3\n")
    Theta_refit = multi.glasso(do.call(rbind, Ehat.list), unlist(Y.indices), gamma.min, Info)
    
    ## return
    return(list(B.refit=B_refit.array, Theta_refit=Theta_refit))
}

#****************************************************#
# Auxiliary functions
#****************************************************#
# cat = function(string, VERBOSE){
#   if(VERBOSE){
#     cat(string)
#   }
# }


#****************************************************#
# EOF
#****************************************************#