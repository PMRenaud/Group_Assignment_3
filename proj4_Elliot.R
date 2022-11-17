ra <- function(th){
  2*th[1]^2 + th[2]^2-4 + 6*th[2] - 7*th[1] 
}

ga <- function(th){
  c(4*th[1]-7, 2*th[2]+6)
}

ha <- function(th){
  h <- matrix(0,2,2)
  h[1,1] <- 4
  h[2,2] <- 2
  h[1,2] <- 0
  h[2,1] <- 0 
  h
}


rc <- function(th){
  2*th[1]^2 + th[2]^2 
}

gc <- function(th){
  c(4*th[1], 2*th[2])
}

hc <- function(th){
  h <- matrix(0,2,2)
  h[1,1] <- 4
  h[2,2] <- 2
  h[1,2] <- 0
  h[2,1] <- 0 
  h
}
rb <- function(th,k=2){
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
} 

gb <- function(th,k=2){
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

hb <- function(th,k=2){
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}




newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                 maxit=100,max.half=20,eps=1e-6){
  
  # f and f.prime are the value of objective function and gradient at initial theta
  ftheta <- func(theta,...)
  f.prime <- grad(theta,...)
  
  # If Hessian is entered as an argument hessian is value of hessian at theta
  if (!is.null(hess)){ 
    hessian <- hess(theta,...)
  }else{ # Hessian is not entered as argument, use finite differencing.
    hessian <- finite_differencing(ftheta,theta,grad,...) 
  }
  # Check that Hessian is positive definite by calling check function
  if (!check_positive_definite(hessian)){
    # If hessian is not positive definite, the function positive definite makes it so
    hessian <- positive_definite(hessian)
  }
  
  # Initialise counter
  count <- 0
  
  ## CAN I TURN THE TWO WHILE LOOPS INTO TH SAME THING BASICALLY???
  
  # Looping until convergence
  while(max(abs(f.prime)) > (tol * (abs(ftheta) + fscale)) && count < maxit){
    f0 <- ftheta
    
    ## Check chol2 inv is working? is it chol2inv(chol())
    theta_new <- theta - chol2inv(chol(hessian)) %*% f.prime
    
    # Update ftheta and f.prime with new values.
    ftheta <- func(theta_new)
    f.prime <- grad(theta_new)
    
    # Can this part be turned into its own function to avoid repetition of code.
    if (!is.null(hess)){
      hessian <- hess(theta_new)
    }else{
      hessian<- finite_differencing(f,theta_new,grad) 
    }
    half_count <- 0 
    # Don't understand this part. Should do the halving part I think.
    while(ftheta >= f0 && half_count < max.half){
      gamma <- gamma / 2
      theta_new <- theta - gamma*chol2inv(chol(hessian)) %*% f.prime
      ftheta <- func(theta_new)
      half_count <- half_count + 1
      if(half_count == max.half){
        cat('Max half')
      }
    }
    theta <- theta_new
    count <- count + 1 # Append counter
  }
  
  ## Output the list of... 
  output <- list(f = func(theta), theta = theta, iter = count, g = f.prime, Hi = chol2inv(chol(hessian)))
  return(output)
}

finite_differencing <- function(f,theta,grad){
  ## like pg 72 in notes
  dim <- length(theta)
  Hfd <- matrix(0,dim,dim) ## finite difference Hessian
  for (i in 1:dim) { ## Looping over parameters
    th1 <- theta;
    th1[i] <- th1[i] + eps        ## Increase th1[i] by eps 
    grad1 <- grad(th1,...) 
    H[i,] <- (grad1 - ftheta)/eps ## Approximate second derivatives, What is ftheta?
    hessian <- H # Is this symmetric, or do we need to do (t(H) * H) / 2?
  }
}

check_positive_definite <- function(A){
  result <- try(chol(A))
  sf <- inherits(result,"matrix")
  if (sf==TRUE){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
positive_definite <- function(A){ # This function can be tidied up a bit?
  if (check_positive_definite(A)==TRUE){
    return(A)
  }else{
    l < 1e-06
    temp < -A
    A <- A +l*diag(1,dim(A)[1])
    result <- try(chol(A)) # Are these two lines needed? - surely they are the same as check_postive_definite?
    sf <- inherits(result,"matrix")
    while (check_positive_definite(A)==TRUE){
      l=10*l
      A <- A +l*diag(1,dim(A)[1])
    }
    return(A)
  }
}
theta = c(2,-2)
newt(theta, rc, gc, hc)

