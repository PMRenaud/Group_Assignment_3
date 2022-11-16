newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                 maxit=100,max.half=20,eps=1e-6){
  
  f <- func(theta)
  grad <- grad(theta)
  if (!is.null(hess)){
  hessian <- hess(theta)
  }else{ 
    hessian <- finite_differencing(f,theta,grad) 
  }
  if (!check_positive_definite(hessian)){
    hessian <- positive_definite(hessian)
    }
  count <- 0
  while(min(abs(grad)) > (tol * (abs(f) + fscale)) && count < maxit){
    f0 <- f
    theta <- theta - chol2inv(hessian) %*% grad
    f <- func(theta)
    grad <- grad(theta)
    if (!is.null(hess)){
      hessian <- hess(theta)
    }else{
      hessian<- finite_differencing(f,theta,grad) 
    }
    hessian <- positive_definite(hessian)
    while(f>=f0){
      f0 <- f
      theta <- theta - 0.5*chol2inv(hessian) %*% grad
      f <- func(theta)
      grad <- grad(theta)
    count <- count + 1
  }
  }
}

finite_differencing <- function(f,theta,grad){
  ## like pg 72 in notes
  dim = length(theta)
  Hfd <- matrix(0,dim,dim) ## finite difference Hessian
  for (i in 1:dim) { ## Looping over parameters
    th1 <- theta;
    th1[i] <- th1[i] + eps        ## Increase th1[i] by eps 
    grad1 <- grad(th1,...) 
    H[i,] <- (grad1 - ftheta)/eps ## Approximate second derivatives
    hessian <- H
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
positive_definite <- function(A){
  if (check_positive_definite(A)==TRUE){
    return(A)
  }else{
    l < 1e-06
    temp < -A
    A <- A +l*diag(1,dim(A)[1])
    result <- try(chol(A))
    sf <- inherits(result,"matrix")
    while (check_positive_definite(A)==TRUE){
      l=10*l
      A <- A +l*diag(1,dim(A)[1])
    }
    return(A)
  }
}

