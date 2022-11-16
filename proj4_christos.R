## Statistical Programming, Practical 4
## Github repo: https://github.com/
## Christos Giannikos (s2436019), Elliot Leishman (s1808740), 
## Patrick Renaud (s2462989)

## The purpose of this practical is to construct an implementation of Newton's 
## method for mininizing functions. This is an iterative method that constructs
## a sequence of values (theta), starting from an initial guess theta(0) that 
## converges towards a minimizer  by using second-order Taylor approximations 
## of the objective function f around the constructed sequence values. The 
## iteration stops when the method determines a value where all elements of the
## gradient vector have an absolute value which is lower than a certain limit.

## In order to implement Newton's method, we have constructed a function named
## newt, which constructs the desired sequence. We have also created some other
## functions which are called by newt. These are finite_differencing, 
## check_positive_definite and positive_definite respectively.

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                 maxit=100,max.half=20,eps=1e-6){
  ## 
  
  ## iteration stops when the method determines a value where all elements of the
  ## gradient vector have an absolute value which is lower than tol times the 
  ## absolute value of the objective function plus fscale
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
    count <- count + 1
    halving_iterations <- 0
    while(f>=f0 && halving_iterations <= max.half){
      f0 <- f
      theta <- theta - 0.5*chol2inv(hessian) %*% grad
      f <- func(theta)
      grad <- grad(theta)
      halving_iterations <- halving_iterations + 1
      count <- count + halving_iterations
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
  ## This function checks whether the input matrix A, which is given as an 
  ## argument, is positive definite. When the input matrix corresponds to a 
  ## Hessian matrix of a function f at some point theta, then this function 
  ## actually returns whether or not the approximated quadratic of f estimated 
  ## at theta, as used in the newt function, has a proper minimum.
  ## Recall that positive definiteness is easily tested by seeing if a Cholesky 
  ## decomposition is possible.
  result <- try(chol(A))
  sf <- inherits(result,"matrix")
  if (sf==TRUE){
    return(TRUE)
  }else{
    return(FALSE)
  }
}
positive_definite <- function(A){
  ## This function accepts as an argument a square matrix A and, if A is not 
  ## a positive definite matrix, it perturbes it to be so. 
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

