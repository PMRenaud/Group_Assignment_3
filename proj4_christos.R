## Statistical Programming, Practical 4
## Github repo: https://github.com/
## Christos Giannikos (s2436019), Elliot Leishman (s1808740), 
## Patrick Renaud (s2462989)

## The purpose of this practical is to construct an implementation of Newton's 
## method for mininizing functions. This is an iterative method that constructs
## a sequence of values (theta), starting from an initial guess theta(0) that 
## converges towards a minimizer. It does so by using second-order Taylor
## approximations of the objective function f around the constructed sequence 
## values. The iteration stops when the method determines a value where all 
## elements of the gradient vector have an absolute value which is lower than a 
## certain limit.


## General outline:
## In order to implement Newton's method, we have constructed a function named
## newt, which constructs the desired sequence. We have also created some other
## functions which are called by newt. These are finite_differencing, 
## check_positive_definite and positive_definite respectively.

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                 maxit=100,max.half=20,eps=1e-6){
  ## This function receives the following arguments:
  
  ## theta, which  is a vector of initial values for the optimization parameters.
  
  ## func, which is the objective function to minimize. Its first argument is
  ## the vector of optimization parameters. Remaining arguments will be passed 
  ## from newt using ‘...’.
  
  ## grad, which is the gradient function. It has the same arguments as func but
  ## returns the gradient vector of the objective w.r.t. the elements of the 
  ## parameter vector.
  
  ## hess, which is the Hessian matrix function. If not supplied then newt 
  ## obtains an approximation to the Hessian by calling the finite_differencing 
  ## we have defined below.
  
  ## tol, which is the convergence tolerance
  
  ## fscale, which provides a rough estimate of the magnitude of func near the 
  ## optimum.
  
  ## maxit, the maximum number of Newton iterations to try before stopping the 
  ## method.
  
  ## max.half, which is the maximum number of times a step should be halved
  ## before concluding that the step has failed to improve the objective function.
  
  ## eps, the finite difference intervals that will be used when the Hessian 
  ## function is not provided.
  
  
  ## iteration stops when the method determines a value where all elements of the
  ## gradient vector have an absolute value which is lower than tol times the 
  ## absolute value of the objective function plus fscale
  
  ## f and gradient are the value of objective function and gradient at initial 
  ## theta
  f <- func(theta)
  gradient <- grad(theta)
  # If the Hessian matrix function is entered as an argument, then we just need 
  ## to calculate the value of hessian at theta
  if (!is.null(hess)){
    hessian <- hess(theta)
  }else{ # If the Hessian function is not entered as argument,in order to 
    ## calculate it a certain point, we have to use finite differencing
    hessian <- finite_differencing(f,theta,gradient) 
  }
  # Check that Hessian is positive definite by calling check function
  if (!check_positive_definite(hessian)){
  # If hessian is not positive definite, the function positive definite 
  # makes it so
    hessian <- positive_definite(hessian)
    }
  count <- 0 # Initialise counter of newton method iterations
  
  # Looping until convergence requirements are met
  while(min(abs(grad)) > (tol * (abs(f) + fscale)) && count < maxit){
    f0 <- f
    theta_new <- theta - chol2inv(hessian) %*% gradient
    f <- func(theta) # Update f and gradient with new values
    gradient <- grad(theta)
    
    if (!is.null(hess)){
      hessian <- hess(theta)
    }else{
      hessian<- finite_differencing(theta,gradient) 
    }
    
    ## if the induced hessian matrix is not positive definite, we make it so
    hessian <- positive_definite(hessian)
    
    count <- count + 1 ## increase the count of method iterations by 1
    
    
    
    halving_iterations <- 0
    while(f>=f0 && halving_iterations <= max.half){
      f0 <- f
      theta <- theta - 0.5*chol2inv(hessian) %*% gradient
      f <- func(theta)
      gradient <- grad(theta)
      halving_iterations <- halving_iterations + 1
      count <- count + halving_iterations
  }
  }
}

finite_differencing <- function(theta,grad){
  ## like pg 72 in notes
  dim = length(theta)
  eps <- 1e-7 ## finite difference interval
  Hfd <- matrix(0,dim,dim) ## initializing finite difference Hessian
  for (i in 1:dim) { ## Looping over parameters
    th1 <- theta;
    th1[i] <- th1[i] + eps  ## Increase th1[i] by eps 
    grad1 <- grad(th1,...)  ## compute resulting gradient 
    Hfd[i,] <- (grad1 - grad(theta))/eps ## Approximate second derivatives
  }
  return(Hfd)
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
    while (check_positive_definite(A)==TRUE){
      l=10*l
      A <- A +l*diag(1,dim(A)[1])
    }
    return(A)
  }
}

