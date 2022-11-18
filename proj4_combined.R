## Statistical Programming, Practical 4
## Github repo: https://github.com/
## Christos Giannikos (s2436019), Elliot Leishman (s1808740), 
## Patrick Renaud (s2462989)

## The purpose of this practical is to construct an implementation of Newton's 
## method for mininizing functions. This is an iterative method that constructs
## a sequence of values (theta), starting from an initial guess theta(0) and 
## converges towards a minimizer. It does so by using second-order Taylor
## approximations of the objective function f around the constructed sequence 
## values. The iteration stops when the method determines a value where all 
## elements of the gradient vector have an absolute value which is lower than a 
## certain limit.


## General outline:
## In order to implement Newton's method, we have constructed a function named
## newt, which constructs the desired sequence. We have also created some other
## functions which are called by newt. These are finite_differencing, 
## pos_def and hessian respectively.



finite_differencing <- function(theta,grad,eps,...){
  ## This function is called when the hessian function is not provided as an 
  ## argument. It uses finite difference approximation in order to calculate 
  ## the corresponding hessian.
  ## It receives as arguments a value for theta, the gradient function 
  ## corresponding to the function that we want to optimize and eps, which give
  ## the finite difference intervals we are going to take.
  
  dim = length(theta)
  Hfd <- matrix(0,dim,dim) ## initializing finite difference Hessian
  for (i in 1:dim) {       ## Looping over parameters
    th1 <- theta;
    th1[i] <- th1[i] + eps  ## Increase th1[i] by eps 
    grad1 <- grad(th1,...)  ## compute resulting gradient 
    Hfd[i,] <- (grad1 - grad(theta,...))/eps ## Approximate second derivatives
  }
  ## Ensure that the induced Hessian is symmetric by taking (t(A)+A)/2 
  Hfd <- (Hfd+t(Hfd))/2
  return(Hfd)
}


pos_def <- function(A,eps){
  ## This function accepts as an argument a square matrix A and, if A is not 
  ## a positive definite matrix, it perturbes it to be so.
  ## This is accomplished by adding a multiple of the identity matrix to it, 
  ## large enough to force positive definiteness. 
  ## Note that if the input matrix A is positive definite, then pos_def simply
  ## returns it. 
  l <- eps
  while(inherits(try(chol(A), silent = TRUE),"matrix") == FALSE){
    ## Enter the loop while A is not positive definite, i.e. as long as Cholesky
    ## decomposition is not possible.
    ## A <- A + l * norm(A)*diag(1, dim(A)[1])
    A <- A + l * diag(1, dim(A)[1])
    l <- 10 * l
  }
  return(A)
}

hessian <- function(theta, grad, eps, hess = NULL, ...){
  ## This function computes the hessian matrix for a function whose gradient is 
  ## given by grad, evaluated at a point theta.
  ## It receives as inputs the value theta that the hessian should be evaluated
  ## at, the gradient grad of the function whose hessian matrix we want to 
  ## calculate and the finite difference intervals eps.
  
  if(is.null(hess)){
    ## If in our initial newt function no hessian is provided as input, then we
    ## perform finite differencing
    f2prime <- finite_differencing(theta, grad, eps, ...)
  }
  
  else{
    ## Since the hessian function is provided, we only need to evaluate it at 
    ## theta.
    f2prime <- hess(theta)
  }
  ## guarantee that the induced hessian is positive definite
  f2prime <- pos_def(f2prime) 
  return(f2prime)
}
newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                 maxit=100,max.half=20,eps=1e-6){
  
  ## This function receives the following arguments:
  
  ## theta, which  is a vector of initial values for the optimization parameters,
  ## func, which is the objective function we aim to minimize, 
  ## grad, which is the gradient function for func,
  ## hess, which is the Hessian matrix function for func. If not supplied then 
  ## newt obtains an approximation to the Hessian by calling the 
  ## finite_differencing function we have defined below,
  ## tol, which is the convergence tolerance,
  ## fscale, which provides a rough estimate of the magnitude of func near the 
  ## optimum,
  ## maxit, which corresponds to the maximum number of Newton iterations to try 
  ## before stopping the method,
  ## max.half, which is the maximum number of times a step should be halved
  ## before concluding that the step has failed to improve the objective 
  ## function and
  ## eps, which give the finite difference intervals that will be used when the 
  ## Hessian function is not provided.
  
  ## After implementing newton's optimization method for the input arguments it
  ## returns the elements:
  ## f, corresponding to the value of the objective function at the minimum.
  ## theta, which is the value of the parameters at the minimum.
  ## iter, denoting the number of iterations taken to reach the minimum.
  ## g, which is the gradient vector at the minimum (which 
  ## should be close to numerical zero).
  ## Hi,i.e. the inverse of the Hessian matrix at the minimum.
  
  # Calculate the value of the function func at theta and its gradient, denoted 
  # fprime at the same point 
  ftheta <- func(theta)
  fprime <- grad(theta)
  
  
  
  if(all(is.finite(ftheta))==FALSE){
    stop('ftheta')
  }
  if(all(is.finite(fprime))==FALSE){
    stop('fprime')
  }
  
  
  counter <-0              # Initialise count
  
  while(max(abs(fprime)) > tol*(abs(ftheta)+fscale) && counter < maxit){
    ## we keep iterating as long as all elements of the gradient vector have 
    # an absolute value lower than tol times the absolute value of the objective 
    ## function plus fscale.
    
    
    f0 <- ftheta           ## Save most recent value of ftheta
    
    ## We form the induced hessian matrix evaluated at the new value for theta. 
    ## By the construction of our hessian function, we know that it will be
    ## positive definite.  
    f2prime <- hessian(theta, grad, eps, hess,...)
    
    ## Now we perform a Cholesky decomposition for the hessian and we calculate 
    ## the new proposal for the optimal theta, which is given by the old value 
    ## for theta plus (Hessian_f(theta))^(-1)*gradient_f(theta)
    
    R <- chol(f2prime) # Cholesky Decomposition of hessian
    theta_new <- theta - backsolve(R, forwardsolve(t(R), fprime))# Update theta
    ftheta <- func(theta_new) # Evaluate f at the new value for theta, theta_new
    
    
     ## To avoid any possibility of divergence, we have to check whether each new 
     ## proposed value for theta actually reduces the objective function. If it 
     ## does not, we backtrack towards the previous value until it does.
     half_count <- 1 # initializing halving counters
     stepsize <- 0.5
     
     ## Loop for halving
     while(ftheta >= f0 && half_count < max.half){
       ## the function exits the loop when the new value for theta has finally 
       ## reduced the objective function, or if it has reached the maximum limit
       ## for halving iterations given by max.half
       
       
       # Update theta_new
       theta_new <- theta - backsolve(R, forwardsolve(t(R), fprime))
       
       # Update stepsize and half_count
       stepsize <- stepsize / 2
       half_count <- half_count+1
       
       # Reevaluate the objective function at the new value of theta
       ftheta <- func(theta_new)
       
       }
     
     if(half_count==max.half){
       ## We issue a warning if the method has failed to reduce the objective 
       ## function after  trying max.half step halvings
       warning('No longer reducing step')
     }
     
     # Calculate the gradient at the induced value of theta
     fprime <- grad(theta_new)
     theta <- theta_new ## Is this line needed?
     counter <- counter+1
     
  }
  
  if(count == maxit){
    warning('Maximum iterations reached without convergence')
    
  }
  ## We calculate the inverse of the hessian matrix at the point of convergence
  invhess <- chol2inv(chol(hessian(theta_new, grad, eps, hess, ...)))
  
  
  
  ## Output our results
  output <- list(f = ftheta, theta = theta, iter = counter, g = fprime, Hi = invhess)
  return(output)
}


   

