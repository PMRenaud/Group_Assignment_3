## Things to check:
## Holds for dimension 1
## Finiteness stuff?
## Warnings and stops
## ... arguments

## Run everything in the debugger before submission


finite_differencing <- function(theta,grad,eps,...){

  dim = length(theta)
  Hfd <- matrix(0,dim,dim) ## initializing finite difference Hessian
  for (i in 1:dim) { ## Looping over parameters
    th1 <- theta;
    th1[i] <- th1[i] + eps  ## Increase th1[i] by eps 
    grad1 <- grad(th1,...)  ## compute resulting gradient 
    Hfd[i,] <- (grad1 - grad(theta,...))/eps ## Approximate second derivatives
  }
  Hfd <- (Hfd+t(Hfd))/2
  return(Hfd)
}

pos_def <- function(A){
  n <- 1e-6
  while(inherits(try(chol(A), silent = TRUE),"matrix") == FALSE){
    A <- A + n * diag(1, dim(A)[1])
    n <- 10 * n
  }
  return(A)
}


## calculates hessian
hessian <- function(theta, grad, eps, hess = NULL, ...){
  if(is.null(hess)){
    f2prime <- finite_differencing(theta, grad, eps, ...)
  }
  
  else{
    ## Update Hessian
    f2prime <- hess(theta)
  }
  
  f2prime <- pos_def(f2prime)
  return(f2prime)
}

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,
                 maxit=100,max.half=20,eps=1e-6){
  
  # Calculate ftheta and fprime
  ftheta <- func(theta)
  fprime <- grad(theta)
  
  # Initialise count
  counter <- 0
  
  while(max(abs(fprime)) > tol*(abs(ftheta)+fscale) && counter < maxit){
    
    
    ## Save most recent value of ftheta
    f0 <- ftheta
    
    
    f2prime <- hessian(theta, grad, eps, hess, ...)
    
    R <- chol(f2prime) # Cholesky Decomposition of hessian
    
    # Update theta
    theta_new <- theta - backsolve(R, forwardsolve(t(R), fprime))
    
    # update ftheta
    ftheta <- func(theta_new)
    
    half_count <- 1 # Think about this?
    stepsize <- 0.5
    
    
    ## Loop for halfing
    while(ftheta >= f0 && half_count < max.half){
      
      # Update theta_new
      theta_new <- theta - backsolve(R, forwardsolve(t(R), fprime))
      
      # Update stepsize and half_count
      stepsize <- stepsize / 2
      half_count <- half_count+1
      
      # update ftheta
      ftheta <- func(theta_new)
      
    }
    
    # Calculate fprime
    fprime <- grad(theta_new)
    
    theta <- theta_new
    
    counter <- counter+1
    
  }
  
  

  invhess <- chol2inv(chol(hessian(theta_new, grad, eps, hess, ...)))
  
  
  
  ## Output the list of stuff - Check this
  output <- list(f = ftheta, theta = theta, iter = counter, g = fprime, Hi = invhess)
  return(output)
}



##-----------------------------------------------------------------------------
## Tests

rb <- function(th,k=2){
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
} 

gb <- function(th,k=2){
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}

hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

theta = c(12,12)
newt(theta, rb, gb, hb)

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
theta = c(12,12)
newt(theta, rd, gd, hd)

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

rd <- function(th){
  th[1]^4 + 3*th[1]^2*th[2]^2 + 2*th[2]^2 
}

gd <- function(th){
  c(4*th[1]^3 + 6*th[1]*th[2]^2, 6*th[1]^2*th[2]+4*th[2])
}

hd <- function(th){
  h <- matrix(0,2,2)
  h[1,1] <- 12*th[1]^2 + 6*th[2]^2
  h[2,2] <- 6*th[1]^2+4
  h[1,2] <- 12*th[1]*th[2]
  h[2,1] <- 12*th[1]*th[2]
  h
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



pos_def(H)




