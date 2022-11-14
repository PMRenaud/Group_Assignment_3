ra <- function(th, k = 3){
  4*th[1]^2 + 2*th[2]^2-4 + 6*th[2] - 7*th[1] 
}

ga <- function(th, k = 3){
  c(8*th[1]-7, 4*th[2]+6)
}

ha <- function(th,k=3){
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




newt <-function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  
  ## Calculating next iteration using 
  ## x_k+1 = x_k - [f''(x_k)]^{-1}/ f'(x_k)
  
  f <- func(theta);f
  f_prime <- grad(theta);f_prime
  f_2prime<- hess(theta);f_2prime
  count <- 1
  
  while(min(abs(f_prime)) > (tol * (abs(f) + fscale)) && count < maxit){
    # Only works if we replace 1 with a half need to check the max.half stuff
    theta <- theta - 0.5 * chol2inv(f_2prime) %*% f_prime
    f <- func(theta)#;cat(f)
    f_prime <- grad(theta)#;cat(f_prime)
    f_2prime<- hess(theta)#;cat(f_2prime)
    count <- count + 1
  
  }
  cat(count)
  return(theta)
}

theta = c(2,-2)
newt(theta, ra, ga, ha)
theta <- newt(theta, ra, ga, ha); theta


