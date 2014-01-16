require("gmp")
require("Rcpp")

## find a prime either a or less than a
## as frequency of primes is generally frequent it should be fairly
## fast even for quite large a
previousprime <- function(a) {
  a <- sub.bigz(a,1)
  while(!isprime(a)) {
    a <- sub.bigz(a,1)
  }
  return(a)
}

## Estimates number of primes as val/ln(val). This is an underestimate by about 10%
est.number.primes <- function(val) {
  return(floor(as.double(val)/log(as.double(val))))
  
}

## Generates a random prime in range 1 - somewhere less than max
## do not use for any kind of good distribution, it will tend not
## to return primes near max (generally at least up to max/2)
## just a rough and ready routine
random.prime <- function(max,base=NULL) {
  
  index <- floor(log(as.double(max))/log(2))
  val <- urand.bigz(1,index)
  if(!is.null(base)) modulus(val) <- base
  return(val)
}

## Generate the multiplicative inverse of val using Bezoult's identity from
## Extended Euclidian Algorithm
mul.inverse <- function(val){
  return(as.bigz(gcdex(val,modulus(val))[2],base))
}


