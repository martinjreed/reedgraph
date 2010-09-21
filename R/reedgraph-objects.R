### Just a simple example showing experiments with S4 classes.
### NOTE methods still PASS BY VALUE, so you generally cannot
### changes slots (class members)
### However setReplaceMethod allows the objectmethod() <- new value
### operation to change slots. If new value is a list
### then it is possible to do as per normal with a C++ or Java member
### method. BUT it is not very worthwhile as setReplaceMethod
### works by creating a copy, changing the slot in the copy then
### replacing original by the copy so
### object <- method(object,value) is just as efficient and
### more or less the same as
### method(object) <- value


setClass("rgtest",representation(y="numeric"))

setGeneric("getY", function(this) standardGeneric("getY"))        
setMethod("getY",
          signature("rgtest"),
          function(this) {
            return(this@y)
          })

### This is how to set an internal slot
### Note it is actually implemented pass-by-value, at the return
### the (changed) copied object replaces the original
setGeneric("setY<-", function(this,value) standardGeneric("setY<-"))        
setReplaceMethod("setY",
          signature("rgtest","numeric"),
          function(this,value) {
            this@y <- value
            this
          })

