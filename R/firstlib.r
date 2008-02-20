.First.lib <- function(lib,pkg)
{
   library.dynam("BPHO",pkg,lib)
   cat("Bayesian Prediction with High-order Models(BPHO) loaded\n", 
       "COPY RIGHT (c) Longhai Li, http://math.usask.ca/~longhai\n",
	"Type ?begin.BPHO for help\n",sep="")

}

