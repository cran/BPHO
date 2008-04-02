.First.lib <- function(lib,pkg)
{
   library.dynam("BPHO",pkg,lib)
   cat("Bayesian Prediction with High-order Models(BPHO) loaded\n", 
       "COPYRIGHT 2007-2008 (c) Longhai Li (http://math.usask.ca/~longhai)\n",
	"Type ?begin.BPHO for help to this package\n",sep="")

}

