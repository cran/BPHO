
compress <- function (
    features, is_sequence=1, order=ncol(features), ptn_file, 
    nos_features=c(), comp=1, quiet=1, no_cases_ign=0 )

{
  
   if (is.vector(features)) features <- matrix(features,1,)
  
   if (is_sequence == 1 & order < ncol(features)) 
      features <- features [, 1:order, drop = FALSE]

   if (file.exists (ptn_file) ) file.remove (ptn_file)

   if (is_sequence == 0 & order == ncol (features) )
   cat( "Note: 'order' is set to number of features. Training may run long.\n")

   if ( order > ncol(features) ) 
   {
      order <- ncol (features)
      warning ("'order' is bigger than p (number of features). I will use p.")
   }

   if (length(nos_features) == 0) nos_features <- apply (features, 2, max)

   .C ("R_compress", as.integer(is_sequence), as.integer(order),
       nrow(features),ncol(features),
       as.integer(t(features)),as.integer(nos_features),
       as.integer(no_cases_ign),ptn_file,as.integer(quiet),
       as.integer(comp),PACKAGE="BPHO") -> tmp #return nothing

}

display_ptn <- function(ptn_file, ix_g=c())
{
  info <- integer(6)

  names(info) <- 
    c("is_sequence","order","#groups", "#patterns","#cases","#features")
  
  ix_g <- as.integer(ix_g - 1) ## minus 1 to be compatible with C codes
  
  if(file.exists(ptn_file))
  {
    info <- .C("R_display_compress",ptn_file,length(ix_g),ix_g,info,
                PACKAGE="BPHO") [[4]]
    info
  }
  else
      cat("File", ptn_file, "does not exist","\n")
}

