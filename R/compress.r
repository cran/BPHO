
#this is a R wrapper function that calls C function "R_compress" to
#find the compressed representation of the parameters for interactions
# Arguments:
# features     the data. The rows are observation,
#              and the columns are features
# nos_fth      a vector indicating the number of possibilities of
#              features
# depth        the deepest order of interactions
# patterns.rda the filename which saves output "Values"

# Values:
# starts_cases  a vector containing the starting index in "values_cases"
#               of indice of cases for each group
# values_cases  a vector of all indice of cases in all groups in order
# starts_ptns   a vector containing the starting index in "values_ptns"
#               of indice of patterns for each group
# values_ptns   a vectors of all indice of patterns in all groups in order

compress <- function(features,nos_fth=c(),no_cases_ign=0,
                     ptn_file=".ptn_file.log",quiet=1,
		     do_comp=1,sequence=1,order=ncol(features))
{
   if(is.vector(features))
      features <- matrix(features,1,)
   if(length(nos_fth)==0)
      for(i in 1:ncol(features))
         nos_fth[i] = max(features[,i])

   .C("R_compress", as.integer(sequence),as.integer(order),
        nrow(features),ncol(features),
        as.integer(t(features)),as.integer(nos_fth),
        as.integer(no_cases_ign),ptn_file,as.integer(quiet),
	as.integer(do_comp),PACKAGE="BPHO")->tmp

}

display_ptn <- function(ptn_file, gids=c()){
    info <- integer(6)
    names(info) <- c("is.sequence","order","#groups",
                     "#patterns","#cases","#features")
    gids <- as.integer(gids)
    if(file.exists(ptn_file))
       return(.C("R_display_compress",ptn_file,length(gids),gids,info,
                 PACKAGE="BPHO")[[4]])
    else
       cat("File", ptn_file, "does not exist","\n")
}

