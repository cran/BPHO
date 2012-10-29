#y.tr is assumed to be coded by 1,2,...,no.cls
training <- function (
    ################## specify for data  ###############################
    train_y, no_cls, mc_file, ptn_file, 
    ################## specify for priors  #############################
    alpha=1, log_sigma_modes=c(), log_sigma_widths=c(), 
    ################## specify for slice sampling ######################
    iters_mc=200, iters_bt=10,  iters_sgm=50,
    w_bt=5, m_bt=20, w_sgm=1, m_sgm=20, ini_log_sigmas=c(), 
    start_over=TRUE )
{   
    order <- display_ptn (ptn_file) ["order"]
       
    if (any (train_y < 1) || any (train_y > no_cls) ) 
        stop ("response cannot smaller than 1 or larger than 'no_cls'")

    ## specifying priors separately for Gaussian and Cauchy if they are empty
    if (alpha == 1) 
    {  if (is.null (log_sigma_modes) )
          log_sigma_modes <- c (1, seq (-1, by = log (0.7), length = order) )
       if (is.null (log_sigma_widths) )
          log_sigma_widths <- c (1e-10, rep (log(10)/2, order) )
    }
    
    if (alpha == 2)
    {  if (is.null (log_sigma_modes) )
          log_sigma_modes <- c (10, seq (5, by = log (0.7), length = order) )
       if (is.null (log_sigma_widths) )
          log_sigma_widths <- c (1e-10, rep (log(10)/2, order) )
    }
    
    ## specifying initial hyperparameters if they are empty
    if (is.null (ini_log_sigmas) )
    {  ini_log_sigmas <- log_sigma_modes
    }
    
    
    ## integerize some variables
    no_cls <- as.integer(no_cls)
    alpha <- as.integer(alpha)
    iters_mc <- as.integer(iters_mc)
    iters_bt <- as.integer(iters_bt)
    iters_sgm <- as.integer(iters_sgm)
    m_bt <- as.integer(m_bt)
    m_sgm <- as.integer(m_sgm)

    if (start_over == TRUE & file.exists (mc_file) ) file.remove (mc_file)
    
    ## transforming response by minus 1
    train_y <- as.integer(train_y-1) ## in C code, response starts from 0
    
    .C("R_training", mc_file,ptn_file,train_y,no_cls,
       iters_mc,iters_bt,iters_sgm,w_bt,w_sgm,m_bt,m_sgm,alpha,
       log_sigma_widths, log_sigma_modes, ini_log_sigmas, PACKAGE="BPHO") -> tmp
}


summary_mc <- function(mc_file)
{
   if(!file.exists(mc_file))
      stop(paste(mc_file,"does not exist"))

   out <- integer(5)
   names(out) <- c("#iters","#class","#groups","order","alpha")

   .C("summary_mc",mc_file,out,PACKAGE="BPHO")[[2]]
}


read_mc <- function(
  group,ix, mc_file, iter_b=1,forward=1,samplesize=c(),quiet=1)
{
   if(!file.exists(mc_file))
      stop(paste(mc_file,"does not exist"))

   info_mc <- summary_mc(mc_file)
   iters <- info_mc["#iters"]


   if (length (samplesize) == 0 || 
       any (c (iter_b < 1, iter_b > iters, 
               iter_b + forward * (samplesize-1) > iters) ) )
   {   
       cat("Use all",iters, "iterations in",mc_file,"\n")
       iter_b <- 1
       forward <- 1
       samplesize <- iters
   }
   
   if(samplesize == 0) {
      cat("The number of iterations is 0, summary information is displayed:\n")
      return(info_mc)
   }
   .C("R_read_mc",group,as.integer(ix-1),mc_file,as.integer(iter_b-1),
      as.integer(forward),as.integer(samplesize), 
      out=rep(0,samplesize),PACKAGE="BPHO")$out
}


read_betas <- function(
   ix_g, ix_cls, mc_file, iter_b=1,forward=1,samplesize=c(),quiet=1)
{
   if(!file.exists(mc_file)) stop(paste(mc_file,"does not exist"))

   no_cls <- summary_mc(mc_file)["#class"]

   read_mc("betas", (ix_g - 1) * no_cls + ix_cls, mc_file, 
           iter_b, forward, samplesize, quiet)
}


medians_betas <- function (mc_file, iter_b=1, forward=1, samplesize=c() )
{
   info_mc <- summary_mc(mc_file)
   
   medians_betas <- matrix (0, info_mc["#groups"], info_mc["#class"])
   colnames (medians_betas) <- paste ("cls", seq (1, info_mc["#class"]), 
             sep ="_" )
   for (i_g in seq (1, info_mc["#groups"]) )
       for (i_cls in seq (1, info_mc["#class"]) )
       {
           medians_betas [i_g, i_cls] <- 
	   median (
	     read_betas (i_g, i_cls, mc_file, iter_b, forward, samplesize) )
       }
   max_medians_betas <- apply (abs (medians_betas), 1, max)
   order_g <- order (max_medians_betas, decreasing = TRUE)
   medians_betas <- medians_betas[order_g,]
   cbind(groupid = order_g, medians_betas)
}

