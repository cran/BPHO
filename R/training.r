#y.tr is assumed to be coded by 1,2,...,no.cls
training <- function(
    mc_file,ptn_file,
    train_y,no_cls,alpha,
    log_sigma_widths,log_sigma_modes,ini_log_sigmas,
    iters_mc,iters_bt,iters_sgm,
    w_bt,w_sgm,m_bt,m_sgm)
{


    train_y <- as.integer(train_y-1);
    no_cls <- as.integer(no_cls);
    iters_mc <- as.integer(iters_mc);
    iters_bt <- as.integer(iters_bt);
    iters_sgm <- as.integer(iters_sgm);
    m_bt <- as.integer(m_bt);
    alpha <- as.integer(alpha);
    m_sgm <- as.integer(m_sgm);

    .C("R_training", mc_file,ptn_file,train_y,no_cls,
       iters_mc,iters_bt,iters_sgm,w_bt,w_sgm,m_bt,m_sgm,alpha,
       log_sigma_widths, log_sigma_modes, ini_log_sigmas, PACKAGE="BPHO")[c()]
}


display_mc <- function(mc_file)
{
   if(!file.exists(mc_file))
      stop(paste(mc_file,"does not exist"))

   out <- integer(5)
   names(out) <- c("#iters","#class","#groups","order","alpha")

   .C("display_mc",mc_file,out,PACKAGE="BPHO")[[2]]
}


read_mc <- function(mc_file,group,ix, iter_b=0,forward=1,n=c(),quiet=1)
{
   if(!file.exists(mc_file))
      stop(paste(mc_file,"does not exist"))

   info_mc <- display_mc(mc_file)
   iters <- info_mc["#iters"]


   if(length(n)==0 || 
      any(c(iter_b < 0, iter_b > iters - 1,iter_b + forward * (n-1) > iters-1)))
   {   
       cat("Use all",iters, "iterations in",mc_file,"\n")
       iter_b =0
       forward = 1
       n = iters
   }
   if(n == 0) {
      cat("The number of iterations is 0, summary information is displayed:\n")
      return(info_mc)
   }
   .C("R_read_mc",group,as.integer(ix),mc_file,as.integer(iter_b),
      as.integer(forward),as.integer(n),out=rep(0,n),PACKAGE="BPHO")$out
}


read_betas <- function(mc_file,ix_g,ix_cls,iter_b=0,forward=1,n=c(),quiet=1)
{
   if(!file.exists(mc_file))
      stop(paste(mc_file,"does not exist"))

   info_mc <- display_mc(mc_file)
   no_cls <- info_mc["#class"]
   iters <- info_mc["#iters"]


   if(length(n)==0 || 
      any(c(iter_b < 0, iter_b > iters - 1,iter_b + forward * (n-1) > iters-1)))
   {   cat("Use all iterations in",mc_file,"\n")
       iter_b =0
       forward = 1
       n = iters
   }
   if(n == 0) {
      cat("The number of iterations is 0, summary information is displayed:\n")
      return(info_mc)
   }

   if(!file.exists(mc_file))
      stop(paste(mc_file,"does not exist"))

   read_mc(mc_file,"betas",ix_g*no_cls+ix_cls - 1,iter_b,forward,n,quiet)
}

display_a_beta <- function(id_beta,mc_file,ptn_file)
{
   info_mc <- display_mc(mc_file)

   no_cls <- info_mc["#class"]
   names(no_cls) <- c()
   i.group <- floor(id_beta/no_cls)
   i.cls <- id_beta - i.group*no_cls + 1

   cat("\nThe information on the pattern related to the beta #",
       id_beta,":\n",sep="")

   display_ptn(ptn_file,i.group)

   cat("\nThe Markov chain samples of this beta for group #",i.group,
       " Class #",i.cls,":\n\n",sep="")

   beta <- read_mc(mc_file,"betas",id_beta)

   list(beta=beta,i.group=i.group,i.cls=i.cls)
}

calc_medians_betas <- function(mc_file,iter_b=0,forward=1,n=c())
{
   info_mc <- display_mc(mc_file)
   no_beta <- info_mc["#groups"] * info_mc["#class"]
   
   medians_betas <- rep(0,no_beta)

   for(i_bt in seq(0,no_beta-1)){
      medians_betas[i_bt+1] <- 
         median(read_mc(mc_file,"betas",i_bt, iter_b,forward,n))
   }
   rbind(ids_betas=seq(0,no_beta-1),medians_betas=medians_betas)   
}

