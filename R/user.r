comp_train_pred <- function(
   ################## specify data information  #####################
   test_x,train_x,train_y,no_cls=c(),nos_fth=c(),
   ################## specify for compression #######################
   is_sequence=1,order,ptn_file=".ptn.log",new_compression=1,do_comp=1,
   ###################### specify for priors  #######################
   alpha=1,log_sigma_widths=c(),log_sigma_modes=c(),
   ################# specify for mc sampling ########################
   mc_file=".mc.log",start_over=FALSE,iters_mc=200,iters_bt=10,
   iters_sgm=50,w_bt=5,w_sgm=1,m_bt=20,m_sgm=20,ini_log_sigmas=c(),
   ################### specify for prediction #######################
   pred_file=c(),iter_b = 100,forward = 1,iters_pred = 100)

{   if(length(log_sigma_modes) != order + 1){
         log_sigma_modes = c(5, seq(0,-5,length=order))
    }
    if(length(log_sigma_widths) != order + 1)
        log_sigma_widths = c(1e-10,rep(2,order))

    #set the initial sigmas by default
    if(length(ini_log_sigmas)!= order + 1)
         ini_log_sigmas = log_sigma_modes

    #set the prior by default
    times <- rep(0,3)
    names(times) <- c("compressing","training","prediction")

    #compressing parameters
    if(new_compression  == 1 || !file.exists(ptn_file)){
       order <- min(order,ncol(train_x))
       if(is_sequence == 1 & order < ncol(train_x)) {
          train_x <- train_x[,1:order,drop=FALSE]
	  nos_fth <- nos_fth[1:order]
       }
       times[1] <- system.time(
       compress(train_x,nos_fth,0,ptn_file,1,do_comp,is_sequence,order))[1]
    }

    #running markov chain
    if(iters_mc > 0){
       if(start_over)  file.remove(mc_file)
       times[2] <- system.time(
       training(mc_file,ptn_file,train_y,no_cls,
                alpha,log_sigma_widths,log_sigma_modes,
                ini_log_sigmas,iters_mc,iters_bt,iters_sgm,
                w_bt,w_sgm,m_bt,m_sgm))[1]
    }

    #making prediction
    pred_result <- c()
    if(iters_pred > 0){
       order <- min(order,ncol(test_x))
       if(is_sequence == 1 & order < ncol(test_x)) {
          test_x <- test_x[,1:order,drop=FALSE]
       }
       times[3] <- system.time(
       pred_result <- predict_bpho(test_x,no_cls,mc_file,ptn_file,
				iter_b,forward,iters_pred))[1]
    if(!is.null(pred_file))
        write.table(pred_result, file = pred_file, row.names = FALSE, sep=",")
    }
    list(pred_result=pred_result,
         files=c(ptn_file,mc_file,pred_file),times=times)
}


cv_comp_train_pred<- function(
   ###################### Specify data,order,no_fold #################
   no_fold=10,train_x,train_y,no_cls=c(),nos_fth=c(),
   #################### specify for compressing#######################
   is_sequence=1,order,ptn_file=".ptn.log",new_compression=1,do_comp=1,
   ###################### specify for priors  ########################
   alpha=1,log_sigma_widths=c(),log_sigma_modes=c(),
   ################# specify for mc sampling #########################
   mc_file=".mc.log",iters_mc=200,iters_bt=10,iters_sgm=50,
   w_bt=5,w_sgm=1,m_bt=20,m_sgm=20,ini_log_sigmas=c(),
   ################### specify for prediction ########################
   pred_file = c(),iter_b = 100,forward = 1,iters_pred = 100)
{
    #set nos_fth and no_cls by default
    if(length(nos_fth) != ncol(train_x))
       nos_fth <- apply(train_x,2,max)

    if(length(no_cls) == 0)
       no_cls = max(train_y)

    n <- nrow(train_x)
    m <- floor(n / no_fold)
    rm <- n - m * no_fold
    pred_result <- data.frame(pred_probs=matrix(0,n,no_cls),y_pred=rep(0,n))

    #start cross-validation
    for(i in 1:no_fold) {

	#prepare test cases
        if( i <= rm)
           testcases <- i + seq(0,m) * no_fold
        else testcases <- i + seq(0, m-1) * no_fold

	out <- comp_train_pred(
	   ################## specify data information  #########
           train_x[testcases,,drop=FALSE],train_x[-testcases,,drop=FALSE],
	   train_y[-testcases],no_cls,nos_fth,
	   ################## specify for compression ###########
           is_sequence,order,ptn_file,new_compression,do_comp,
           ###################### specify for priors  ###########
	   alpha,log_sigma_widths,log_sigma_modes,
           ################# specify for mc sampling ############
	   mc_file,start_over=TRUE,iters_mc,iters_bt,
	   iters_sgm,w_bt,w_sgm,m_bt,m_sgm,ini_log_sigmas,
           ################## specify for prediction ############
	   pred_file=c(),iter_b,forward,iters_pred)

	pred_result[testcases,] <- out$pred_result
	times <- out$times
	files <- out$files
    }
    if( !is.null(pred_file) )
        write.csv(pred_result,file = pred_file, row.names = FALSE)

    list(times=times,pred_result=pred_result,
         files=c(ptn_file,mc_file,pred_file))
}
