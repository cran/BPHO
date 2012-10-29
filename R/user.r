comp_train_pred <- function(
    ################## specify data information  #######################
    test_x = c(),train_x,train_y,no_cls,nos_features=c(), 
    ################## specify for compression #########################
    is_sequence=1,order,ptn_file=c(), new_comp = 1, comp = 1, 
    no_cases_ign = 0,
    ################## specify for priors  #############################
    alpha=1,log_sigma_widths=c(),log_sigma_modes=c(),
    ################## specify for mc sampling #########################
    mc_file=c(),start_over=TRUE,iters_mc=200,iters_bt=10,
    iters_sgm=50,w_bt=5,m_bt=20,w_sgm=1,m_sgm=20,ini_log_sigmas=c(),
    ################## specify for prediction ##########################
    pred_file=c(),iter_b = 100,forward = 2,samplesize = 50)

{   
  times <- rep(0,3)
  names(times) <- c("compressing","training","prediction")

  if (is.null (ptn_file) ) 
      ptn_file <- paste(".ptn",".o",order, ".log",sep="")

  if (is.null (mc_file) ) 
      mc_file <- paste(".mc",".o",order, ".alpha", alpha,".log",sep="")

  if (is.null (pred_file)) 
      pred_file <- paste(".pred",".o",order,".alpha",alpha,".csv",sep="")

  #compressing parameters (transforming features)
  if(new_comp  == 1 || !file.exists(ptn_file))
  {
    times[1] <- system.time(
    compress(features = train_x, nos_features = nos_features,
             ptn_file = ptn_file, comp = comp, is_sequence = is_sequence, 
             order = order, no_cases_ign = no_cases_ign, quiet = 1)) [1]
  }

  #running markov chain
  if(iters_mc > 0)
  {

    times[2] <- system.time(
    training(mc_file = mc_file, ptn_file = ptn_file, 
             train_y = train_y, no_cls = no_cls,
             alpha = alpha, log_sigma_widths = log_sigma_widths,
             log_sigma_modes = log_sigma_modes,
             ini_log_sigmas = ini_log_sigmas, start_over = start_over,
             iters_mc = iters_mc, iters_bt = iters_bt, iters_sgm = iters_sgm,
             w_bt = w_bt, m_bt = m_bt, w_sgm = w_sgm,  m_sgm = m_sgm)) [1]
  }

  #making prediction
  pred_result <- c()
  
  if (samplesize > 0 & !is.null (test_x) )
  {
    times[3] <- system.time (
    pred_result <- predict_bpho ( 
         test_x = test_x, mc_file = mc_file, ptn_file = ptn_file,
         iter_b = iter_b, forward = forward, samplesize = samplesize) )[1]
    write.table (pred_result, file = pred_file, row.names = FALSE, sep=",")
  }
  
  list (pred_result = pred_result,files = c(ptn_file,mc_file,pred_file),
        times = times)
}

cv_comp_train_pred <- function(
    ###################### Specify data,order,no_fold ##################
    no_fold=10,train_x,train_y,no_cls=c(),nos_features=c(),
    ###################### specify for compressing #####################
    is_sequence=1,order,ptn_file=c(),comp=1, no_cases_ign = 0,
    ###################### specify for priors  #########################
    alpha=1,log_sigma_widths=c(),log_sigma_modes=c(),
    ###################### specify for mc sampling #####################
    mc_file=c(),iters_mc=200,iters_bt=10,iters_sgm=50,
    w_bt=5,w_sgm=1,m_bt=20,m_sgm=20,ini_log_sigmas=c(),
    ###################### specify for prediction ######################
    pred_file=c(),iter_b=100,forward=2,samplesize=50)
{
  #set nos_features and no_cls by default
  n <- nrow(train_x)
  if( length(nos_features) != ncol(train_x) ) 
      nos_features <- apply(train_x,2,max)
  if( length(no_cls) == 0 ) no_cls <- max(train_y)
  pred_result <- data.frame(pred_probs=matrix(0,n,no_cls),y_pred=rep(0,n))

  if (is.null (ptn_file)) 
      ptn_file <- paste(".ptn",".o",order, ".log",sep="")

  if (is.null (mc_file)) 
      mc_file <- paste(".mc",".o",order, ".alpha", alpha,".log",sep="")

  if (is.null (pred_file)) 
      pred_file <- paste(".pred",".o",order,".alpha",alpha,".csv",sep="")

  m <- floor(n / no_fold)
  rm <- n - m * no_fold
  #start cross-validation
  times <- rep (0,3)
  for( i in 1:no_fold ) 
  {
    #prepare test cases
    if( i <= rm ) testcases <- i + seq(0,m) * no_fold
      else testcases <- i + seq(0, m-1) * no_fold

    out <- comp_train_pred(
            ################## specify data information  #####################
            test_x = train_x[testcases,,drop=FALSE],
            train_x = train_x[-testcases,,drop=FALSE],
            train_y = train_y[-testcases], 
            no_cls = no_cls, nos_features = nos_features,
            ################## specify for compression #######################
            is_sequence = is_sequence,
            order = order,ptn_file = ptn_file, new_comp = 1,comp = comp,
            no_cases_ign = no_cases_ign,
            ###################### specify for priors  #######################
            alpha = alpha,log_sigma_widths = log_sigma_widths,
            log_sigma_modes = log_sigma_modes,
            ################# specify for mc sampling ########################
            mc_file = mc_file,start_over=TRUE,
            iters_mc = iters_mc, iters_bt = iters_bt, iters_sgm = iters_sgm,
            w_bt = w_bt, w_sgm = w_sgm, m_bt = m_bt, m_sgm = m_sgm, 
            ini_log_sigmas = ini_log_sigmas,
            ################## specify for prediction ########################
            pred_file = c(),iter_b = iter_b, forward = forward,
            samplesize = samplesize)

    pred_result[testcases,] <- out$pred_result
    times <- times + out$times
  }
  
  ## evaluate predictions
  eval_result <- evaluate_prediction (test_y = train_y, 
                 pred_result = pred_result, file_eval_details = pred_file)
  
  list (eval_details = eval_result$eval_details, 
        error_rate = eval_result$error_rate, 
        amll = eval_result$amll, 
        times = times, files = c(ptn_file,mc_file,pred_file) 
       )
}
