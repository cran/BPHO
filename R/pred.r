predict_bpho <- function(test_x, mc_file, ptn_file, iter_b, forward, samplesize)
{
   if (is.vector(test_x)) test_x <- matrix(test_x,1,length(test_x))

   info_mc <- summary_mc (mc_file)
   info_ptn <- display_ptn (ptn_file)
   no_cls <- info_mc ["#class"]
   order <- info_ptn ["order"]
   is_sequence <- info_ptn ["is_sequence"]
   if (info_ptn["#groups"] != info_mc["#groups"])
      stop(paste(mc_file,"and",ptn_file,"does not match"))

   if (any( c(iter_b ,iter_b + forward * (samplesize-1)) > info_mc["#iters"]) )
      stop("The specified MC indice out of range")

   if (is_sequence == 1 & order < ncol(test_x)) 
   {
      test_x <- test_x [, 1:order, drop = FALSE]
   }
   
   n <- nrow (test_x)
   p <- ncol (test_x)
   
   pred_probs <- matrix(
     .C("R_pred",n,p,as.integer(t(test_x)),
        as.integer(no_cls),prediction=rep(0,n*no_cls),
        mc_file,ptn_file,as.integer(iter_b),as.integer(forward),
        as.integer(samplesize),PACKAGE="BPHO")$prediction,
        n,no_cls,byrow=TRUE)

   y_pred <- apply(pred_probs,1,which.max)

   data.frame(pred_probs=pred_probs,y_pred=y_pred)
}

evaluate_prediction <- function(test_y,pred_result,file_eval_details=c())
{
   if(is.vector(pred_result))
     stop ("There is only a test case. Cannot make evaluation.")
    
   if( nrow (pred_result) != length (test_y) ) 
     stop ("True values of response and prediction result NOT match")
   wrong <- 1*(pred_result[,"y_pred"] != test_y)
   error_rate <- mean(wrong)
   probs_true <- test_y
   for( i in 1:length(test_y) )
      probs_true[i] <- pred_result[i,test_y[i]]
   amll <- mean( -log(probs_true) )
   eval_details <- data.frame(y_true=test_y,y_pred=pred_result[,"y_pred"],
                       wrong=wrong,probs_at_true=probs_true)
   if( !is.null(file_eval_details) )
   {
      sink(file_eval_details)
      print(eval_details)
      sink()
   }
   list(eval_details=eval_details, error_rate=error_rate, amll=amll)
}

split_cauchy <- function(n,s, sigma1,sigmasum,debug=1)
{
   x <- rep(0,n);
   for(i in 1:n)
   {  x[i] <- .C("split_cauchy",0,s,sigma1,sigmasum,as.integer(debug),
                 PACKAGE="BPHO")[[1]]
   }
   x
}
