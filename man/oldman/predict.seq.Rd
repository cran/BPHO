\name{comp_train_pred_seq}

\alias{comp_train_pred_seq}

\title{Bayesian prediction for a sequence with high order
 interactions of previous states}

\description{This function predicts the next state of a sequence with  its
  previous states by learning from training data. Note that the sequence is
  assumed to be placed in anti-time order, i.e. the first state is  closest to
  the one one want to predict. We take the high-order interactions into account.
  We use a compression technique to reduce the number  of parameters associated
  with high-order interactions. The model is trained  with a sort of Markov
  chain sampling method, called slice sampling, which is  developed by Professor
  Radford Neal in University of Toronto. 

  The compression information and Markov chain iterations are  stored in binary
  log files (Note that they can not be seen with human eyes). The information
  can be displayed with function  `display_compress_seq' or `read_mc_g'. 

  This function is a compact function, and can be used to compress parameters,
  or train the model with Markov chain, or predict for test cases, or
  combinations them, by setting other relevant parameters to be 0.}

\usage{
comp_train_pred_seq (
        #Specify data information
        test_x,train_x,train_y,no_cls=c(),nos_fth=c(),
	#Specify compression information        
        order,ptn_file="ptn.log",new_compression=1,do_comp=1,  
	#Specify prior information     
        alpha=1,sigma_precisions=c(),sigma_means=c(),
	#Specify slice sampling information
        mc_file="mc.log",start_over=FALSE,iters_mc=200,iters_bt=10,
	iters_sgm=10,w_bt=50,w_sgm=50,m_bt=50,m_sgm=50,ini_sigmas=c(),
	#Specify prediction 
        pred_file=c(),iter_b = 100,forward = 1,iters_pred = 200)
}

\arguments{

\item{test_x}{Categorical inputs (also called features,covariates,independent
variables,explanary variables,predictors) of test data on which your  prediction
is based.  The row is subject and the column is input. Must be coded with 1,2,...
(0 is reserved for other  use). }

\item{train_x}{Categorical inputs of training data of the same 
               format as `test_x'.}
	       
\item{train_y}{Response of training data, a vector of length equal to the row 
	       of `train_x'. Must be coded with 1,2,... `no_cls' .}

\item{no_cls}{a number indicating the number of possibilities (classes) of
               response, by default the maximum value in `train_y'.}
	              	        	      
\item{nos_fth}{a vector indicating the number of possibilites (classes) of 
                each feature, by default the maximum value of each feature.}

\item{order}{the order of interactions one want to consider, by default the 
              total number of features, i.e. ncol(train_x).}

\item{ptn_file}{a character sequence specifying the name of the file to which 
		 the compression result is sent. The information of the
		 resulting file can be read using function
		 `display_compress_seq'.}

\item{new_compression}{`new_compression=1' indicates doing the compression
again.}


\item{do_comp}{`do_comp=1' indicates doing compression, otherwise using 
                original parameterization. This is used only to test program.
	       }

\item{alpha}{`alpha=1' means that Cauchy prior is used, `alpha=2' means 
               that Gaussian prior is used.}
	       			       		 	      		
\item{sigmas_precisions}{A vector of precision parameters of the 
			 Inverse-Chisquare prior for `sigmas', starting with
			 order 0 (intercept),  having length of `order'+1.By
			 default, 1e20 for intercept and 1.5 for others.}
			 
\item{sigmas_means}{A vector of mean parameters of the  Inverse-Chisquare prior
for `sigmas', of the same  format as `sigmas_precisions'. By default, if
`Cauchy' prior is used,  0.5 for intercept and 0.005/order for others,  if
`Gaussian' prior is used, 100 for intercept, 0.05/order for
others}			  

\item{mc_file}{A character sequence specifying the name of the file to which
Markov chain is sent} 

\item{start_over}{`start_over=TRUE' means that the old Markov chain file is
deleted before a Markov chain sampling starts, otherwise the Markov chain will  
continue with the last iteration stored in it.}

\item{iters_mc}{Number of iterations of Markov chain that will be run. It could 
  be 0, meaning no Markov chain sampling will be run.}

\item{iters_bt}{Number of iterations of slice sampling for `betas' for each
iteration of Markov chain sampling}

\item{iters_sgm}{Number of iterations of slice sampling for `sigmas' for each
iteration of slice sampling for `betas'}

\item{w_bt}{Width parameter for slice sampling for `betas', i.e. the length of
itervals of every stepping out.}

\item{w_sgm}{Width parameter for slice sampling for `sigmas'.}

\item{m_bt}{Maximum number of stepping out for slice sampling for `betas'.}

\item{m_sgm}{Maximum number of stepping out for slice sampling for `sigmas'.}

\item{ini_sigmas}{Initial values for `sigmas', by default `sigmas_means'.}

\item{pred_file}{A character sequence specifying the file to which the
prediction result is sent. If `pred_file=c()', the prediction result is printed 
out on screen (sent to standard output).}

\item{iter_b}{The first iteration of Markov chain used to make prediction.}

\item{forward}{1 of `forward' Markov chain starting from `iter_b' is used to 
               make prediction.}
\item{iters_pred}{The number of iterations of Markov chain used to 
                   make prediction.}	       

}

\value{

\item{times}{The time in second for, as this order, compressing parameters,
training the model with compressed parameters, predicting on test cases}

\item{pred_result}{a data frame with first 'no_cls' columns being the predictive
probability and the next column being the predicted response value is returned.}

\item{files}{Three character sequences: the 1st is the file storing compression
information, the 2nd is the file storing Markov chain, and the 3rd one is the
file containing the detailed prediction result, i.e., pred_result}

}

\seealso{display_compress_seq, read_mc_g, gendata_seq,evaluate_prediction}


\examples{
#save the following R script and run it after modifying for your problem
#library the pre-complied R package, `lib.loc' specify where 
#the library tree is built
#    library("SeqPredHighOrderInteract",lib.loc="~/rlib")

#generate data with a built-in scheme

#compress parameters, training model, making prediction
mc_file="mc.log"
ptn_file="ptn.log"
pred_file="pred.csv"
order=3
alpha=1

\dontrun{comp_train_pred_seq(
################## specify data information  ###############
test_x=d1.ts$X,train_x=d1.tr$X,train_y=d1.tr$y,
no_cls=2,nos_fth=rep(2,10),
################## specify for compression #################
order=order,ptn_file=ptn_file,new_compression=1,do_comp=1,
###################### specify for priors  #################
alpha=alpha,sigma_precisions=c(),sigma_means=c(),
################# specify for mc sampling ##################
mc_file=mc_file,start_over=TRUE,iters_mc=200,
iters_bt=5,iters_sgm=10,w_bt=50,w_sgm=50,
m_bt=50,m_sgm=50,ini_sigmas=c(),
################## specify for prediction ##################
pred_file=pred_file,iter_b = 100,forward = 2,iters_pred = 50)
}


}


