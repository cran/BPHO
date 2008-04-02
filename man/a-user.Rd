\name{comp_train_pred}

\alias{begin.BPHO}

\alias{comp_train_pred}

\alias{cv_comp_train_pred}

\title{User-level functions for compressing parameters, training the models with MCMC, and making predictions for test cases}

\description{The function \code{comp_train_pred} can be used for three tasks: compressing parameter, training the models with MCMC, and making prediction for test cases. When \bold{\code{new_compression=1}}, it compresses parameters based on training cases  and the information about parameter compression is written to the binary file \code{ptn_file}. When \bold{\code{new_compression=0}}, it uses the existing \code{ptn_file}. When \bold{\code{iters_mc > 0}}, it trains the models with Markov chain Monte Carlo  and the Markov chain iterations are written to the binary file \code{mc_file}. The methods of writing to and reading from the files \code{ptn_file} and \code{mc_file} can be found from the documentations \code{\link{compression}} and \code{\link{training}}. When \bold{\code{iters_pred > 0}}, it predicts the responses of test cases and the result is written to the file \code{pred_file} and also returned as a value of this function.

The function \code{cv_comp_train_pred} is a short-cut function for performing cross-validation with the function \code{comp_train_pred}.

The argument \code{is_sequence=1} indicates that a sequence prediction model is fitted to the data, and \code{is_sequence=1} indicates that a general classification model based on discrete predictor variables is fitted.
}

\usage{
comp_train_pred(
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

cv_comp_train_pred(
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
}

\arguments{

\item{test_x}{Discrete features (also called inputs,covariates,independent variables, explanatory variables, predictor variables) of test data on which the  predictions are based.  The row is subject and the columns are inputs, which are coded with 1,2,..., with 0 reserved to represent that this input is not considered in a pattern. When the sequence prediction models are fitted, it is assumed that the first column is the state closest to the response. For example, a sequence `x1,x2,x3,x4' is saved in \code{test_x} as `x4,x3,x2,x1', for predicting the response `x5'.}

\item{train_x}{Discrete features of training data of the same format as \code{test_x}.}

\item{train_y}{Discrete response of training data, a vector with length equal to the row of \code{train_x}. Assumed to be coded with 1,2,... \code{no_cls} .}

\item{no_cls}{the number of possibilities (classes) of the response, default to the maximum value in \code{train_y}.}

\item{nos_fth}{a vector, with each element storing the number of possibilities (classes) of each feature, default to the maximum value of each feature.}

\item{is_sequence}{\code{is_sequence=1} indicates that sequence prediction models are fitted to the data, and \code{is_sequence=0} indicates that general classification models based on discrete predictor variables are fitted.
}

\item{no_fold}{Number of folders in cross-validation.}

\item{order}{the order of interactions considered, default to the total number of features, i.e. \code{ncol(train_x)}.}

\item{ptn_file}{a character string, the name of the binary file to which the compression result is written. The method of writing to and reading from \code{ptn_file} can be found from the documentation for \code{\link{compression}}.}

\item{new_compression}{\code{new_compression=1} indicates removing the old file \code{ptn_file} if it exists and doing the compression once again. \code{new_compression=0} indicates using the old file \code{ptn_file} without doing compression once again. Note that when \code{new_compression=0}, the specification related to training cases does not take effect.}

\item{do_comp}{\code{do_comp=1} indicates doing compression, and \code{do_comp=0} indicates using original parametrization. This is used only to make comparison. In practice, we definitely recommend using our compression technique to reduce the number of parameters.}

\item{alpha}{\code{alpha=1} indicates that Cauchy prior is used, \code{alpha=2} indicates that Gaussian prior is used.}

\item{log_sigma_widths, log_sigma_modes}{two vectors of length \code{order+1}, which are interpreted as follows: the Gaussian distribution with location \code{log_sigma_modes[o]} and standard deviation \code{log_sigma_widths[o]} is the prior for `log(sigmas[o])', which is the hyperparameter (width parameter of Gaussian distribution or Cauchy distribution) for the regression coefficients (i.e. `beta's) associated with the interactions of order `o'.}

\item{mc_file}{A character string, the name of the binary file to which Markov chain is written. The method of writing to and reading from \code{mc_file} can be found from the documentation for \code{\link{training}}.}

\item{start_over}{\code{start_over=TRUE} indicates that the existing file \code{mc_file} is deleted before a Markov chain sampling starts, otherwise the Markov chain will continue from the last iteration stored in \code{mc_file}. }

\item{iters_mc,iters_bt,iters_sgm}{\code{iters_mc} iterations of super-transition will be run. Each super-transition consists of \code{iters_bt} iterations of updating `beta's, and for each updating of `beta's, the hyperparameters `log(sigma)'s are updated \code{iters_sgm} times. When \code{iters_mc=0}, no Markov chain sampling will be run and other arguments related to Markov chain sampling take no effect.}

\item{w_bt,w_sgm, m_bt,m_sgm}{\code{w_bt} is the amount of stepping-out in updating `beta' with slice sampling,  \code{m_bt} is the maximum number of stepping-out in slice sampling for updating `beta'. \code{w_sgm} and \code{m_sgm} are intepreted similarly for sampling for `log(sigma)'.}

\item{ini_log_sigmas}{Initial values of `log(sigma)', default to \code{log_sigma_modes}.}

\item{pred_file}{A character string, the name of the file to which the prediction result is written. If \code{pred_file=c()}, the prediction result is printed out on screen (or sent to standard output).}

\item{iter_b, forward, iters_pred}{Starting from \code{iter_b}, one of every \code{forward} Markov chain samples, with the number of total samples being <= \code{iters_pred} and the maximum usable in the file \code{mc_file}, is used to make prediction. }

}

\value{

\item{times}{The time in second for, as this order, compressing parameters,
training the model, predicting for test cases}

\item{pred_result}{a data frame with first \code{no_cls} columns being the predictive probability and the next column being the predicted response value is returned.}

\item{files}{Three character strings: the 1st is the name of the file storing compression information, the 2nd is the name of the file storing Markov chain, and the 3rd one is the name of the file containing the detailed prediction result, i.e., \code{pred_result}}

}

\seealso{\code{\link{gendata}},\code{\link{compression}},\code{\link{training}},\code{\link{prediction}}}

\author{Longhai Li, \url{http://math.usask.ca/~longhai}}

\references{\url{http://math.usask.ca/~longhai/doc/seqpred/seqpred.abstract.html}}

\examples{
##  loading package
##  library("BPHO",lib.loc="~/rlib")

#####################################################################
########The following are demonstrations of using the whole package
#####################################################################

## generate data from a hidden Markov model
data_hmm <- gen_hmm(n=200,p=10,no_h=8,no_o=2,
                    prob_h_stay=0.8,prob_o_stay=0.8)

## compressing parameters, training model, making prediction
comp_train_pred(
	################## specify data information  ################
	test_x=data_hmm$X[1:100,],train_x=data_hmm$X[-(1:100),],
	train_y=data_hmm$y[-(1:100)],no_cls=2,nos_fth=rep(2,10),
	################## specify for compression ##################
	is_sequence=1,order=4,ptn_file=".ptn_file.log",
	new_compression=1,do_comp=1,
	###################### specify for priors  ##################
	alpha=1,log_sigma_widths=c(),log_sigma_modes=c(),
	################# specify for mc sampling ###################
	mc_file=".mc_file.log",start_over=TRUE,iters_mc=100,
	iters_bt=1,iters_sgm=2,w_bt=5,w_sgm=1,
	m_bt=20,m_sgm=20,ini_log_sigmas=c(),
	################## specify for prediction ###################
	pred_file=".pred_file.csv",iter_b = 10,forward = 1,
	iters_pred = 90)

## display summary information about compression
display_ptn(ptn_file=".ptn_file.log")

## display the pattern information for group 1 and group 2
display_ptn(ptn_file=".ptn_file.log",gid=c(1,2))

## display the general information of Markov chain sampling
display_mc(mc_file=".mc_file.log")

## read Markov chain values of log-likelihood from  ".mc_file.log"
read_mc(mc_file=".mc_file.log",group="lprobs",ix=0,
        iter_b=0,forward=1,n=100)

## particularly read `betas' by specifying the group and class id
read_betas(mc_file=".mc_file.log",ix_g=5,ix_cls=2,
           iter_b=0,forward=1,n=100)

## display the information on the pattern related to a `beta'
display_a_beta(mc_file=".mc_file.log",
               ptn_file=".ptn_file.log",id_beta=5)

## calculate the medians of samples of each 'beta'
calc_medians_betas(mc_file=".mc_file.log",iter_b=10,forward=1,n=90)

## evaluate prediction with true values of the response
evaluate_prediction(
       test_y=data_hmm$y[1:100],
       pred_result=read.csv(".pred_file.csv"),
       file_eval_details="eval_details")

#perform cross-validation with training data only
cv_comp_train_pred(
	################## specify data information  ################
	no_fold=2,train_x=data_hmm$X[-(1:100),],
	train_y=data_hmm$y[-(1:100)],no_cls=2,nos_fth=rep(2,10),
	################## specify for compression ##################
	is_sequence=1,order=4,ptn_file=".ptn_file.log",
	new_compression=1,do_comp=1,
	###################### specify for priors  ##################
	alpha=1,log_sigma_widths=c(),log_sigma_modes=c(),
	################# specify for mc sampling ###################
	mc_file=".mc_file.log",iters_mc=100,
	iters_bt=1,iters_sgm=2,w_bt=5,w_sgm=1,
	m_bt=20,m_sgm=20,ini_log_sigmas=c(),
	################## specify for prediction ###################
	pred_file=".pred_file.csv",iter_b = 10,forward = 1,
	iters_pred = 90)

#####################################################################
#####################################################################

## generating a classification data
data_class <- gen_bin_ho(n=400,p=3,order=3,alpha=1,
                 sigmas=c(0.3,0.2,0.1),nos_features=c(4,4,4),beta0=0)

## compressing parameters, training model, making prediction
comp_train_pred(
	################## specify data information  ################
	test_x=data_class$X[1:100,],train_x=data_class$X[-(1:100),],
	train_y=data_class$y[-(1:100)],no_cls=2,nos_fth=rep(4,3),
	################## specify for compression ##################
	is_sequence=0,order=3,ptn_file=".ptn_file.log",
	new_compression=1,do_comp=1,
	###################### specify for priors  ##################
	alpha=1,log_sigma_widths=c(),log_sigma_modes=c(),
	################# specify for mc sampling ###################
	mc_file=".mc_file.log",start_over=TRUE,iters_mc=500,
	iters_bt=1,iters_sgm=5,w_bt=5,w_sgm=0.5,
	m_bt=20,m_sgm=20,ini_log_sigmas=c(),
	################## specify for prediction ###################
	pred_file=".pred_file.csv",iter_b = 100,forward = 1,
	iters_pred = 400)

## display summary information about compression
display_ptn(ptn_file=".ptn_file.log")

## display the pattern information for group 1 and group 2
display_ptn(ptn_file=".ptn_file.log",gid=c(1,2))

## display the general information of Markov chain sampling
display_mc(mc_file=".mc_file.log")

## read Markov chain values of log-likelihood from ".mc_file.log"
read_mc(mc_file=".mc_file.log",group="lprobs",ix=0,
        iter_b=0,forward=1,n=100)

## particularly read `betas' by specifying the group and class id
read_betas(mc_file=".mc_file.log",ix_g=5,ix_cls=2,
           iter_b=0,forward=1,n=500)

## display the information on the pattern related to a `beta'
display_a_beta(mc_file=".mc_file.log",ptn_file=".ptn_file.log",
               id_beta=5)

## calculate the medians of samples of each 'beta'
calc_medians_betas(mc_file=".mc_file.log",iter_b=100,forward=1,n=400)

## evaluate prediction with true values of the response
evaluate_prediction(
       test_y=data_class$y[1:100],
       pred_result=read.csv(".pred_file.csv"),
       file_eval_details="eval_details")
}

\keyword{classif}
