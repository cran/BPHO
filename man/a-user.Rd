\name{comp_train_pred}

\alias{begin.BPHO}

\alias{comp_train_pred}

\alias{cv_comp_train_pred}

\title{User-level functions for compressing parameters, training the models with MCMC, and making predictions for test cases}

\description{The function \code{comp_train_pred} can be used for all of three tasks: compressing parameter, training the models with MCMC, and making prediction for test cases. When \bold{\code{new_comp=1}}, it compresses parameters based on training cases  and the information about parameter compression is written to the binary file \code{ptn_file}. When \bold{\code{new_comp=0}}, it uses the existing \code{ptn_file}. When \bold{\code{iters_mc > 0}}, it trains the models with Markov chain Monte Carlo  and the Markov chain iterations are written to the binary file \code{mc_file}. The methods of writing to and reading from the files \code{ptn_file} and \code{mc_file} can be found from the documentations \code{\link{compression}} and \code{\link{training}}. When \bold{\code{iters_pred > 0}} and \bold{\code{test_x}} isn't empty, it predicts the responses of test cases and the result is written to the file \code{pred_file} and also returned as a value of this function.

The function \code{cv_comp_train_pred} is a short-cut function for performing cross-validation with the function \code{comp_train_pred}.

The argument \code{is_sequence=1} indicates that a sequence prediction model is fitted to the data, and \code{is_sequence=1} indicates that a general classification model based on discrete predictor variables is fitted.
}

\usage{
comp_train_pred(
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

cv_comp_train_pred(
    ###################### Specify data,order,no_fold ##################
    no_fold=10,train_x,train_y,no_cls=c(),nos_features=c(),
    ###################### specify for compressing #####################
    is_sequence=1,order,ptn_file=c(),comp=1,no_cases_ign = 0,
    ###################### specify for priors  #########################
    alpha=1,log_sigma_widths=c(),log_sigma_modes=c(),
    ###################### specify for mc sampling #####################
    mc_file=c(),iters_mc=200,iters_bt=10,iters_sgm=50,
    w_bt=5,w_sgm=1,m_bt=20,m_sgm=20,ini_log_sigmas=c(),
    ###################### specify for prediction ######################
    pred_file=c(),iter_b=100,forward=2,samplesize=50)
}

\arguments{

\item{test_x}{Discrete features (also called inputs,covariates,independent variables, explanatory variables, predictor variables) of test data on which the  predictions are based.  The row is subject and the columns are inputs, which are coded with 1,2,..., with 0 reserved to represent that this input is not considered in a pattern. When the sequence prediction models are fitted, it is assumed that the first column is the state closest to the response. For example, a sequence `x1,x2,x3,x4' is saved in \code{test_x} as `x4,x3,x2,x1', for predicting the response `x5'. It could be empty if there isn't prediction needed.}

\item{train_x}{Discrete features of training data of the same format as \code{test_x}.}

\item{train_y}{Discrete response of training data, a vector with length equal to the row of \code{train_x}. Assumed to be coded with 1,2,... \code{no_cls} .}

\item{no_cls}{the number of possibilities (classes) of the response, default to the maximum value in \code{train_y}.}

\item{nos_features}{a vector, with each element storing the number of possibilities (classes) of each feature, default to the maximum value of each feature.}

\item{is_sequence}{\code{is_sequence=1} indicates that sequence prediction models are fitted to the data, and \code{is_sequence=0} indicates that general classification models based on discrete predictor variables are fitted.
}

\item{no_fold}{Number of folders in cross-validation.}

\item{order}{the order of interactions considered, default to the total number of features, i.e. \code{ncol(train_x)}.}

\item{ptn_file}{a character string, the name of the binary file to which the compression result is written. The method of writing to and reading from \code{ptn_file} can be found from the documentation for \code{\link{compression}}.}

\item{new_comp}{\code{new_comp=1} indicates removing the old file \code{ptn_file} if it exists and doing the compression once again. \code{new_comp=0} indicates using the old file \code{ptn_file} without doing compression once again. Note that when \code{new_comp=0}, the specification related to training cases does not take effect.}

\item{no_cases_ign}{When the number of training cases for a pattern is no more than \code{no_cases_ign}, this pattern will be ignored, default to 0.}


\item{comp}{\code{comp=1} indicates doing compression, and \code{comp=0} indicates using original parametrization. This is used only to make comparison. In practice, we definitely recommend using our compression technique to reduce the number of parameters.}

\item{alpha}{\code{alpha=1} indicates that Cauchy prior is used, \code{alpha=2} indicates that Gaussian prior is used.}

\item{log_sigma_widths, log_sigma_modes}{two vectors of length \code{order+1}, which are interpreted as follows: the Gaussian distribution with location \code{log_sigma_modes[o]} and standard deviation \code{log_sigma_widths[o]} is the prior for `log(sigmas[o])', which is the hyperparameter (width parameter of Gaussian distribution or Cauchy distribution) for the regression coefficients (i.e. `beta's) associated with the interactions of order `o'. If they are set to empty, the program will specify them automatically. By default, \code{log_sigma_widths} are all equal to 1, and \code{log_sigma_modes} starts at 2 for order 1 if Gaussian prior is used, and -1 if Cauchy prior is used, then decreases by 0.5 per order.}

\item{mc_file}{A character string, the name of the binary file to which Markov chain is written. The method of writing to and reading from \code{mc_file} can be found from the documentation for \code{\link{training}}.}

\item{start_over}{\code{start_over=TRUE} indicates that the existing file \code{mc_file} is deleted before a Markov chain sampling starts, otherwise the Markov chain will continue from the last iteration stored in \code{mc_file}. }

\item{iters_mc,iters_bt,iters_sgm}{\code{iters_mc} iterations of super-transition will be run. Each super-transition consists of \code{iters_bt} iterations of updating `beta's, and for each updating of `beta's, the hyperparameters `log(sigma)'s are updated \code{iters_sgm} times. When \code{iters_mc=0}, no Markov chain sampling will be run and other arguments related to Markov chain sampling take no effect.}

\item{w_bt, m_bt, w_sgm, m_sgm}{\code{w_bt} is the amount of stepping-out in updating `beta' with slice sampling,  \code{m_bt} is the maximum number of stepping-out in slice sampling for updating `beta'. \code{w_sgm} and \code{m_sgm} are intepreted similarly for sampling for `log(sigma)'.}

\item{ini_log_sigmas}{Initial values of `log(sigma)', default to \code{log_sigma_modes}.}

\item{pred_file}{A character string, the name of the file to which the prediction result is written. If \code{pred_file=c()}, the prediction result is printed out on screen (or sent to standard output).}

\item{iter_b, forward, samplesize}{Starting from \code{iter_b}, one of every \code{forward} Markov chain samples, with the number of total samples being <= \code{samplesize} and the maximum usable in the file \code{mc_file}, is used to make prediction. }

}

\value{

\item{}{The function \code{comp_train_pred} returns the following values:}

\item{times}{The time in second for, as this order, compressing parameters,
training the model, predicting for test cases}

\item{pred_result}{a data frame with first \code{no_cls} columns being the predictive probability and the next column being the predicted response value is returned.}

\item{files}{Three character strings: the 1st is the name of the file storing compression information, the 2nd is the name of the file storing Markov chain, and the 3rd one is the name of the file containing the detailed prediction result, i.e., \code{pred_result}.}

\item{}{The function \code{cv_comp_train_pred} returns the following additional values:}

\item{eval_details}{a data frame. The first column is the true response, the second is the guessed value by taking the label of class with largest predictive probability, the third is indicator whether a wrong decision is make, the last column is the predictive probability at the true class.}

\item{error_rate}{the proportion of wrong prediction.}

\item{amll}{the average of minus log probabilities at true class, i.e. the average of the logarithms of the last column of \code{eval_details}.}


}

\seealso{\code{\link{gendata}},\code{\link{compression}},\code{\link{training}},\code{\link{prediction}}}

\author{Longhai Li, \url{http://math.usask.ca/~longhai}}

\references{\url{http://math.usask.ca/~longhai/doc/seqpred/seqpred.abstract.html}}

\examples{
#####################################################################
########  These are demonstrations of using the whole package
#####################################################################

#####################################################################
########  Apply to Sequence Prediction Models
#####################################################################


## generate data from a hidden Markov model
data_hmm <- gen_hmm(n = 200, p = 10, no_h = 4, no_o = 2, 
                    prob_h_stay = 0.8, prob_o_stay = 0.8)

## compress parameters (transforming features in training data.)
compress (features = data_hmm$X[-(1:100),], is_sequence = 1, order = 4,
          ptn_file = ".ptn_file.log")

## print summary information of '.ptn_file.log'
display_ptn (ptn_file = ".ptn_file.log")

## draw 50 samples with slice sampling, saved in '.mc_file.log'
training (train_y = data_hmm$y[-(1:100)], no_cls = 2, 
          mc_file = ".mc_file.log", ptn_file = ".ptn_file.log", 
          iters_mc = 50)

## draw 50 more samples from the last iteration in '.mc_file.log'
training (train_y = data_hmm$y[-(1:100)], no_cls = 2, 
          mc_file = ".mc_file.log", ptn_file = ".ptn_file.log", 
          iters_mc = 50, start_over = FALSE)

## display general information of Markov chain stored in '.mc_file.log'
summary_mc (mc_file = ".mc_file.log")

## find medians of Markov chain samples of all betas
med_betas <- medians_betas (mc_file = ".mc_file.log", 
             iter_b = 50, forward = 1, samplesize = 50)

## take the id of group with highest absolute median
g_impt <- med_betas [1, "groupid"] 

## particularly read `betas' by specifying the group and class id
plot ( read_betas (ix_g = g_impt, ix_cls = 2, mc_file = ".mc_file.log", 
       iter_b = 50, forward = 1, samplesize = 50) )

## get information about pattern groups
display_ptn (".ptn_file.log", ix_g = g_impt )

## read Markov chain values of log-likelihood from  '.mc_file.log'
plot ( read_mc (group = "lprobs", ix = 1, mc_file = ".mc_file.log", 
       iter_b = 50, forward = 1, samplesize = 50) )

## predict for test cases
pred_probs <- predict_bpho (test_x = data_hmm$X[1:100,],
              mc_file = ".mc_file.log", ptn_file = ".ptn_file.log",
              iter_b = 50, forward = 1, samplesize = 50 )
              
## evaluate predictions with true value of response
evaluate_prediction (data_hmm$y[1:100], pred_probs)

#####################################################################
########  Apply to General Classification Models
#####################################################################

## generating a classification data
data_class <- gen_bin_ho(n = 200, p = 3, order = 3, alpha = 1,
              sigmas = c(0.3,0.2,0.1), nos_features = rep(4,3), beta0 = 0)

## compress parameters (transforming features in training data.)
compress (features = data_class$X[-(1:100),], is_sequence = 0, order = 20,
          ptn_file = ".ptn_file.log")

## print summary information of '.ptn_file.log'
display_ptn (ptn_file = ".ptn_file.log")

## draw 50 samples with slice sampling, saved in '.mc_file.log'
training (train_y = data_class$y[-(1:100)], no_cls = 2, 
          mc_file = ".mc_file.log", ptn_file = ".ptn_file.log", 
          iters_mc = 50)

## draw 50 more samples from the last iteration in '.mc_file.log'
training (train_y = data_class$y[-(1:100)], no_cls = 2, 
          mc_file = ".mc_file.log", ptn_file = ".ptn_file.log", 
          iters_mc = 50, start_over = FALSE)

## display general information of Markov chain stored in '.mc_file.log'
summary_mc (mc_file = ".mc_file.log")

## find medians of Markov chain samples of all betas
med_betas <- medians_betas (mc_file = ".mc_file.log", 
             iter_b = 50, forward = 1, samplesize = 50)

## take the id of group with highest absolute median
g_impt <- med_betas [1, "groupid"] 

## particularly read `betas' by specifying the group and class id
plot ( read_betas (ix_g = g_impt, ix_cls = 2, mc_file = ".mc_file.log", 
       iter_b = 50, forward = 1, samplesize = 50) )

## get information about pattern groups
display_ptn (".ptn_file.log", ix_g = g_impt )

## read Markov chain values of log-likelihood from  '.mc_file.log'
plot ( read_mc (group = "lprobs", ix = 1, mc_file = ".mc_file.log", 
       iter_b = 50, forward = 1, samplesize = 50) )

## predict for test cases
pred_probs <- predict_bpho (test_x = data_class$X[1:100,],
              mc_file = ".mc_file.log", ptn_file = ".ptn_file.log",
              iter_b = 50, forward = 1, samplesize = 50 )
              
## evaluate predictions with true value of response
evaluate_prediction (data_class$y[1:100], pred_probs)

#####################################################################
########  Demonstrations of using a single function 'comp_train_pred'
#####################################################################

## generating a classification data
data_class <- gen_bin_ho(n = 200, p = 3, order = 3, alpha = 1,
              sigmas = c(0.3,0.2,0.1), nos_features = rep(4,3), beta0 = 0)

## carry out compression, training, and prediction with a function 
comp_train_pred (
    ################## specify data information  #######################
    test_x = data_class$X[1:100,], train_x = data_class$X[-(1:100),],
    train_y = data_class$y[1:100], no_cls = 2, nos_features = rep(4,3),
    ################## specify for compression #########################
    is_sequence = 0, order = 3, ptn_file=".ptn.log", new_comp = 1, comp = 1,
    ################## specify for priors  #############################
    alpha = 1, log_sigma_widths = c(), log_sigma_modes = c(),
    ################## specify for mc sampling #########################
    mc_file = ".mc.log", iters_mc = 200, iters_bt = 10,
    iters_sgm = 50, w_bt = 5, m_bt = 20, w_sgm = 1, m_sgm = 20,
    ini_log_sigmas = c(), start_over = TRUE, 
    ################## specify for prediction ##########################
    pred_file = "pred_file.csv", iter_b = 100, forward = 2, samplesize = 50)

## evaluate predictions with true value of response
evaluate_prediction (data_class$y[1:100], read.csv("pred_file.csv") )

#####################################################################
#####  Demonstrations of using a single function 'cv_comp_train_pred'
#####################################################################

## generating a classification data
data_class <- gen_bin_ho(n = 200, p = 3, order = 3, alpha = 1,
              sigmas = c(0.3,0.2,0.1), nos_features = rep(4,3), beta0 = 0)

## carry out cross-validation with data set data_class
cv_comp_train_pred (
    ################## Specify data,order,no_fold ######################
    no_fold = 2, train_x = data_class$X,train_y = data_class$y,
    no_cls = 2, nos_features = rep(4,3),
    ################## specify for compression #########################
    is_sequence = 0, order = 3, ptn_file=".ptn.log",  comp = 1,
    ################## specify for priors  #############################
    alpha = 1,log_sigma_widths = c(),log_sigma_modes = c(),
    ################## specify for mc sampling #########################
    mc_file = ".mc.log", iters_mc = 200, iters_bt = 10,
    iters_sgm = 50,w_bt = 5,m_bt = 20,w_sgm = 1,m_sgm = 20,
    ini_log_sigmas = c(), 
    ################## specify for prediction ##########################
    pred_file = "pred_file.csv", iter_b = 100, forward = 2, samplesize = 50)


}

\keyword{classif}
