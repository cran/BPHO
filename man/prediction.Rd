\name{prediction}

\alias{prediction}
\alias{predict_bpho}
\alias{evaluate_prediction}
\alias{split_cauchy}

\title{Functions related to prediction}

\description{

The function \code{predict_bpho} predicts the response of test cases.

The function \code{evaluate_prediction} evaluates the performance of the prediction in terms of average minus log probabilities and error rate. The function \code{split_cauchy} draws samples from a Cauchy distribution of two variables constraint to that their sum is fixed.

}

\usage{
predict_bpho(test_x,no_cls,mc_file,ptn_file,iter_b,forward,
             iters_pred)
evaluate_prediction(test_y,pred_result,file_eval_details=c())
split_cauchy(n,s, sigma1,sigmasum,debug=1)
}

\arguments{

\item{test_x}{Discrete features (also called inputs,covariates,independent variables, explanatory variables, predictor variables) of test data on which the  predictions are based.  The row is subject and the columns are inputs, which are coded with 1,2,..., with 0 reserved to represent that this input is not considered in a pattern. When the sequence prediction models are fitted, it is assumed that the first column is the state closest to the response. For example, a sequence `x1,x2,x3,x4' is saved in \code{test_x} as `x4,x3,x2,x1', for predicting the response `x5'.}

\item{test_y}{Discrete responses of test data, a vector with length equal to the row of \code{test_x}. Assumed to be coded with 1,2,... \code{no_cls}.}

\item{no_cls}{the number of possibilities (classes) of the response.}

\item{ptn_file}{a character string, the name of the binary file to which the compression result is saved. The method of writing to and reading from \code{ptn_file} can be found from the documentation \code{\link{compression}}.}

\item{mc_file}{A character string, the name of the binary file to which Markov chain is written. The method of writing to and reading from \code{mc_file} can be found from the documentation \code{\link{training}}.}

\item{iter_b, forward, iters_pred}{Starting from \code{iter_b}, one of every \code{forward} Markov chain samples, with the number of total samples being <= \code{iters_pred} and the maximum usable in the file \code{mc_file}, is used to make prediction.}

\item{pred_result}{the value returned from the function \code{predict_bpho}.}

\item{file_eval_details}{the details of evaluation is sent to the file \code{file_eval_details}.}

\item{n}{number of samples one wishes to obtain.}

\item{s}{sum of two Cauchy random variables.}

\item{sigma1}{scale parameter for the first Cauchy random variable.}

\item{sigmasum}{the sum of scale parameters for two Cauchy random variables.}

\item{debug}{indicator whether you are debugging the C program.}

}

\value{

The function \code{predict_bpho} returns a data frame, with the first \code{no_cls} columns storing the predictive probabilities for each class, and the last column is the guess for the response by choosing the label of the class with largest predictive probability.

The function \code{evaluate_prediction} returns the following values:

\item{eval_details}{a data frame. The first column is the true response, the second is the guessed value by taking the label of class with largest predictive probability, the third is indicator whether a wrong decision is make, the last column is the predictive probability at the true class.}

\item{error_rate}{the proportion of wrong prediction.}

\item{amll}{the average of minus log probabilities at true class, i.e. the average of the logarithms of the last column of \code{eval_details}.}

The function \code{split_cauchy} returns a vector of \code{n} random numbers.

}

\examples{
## the function `predict_bpho' is demonstrated with the function
## `\link{comp_train_pred}' which calls `predict_bpho' inside.

## examples of 'evaluate_prediction' can be found from
## the documentation for \link{comp_train_pred}.


## testing the function split_cauchy
split_cauchy(100,10,1,5)

}

\seealso{\code{\link{comp_train_pred}},\code{\link{compression}},\code{\link{training}}}

\keyword{classif}