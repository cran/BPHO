\name{gendata}

\alias{gendata}
\alias{gen_hmm}
\alias{gen_bin_ho}
\alias{text_to_number}
\alias{gen_X}

\title{Functions for generating data sets}

\description{ \code{gen_hmm} generates sequences using hidden Markov models. \code{gen_bin_ho} generates general discrete data using logistic models, with high-order interactions considered; the response is binary. \code{text_to_number} converts an English text file into sequence of 1 (special symbols such as space, symbol),2 (vowl),3 (consonant).
}

\usage{
gen_hmm(n,p,no_h,no_o,prob_h_stay, prob_o_stay)
gen_bin_ho(n,p,order,alpha,sigmas,nos_features,beta0)
text_to_number(p,file)
gen_X(n,p,K)
}

\arguments{

\item{n}{number of cases.}

\item{p}{number of features, or length of sequence.}

\item{K}{number of possibilities for each feature.}

\item{no_h}{number of states of hidden Markov chain.}

\item{no_o}{number of states of output in hidden Markov model.}

\item{prob_h_stay}{In simulating the hidden Markov chain, a chain will stay in its previous state with probability \code{prob_h_stay}, and move to other states with some minor probabilities adding up to \code{1-prob_h_stay}.}

\item{prob_o_stay}{In simulating the output state of hidden Markov model, the "output" is equal to ("hidden state" mod \code{no_o})+1 with probability \code{prob_o_stay} and equally likely other states.}

\item{order}{the order of interactions considered in simulating data from general classification models.}

\item{alpha}{\code{alpha=2} indicates that Gaussian distributions are used to generate the ``beta"s and \code{alpha=1} indicates that Cauchy distributions are used.}

\item{sigmas}{hyperparameters in generating "beta"s, a vector of length \code{order}.}

\item{nos_features}{number of states for each feature, i.e., the number of possibilities for each feature. A vector of length \code{p}.}

\item{beta0}{intercept of linear function in generating classification data.}

\item{file}{name of the file containing text file, a character string.}
}

\value{

\item{X}{values of predictors, a matrix. Each row is a case. For sequence, the data for each case (a row) is placed in the reverse order of time. For example, sequence "x1,x2,x3" is represented with a row of \code{X}: x3,x2,x1. The values of predictor \code{X} are coded by 1,2,3,...,\code{nos_features}. The function \code{gen_X} generates only this matrix.}

\item{y}{values of the response, a vector, coded by 1,2,...}

\item{betas}{a matrix of two columns saving the values of ``betas" used in generating classification data. The first column is the absolute identity of this beta, and the 2nd column is the value. The total number of ``betas" is saved in \code{no_betas}.}

}

\seealso{\code{\link{comp_train_pred}}}

\examples{
data_hmm <- gen_hmm(100,10,8,2,0.8,0.8)
data_bin_ho <- gen_bin_ho(100,3,2,1,c(5,2),c(3,3,3),0)
X <- gen_X(100,5,3)
}

\keyword{datagen}