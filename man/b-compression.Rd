\name{compression}

\alias{compression}

\alias{compress}

\alias{display_ptn}

\title{Functions related to parameter compression}

\description{
The function \code{compress} groups the patterns in a way such that the interaction patterns in a group are expressed by the same training cases. In training the models with MCMC, we need to use only one parameter for each group, which represents the sum of all the parameters in this group. The original parameters are seemly compressed. A large amount of training time is saved by this compression techniques.

The result of this grouping is saved in a binary file in a way such that it can be retrieved as a linked list in C, with each node consisting of a description (an integer vector of fixed length) of the group of patterns and the indice (an integer vector of varying length, with 0 for the first training case) of training cases expressing this group of patterns. This file is needed to train the models with MCMC and to predict the responses of test cases using the function \code{\link{comp_train_pred}}.

The function \code{display_ptn} displays the summary information about this compression, such as the number of groups and total number of patterns expressed by the training cases. When \code{ix_g} is nonempty, it also displays the detailed information about the groups specified by \code{ix_g}, such as the pattern description and the indice of training cases associated with this group.

}

\usage{
compress (
    features, is_sequence=1, order=ncol(features), ptn_file, 
    nos_features=c(), comp=1, quiet=1, no_cases_ign=0 )
display_ptn ( ptn_file, ix_g=c() )
}

\arguments{

\item{features}{Discrete features (also called features,covariates,independent variables, explanatory variables, predictor variables) of training data on which the  predictions are based.  The row is subject and the columns are inputs, which are coded with 1,2,..., with 0 reserved to represent that this input is not considered in a pattern. When the sequence prediction models are fitted, it is assumed that the first column is the state closest to the response. For example, a sequence `x1,x2,x3,x4' is saved in \code{test_x} as `x4,x3,x2,x1', for predicting the response `x5'.}

\item{is_sequence}{\code{is_sequence=1} indicates that sequence prediction models are fitted to the data, and \code{is_sequence=0} indicates that general classification models based on discrete predictor variables are fitted.}

\item{order}{the order of interactions considered, default to total number of features, i.e. \code{ncol(features)}.}

\item{ptn_file}{a character string, the name of the binary file to which the compression result is written.}

\item{nos_features}{a vector, with each element storing the number of possibilities (classes) of each feature, default to the maximum value of each feature.}

\item{comp}{\code{comp=1} indicates doing compression, and \code{comp=0} indicates using original parametrization. This is used only to make comparison. In practice, I of course recommend using our compression technique to reduce the number of parameters.}

\item{quiet}{If \code{quiet=0}, some messages during compression are printed on screen for monitor the compression, if \code{quiet=1} the function works silently.}

\item{no_cases_ign}{When the number of training cases for a pattern is no more than \code{no_cases_ign}, this pattern will be ignored, default to 0.}


\item{ix_g}{an integer vector, containing the indice of groups whose information you want to display, with 1 for the first group.}


}

\value{

\item{}{The function \code{compress} returns no value. Instead, it saves the result of compression in the file \code{ptn_file}.}

\item{}{The function \code{display_ptn} returns a vector of 6 numbers.  Their meanings are as follows: \code{is.sequence} -- indictor whether a sequence model is fitted,\code{order} -- the maximum order of interactions considered, \code{#groups} -- the number of groups found, \code{#patterns} -- the number of interaction patterns expressed by the training cases,\code{#cases} -- the number of training cases,\cr \code{#features} -- the number of features.}

\item{}{When \code{ix_g} is nonempty, it also displays the details about the queried groups. The information printed on screen for each group is read as follows. Under \bold{\code{superpatterns}}, it displays a compact description of the pattern group, which is in a special format defined in the references associated with this software. Under \bold{\code{expression}}, it displays the indice of training cases that express this group of patterns. Under \bold{\code{sigmas}}, it displays the number of patterns with a certain order, starting from order 0. This information is needed to compute the width parameter of the regression coeficient associated with this group from the values of hyperparameters `sigma's.}

}

\examples{

## generate features
features <- gen_X(50,5,2)

## compressing the parameter based on 'features'
compress ( features, is_sequence=1, ptn_file=".ptn_file.log")

## display the summary information in the file ".ptn_file.log"
display_ptn(".ptn_file.log")

## display the information for group #2 and #3
display_ptn(".ptn_file.log",ix_g=c(2,3))

}

\seealso{\code{\link{comp_train_pred}},  \code{\link{training}},  \code{\link{prediction}}}

\keyword{classif}
