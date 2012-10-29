\name{training}

\alias{training}
\alias{summary_mc}
\alias{read_mc}
\alias{read_betas}
\alias{medians_betas}

\title{Functions related to Markov chain sampling}

\description{
The models are trained with Markov chain Monte Carlo (MCMC) methods. Slice sampling is used to update `beta's, the regression coefficients for groups, and  `log(sigma)', where `sigma' is the width parameter of the prior for `beta'.

The function \code{training} carries out the Markov chain sampling, saving the Markov chain samples in a binary file \code{mc_file}.

The function \code{summary_mc} displays the summary information in the file \code{mc_file}.

The function \code{read_mc} reads the Markov chain samples from the file \code{mc_file} at given iterations.

The function \code{read_betas} is based on the function \code{read_mc}. It specifically reads the `beta' for given group and class identities.

The function \code{medians_betas} returns the medians of the Markov chain samples for all `beta's at specified iterations. This function is for discovering important interaction patterns. An interaction pattern with large absolute medians is suspected to be an important pattern for predicting the response.
}

\usage{
training (
    ################## specify for data  ###############################
    train_y, no_cls, mc_file, ptn_file, 
    ################## specify for priors  #############################
    alpha=1, log_sigma_modes=c(), log_sigma_widths=c(), 
    ################## specify for slice sampling ######################
    iters_mc=200, iters_bt=10,  iters_sgm=50,
    w_bt=5, m_bt=20, w_sgm=1, m_sgm=20, ini_log_sigmas=c(), 
    start_over=TRUE )
summary_mc(mc_file)
read_mc(group,ix, mc_file, iter_b=1,forward=1,samplesize=c(),quiet=1)
read_betas(ix_g, ix_cls, mc_file, iter_b=1,forward=1,samplesize=c(),quiet=1)
medians_betas(mc_file, iter_b=1, forward=1, samplesize=c() )
}

\arguments{

\item{mc_file}{A character string, the name of the binary file to which Markov chain is written.}

\item{group}{
A character string giving the group name of values.

It can be one of 'lprobs',lsigmas','betas', 'evals'.

Group 'lprobs' contains: the values of log probabilities of data given the values of `beta's (identified by \code{ix=0}), the value of log prior of 'beta's given `sigma's (identified by \code{ix=1}), the value of log prior of 'log(sigma)'s (identified by \code{ix=2}), and the value of log posterior (identified by \code{ix=3}), which is the sum of the previous three values.

Group 'lsigmas' contains: the values of hyperparameters 'log(sigma)', with \code{ix} indicating the order, starting from 0.

Group 'betas' contains: the values of 'betas', with \code{ix} indicating the index of `beta'. The `beta's in each iteration is placed as that the \code{no_cls} values of `beta's for pattern group `i' are followed by the next \code{no_cls} values for pattern group 'i+1'. The smallest index is 0.

Group `evals' contains: the average times of evaluating the posterior distribution in updating each `beta' using slice sampling (identified by \code{ix=0}), and the average rejection rate of updating each `log(sigma)' with Metropolis sampling (identified by \code{ix=1}).
}

\item{ix}{index of parameters inside each group, as discussed for \code{group} above.}

\item{ix_g}{index of pattern group, starting from 1.}

\item{ix_cls}{index of class, ranging from \bold{1} to \bold{\code{no_cls}}.}

\item{iter_b, forward, samplesize}{Starting from \code{iter_b}, one of every \code{forward} Markov chain samples, with the number of total samples being <= \code{samplesize} and the maximum usable in the file \code{mc_file}, is read. }

\item{train_y}{Discrete response of training data. Assumed to be coded with 1,2,... \code{no_cls}.}

\item{no_cls}{the number of possibilities (classes) of the response, default to the maximum value in \code{train_y}.}

\item{alpha}{\code{alpha=1} indicates that Cauchy prior is used, \code{alpha=2} indicates that Gaussian prior is used.}

\item{log_sigma_widths, log_sigma_modes}{two vectors of length \code{order+1}, which are interpreted as follows: the Gaussian distribution with location \code{log_sigma_modes[o]} and standard deviation \code{log_sigma_widths[o]} is the prior for `log(sigmas[o])', which is the hyperparameter (width parameter of Gaussian distribution or Cauchy distribution) for the regression coefficients (i.e. `beta's) associated with the interactions of order `o'. If they are set to empty, the program will specify them automatically. By default, \code{log_sigma_widths} are all equal to 1, and \code{log_sigma_modes} starts at 2 for order 1 if Gaussian prior is used, and -1 if Cauchy prior is used, then decreases by 0.5 per order.}

\item{ptn_file}{a character string, the name of the binary file where the compression result is saved. The method of writing to and reading from \code{ptn_file} can be found from the documentation for \code{\link{compression}}.}

\item{iters_mc,iters_bt,iters_sgm}{\code{iters_mc} iterations of super-transition will be run. Each super-transition consists of \code{iters_bt} iterations of updating `beta's, and for each updating of `beta's, the hyperparameters `log(sigma)'s are updated \code{iters_sgm} times. When \code{iters_mc=0}, no Markov chain sampling will be run and other arguments related to Markov chain sampling take no effect.}

\item{w_bt, m_bt, w_sgm, m_sgm}{\code{w_bt} is the amount of stepping-out in updating `beta' with slice sampling,  \code{m_bt} is the maximum number of stepping-out in slice sampling for updating `beta'. \code{w_sgm} and \code{m_sgm} are intepreted similarly for sampling for `log(sigma)'.}

\item{ini_log_sigmas}{Initial values of `log(sigma)', default to \code{log_sigma_mode}.}

\item{quiet}{\code{quiet=1} suppresses the messages printed during reading the file \code{mc_file}.}

\item{start_over}{\code{start_over=TRUE} indicates that the existing file \code{mc_file} is deleted before a Markov chain sampling starts, otherwise the Markov chain will continue from the last iteration stored in \code{mc_file}. }


}

\value{
\item{}{The function \code{summary_mc} returns a vector with names as \cr \code{#iters},\code{#class},\code{#groups},\code{order},\code{alpha}.}

\item{}{The function \code{read_mc} returns the Markov chain samples for a variable at specified iterations.}

\item{}{The function \code{read_betas} returns the Markov chain samples for a `beta' at specified iterations.}


\item{}{The function \code{medians_betas} returns the medians of Markov chain samples of all `beta's at given iterations, displayed as a matrix, with rows for different pattern group, and columns for different value of response.}

\item{}{The function \code{training} returns no value. Instead, the Markov chain samples are written to the binary file \code{mc_file}.}
}

\examples{

## these functions are demonstrated in the section `\link{comp_train_pred}'.

}

\seealso{\code{\link{comp_train_pred}},\code{\link{compression}},\code{\link{prediction}}}

\keyword{classif}