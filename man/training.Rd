\name{training}

\alias{training}
\alias{display_mc}
\alias{read_mc}
\alias{read_betas}
\alias{display_a_beta}
\alias{calc_medians_betas}
\title{Functions related to training the models with MCMC methods}

\description{
The models are trained with Markov chain Monte Carlo (MCMC) methods. Slice sampling is used to update `beta's, the regression coefficients for groups, Metropolis sampling with Normal proposal distribution is used to update `log(sigma)'s, where `sigma' is the scale parameter of the prior for `beta's.

The function \code{training} carries out the Markov chain sampling, saving the Markov chain in a binary file \code{mc_file}.

The function \code{display_mc} displays the summary information in the file \code{mc_file}.

The function \code{read_mc} reads the Markov chain samples from the file \code{mc_file} at given iterations.

The function \code{read_betas} is based on the function \code{read_mc}. It specifically reads the `beta' for given group and class identities.

The function \code{display_a_beta} displays both the pattern information for the group associated with the `beta' specified by \code{id_beta}, and also return the full Markov chain samples of this `beta'.

The function \code{calc_medians_betas} returns the medians of the Markov chain samples for all `beta's at specified iterations. This function is used to investigate useful interaction patterns. An interaction pattern with large absolute medians is highly suspected to be an important pattern for the response. 
}

\usage{
display_mc(mc_file)
read_mc(group,ix, mc_file,iter_b,forward,n,quiet=1)
read_betas(mc_file,ix_g,ix_cls,iter_b,forward,n,quiet=1)
display_a_beta(mc_file,ptn_file, id_beta)
training(mc_file,ptn_file, train_y,no_cls,
         alpha,sigma_precisions,log_sigma_modes,ini_log_sigmas,
         iters_mc,iters_bt,iters_sgm,
         w_bt,w_sgm,m_bt)
}

\arguments{

\item{mc_file}{A character string, the name of the binary file to which Markov chain is written.}

\item{group}{A character string giving the group name of values. It can be one of 'lprobs',lsigmas','betas', 'eval'. 

Group 'lprobs' contains the values of log likelihood (for \code{ix=0}), log prior of 'betas' (\code{ix=1}), log prior of 'lsigmas' (\code{ix=2}), log posterior (\code{ix=3}).  

Group 'lsigmas' contains the values of hyperparameters 'log(sigma)' for each order, starting from order 0 (\code{ix=0}).  

Group 'betas' contains the values of 'betas'. \code{ix} indicates the index of `beta'. The `beta's in each iteration is placed as follows: the \code{no_cls} for group `i' is followed by the next \code{no_cls} for group 'i+1'. The smallest index is 0.

Group `evals' contains the average times of evaluating the posterior distribution in updating a `beta' using slice sampling (\code{ix=0}), and the average rejection rate of updating `log(sigma)' with Metropolis sampling. 
}

\item{ix}{index of parameters inside each group, as discussed for \code{group}.}

\item{ix_g}{index of group, starting from 0.}

\item{ix_cls}{index of class, \bold{starting from 1}.}

\item{id_beta}{index of `beta', starting from 0.}

\item{iter_b, forward, n}{Starting from \code{iter_b}, one of every \code{forward} Markov chain samples, with the number of total samples being <= \code{n} and the maximum usable in the file \code{mc_file}, is read. }

\item{quiet}{\code{quiet=1} suppresses the messages printed during reading the file \code{mf_file}.}

\item{no_cls}{the number of possibilities (classes) of the response, default to the maximum value in \code{train_y}.}

\item{alpha}{\code{alpha=1} indicates that Cauchy prior is used, \code{alpha=2} indicates that Gaussian prior is used.}

\item{sigma_precisions, log_sigma_modes}{two vectors of length \code{order+1}, whose meanings are interpreted as follows: the Gaussian distribution with location \code{log_sigma_modes[o]} and standard deviation \code{1/sigmas_precisions[o]} is the prior distribution for logarithm of `sigmas[o]', which is the hyperparameter (scale parameter of Gaussian distribution or Cauchy distribution) for the regression coefficients (i.e. `beta's) associated with the interactions of order `o'.}

\item{ptn_file}{a character string, the name of the binary file to which the compression result is written. The method of writing to and reading from \code{ptn_file} can be found from the documentation for \link{compression}.}

\item{iters_mc,iters_bt,iters_sgm}{\code{iters_mc} iterations of super-transition will be run. Each super-transition consists of \code{iters_bt} iterations of updating `beta's, and for each updating of `beta's, the hyperparameters `log(sigma)'s are updated \code{iters_sgm} times. When \code{iters_mc=0}, no Markov chain sampling will be run and other arguments related to Markov chain sampling take no effect.}

\item{w_bt,w_sgm, m_bt}{\code{w_bt} is the amount of stepping-out in updating `beta' with slice sampling. \code{w_sgm} is the standard deviation of Gaussian distribution, which is used as the proposal distribution for updating `log(sigma)' with Metropolis sampling. \code{m_bt} is the maximum number of stepping-out in slice sampling for updating `beta'.}

\item{ini_log_sigmas}{Initial values of `log(sigma)', default to \code{log_sigma_mode}.}

\item{train_y}{Discrete response of training data. Assumed to be coded with 1,2,... \code{no_cls}.}

}

\value{
\item{}{The function \code{display_mc} returns a vector with names \code{#iters},\code{#class},\code{#groups},\code{order},\code{alpha}.}

\item{}{The function \code{read_mc} returns the Markov chain samples for a variable at specified iterations.}

\item{}{The function \code{read_betas} returns the Markov chain samples for a `beta' at specified iterations.}

\item{}{The function \code{display_a_beta} displays the pattern information for the specified pattern, and also returns the Markov chain samples of this `beta'.}

\item{}{The function \code{calc_medians_betas} returns the medians of Markov chain samples of all `beta's at given iterations.}  

\item{}{The function \code{training} returns no value. Instead, the Markov chain is written to the file \code{mc_file}.}
}

\examples{
## examples are given in \link{begin.BPHO}.
}

\seealso{\link{begin.BPHO},\link{compression},\link{prediction}}

\keyword{classif}