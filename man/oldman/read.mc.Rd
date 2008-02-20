\name{read_mc_g}

\alias{read_mc_g}

\title{Read Markov chain values from log file saved by `comp_train_pred_seq'}

\description{
This function reads Markov chain values saved by `comp_train_pred_seq'
}

\usage{
read_mc_g(group=c("lprobs","sigmas","betas"), 
          ix, mc_file,iter_b,forward,n,quiet)
}

\arguments{

\item{group}{A character sequence giving the group name of values. It can be 
one of ``lprobs'',``sigmas'',``betas''. 

Group ``lprobs'' includes, in order, log likelihood (ix=0), log prior
distribution of ``betas'' (ix=1), log prior distribution  of ``sigmas'' (ix=2),
log posterior distribution (i.e. sum of former three term)(ix=3).  

Group ``sigmas'' includes the values of hyperparameters associated with each
order, starting from order 0 (ix=0).  

Group ``betas'' includes the values of ``betas''. To read value of
betas, please use function `read_betas_g', which is a wrapper of this function
for only ``betas''.}

\item{ix}{index of parameters inside group. Starting from 0.}

\item{mc_file}{Name of the file with Markov chain iterations}

\item{iter_b}{The first iteration one like to read}

\item{forward}{The number of iteractions skipped}

\item{n}{The number of iterations one like to read}

\item{quiet}{If FALSE print out some relevant information}
}

\value{ 

It returns the values of parameters one requests with a vector of length `n'.
}

\seealso{comp_train_pred_seq,read_betas_g}

\examples{
\dontrun{ read_mc_g("lprobs",0,"mc.log",10,5,400,1) }

}


\author{

Longhai Li, Department of Statistics, University of Toronto,
\email{longhai@stat.utoronto.ca},
\url{http://www.utstat.toronto.edu/longhai}

}
