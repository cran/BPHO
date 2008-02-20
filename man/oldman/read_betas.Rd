\name{read_betas_g}

\alias{read_betas_g}

\title{Read `beta' values from Markov chain log file saved 
by `comp_train_pred_seq'}

\description{
This function reads `beta' values saved by `comp_train_pred_seq'
}

\usage{
read_betas_g(mc_file,ix_g,ix_cls,no_cls,iter_b,forward,n,quiet=1)
}

\arguments{

\item{mc_file}{Name of the file with Markov chain iterations}

\item{ix_g}{Index of group}

\item{ix_cls}{Index of Class}

\item{no_cls}{Total number of Classes for response}

\item{iter_b}{The first iteration one like to read}

\item{forward}{The number of iteractions skipped}

\item{n}{The number of iterations one like to read}

\item{quiet}{If FALSE print out some relevant information}
}

\value{ 

It returns the values of `beta' one requests with a vector of length `n'.
}

\seealso{comp_train_pred_seq,read_mc_g}

\author{

Longhai Li, Department of Statistics, University of Toronto,
\email{longhai@stat.utoronto.ca},
\url{http://www.utstat.toronto.edu/longhai}

}
