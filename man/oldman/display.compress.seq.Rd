\name{display_compress_seq}

\alias{display_compress_seq}

\title{Display Compression Information of `comp_train_pred_seq'}

\description{
This function displays the compression information saved 
by `comp_train_pred_seq'
}

\usage{ display_compress_seq(file, gids) }

\arguments{

\item{file}{A character sequence specifying the name of the file containing  the
	    compression information.}

\item{gids}{A vector giving the ID of the groups whose information 
            you want to disaplay. It can be NULL, meaning you want to display 
	    only the general information.}
}

\value{ 

It returns a vector of 4 elements: the 1st one is the number of training cases,
the 2nd one is `order' considered, which is also length of predictor sequence,
the 3rd one is the number of groups of parameters after compression, the 4th 
one is the number of original patterns expressed in training data

}

\seealso{comp_train_pred_seq}

\examples{
#display only general inforamtion
\dontrun{ display_compress_seq("ptn.log",c()) }
#dispaly information for group 0,1,2
\dontrun{ display_compress_seq("ptn.log",c(0,1,2))}


}


\author{

Longhai Li, Department of Statistics, University of Toronto,
\email{longhai@stat.utoronto.ca},
\url{http://www.utstat.toronto.edu/longhai}

}
