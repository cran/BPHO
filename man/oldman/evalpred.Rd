\name{evaluate_prediction}

\alias{evaluate_prediction}

\title{Evaluate the prediction result given by `comp_train_pred_seq'}

\description{
Evaluate the prediction result given by `comp_train_pred_seq'
}

\usage{
evaluate_prediction(y_true,pred_result,file_detail=c())
}

\arguments{

\item{y_true}{True values of response}

\item{pred_result}{The returned value from `comp_train_pred_seq'}

\item{file_detail}{The name of file to which the detailed evaluation result is
		   sent} 

}

\value{ 

\item{error.rate}{Fraction of wrong prediction}

\item{amll}{Average Minus Log Probabilities (Likelihood)}

}

\seealso{comp_train_pred_seq}

\author{

Longhai Li, Department of Statistics, University of Toronto,
\email{longhai@stat.utoronto.ca},
\url{http://www.utstat.toronto.edu/longhai}

}
