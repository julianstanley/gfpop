% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gfpop.R
\name{itergfpop}
\alias{itergfpop}
\title{Graph-constrained functional pruning optimal partitioning iterated}
\usage{
itergfpop(data, mygraph, type = "mean", weights = NULL,
  iter.max = 100, D.init = 1)
}
\arguments{
\item{data}{vector of data to segment}

\item{mygraph}{dataframe of class graph to constrain the changepoint inference}

\item{type}{a string defining the cost model to use: "gauss", "variance", "poisson", "exp", "negbin"}

\item{weights}{vector of weights (positive numbers), same size as data}

\item{iter.max}{maximal number of iteration of the gfpop function}

\item{D.init}{initialisation of the number of segments}
}
\value{
a gfpop object = (changepoints, states, forced, parameters, globalCost, Dvect)
\describe{
\item{\code{changepoints}}{is the vector of changepoints (we give the last element of each segment)}
\item{\code{states}}{is the vector giving the state of each segment}
\item{\code{forced}}{is the vector specifying whether the constraints of the graph are active (=1) or not (=0)}
\item{\code{parameters}}{is the vector of successive parameters of each segment}
\item{\code{globalCost}}{is a number equal to the total loss: the minimal cost for the optimization problem with all penalty values excluded}
\item{\code{Dvect}}{is a vector of integers. The successive tested D in the Birgé Massart penalty until convergence}
 }
}
\description{
Functional pruning optimal partitioning with a graph structure to take into account constraints on consecutive segment parameters. This is an iterated version of the main gfpop function using a Birgé Massart like penalty
}
