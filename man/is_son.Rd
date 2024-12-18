% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{is_son}
\alias{is_son}
\title{Check if a node is son}
\usage{
is_son(tree, parents, son)
}
\arguments{
\item{tree}{Phylogenetic tree.}

\item{parents}{A vector with a set of parent nodes where to check whether the node belongs to.}

\item{son}{Node to be tested.}
}
\value{
A boolean value about whether the son is in parents.
}
\description{
A function to check wether a tip belongs to a specific node.
}
\author{
Moisès Bernabeu
}
