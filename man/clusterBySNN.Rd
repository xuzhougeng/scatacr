% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{clusterBySNN}
\alias{clusterBySNN}
\title{Find Cluster using SNN graph clustering methods}
\usage{
clusterBySNN(obj, min.group.size = 200, dims.use = 1:50,
  init.resolution = 0.8, n.start = 10, k.param = 20,
  prune.SNN = 1/15, annoy.metric = "euclidean", ...)
}
\arguments{
\item{obj}{Matrix  cell-by-feature matrix}

\item{min.group.size}{the mininum cell required to be cluster,
200 is required according to 2019, NBT.
Set it to NULL if the peaks are called}

\item{dims.use}{Dimensions of reduction to use as input}

\item{init.resolution}{Value of the resolution parameter,
use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities}

\item{n.start}{Number of random starts}

\item{k.param}{Defines k for the k-nearest neighbor algorithm}

\item{prune.SNN}{Sets the cutoff for acceptable Jaccard index
when computing the neighborhood overlap for the
SNN construction. Any edges with values less than
or equal to this will be set to 0 and removed from
the SNN graph. Essentially sets the strigency of pruning
(0 — no pruning, 1 — prune everything).}

\item{annoy.metric}{Distance metric for annoy. Options include:
euclidean, cosine, manhattan, and hamming}

\item{...}{Arguments passed to Seurat's
FindNeighbors and FindClusters}
}
\value{
named vector
}
\description{
Find Cluster using SNN graph clustering methods
}
