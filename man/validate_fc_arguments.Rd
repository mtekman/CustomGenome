% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/custom_genome.R
\name{validate_fc_arguments}
\alias{validate_fc_arguments}
\title{Validate the featureCounts Ellipsis Arguments for RSubread's
featureCounts function.}
\usage{
validate_fc_arguments(...)
}
\arguments{
\item{...}{An ellipsis argument that is passed in Rsubreads's
\verb{featureCounts' function. Typical values are }GTF.featuretype'
(default: "exon"), `isPairedEnd' (default: TRUE).}
}
\value{
A list containing validated and/or modified ellipsis
arguments.
}
\description{
Convert ellipsis arguments into a list and validate
their value.
}
