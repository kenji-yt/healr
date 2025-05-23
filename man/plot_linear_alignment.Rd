% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_utils.R
\name{plot_linear_alignment}
\alias{plot_linear_alignment}
\title{Plot linear regression of copy number at anchors between each pair of progenitors}
\usage{
plot_linear_alignment(
  alignment,
  view_samples = FALSE,
  output_dir = FALSE,
  color = "blue4",
  alpha = 0.1,
  width = 0.2,
  height = 0.2
)
}
\arguments{
\item{alignment}{A heal alignment object created with get_heal_alignment().}

\item{view_samples}{A vector of sample names to plot (as character)('FALSE' by default).}

\item{output_dir}{The name of a directory to write all plots to. Will create one if nonexistent.}

\item{color}{The color of the points ("blue4" by default).}

\item{alpha}{The transparency of the points (0.1 by default).}

\item{width}{The jitter of points along the x axis.}

\item{height}{The jitter of points along the x axis.}
}
\value{
Nothing. Plots are shown and/or saved to output_dir.
}
\description{
Plot linear regression of copy number at anchors between each pair of progenitors
}
