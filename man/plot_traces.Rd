% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessors.R
\name{plot_traces}
\alias{plot_traces}
\title{Plot sample traces}
\usage{
plot_traces(
  fragments_list,
  show_peaks = TRUE,
  n_facet_col = 2,
  sample_subset = NULL,
  xlim = NULL,
  ylim = NULL
)
}
\arguments{
\item{fragments_list}{A list of fragments_repeats or fragments_trace objects containing fragment data.}

\item{show_peaks}{If peak data are available, TRUE will plot the peaks on top of the trace as dots.}

\item{n_facet_col}{A numeric value indicating the number of columns for faceting in the plot.}

\item{sample_subset}{A character vector of unique ids for a subset of samples to plot}

\item{xlim}{the x limits of the plot. A numeric vector of length two.}

\item{ylim}{the y limits of the plot. A numeric vector of length two.}
}
\value{
plot traces from fragments object
}
\description{
Plot the raw trace data
}
\details{
A plot of the raw signal by bp size. Red vertical line indicates the scan was
flagged as off-scale. This is in any channel, so use your best judgment to determine
if it's from the sample or ladder channel.
}
\examples{

file_list <- instability::cell_line_fsa_list

  test_ladders <- find_ladders(file_list,
                               ladder_channel = "DATA.105",
                               signal_channel = "DATA.1",
                               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                               hq_ladder = FALSE,
                               max_combinations = 2500000,
                               ladder_selection_window = 8)

 # Manually inspect the ladders
 plot_traces(test_ladders[1], n_facet_col = 1)


}