% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessors.R
\name{extract_trace_table}
\alias{extract_trace_table}
\title{Extract traces}
\usage{
extract_trace_table(fragments_trace_list)
}
\arguments{
\item{fragments_trace_list}{a list of fragments objects}
}
\value{
returns a dataframe of the raw trace data. Each row representing a single scan.
}
\description{
Extract the raw trace from a list of fragments objects
}
\examples{

test_ladders <- find_ladders(instability::cell_line_fsa_list[1])

extracted_traces <- extract_trace_table(test_ladders)

}
