% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessors.R
\name{find_alleles}
\alias{find_alleles}
\title{Find Alleles in Fragments List}
\usage{
find_alleles(
  fragments_list,
  number_of_peaks_to_return = 1,
  peak_region_size_gap_threshold = 6,
  peak_region_height_threshold_multiplier = 1
)
}
\arguments{
\item{fragments_list}{A list of fragment objects containing peak data.}

\item{number_of_peaks_to_return}{Number of main peaks to be returned for each fragment. Must either be 1 or 2 if a normal sized allele is also to be found.}

\item{peak_region_size_gap_threshold}{Gap threshold for identifying peak regions. The peak_region_size_gap_threshold is a parameter used to determine the maximum allowed gap between peak sizes within a peak region. Adjusting this parameter affects the size range of peaks that can be grouped together in a region. A smaller value makes it more stringent, while a larger value groups peaks with greater size differences, leading to broader peak regions that may encompass wider size ranges.}

\item{peak_region_height_threshold_multiplier}{Multiplier for the peak height threshold. The peak_region_height_threshold_multiplier parameter allows adjusting the threshold for identifying peak regions based on peak heights. Increasing this multiplier value will result in higher thresholds, making it more stringent to consider peaks as part of a peak region. Conversely, reducing the multiplier value will make the criteria less strict, potentially leading to more peaks being grouped into peak regions. It's important to note that this parameter's optimal value depends on the characteristics of the data and the specific analysis goals. Choosing an appropriate value for this parameter can help in accurately identifying meaningful peak regions in the data.}
}
\value{
A list of fragments with identified main alleles.
}
\description{
This function identifies main alleles within each fragment in a list of fragments.
}
\details{
This function finds the main alleles for each fragment in the list by identifying clusters of peaks ("peak regions")
with the highest signal intensities. This is based on the idea that PCR amplicons of repeats have broad peaks and PCR artififacts that help identifying the alleles.
The number of peaks to be returned, and the parameters for identifying peak regions can be customized.
It's important to note that both peak_region_height_threshold_multiplier and peak_region_size_gap_threshold influence the criteria for identifying peak regions, and finding the right balance between them is crucial.
}
\examples{
file_list <- instability::cell_line_fsa_list

test_ladders <- find_ladders(file_list[1])

fragments_list <- find_fragments(test_ladders,
  min_bp_size = 300
)


test_alleles <- find_alleles(
  fragments_list = fragments_list,
  number_of_peaks_to_return = 1,
  peak_region_size_gap_threshold = 6,
  peak_region_height_threshold_multiplier = 1
)
}
