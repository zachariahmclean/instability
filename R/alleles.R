################## general helper functions for finding peaks ###################
# maxima -----------------------------------------------------------------------
find_maxima <- function(x) {
  maxima <- logical(length(x))

  for (i in 2:(length(x) - 1)) {
    maxima[i] <- (x[i] > x[i - 1] && x[i] > x[i + 1])
  }

  maxima[1] <- (length(x) > 1) && (x[1] > x[2])
  maxima[length(x)] <- (length(x) > 1) && (x[length(x)] > x[length(x) - 1])

  return(maxima)
}



# find peak regions -------------------------------------------------------

find_peak_regions <- function(height, size, size_gap_threshold, region_height_threshold_multiplier) {
  peak_regions <- rep(NA_real_, length(height))
  mean_heights <- mean(height) * region_height_threshold_multiplier

  for (i in seq_along(height)) {
    if (height[i] < mean_heights || i == 1 || i == length(height)) {
      peak_regions[i] <- NA_real_
    } else if (height[i - 1] < mean_heights && height[i + 1] < mean_heights) {
      peak_regions[i] <- NA_real_
    } else {
      # check to see if peaks before it are within the size threshold
      current_size <- size[i]
      valid_lower_peaks <- which(size < current_size & size > current_size - size_gap_threshold)
      unique_regions <- unique(na.omit(peak_regions))

      #broke code around here
      #need to check if valid lower peaks are already in a region. How did I used to do that?



      if (length(valid_lower_peaks) > 0) {
        browser()
        if (length(unique_regions) > 0) {
          peak_regions[i] <- unique_regions[length(unique_regions)]
        } else {
          peak_regions[i] <- 1
        }
      } else {
        if (length(unique_regions) > 0) {
          peak_regions[i] <- unique_regions[length(unique_regions)] + 1
        } else {
          peak_regions[i] <- 1
        }
      }
    }
  }



  # Add a buffer to each region
  unique_regions <- unique(na.omit(peak_regions))

  for (i in seq_along(unique_regions)) {
    region_positions <- which(peak_regions == i)

    if (length(region_positions) > 0) {
      # Handle boundary cases
      if (region_positions[1] == 1) {
        region_positions <- c(region_positions, region_positions[length(region_positions)] + 1)
      } else if (region_positions[length(region_positions)] == length(smoothed_height)) {
        region_positions <- c(region_positions[1] - 1, region_positions)
      } else {
        region_positions <- c(region_positions[1] - 1, region_positions, region_positions[length(region_positions)] + 1)
      }

      # Ensure we are not accessing out of bounds
      region_positions <- region_positions[region_positions > 0 & region_positions <= length(smoothed_height)]
      peak_regions[region_positions] <- i
    }
  }

  return(peak_regions)
}


# find peaks ------------------------------------------------------------------
find_peaks <- function(heights, size, peak_region_size_gap_threshold, peak_region_height_threshold_multiplier) {
  # Function to find positions of local maxima
  find_maxima_positions <- function(heights, positions) {
    maxima_positions <- which(find_maxima(heights))

    # If finding maxima fails, pick the tallest peak
    if (length(maxima_positions) == 0) {
      maxima_positions <- which(heights == max(heights))
    }

    return(positions[maxima_positions])
  }

  # Function to identify peak regions based on height and size


  # Function to find position of the tallest peak
  find_tallest_peak_position <- function(heights, positions) {
    tallest_peak_position <- which(heights == max(heights))

    # If there are multiple peaks with the same height, choose the largest
    if (length(tallest_peak_position) > 1) {
      tallest_peak_position <- tallest_peak_position[length(tallest_peak_position)]
    }
    return(positions[tallest_peak_position])
  }

  #### Find all top peaks

  # Find peak regions
  peak_regions <- find_peak_regions(heights, size, peak_region_size_gap_threshold, peak_region_height_threshold_multiplier)

  # Find unique peak regions
  unique_regions <- unique(na.omit(peak_regions))
  top_peaks <- numeric(length(unique_regions))

  # Find the tallest peak within each peak region
  for (i in seq_along(unique_regions)) {
    region_positions <- which(peak_regions == i)

    # # Add buffer positions to the region for finding maxima
    # region_positions <- add_position_buffer(region_positions, length(peak_regions))

    # Find maxima positions within the region
    maxima_positions <- find_maxima_positions(heights[region_positions], region_positions)

    # Find the position of the tallest peak within the maxima positions
    tallest_peak_position <- find_tallest_peak_position(heights[maxima_positions], maxima_positions)

    top_peaks[i] <- tallest_peak_position
  }

  # Return the peak regions and top regional peaks
  return(list(peak_regions = peak_regions, top_regional_peaks = top_peaks))
}

# candidate alleles ------------------------------------------------------------
find_candidate_peaks <- function(heights,
                                 size,
                                 peak_region_size_gap_threshold,
                                 peak_region_height_threshold_multiplier,
                                 number_of_peaks_to_return) {


  # Find peak regions
  peak_regions <- find_peak_regions(heights, size, peak_region_size_gap_threshold,
                                    peak_region_height_threshold_multiplier)

  #find all possible peaks
  all_peaks <- pracma::findpeaks(heights, peakpat = "[+]{1,}[0]*[-]{1,}")


  # Find unique peak regions
  unique_regions <- unique(na.omit(peak_regions))
  top_peaks <- numeric(length(unique_regions))

  # Find the tallest peak within each peak region
  for (i in seq_along(unique_regions)) {
    region_positions <- which(peak_regions == i)

    # Find the position of the tallest peak within the maxima positions
    peak_subset <- all_peaks[which(all_peaks[,2] %in% region_positions), ]
    tallest_peak_position <- all_peaks[which.max(peak_subset[,1]), 2]

    top_peaks[i] <- tallest_peak_position
  }



  # do a first pass and if only one significant peak region found when we expect two,
  # see if there are two significant maxima in the region
  # this is to identify alleles close in size and homozygous


  #Too do
  ##if using findpeaks, need to remove find maxima call down here. might help simplify matters

  if (length(output$top_regional_peaks) == 1 &
    number_of_peaks_to_return == 2) {
    region_positions <- which(output$peak_regions == 1)
    region_heights <- heights[region_positions]
    region_maxima <- find_maxima(region_heights)
    significant_maxima <-
      which(region_heights > max(region_heights) * 0.5)
    # chose the two tallest maxima if more than one peak has now been found
    if (length(significant_maxima) > 1) {
      sig_maxima_heights <- region_heights[significant_maxima]
      second_tallest_height <-
        sig_maxima_heights[order(sig_maxima_heights, decreasing = TRUE)][2]
      top_peaks <-
        region_positions[which(region_heights[significant_maxima] >= second_tallest_height)][1:2]
    }# deal with case where the peak after the wt peak is pretty high, perhaps indicating heterozygous +1
    else if (heights[output$top_regional_peaks[[1]] + 1] / heights[output$top_regional_peaks[[1]]] > 0.5) {
      top_peaks <- c(output$top_regional_peaks[[1]], output$top_regional_peaks[[1]] + 1)
    } # homozygous
    else {
      top_peaks <- c(output$top_regional_peaks[[1]], output$top_regional_peaks[[1]])
    }
  } else {
    top_peaks <- output$top_regional_peaks

    top_peaks <-
      top_peaks[order(heights[top_peaks], decreasing = TRUE)][1:number_of_peaks_to_return]
  }


  return(list(
    top_peaks = top_peaks,
    peak_regions = output$peak_regions
  ))
}

################### R6 Class Method Helpers #################################

find_main_peaks_helper <- function(fragments_repeats_class,
                                   number_of_peaks_to_return,
                                   peak_region_size_gap_threshold,
                                   peak_region_height_threshold_multiplier) {
  # do some checks
  if (!number_of_peaks_to_return %in% c(1, 2)) {
    stop("number_of_peaks_to_return must be 1 or 2",
      call. = FALSE
    )
  }

  # this abstraction doesn't work because the peak data is hard coded below
  # change it so the peak data is passed as an argument

  #first select if working off repeat size or bp size
  fragment_heights <- if (is.null(fragments_repeats_class$repeat_table_df)) fragments_repeats_class$peak_table_df$height else fragments_repeats_class$repeat_table_df$height
  fragment_sizes <- if (is.null(fragments_repeats_class$repeat_table_df)) fragments_repeats_class$peak_table_df$size else fragments_repeats_class$repeat_table_df$repeats

  top_peaks <- find_candidate_peaks(
    fragment_heights,
    fragment_sizes,
    number_of_peaks_to_return = number_of_peaks_to_return,
    peak_region_size_gap_threshold = peak_region_size_gap_threshold,
    peak_region_height_threshold_multiplier = peak_region_height_threshold_multiplier
  )

  if (length(top_peaks$top_peaks) == 0) {
    warning(paste0(fragments_repeats_class$unique_id, ": No main alleles identified"))
  }

  # change this so that it populates either repeat, size

  if (number_of_peaks_to_return == 2) {
    # shorter repeat allele
    if (is.null(fragments_repeats_class$repeat_table_df)) {
      fragments_repeats_class$allele_2_size <- fragment_sizes[top_peaks$top_peaks[1]]
    } else {
      fragments_repeats_class$allele_2_repeat <- fragment_sizes[top_peaks$top_peaks[1]]
    }
    fragments_repeats_class$allele_2_height <- fragment_heights[top_peaks$top_peaks[1]]

    # longer repeat allele
    if (is.null(fragments_repeats_class$repeat_table_df)) {
      fragments_repeats_class$allele_1_size <- fragment_sizes[top_peaks$top_peaks[2]]
    } else {
      fragments_repeats_class$allele_1_repeat <- fragment_sizes[top_peaks$top_peaks[2]]
    }
    fragments_repeats_class$allele_1_height <- fragment_heights[top_peaks$top_peaks[2]]
  } else if (number_of_peaks_to_return == 1) {
    # shorter repeat allele doesn't exist
    fragments_repeats_class$allele_2_size <- NA_real_
    fragments_repeats_class$allele_2_repeat <- NA_real_
    fragments_repeats_class$allele_2_height <- NA_real_

    # longer repeat allele
    if (is.null(fragments_repeats_class$repeat_table_df)) {
      fragments_repeats_class$allele_1_size <- fragment_sizes[top_peaks$top_peaks[1]]
    } else {
      fragments_repeats_class$allele_1_repeat <- fragment_sizes[top_peaks$top_peaks[1]]
    }
    fragments_repeats_class$allele_1_height <- fragment_heights[top_peaks$top_peaks[1]]
  }

  # peak_regions
  fragments_repeats_class$.__enclos_env__$private$peak_regions <-
    top_peaks$peak_regions

  return(fragments_repeats_class)
}
