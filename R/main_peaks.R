# general helper functions for finding peaks
##
find_maxima = function(x) {
  maxima <- logical(length(x))

  for (i in 2:(length(x) - 1)) {
    maxima[i] <- (x[i] > x[i - 1] && x[i] > x[i + 1])
  }

  maxima[1] <- (length(x) > 1) && (x[1] > x[2])
  maxima[length(x)] <- (length(x) > 1) && (x[length(x)] > x[length(x) - 1])

  return(maxima)
}

##
find_peaks = function(heights,
                      size,
                      peak_region_size_gap_threshold,
                      peak_region_height_threshold_multiplier){
  find_maxima_positions = function(heights,
                                   positions) {

    maxima_positions <- which(find_maxima(heights))
    #if finding maxima fails, just pick the tallest peak
    if (length(maxima_positions) == 0) {
      maxima_positions <- which(heights == max(heights))
    }

    return(positions[maxima_positions])
  }
##
find_peak_regions = function(height,
                             size,
                             size_gap_threshold = peak_region_size_gap_threshold,
                             region_height_threshold_multiplier = peak_region_height_threshold_multiplier) {
    smoothed_height <- as.vector(smooth(height))
    peak_region <- rep(NA_real_, length(smoothed_height)) # make empty vector
    mean_heights <- mean(smoothed_height) * region_height_threshold_multiplier
    for (i in seq_along(smoothed_height)) {
      #figure out the distance to the peak after
      if(i == length(size)){ #if its the final peak, there's no gap so make it zero
        lead_gap <- 0
      }
      else{
        lead_gap <- size[[i + 1]] - size[[i]]
      }

      #find the peak regions
      #if it's below the average peak height, it gets no region (NA)
      if (smoothed_height[[i]] < mean_heights) {
        peak_region[[i]] <- NA_real_
      }   #If if the first or last, it also gets no region
      else if(i == 1 | i == length(smoothed_height)){
        peak_region[[i]] <- NA_real_
      }  #peak either side below threshold, this could still be peak region?
      else if(smoothed_height[[i - 1]] < mean_heights & smoothed_height[[i + 1]] < mean_heights){
        peak_region[[i]] <- NA_real_
      } # if it's the first element above threshold, can start the first region
      else if (all(is.na(peak_region)) == TRUE) {
        peak_region[[i]] <- 1
      } #if there's a peak in the size gap, then make it part of the same region
      else if (any(!is.na(peak_region[which(size < size[[i]] & size > size[[i]] - size_gap_threshold)]))) {
        unique_regions <- unique(na.omit(peak_region[which(size < size[[i]] & size > size[[i]] - size_gap_threshold)]))
        peak_region[[i]] <- unique_regions[length(unique_regions)]
      } #no peak before, then start a new region
      else if (any(is.na(peak_region[which(size < size[[i]] & size > size[[i]] - size_gap_threshold)]))) {
        above_threshold <- na.omit(peak_region)
        peak_region[[i]] <- above_threshold[length(above_threshold)] + 1
      } #if there's anything left make it NA
      else {
        peak_region[[i]] <- NA_real_
      }
    }
    return(peak_region)
  }
##
  add_position_buffer = function(region_positions, full_vector_length) {
      if (region_positions[1] == 1) {
        region_positions <- c(region_positions, region_positions[length(region_positions)] + 1)
      } #if the peaks are the last ones in the vector
      else if (region_positions[length(region_positions)] == full_vector_length) {
        region_positions <- c(region_positions[1] - 1 , region_positions)
      }
      else{
        region_positions <-
          c(region_positions[1] - 1 ,
            region_positions,
            region_positions[length(region_positions)] + 1)
      }
    }
##
  find_tallest_peak_position = function(heights,
                                          positions) {
      tallest_peak_position <- which(heights == max(heights))
      #if there are multiple peaks with the same height, choose the largest
      if (length(tallest_peak_position) > 1) {
        tallest_peak_position <- tallest_peak_position[length(tallest_peak_position)]
      }
      return(positions[tallest_peak_position])
    }

  #### find all top peaks

  peak_regions <- find_peak_regions(heights, size)

  unique_regions <- unique(na.omit(peak_regions))
  top_peaks <- vector("numeric", length(unique_regions))
  #find the top peak within each peak region
  for (i in seq_along(unique_regions)) {
    region_positions <- which(peak_regions == i)
    #add a buffer for finding maxima
    region_positions <- add_position_buffer(region_positions, length(peak_regions))
    #Find the maxima within the region
    maxima_positions <- find_maxima_positions(heights[region_positions],
                                                   region_positions)
    # pick tallest and deal with cases where there's equally tall peaks
    tallest_peak_position <- find_tallest_peak_position(heights[maxima_positions],
                                                             maxima_positions)
    top_peaks[[i]] <- tallest_peak_position

  }

  return(list(peak_region = peak_regions,
              top_regional_peaks = top_peaks))

}
##
find_candidate_peaks = function(heights,
                        size,
                        peak_region_size_gap_threshold,
                        peak_region_height_threshold_multiplier,
                        number_of_peaks_to_return,
                        return_peak_regions ){
  #find peaks

  output<- find_peaks(heights = heights,
                      size = size,
                      peak_region_size_gap_threshold = peak_region_size_gap_threshold,
                      peak_region_height_threshold_multiplier = peak_region_height_threshold_multiplier)

  #do a first pass and if only one significant peak region found when we expect two,
  #see if there are two significant maxima in the region
  #this is to identify alleles close in size and homozygous


  if(length(output$top_regional_peaks) == 1 & number_of_peaks_to_return == 2){
    region_positions <- which(output$peak_region == 1)
    region_heights <- heights[region_positions]
    region_maxima <- find_maxima(region_heights)
    significant_maxima <- which(region_heights > max(region_heights) * 0.5)
    #chose the two tallest maxima if more than one peak has now been found
    if(length(significant_maxima) > 1){
      sig_maxima_heights <- region_heights[significant_maxima]
      second_tallest_height <- sig_maxima_heights[order(sig_maxima_heights, decreasing = TRUE)][2]
      top_peaks <- region_positions[which(region_heights[significant_maxima] >= second_tallest_height)][1:2]
    } #deal with case where the peak after the wt peak is pretty high, perhaps indicating heterozygous +1
    else if(region_heights[significant_maxima + 1] / region_heights[significant_maxima] > 0.5){
      top_peaks <- region_positions[c(significant_maxima, significant_maxima + 1)]
    } #homozygous
    else {
      top_peaks <- c(region_positions[significant_maxima], region_positions[significant_maxima])
    }
  }
  else {
    top_peaks <- output$top_regional_peaks

    top_peaks <- top_peaks[order(heights[top_peaks], decreasing = TRUE)][1:number_of_peaks_to_return]
  }

  if(return_peak_regions == FALSE){
    return(top_peaks)
  }
  else{
    return(list(top_peaks = top_peaks,
                peak_region = output$peak_region))
  }



}


# class abstraction -------------------------------------------------------

find_main_peaks_helper <- function(
    fragments_class,
    fragment_sizes,
    fragment_heights,
    data_type,
    number_of_peaks_to_return,
    peak_region_size_gap_threshold,
    peak_region_height_threshold_multiplier){

 # do some checks
  if(!number_of_peaks_to_return %in% c(1,2)){
    stop("number_of_peaks_to_return must be 1 or 2",
         call. = FALSE)
  }

#this abstraction doesn't work because the peak data is hard coded below
  #change it so the peak data is passed as an argument

  top_peaks <- find_candidate_peaks(
    fragment_heights,
    fragment_sizes,
    number_of_peaks_to_return = number_of_peaks_to_return,
    peak_region_size_gap_threshold = peak_region_size_gap_threshold,
    peak_region_height_threshold_multiplier = peak_region_height_threshold_multiplier,
    return_peak_regions = TRUE
  )

  if(length(top_peaks$top_peaks) == 0){
    warning(paste0(self2$unique_id, ": No main alleles identified"))
  }

  # data type
  # need to do this eval/parse below so that this helper function can be used in either bp_fragments or repeats_fragments
  type <- ifelse(data_type == "bp_size", "size", "repeat")

  if(number_of_peaks_to_return == 2){
    #shorter repeat allele
    eval(parse(text = paste0("fragments_class$allele_2_", type, " <- fragment_sizes[top_peaks$top_peaks[1]]")))
    fragments_class$allele_2_height <- fragment_heights[top_peaks$top_peaks[1]]
    #longer repeat allele
    eval(parse(text = paste0("fragments_class$allele_1_", type, "  <- fragment_sizes[top_peaks$top_peaks[2]]")))
    fragments_class$allele_1_height <- fragment_heights[top_peaks$top_peaks[2]]

  }else if(number_of_peaks_to_return == 1){
    #shorter repeat allele doesn't exist
    eval(parse(text = paste0("fragments_class$allele_2_", type, " <- NA_real_")))
    fragments_class$allele_2_height <- NA_real_
    #longer repeat allele
    eval(parse(text = paste0("fragments_class$allele_1_", type, "  <- fragment_sizes[top_peaks$top_peaks[1]]")))
    fragments_class$allele_1_height <- fragment_heights[top_peaks$top_peaks[1]]
  }

  #peak_regions
  fragments_class$.__enclos_env__$private$peak_regions <- top_peaks$peak_regions

  return(fragments_class)
}



