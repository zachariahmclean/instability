# General helper functions --------------------------------------------------------

nCr <- function(n, r){
  c = factorial(n) / (factorial(r) * factorial(n - r))
  return(c)
}

detrend_signal = function(x, bins = 50){
  x_split <- split(x, ceiling(seq_along(x)/(length(x) / bins)))
  median_signal <- sapply(x_split, median)

  x_detrended_list <- vector('list', length = length(x_split))
  for (i in seq_along(x_detrended_list)) {
    x_detrended_list[[i]] <- x_split[[i]] - median_signal[i]
  }
  x_detrended <- unlist(x_detrended_list, use.names = FALSE)

  return(x_detrended)
}


# method specific

process_ladder_signal = function(ladder,
                                 scans,
                                 spike_location,
                                 smoothing_window){
  ladder_df <- data.frame(signal = ladder, scan = scans)
  ladder_df <- ladder_df[which(ladder_df$scan >= spike_location),]
  ladder_df$detrended_signal <- detrend_signal(ladder_df$signal)
  ladder_df$smoothed_signal <- pracma::savgol(ladder_df$detrended_signal,
                                              smoothing_window)
  return(ladder_df)
}

find_ladder_peaks = function(ladder_df,
                             n_reference_sizes){
  median_signal <- median(ladder_df$smoothed_signal, na.rm = TRUE)
  sd_signal <- stats::sd(ladder_df$smoothed_signal, na.rm = TRUE)

  ladder_peaks = vector("numeric")
  ladder_peak_threshold = 1

  while (length(ladder_peaks) < n_reference_sizes) {

    peaks <- pracma::findpeaks(ladder_df$smoothed_signal,
                               peakpat = "[+]{6,}[0]*[-]{6,}", #see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
                               minpeakheight = median_signal + sd_signal * ladder_peak_threshold)

    ladder_peaks <- ladder_df$scan[peaks[,2]]

    # lower the threshold for the next cycle
    ladder_peak_threshold = ladder_peak_threshold - 0.01

    # provide an exit if there are not enough peaks found
    if(sd_signal * ladder_peak_threshold <= 0){
      break
    }
  }

  #go through raw signal and make sure that the identified scan in the smoothed signal is still the highest
  #it will also deal with cases where the scans have the same height (which.max will chose first)
  n_scans <- length(ladder_df$scan)
  window_width <- 3
  peak_position <- numeric(length(ladder_peaks))
  for (i in seq_along(peak_position)) {


    if(ladder_peaks[i] + window_width > 1 & ladder_peaks[i] + window_width < n_scans ){ # make sure that the subsetting would be in bounds when taking window into account
      max_peak <- which.max(ladder_df$signal[(ladder_peaks[i] - window_width):(ladder_peaks[i] + window_width)])

      peak_position[i] <- ladder_peaks[i] - window_width -1 + max_peak
    }
    else{
      peak_position[i] <- ladder_peaks[i]
    }
  }

  return(ladder_peaks)
}



ladder_iteration = function(reference_sizes,
                            observed_sizes,
                            choose,
                            max_combinations){

  find_best_combination <- function(recombinations,
                                    reference_sizes){
    rsq_vector <- vector("numeric", ncol(recombinations))
    for (i in 1:ncol(recombinations)) {
      rsq_vector[[i]] <-  stats::cor(reference_sizes, recombinations[,i]) ^ 2
    }
    #select the values from the original recombinations matrix not the one with the already selected values added for the regression
    selected_recombinations <- recombinations[, which.max(rsq_vector)]

    return(selected_recombinations)
  }

  assigned_reference = vector('numeric')
  assigned_observed = vector('numeric')
  max_iterations <- 1000
  iteration_count <- 1

  #Keep going until the the number reference peaks left is small
  while (length(reference_sizes) > 0 && iteration_count < max_iterations) {

    # observed_sizes and reference_sizes reset each loop

    # calculate the window to look at for this loop
    n_observed = length(observed_sizes)
    n_refernce = length(reference_sizes)
    remainder = n_refernce - choose #how many peaks will be left to be assigned after this selection
    start_window = n_observed - remainder


    if(n_observed == n_refernce){
      assigned_observed = append(assigned_observed, observed_sizes)
      assigned_reference = append(assigned_reference, reference_sizes)
      break
    }
    else{
      n_recombinations = nCr(length(observed_sizes[1:start_window]), choose)
      if(n_recombinations > max_combinations){
        stop(call. = FALSE,
             paste0("Too many combinations to test (", n_recombinations,")"))
      }

      recombinations = utils::combn(observed_sizes[1:start_window] , choose)

      selected_recombinations = find_best_combination(
        recombinations = recombinations,
        reference_sizes = reference_sizes[1:choose]
      )
    }

    #assign the selections
    assigned_observed = append(assigned_observed, selected_recombinations)
    assigned_reference = append(assigned_reference, reference_sizes[1:choose])

    #set up up the next iteration
    last_selected_scan = selected_recombinations[length(selected_recombinations)]
    last_selected_scan_position = which(observed_sizes == last_selected_scan)

    last_selected_reference = reference_sizes[choose]
    last_selected_reference_position = which(reference_sizes == last_selected_reference)

    reference_sizes_left = length(reference_sizes) - last_selected_reference_position

    # fix introduction of NAs by choose being longer than remaining sizes
    if(reference_sizes_left < choose){
      choose = reference_sizes_left
    }

    #break for final iteration
    if(choose == 0 ){
      break
    }

    #deal with situation of just a couple of reference peaks left over, but you need a good amount to make correlation
    # therefore increase chose to the max
    if(reference_sizes_left - choose < 4){
      choose = reference_sizes_left
    }
    else if(reference_sizes_left < choose){ # fix introduction of NAs by choose being longer than remaining sizes
      choose = reference_sizes_left
    }

    reference_sizes = reference_sizes[(last_selected_reference_position + 1):length(reference_sizes)]
    observed_sizes = observed_sizes[(last_selected_scan_position + 1):length(observed_sizes)]

    iteration_count <- iteration_count + 1

  }

  return(data.frame(scan = assigned_observed, size = assigned_reference))

}







fit_ladder <- function(ladder,
                       scans,
                       ladder_sizes,
                       spike_location,
                       smoothing_window,
                       max_combinations,
                       ladder_selection_window){



  if (is.null(scans)){
    scans <- 0:(length(ladder) - 1)
  }

  if(is.null(spike_location)){
    spike_location <- which.max(ladder) + 50
  }

  # print(paste("Spike scan number:", spike_location))

  ladder_df <- process_ladder_signal(
    ladder = ladder,
    scans = scans,
    spike_location = spike_location,
    smoothing_window = smoothing_window)

  ladder_peaks <- find_ladder_peaks(
    ladder_df = ladder_df,
    n_reference_sizes = length(ladder_sizes))

  peaks_fit_df <- ladder_iteration(
      reference_sizes = rev(ladder_sizes), #start away from the spike going backwards
      observed_sizes = rev(ladder_peaks),
      choose = ladder_selection_window,
      max_combinations = max_combinations
    )

  peaks_not_fit <- ladder_peaks[which(!ladder_peaks %in% peaks_fit_df$scan)]

  peaks_not_fit_df <- data.frame(scan = peaks_not_fit,
                                 size = rep(NA_real_, length(peaks_not_fit)))


  combined_ladder_peaks <- rbind(peaks_fit_df, peaks_not_fit_df)
  combined_ladder_peaks <- combined_ladder_peaks[order(combined_ladder_peaks$scan), ]

  return(combined_ladder_peaks)

}


ladder_rsq_warning_helper <- function(framents_trace,
                                      rsq_threshold){


  rsq <- sapply(framents_trace$mod_parameters, function(x) suppressWarnings(summary(x$mod)$r.squared))
  if(any(rsq < rsq_threshold)){
    size_ranges <- sapply(framents_trace$mod_parameters, function(x) x$mod$model$yi)
    size_ranges <- size_ranges[ , which(rsq < rsq_threshold), drop = FALSE]
    size_ranges_vector <- vector('numeric', ncol(size_ranges))
    for (j in seq_along(size_ranges_vector)) {
      size_ranges_vector[j] <- paste0(size_ranges[1,j],"-",size_ranges[3,j])

    }
    warning(call. = FALSE,
            paste("sample", framents_trace$unique_id, "has badly fitting ladder for bp sizes:",
                  paste0(size_ranges_vector, collapse = ", ")))
  }
}


# bp sizing ---------------------------------------------------------------


local_southern_fit <- function(x, y) {

  #do some quality control. There should be no missing values and vectors should be same length
  if(length(x) != length(y)){
    stop(call. = FALSE,
         "local_southern_fit error: ladder scan and size vectors different lengths")
  }
  else if(any(is.na(x)) | any(is.na(y))){
    stop(call. = FALSE,
         "local_southern_fit error: missing values in ladder scan or size")
  }


  # Sort the data points by x values
  sorted_indices <- order(x)
  x_sorted <- x[sorted_indices]
  y_sorted <- y[sorted_indices]

  # Function to calculate the fitting constants for each group of three neighboring points
  mod_list <- vector("list", length = length(x_sorted) - 2)

  for (i in 1:(length(x_sorted) - 2)) {
    xi <- x_sorted[i:(i+2)]
    yi <- y_sorted[i:(i+2)]
    mod_list[[i]] <- list(mod = lm(yi ~ xi),
                          first = xi[1],
                          last = xi[3])
  }

  return(mod_list)
}

local_southern_predict <- function(local_southern_fit, scans) {

  #total number of groups to brake the scans into:
  ladder_scan_pos <- sapply(local_southern_fit, function(fit) fit$first)

  # Find the nearest ladder position for each scan position
  nearest_ladder_index <- sapply(scans, function(scan) which.min(abs(scan - ladder_scan_pos)))

  # Assign the scan positions to corresponding groups based on nearest ladder position
  group_assignments <- rep(NA, length(scans))
  for (i in seq_along(nearest_ladder_index)) {
    group_assignments[i] <- nearest_ladder_index[i]
  }

  scan_split <- split(scans, group_assignments)

  size_split <- vector("list", length = length(scan_split))
  for(i in seq_along(scan_split)){

    if(i == 1 | i == length(scan_split)){
      size_split[[i]] = stats::predict(local_southern_fit[[i]]$mod, data.frame(xi = scan_split[[i]]))
    }
    else{
      lower_prediction <- stats::predict(local_southern_fit[[i-1]]$mod, data.frame(xi = scan_split[[i]]))
      upper_prediction <- stats::predict(local_southern_fit[[i]]$mod, data.frame(xi = scan_split[[i]]))
      size_split[[i]] <- (lower_prediction + upper_prediction) / 2
    }
  }

  size <- unlist(size_split)

  return(size)
}



find_ladder_helper <- function(fragments_trace,
                               fsa,
                               ladder_channel,
                               signal_channel,
                               ladder_sizes,
                               spike_location,
                               scan_subset,
                               smoothing_window,
                               max_combinations,
                               ladder_selection_window,
                               show_progress_bar){

  if(show_progress_bar){
    pb <- utils::txtProgressBar(min = 0, max = length(ladder_list), style = 3)
  }

  fragments_trace$raw_ladder <- fsa$Data[[ladder_channel]]
  fragments_trace$raw_data <- fsa$Data[[signal_channel]]
  fragments_trace$scan <- 0:(length(fsa$Data[[signal_channel]])- 1)
  fragments_trace$off_scale_scans <- fsa$Data$OfSc.1

  #allow user to subset to particular scans
  if(!is.null(scan_subset)){

    fragments_trace$raw_ladder <- fragments_trace$raw_ladder[scan_subset[1]:scan_subset[2]]
    fragments_trace$raw_data <- fragments_trace$raw_data[scan_subset[1]:scan_subset[2]]
    fragments_trace$scan <- fragments_trace$scan[scan_subset[1]:scan_subset[2]]

    # set spike location since it's automatically set usually, and user may select scans to start after
    spike_location <- scan_subset[1]
  }

  #ladder
  ladder_df <- fit_ladder(
    ladder = fragments_trace$raw_ladder,
    scans = fragments_trace$scan,
    ladder_sizes = ladder_sizes,
    spike_location = spike_location,
    smoothing_window = smoothing_window,
    max_combinations = max_combinations,
    ladder_selection_window = ladder_selection_window)

  fragments_trace$ladder_df <- ladder_df

  # predict bp size
  ladder_df <- ladder_df[which(!is.na(ladder_df$size)), ]
  ladder_df <- ladder_df[which(!is.na(ladder_df$scan)), ]
  fragments_trace$mod_parameters <- local_southern_fit(ladder_df$scan, ladder_df$size)

  predicted_size <- local_southern_predict(local_southern_fit = fragments_trace$mod_parameters , scans = fragments_trace$scan)

  fragments_trace$trace_bp_df <- data.frame(
    unique_id = rep(fragments_trace$unique_id, length(fragments_trace$scan)),
    scan = fragments_trace$scan,
    size = predicted_size,
    signal = fragments_trace$raw_data,
    ladder_signal = fragments_trace$raw_ladder,
    off_scale = fragments_trace$scan %in% fragments_trace$off_scale_scans
  )

  #make a warning if one of the ladder modes is bad
  ladder_rsq_warning_helper(fragments_trace,
                            rsq_threshold = 0.998)

  if(show_progress_bar){
      utils::setTxtProgressBar(pb, i)

  }

  return(fragments_trace)


}







# peak calling ------------------------------------------------------------


find_fragment_peaks <- function(trace_bp_df,
                                smoothing_window,
                                minimum_peak_signal,
                                ...){

  # smoothed_signal <- moving_average(trace_bp_df$signal,
  #                                   n = smoothing_window)
  #
  smoothed_signal <- pracma::savgol(trace_bp_df$signal,
                                    smoothing_window)

  #deals with cases of user overriding values
  if("peakpat" %in% ...names()){
    peaks <- pracma::findpeaks(smoothed_signal,
                               minpeakheight = minimum_peak_signal,
                               ...)
  }
  else if("minpeakheight" %in% ...names()) {
    peaks <- pracma::findpeaks(smoothed_signal,
                               peakpat = "[+]{6,}[0]*[-]{6,}", #see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
                               ...)
  }
  else{
    peaks <- pracma::findpeaks(smoothed_signal,
                               peakpat = "[+]{6,}[0]*[-]{6,}", #see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
                               minpeakheight = minimum_peak_signal,
                               ...)
  }

  n_scans <- length(trace_bp_df$signal)
  window_width <- 3

  #go through raw signal and make sure that the identified scan in the smoothed signal is still the highest
  #it will also deal with cases where the scans have the same height (which.max will chose first)
  peak_position <- numeric(nrow(peaks))
  for (i in seq_along(peak_position)) {


    if(peaks[i,2] + window_width > 1 & peaks[i,2] + window_width < n_scans ){ # make sure that the subsetting would be in bounds when taking window into account
      max_peak <- which.max(trace_bp_df$signal[(peaks[i,2] - window_width):(peaks[i,2] + window_width)])

      peak_position[i] <- peaks[i,2] - window_width -1 + max_peak
    }
    else{
      peak_position[i] <- peaks[i,2]
    }
  }

  df <- trace_bp_df[peak_position, c("scan", "size", "signal", "off_scale")]
  colnames(df) <- c("scan", "size", "height", "off_scale")

  # remove shoulder peaks
  df2 <- deshoulder(df, shoulder_window = 1.5)

  return(df2)
}



# ladder fixing -----------------------------------------------------------


ladder_fix_helper <- function(fragments_trace,
                              replacement_ladder_df){
  fragments_trace_copy <- fragments_trace$clone()

  fragments_trace_copy$ladder_df <- replacement_ladder_df
  ladder_df <- fragments_trace_copy$ladder_df[which(!is.na(fragments_trace_copy$ladder_df$size)), ]
  ladder_df <- ladder_df[which(!is.na(ladder_df$scan)), ]
  fragments_trace_copy$mod_parameters <- local_southern_fit(ladder_df$scan, ladder_df$size)

  predicted_size <- local_southern_predict(local_southern_fit = fragments_trace_copy$mod_parameters , scans = fragments_trace_copy$scan)

  fragments_trace_copy$trace_bp_df <- data.frame(
    unique_id = rep(fragments_trace_copy$unique_id, length(fragments_trace_copy$scan)),
    scan = fragments_trace_copy$scan,
    size = predicted_size,
    signal = fragments_trace_copy$raw_data,
    ladder_signal = fragments_trace_copy$raw_ladder
  )

  #make a warning if one of the ladder modes is bad
  ladder_rsq_warning_helper(fragments_trace_copy,
                            rsq_threshold = 0.998)

  return(fragments_trace_copy)

}


ladder_self_mod_predict <- function(fragments_trace,
                              size_threshold,
                              size_tolerance,
                              rsq_threshold){
  fragments_trace_copy <- fragments_trace$clone()

  ladder_sizes <- fragments_trace_copy$ladder_df[which(!is.na(fragments_trace_copy$ladder_df$size)), "size"]
  ladder_peaks <- fragments_trace_copy$ladder_df$scan

  mod_validations <- vector("list", length(fragments_trace_copy$mod_parameters))
  for (i in seq_along(fragments_trace_copy$mod_parameters)) {

    predictions <- stats::predict(fragments_trace_copy$mod_parameters[[i]]$mod, newdata = data.frame(xi = ladder_peaks))
    low_size_threshold <- fragments_trace_copy$mod_parameters[[i]]$mod$model$yi[1] - size_threshold
    high_size_threshold <- fragments_trace_copy$mod_parameters[[i]]$mod$model$yi[3] + size_threshold

    predictions_close <- predictions[which(predictions > low_size_threshold & predictions < high_size_threshold)]
    mod_sizes <- fragments_trace_copy$mod_parameters[[i]]$mod$model$yi

    ladder_hits <- sapply(predictions_close, function(x){
      diff = ladder_sizes - x
      ladder_hit <- ladder_sizes[which(diff > -size_tolerance & diff < size_tolerance)]
      if(length(ladder_hit) == 0){
        return(NA_real_)
      }
      else{
        return(ladder_hit)
      }
    })


    ladder_hits <- ladder_hits[!ladder_hits %in% mod_sizes & !is.na(ladder_hits)]

    mod_validations[[i]]$predictions <- predictions
    mod_validations[[i]]$predictions_close <- predictions_close
    mod_validations[[i]]$mod_sizes <- mod_sizes
    mod_validations[[i]]$rsq <- summary(fragments_trace_copy$mod_parameters[[i]]$mod)$r.squared
    mod_validations[[i]]$ladder_hits <- ladder_hits
    mod_validations[[i]]$mod <- fragments_trace_copy$mod_parameters[[i]]$mod
  }

  #predicted to be a good model if it has some ladder hits and good rsq
  #use the good models to predict the rest of the peaks
  valid_models_tf <- sapply(mod_validations, function(x) ifelse(length(x$ladder_hits) > 0 & x$rsq > rsq_threshold, TRUE, FALSE))
  valid_models <- mod_validations[which(valid_models_tf)]
  predicted_sizes_list <- lapply(valid_models, function(x){
    sizes <- vector("numeric", length(x$predictions))
    for (i in seq_along(x$predictions)) {
      if(x$predictions[i] %in% x$predictions_close){
        sizes[i] <- x$predictions[i]
      }
      else{
        sizes[i] <- NA_real_
      }
    }

    return(sizes)

  })


  predicted_sizes_matrix <- do.call(cbind, predicted_sizes_list)
  predicted_sizes_avg <- numeric(nrow(predicted_sizes_matrix))
  for (i in seq_along(predicted_sizes_avg)) {
    predicted_sizes_avg[i] <- median(predicted_sizes_matrix[i,],
                                     na.rm = TRUE)
  }

  confirmed_sizes <- unique(unlist(lapply(valid_models, function(x) x$mod_sizes)))
  unconfirmed_sizes <- ladder_sizes[which(!ladder_sizes %in% confirmed_sizes)]

  assigned_size <- sapply(predicted_sizes_avg, function(x){
    diff = ladder_sizes - x
    ladder_hit <- ladder_sizes[which(diff > -size_tolerance*2 & diff < size_tolerance*2)]
    if(length(ladder_hit) == 0){
      return(NA_real_)
    }
    else{
      return(ladder_hit)
    }
  })

  assigned_df <- data.frame(
    scan = ladder_peaks[which(assigned_size %in% unconfirmed_sizes)],
    size = assigned_size[which(assigned_size %in% unconfirmed_sizes)]
  )

  #bind rows back with the sizes that were infered to be correct

  ladder_df <- rbind(assigned_df,
        fragments_trace_copy$ladder_df[which(fragments_trace_copy$ladder_df$size %in% confirmed_sizes),] )

  # now just rerun the bp sizing
  fixed_fragment_trace <- ladder_fix_helper(
    fragments_trace_copy,
    ladder_df
  )

  return(fixed_fragment_trace)

}



