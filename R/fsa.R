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

moving_average = function(x, n = 10){
  #centered, n =  total scans
  as.vector(stats::filter(x, rep(1 / n, n), sides = 2))
}

# method specific

process_ladder_signal = function(ladder,
                                 scans,
                                 spike_location,
                                 smoothing_window){
  ladder_df <- data.frame(signal = ladder)
  ladder_df$scan <- scans

  if(is.null(spike_location)){
    spike_location <- which.max(ladder_df$signal) + 45
  }

  ladder_df <- ladder_df[spike_location:nrow(ladder_df),]

  ladder_df$detrended_signal <- detrend_signal(ladder_df$signal)
  ladder_df$smoothed_signal <- moving_average(ladder_df$detrended_signal, n = smoothing_window)
  # peaks <- pracma::findpeaks(ladder_df$smoothed_signal,
  #                                       peakpat = "[+]{1,}[0]+[-]{1,}" #see https://stackoverflow.com/questions/46864211/find-sustained-peaks-using-pracmafindpeaks
  #                                       )
  # ladder_df$maxima <- find_maxima(ladder_df$scan %in% peaks[,2])
  ladder_df$maxima <- find_maxima(ladder_df$smoothed_signal)

  return(ladder_df)
}

find_ladder_peaks = function(ladder_df,
                             n_reference_sizes){
  median_signal <- median(ladder_df$smoothed_signal, na.rm = TRUE)
  sd_signal <- sd(ladder_df$smoothed_signal, na.rm = TRUE)

  ladder_peaks = vector("numeric")
  ladder_peak_threshold = 1

  while (length(ladder_peaks) < n_reference_sizes) {
    ladder_df$peak <- ladder_df$smoothed_signal > median_signal + sd_signal * ladder_peak_threshold & ladder_df$maxima
    ladder_peaks <- ladder_df[which(ladder_df$peak), "scan"]

    ladder_peak_threshold = ladder_peak_threshold - 0.01
  }

  return(ladder_peaks)
}



ladder_iteration = function(reference_sizes,
                            observed_sizes,
                            choose = 4,
                            max_combinations = 2500000){

  find_best_combination <- function(recombinations,
                                    reference_sizes){
    rsq_vector <- vector("numeric", ncol(recombinations))
    for (i in 1:ncol(recombinations)) {
      rsq_vector[[i]] <-  cor(reference_sizes, recombinations[,i]) ^ 2
    }
    #select the values from the original recombinations matrix not the one with the already selected values added for the regression
    selected_recombinations <- recombinations[, which.max(rsq_vector)]

    return(selected_recombinations)
  }

  assigned_reference = vector('numeric')
  assigned_observed = vector('numeric')

  #Keep going until the the number reference peaks left is small
  while (length(reference_sizes) > 0) {

    #calculate the window to look at
    n_observed = length(observed_sizes)
    n_refernce = length(reference_sizes)
    remainder = n_refernce - choose #how many peaks will left to be assigned after this selection
    start_window = n_observed - remainder  #What is going on here? something to do with keeping track of how many peaks are left to assign in the observed


    if(start_window == 1){ # this means that the length of the remaining lists is the same so we can line them up
      selected_recombinations = observed_sizes
    }
    else{
      n_recombinations = nCr(length(observed_sizes[1:start_window]), choose)
      if(n_recombinations > max_combinations){
        stop(call. = FALSE,
             paste0("Too many combinations to test (", n_recombinations,")"))
      }

      recombinations = combn(observed_sizes[1:start_window] , choose)
      # print(paste0(length(observed_sizes[1:start_window]), " pick ", choose))
      # print(paste0(ncol(recombinations), " combinations"))

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

  }

  return(data.frame(scan = assigned_observed, size = assigned_reference))

}







fit_ladder <- function(ladder,
                       scans = NULL,
                       ladder_sizes = NULL,
                       hq_ladder=TRUE,
                       spike_location = NULL,
                       smoothing_window = 5,
                       max_combinations = 2500000,
                       ladder_selection_window = 5){



  if (is.null(scans)){
    scans <- 0:(length(ladder) - 1)
  }

  if (is.null(ladder_sizes)){
    ladder_sizes <- c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500)
  }

  if(hq_ladder){
    ladder_sizes <- ladder_sizes[!ladder_sizes %in% c(35, 250, 340)]
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
      size_split[[i]] = predict(local_southern_fit[[i]]$mod, data.frame(xi = scan_split[[i]]))
    }
    else{
      lower_prediction <- predict(local_southern_fit[[i-1]]$mod, data.frame(xi = scan_split[[i]]))
      upper_prediction <- predict(local_southern_fit[[i]]$mod, data.frame(xi = scan_split[[i]]))
      size_split[[i]] <- (lower_prediction + upper_prediction) / 2
    }
  }

  size <- unlist(size_split)

  return(size)
}


# peak calling ------------------------------------------------------------


#' Title
#'
#' @param fragments_trace
#' @param smoothing_window
#' @param minumum_peak_signal
#' @param min_bp_size
#' @param max_bp_size
#'
#' @return
#' @importFrom pracma findpeaks
#'
#' @examples
find_fragment_peaks <- function(fragments_trace,
                                smoothing_window,
                                minumum_peak_signal,
                                min_bp_size,
                                max_bp_size){

  smoothed_signal <- moving_average(fragments_trace$trace_bp_df$signal,
                                    n = smoothing_window)

  peaks <- pracma::findpeaks(smoothed_signal,
                             peakpat = "[+]{3,}[0]*[-]{3,}", #see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
                             minpeakheight = minumum_peak_signal)
  n_scans <- length(fragments_trace$trace_bp_df$signal)
  window_width <- 3

  peak_position <- numeric(nrow(peaks))
  for (i in seq_along(peak_position)) {


    if(peaks[i,2] + window_width > 1 & peaks[i,2] + window_width < n_scans ){
      max_peak <- which.max(fragments_trace$trace_bp_df$signal[(peaks[i,2] - window_width):(peaks[i,2] + window_width)])

      peak_position[i] <- peaks[i,2] - window_width -1 + max_peak
    }
    else{
      peak_position[i] <- peaks[i,2]
    }
  }


  df <- fragments_trace$trace_bp_df[peak_position, c("scan", "size", "signal")]
  colnames(df) <- c("scan", "size", "height")
  df$unique_id <- rep(fragments_trace$unique_id, nrow(df))
  df <- df[which(df$size > min_bp_size & df$size < max_bp_size), ]

  return(df)
}



# ladder fixing -----------------------------------------------------------




ladder_self_mod_predict <- function(fragments_trace,
                              size_threshold,
                              size_tolerance,
                              rsq_threshold){
  fragments_trace_copy <- fragments_trace$clone()

  ladder_sizes <- fragments_trace_copy$ladder_df[which(!is.na(fragments_trace_copy$ladder_df$size)), "size"]
  ladder_peaks <- fragments_trace_copy$ladder_df$scan

  mod_validations <- vector("list", length(fragments_trace_copy$mod_parameters))
  for (i in seq_along(fragments_trace_copy$mod_parameters)) {

    predictions <- predict(fragments_trace_copy$mod_parameters[[i]]$mod, newdata = data.frame(xi = ladder_peaks))
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

  ladder_df <- rbind(assigned_df,
        fragments_trace_copy$ladder_df[which(fragments_trace_copy$ladder_df$size %in% confirmed_sizes),] )

  fragments_trace_copy$ladder_df <- ladder_df

  return(fragments_trace_copy)

}




