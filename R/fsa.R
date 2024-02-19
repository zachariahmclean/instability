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

local_southern_fit <- function(x, y) {
  browser()
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


# method specific

process_ladder_signal = function(ladder,
                                 scans,
                                 spike_location){
  ladder_df <- data.frame(signal = ladder)
  ladder_df$scan <- scans

  if(is.null(spike_location)){
    spike_location <- which.max(ladder_df$signal) + 45
  }

  ladder_df <- ladder_df[spike_location:nrow(ladder_df),]

  ladder_df$detrended_signal <- detrend_signal(ladder_df$signal)
  ladder_df$smoothed_signal <- moving_average(ladder_df$detrended_signal)
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

    ladder_peak_threshold = ladder_peak_threshold - 0.05
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
    spike_location = spike_location)

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









