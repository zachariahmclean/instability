######################## Helper functions ####################################


deshoulder <- function(peak_table_df, shoulder_window) {
  fragment_size <- peak_table_df$size
  heights <- peak_table_df$height
  fragment_size_deshoulder <- numeric()
  for (i in seq_along(fragment_size)) {
    if (i == 1) {
      if (fragment_size[i + 1] - fragment_size[i] > shoulder_window) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else if (heights[i] > heights[i + 1]) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else {
        next
      }
    } else if (i == length(fragment_size)) {
      if (fragment_size[i - 1] - fragment_size[i] > shoulder_window) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else if (heights[i] > heights[i - 1]) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else {
        next
      }
    } else if (fragment_size[i + 1] - fragment_size[i] < shoulder_window | fragment_size[i] - fragment_size[i - 1] < shoulder_window) {
      before <- fragment_size[i + 1] - fragment_size[i] < shoulder_window
      after <- fragment_size[i] - fragment_size[i - 1] < shoulder_window
      before_and_higher <- heights[i] > heights[i + 1]
      after_and_higher <- heights[i] > heights[i - 1]

      if (before & after) {
        if (before_and_higher & after_and_higher) {
          fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
        } else {
          next
        }
      } else if (before && before_and_higher) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else if (after && after_and_higher) {
        fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
      } else {
        next
      }
    } else {
      fragment_size_deshoulder <- append(fragment_size_deshoulder, fragment_size[i])
    }
  }

  new_peak_table_df <- peak_table_df[which(peak_table_df$size %in% fragment_size_deshoulder), ]

  return(new_peak_table_df)
}






# force whole repeat unit algorithm -----------------------------------------------------------
np_repeat <- function(size,
                      main_peak_size,
                      main_peak_repeat,
                      repeat_size) {
  # For loop to get values
  # first find the distance to main peak and the peak before
  size_delta_from_main_peak <- size - main_peak_size
  # create a numeric vector of length of peak to store values
  np_repeat <- vector(mode = "numeric", length = length(size))
  # note that the loop has to start in the middle at the main peak since it's looking for the previous value in the new vector
  for (i in seq_along(size)) {
    # set main peak
    if (i == which(size_delta_from_main_peak == 0)) {
      np_repeat[[i]] <- main_peak_repeat
    }
    # calculate for peaks greater than main peak
    if (i > which(size_delta_from_main_peak == 0)) {
      # calculate size distance to nearest peak in cag length,
      # add that on to previous cag, then round to whole cag
      np_repeat[[i]] <-
        np_repeat[[i - 1]] + round((size[[i]] - size[[i - 1]]) / repeat_size)
    }
  }
  # --- reverse order and find smaller repeats ---
  size_rev <- rev(size)
  np_repeat_rev <- rev(np_repeat)
  size_delta_from_main_peak_rev <- rev(size_delta_from_main_peak)
  # second loop that goes along the vector in the from larger to smaller bp peaks
  for (i in seq_along(size_rev)) {
    # calculate peaks greater than main peak (smaller bp peaks since reversed)
    if (i > which(size_delta_from_main_peak_rev == 0)) {
      # note this is a subtraction since the size_delta_peak_after is a negative value
      np_repeat_rev[[i]] <-
        np_repeat_rev[[i - 1]] + round((size_rev[[i]] - size_rev[[i - 1]]) / repeat_size)
    }
  }
  return(rev(np_repeat_rev))
}


# fft algo ----------------------------------------------------------------


####### fft helpers
fragments_fft = function(df){
  Fs <- 1  # Arbitrary sampling frequency
  Ts <- 1 / Fs  # Arbitrary sampling time

  signal_periodic <- as.vector(pracma::detrend(df$signal))

  n <- length(signal_periodic)
  xf <- seq(0, 1 / (2 * Ts), length.out = n %/% 2)  # int div in R is % / %
  yf <- stats::fft(signal_periodic)
  yf <- (2 / n) * abs(yf[seq(1, n %/% 2)])

  data.frame("xf" = xf,
             "yf" = yf)
}


find_main_freq = function(fft_df, skip_rows = 3){
  #skip rows to avoid noise at start
  df2 <- fft_df[skip_rows+1:nrow(fft_df), ]
  peaks <- pracma::findpeaks(df2$yf)
  fund_freq_position <- peaks[which(peaks[,1] == max(peaks[,1], na.rm = TRUE)), 2][1]
  freq <- df2[fund_freq_position, "xf"] * 3

  return(freq)
}

find_scan_period <- function(df,
                             main_peak_scan){
  fft_df <- fragments_fft(df)

  #detrend signal
  fft_df$yf <- detrend_signal(fft_df$yf)

  main_freq <- find_main_freq(fft_df)
  pure_wave <- cos(2 * main_freq * (df$scan - main_peak_scan))
  cos_max <- pracma::findpeaks(pure_wave, nups = 3)
  scans_diffs <- diff(df[cos_max[,2], "scan"])
  peak_scan_peroid <- round(median(scans_diffs))


  return(peak_scan_peroid)
}

find_peaks_by_scan_period <- function(df,
                                      main_peak_scan,
                                      peak_scan_period,
                                      direction,
                                      window){

  if(direction == 1){
    df_post_main <- df[which(df$scan > main_peak_scan), ]

  } else{
    df_post_main <- df[which(df$scan < main_peak_scan), ]
    df_post_main <- df_post_main[order(df_post_main$scan, decreasing = TRUE), ]
    peak_scan_period = peak_scan_period * -1
  }

  called_peaks <- numeric()
  current_scan_position <- main_peak_scan + peak_scan_period
  while (TRUE) {

    window_range <- (current_scan_position - window):(current_scan_position + window)
    window_df <- df_post_main[df_post_main$scan %in% window_range, ]

    if (nrow(window_df) > 0) {
      tallest_in_window <- window_df[which.max(window_df$signal), "scan"]
      called_peaks <- c(called_peaks, tallest_in_window)

      # Update current scan position
      current_scan_position <- tallest_in_window + peak_scan_period
    } else {
      # If no more data points, terminate
      break
    }
  }

  return(called_peaks)
}


fft_repeat_caller <- function(fragments_repeat,
                              scan_peak_window = 3,
                              fragment_window = 3*5){

  if(is.na(fragments_repeat$allele_1_size)){
   df <- data.frame("unique_id" = character(),
                    "scan" = numeric(),
                    "size" = numeric(),
                    "signal" = numeric())

    return(df)
  }

  fragment_window_positions <- which(fragments_repeat$trace_bp_df$size > fragments_repeat$allele_1_size - fragment_window & fragments_repeat$trace_bp_df$size < fragments_repeat$allele_1_size + fragment_window )
  window_df <- fragments_repeat$trace_bp_df[fragment_window_positions, ]
  main_peak_scan <- window_df[which(window_df$size == fragments_repeat$allele_1_size), "scan"]

  peak_scan_peroid <- find_scan_period(window_df, main_peak_scan)

  pos_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
                                         main_peak_scan,
                                         peak_scan_peroid,
                                         direction = 1,
                                         window = scan_peak_window)

  neg_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
                                         main_peak_scan,
                                         peak_scan_peroid,
                                         direction = -1,
                                         window = scan_peak_window)

  peak_table <- fragments_repeat$trace_bp_df
  peak_table <- peak_table[which(peak_table$scan %in% c(neg_peaks, main_peak_scan, pos_peaks)), ]
  peak_table <- peak_table[which(peak_table$size > fragments_repeat$.__enclos_env__$private$min_bp_size & peak_table$size < fragments_repeat$.__enclos_env__$private$max_bp_size), ]


  return(peak_table)


}






size_period_repeat_caller <- function(fragments_repeat,
                                      size_period,
                                      scan_peak_window = 3,
                                      fragment_window = 3*5){

  if(is.na(fragments_repeat$allele_1_size)){
    df <- data.frame("unique_id" = character(),
                     "scan" = numeric(),
                     "size" = numeric(),
                     "signal" = numeric())

    return(df)
  }

  fragment_window_positions <- which(fragments_repeat$trace_bp_df$size > fragments_repeat$allele_1_size - fragment_window & fragments_repeat$trace_bp_df$size < fragments_repeat$allele_1_size + fragment_window )
  window_df <- fragments_repeat$trace_bp_df[fragment_window_positions, ]
  main_peak_scan <- window_df[which(window_df$size == fragments_repeat$allele_1_size), "scan"]


  # determine period

  peak_scan_peroid <- round(size_period / median(diff(window_df$size)))

  pos_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
                                         main_peak_scan,
                                         peak_scan_peroid,
                                         direction = 1,
                                         window = scan_peak_window)

  neg_peaks <- find_peaks_by_scan_period(fragments_repeat$trace_bp_df,
                                         main_peak_scan,
                                         peak_scan_peroid,
                                         direction = -1,
                                         window = scan_peak_window)

  peak_table <- fragments_repeat$trace_bp_df
  peak_table <- peak_table[which(peak_table$scan %in% c(neg_peaks, main_peak_scan, pos_peaks)), ]
  peak_table <- peak_table[which(peak_table$size > fragments_repeat$.__enclos_env__$private$min_bp_size & peak_table$size < fragments_repeat$.__enclos_env__$private$max_bp_size), ]


  return(peak_table)


}













# repeat length correction -------------------------------------------------------

model_repeat_length <- function(fragments_list,
                                repeat_size,
                                assay_size_without_repeat,
                                repeat_length_correction) {
  calling_close_neighbouring_repeats <- function(controls_fragments) {
    # use np_repeat to accurately call the repeat length of the neighboring peaks
    # extract a dataframe of the called repeats that can then be used to make a model
    controls_fragments_df_list <- lapply(controls_fragments, function(x) {
      df_length <- nrow(x$peak_table_df)
      # identify peaks close to modal peak and at least 20% as high
      main_peak_delta <- x$peak_table_df$size - x$allele_1_size
      height_prop <- x$peak_table_df$height / x$allele_1_height
      peak_cluster <- vector("logical", length = nrow(x$peak_table_df))
      for (i in seq_along(main_peak_delta)) {
        if (abs(main_peak_delta[[i]]) < 30 & height_prop[[i]] > 0.2) {
          peak_cluster[[i]] <- TRUE
        } else {
          peak_cluster[[i]] <- FALSE
        }
      }
      cluster_df <- x$peak_table_df[peak_cluster, ]
      cluster_df_length <- nrow(cluster_df)
      # use np_repeat method to accurately call the neighboring repeats
      data.frame(
        unique_id = rep(x$unique_id, cluster_df_length),
        size = cluster_df$size,
        validated_repeats = np_repeat(
          size = cluster_df$size,
          main_peak_size = x$allele_1_size,
          main_peak_repeat = x$size_standard_repeat_length,
          repeat_size = repeat_size
        ),
        height = cluster_df$height,
        plate_id = rep(x$plate_id, cluster_df_length)
      )
    })

    controls_repeats_df <- do.call(rbind, controls_fragments_df_list)
  }

  # correct repeats
  if (repeat_length_correction == "from_metadata") {
    ## first pull out a dataframe for all samples with a column that indicates if it's a positive control or not
    extracted <- lapply(fragments_list, function(x) {
      data.frame(
        unique_id = x$unique_id,
        size_standard = x$size_standard,
        allele_1_size = x$allele_1_size,
        plate_id = x$plate_id,
        size_standard_sample_id = x$size_standard_sample_id
      )
    })
    extracted_df <- do.call(rbind, extracted)

    # Check to see if there are controls, if there are none, give error
    if (!any(extracted_df$size_standard == TRUE)) {
      stop("No repeat-length control samples were detected. Ensure that the metadata has been added to the samples with 'add_metadata()' and check your metadata to make sure 'TRUE' is indicated in the appropriate column to indicate samples that are to be used for predicting the repeat length",
        call. = FALSE
      )
    }
    # pull out the controls
    controls_df <- extracted_df[which(extracted_df$size_standard == TRUE), , drop = FALSE]
    controls_fragments <- fragments_list[which(names(fragments_list) %in% controls_df$unique_id)]
    controls_repeats_df <- calling_close_neighbouring_repeats(controls_fragments)
  } else if (repeat_length_correction == "from_genemapper") {
    # Do some checks by identify samples that have genemapper called alleles
    controls_samples_list <- lapply(fragments_list, function(x) x$peak_table_df[which(!is.na(x$peak_table_df$allele)), ])
    controls_samples_df <- do.call(rbind, controls_samples_list)
    if (nrow(controls_samples_df) == 0) {
      stop(paste("Correction could not go ahead because no genemapper alleles could be indetified"),
        call. = FALSE
      )
    }

    # pick the closest peak to the main peak size and temporarily make that allele_1
    controls_fragments <- lapply(
      fragments_list[which(names(fragments_list) %in% unique(controls_samples_df$unique_id))],
      function(x) {
        genemapper_alleles <- controls_samples_df[which(controls_samples_df$unique_id == x$unique_id), ]
        allele_1_delta_abs <- abs(genemapper_alleles$size - x$allele_1_size)
        closest_to_allele_1 <- which(allele_1_delta_abs == min(allele_1_delta_abs))
        selected_genemapper_allele <- genemapper_alleles[closest_to_allele_1[1], ]

        # make sure it doesn't modify in place and mess up the selection of the real main peak
        y <- x$clone()
        y$allele_1_size <- selected_genemapper_allele$size
        y$size_standard <- TRUE
        y$size_standard_repeat_length <- selected_genemapper_allele$allele
        return(y)
      }
    )

    controls_repeats_df <- calling_close_neighbouring_repeats(controls_fragments)
  }

  # Check to see if there are controls for each plate, if there are no controls for a plate, give error
  all_plate_ids <- lapply(fragments_list, function(x) x$plate_id)
  control_plate_ids <- unique(controls_repeats_df$plate_id)
  if (length(unique(control_plate_ids)) != length(unique(all_plate_ids))) {
    plates_missing_controls <- paste0(all_plate_ids[which(!all_plate_ids %in% control_plate_ids)], collapse = ", ")
    stop(paste("Plate(s)", plates_missing_controls, "have no repeat-length control samples"),
      call. = FALSE
    )
  }

  # identify size stds with shared id for more quality control
  # can compare to each other to make sure that the same peak in the distribution has been selected as the modal
  standard_sample_ids <-sapply(controls_fragments, function(x) x$size_standard_sample_id)
  if(any(!is.na(standard_sample_ids))){
    unique_standard_sample_ids <- unique(standard_sample_ids)

    for (i in seq_along(unique_standard_sample_ids)) {
      ids_i <- controls_df[which(controls_df$size_standard_sample_id == unique_standard_sample_ids[i]), "unique_id"]
      controls_repeats_df_i <- controls_repeats_df[which(controls_repeats_df$unique_id %in% ids_i), ]
      controls_repeats_df_list_i <- split(controls_repeats_df_i, controls_repeats_df_i$unique_id)
      differences_i <- sapply(controls_repeats_df_list_i, function(x) {
        # this compares the average (weighted mean) and the mode
        # if there's a little wobble at the top with two cloes peaks
        # then the average size shouldn't change much
        # but the difference between the mode and the average changes a whole repeat unit
        # which can indicate to us that one of the stds might be off
        weighted.mean(x$size, x$height) - x$size[which.max(x$height)]
      })
      differences_of_differences_i <- lapply(differences_i, function(x) differences <- x - differences_i)
      if(any(sapply(differences_of_differences_i, function(x) abs(x) > repeat_size * 0.8))){
        # how do you dtermine which of the samples might be off? 
        # perhaps we could try and help the person figure that out
        # but that is complicated, instead give them a warning
        warning(
          call. = FALSE, 
          paste0("Warning! It looks like at least one of the samples in the size standard group '",
          unique_standard_sample_ids[i], "' has a different modal peak than the other samples. ",
          "It's possible that the modal peak has shifted to a different spot in the distribution in at least of one the runs. ",
          "Use plot_size_standard_samples() to visualize and identify this sample, the update metadata with the correct repeat length of the modal peak."
        )
        )
      }
    }
  }
  
  message(paste0("Repeat correction model: ", length(unique(controls_repeats_df$unique_id)), " samples used to build model"))

  # Can now make a model based on the bp size and the known repeat size
  if (length(unique(controls_repeats_df$plate_id)) == 1) {
    # when there's only one plate just set up simple lm
    correction_mods <- stats::lm(validated_repeats ~ size, data = controls_repeats_df)
    repeat_bp_size <- round(1 / correction_mods$coefficients[2], 2)
    message(paste0("Repeat correction model: ", repeat_bp_size, " bp increase per repeat"))
  } else {
    # when there are multiple samples a linear model can be made using the modal peak and the known repeat length of the modal peak
    correction_mods <- lm(validated_repeats ~ size * plate_id, data = controls_repeats_df)
  }

  # check to see if any samples look off

  controls_repeats_df$predicted_repeat <- stats::predict.lm(correction_mods, controls_repeats_df)
  controls_repeats_df$residuals <- correction_mods$residuals
  message(paste0("Repeat correction model: Average repeat residual ", round(mean(controls_repeats_df$residuals), 10)))

  if (any(abs(controls_repeats_df$residuals) > 0.3)) {
    message("Repeat correction model: Warning! The following samples may be off and need investigaion. It's possible that at least one of these samples has the incorrect repeat length indicated in the metadata.")

    samples_all_controls <- unique(controls_repeats_df$unique_id)
    samples_high_diff <- unique(controls_repeats_df[which(abs(controls_repeats_df$residuals) > 0.5), "unique_id"])
    for (i in seq_along(samples_high_diff)) {
      sample_id <- samples_high_diff[i]
      sample_control_df <- controls_repeats_df[which(controls_repeats_df$unique_id == sample_id), ]
      sample_control_peaks_n <- nrow(sample_control_df)
      sample_control_peaks_off_df <- sample_control_df[which(abs(sample_control_df$residuals) > 0.5), ]
      sample_control_peaks_off_df_n <- nrow(sample_control_peaks_off_df)

      message(paste0(sample_id, " has ", sample_control_peaks_off_df_n, "/", sample_control_peaks_n, " peaks used for making model with high residual repeat size (average residual ", round(mean(sample_control_df$residuals), 2), " repeats)"))
    }
  }
  return(list(correction_mods = correction_mods, controls_repeats_df = controls_repeats_df))
}


##################### R6 Class Method Helpers ##################################
# bp_fragments
add_repeats_helper <- function(fragments_repeats,
                               assay_size_without_repeat,
                               repeat_size,
                               repeat_calling_algorithm,
                               repeat_calling_algorithm_size_window_around_allele,
                               repeat_calling_algorithm_peak_assignment_scan_window,
                               repeat_calling_algorithm_size_period,
                               force_whole_repeat_units,
                               correct_repeat_length) {
  ##
  ### in this function, need to do the following in the correct order
  #### 1) calculate repeats by assay_size_without_repeat and repeat size
  #### 2) correct repeat length with positive control
  #### 3) use a method to calculate repeats

  # check to make sure all the required inputs for the function have been given
  if (fragments_repeats$.__enclos_env__$private$find_main_peaks_used == FALSE) {
    stop(paste0(fragments_repeats$unique_id, " requires main alleles to be identified before repeats can be called. Find alleles using 'find_main_peaks()' whitin the class, or use the 'find_alleles()' accesesor to find the main peaks across a list of 'fragments_repeats' objects"),
      call. = FALSE
    )
  } else if (correct_repeat_length == TRUE & is.null(fragments_repeats$.__enclos_env__$private$correction_mod)) {
    stop("Correcting the repeat length requires a model based on positive controls, so 'correct_repeat_length' & 'correction_mod' inputs are not meant for users to directly use. To correct the repeat length, you need to work on the 'fragments_repeats' objects in a list format and use accessor functions. On a list of 'fragments_repeats' objects, i) use 'add_metadata()' to indicate which samples are positive controls, and ii) use 'find_alleles()' accesesor function to call and correct repeat lengths across all samples",
      call. = FALSE
    )
  }

  # only continue from here if main peaks were successfully found, otherwise, don't return repeat data (ie it can be an empty df)
  repeat_class <- fragments_repeats$clone()

  if (is.na(fragments_repeats$allele_1_size) | is.na(fragments_repeats$allele_1_height)) {
    repeat_class$.__enclos_env__$private$repeats_not_called_reason <- "No main peaks"
    warning(paste0(repeat_class$unique_id, ": repeats were not called (no main peaks in sample)"),
      call. = FALSE
    )
    # populate with empty dataframe to help the rest of the pipeline
    repeat_class$repeat_table_df <- data.frame(
      unique_id = character(),
      size = numeric(),
      height = numeric(),
      repeats = numeric(),
      off_scale = logical()
    )

  } else {


    # repeat calling algorithm
    if(repeat_calling_algorithm == "simple"){
      repeat_table_df <- data.frame(
        unique_id = fragments_repeats$peak_table_df$unique_id,
        size = fragments_repeats$peak_table_df$size,
        height = fragments_repeats$peak_table_df$height,
        calculated_repeats = (fragments_repeats$peak_table_df$size - assay_size_without_repeat) / repeat_size,
        repeats = (fragments_repeats$peak_table_df$size - assay_size_without_repeat) / repeat_size,
        off_scale = ifelse(any(colnames(fragments_repeats$peak_table_df) == "off_scale"),
                           fragments_repeats$peak_table_df$off_scale,
                           rep(FALSE, nrow(fragments_repeats$peak_table_df))
        )
      )

    } else if(repeat_calling_algorithm == "fft"){

      #check to see that fragments repeats has trace data since that is required.
      if(is.null(fragments_repeats$trace_bp_df)){
        stop("fft algorithim requires trace data. Use fsa samples rather than peak table is inputs into the pipeline.",
             call. = FALSE)
      }

      fft_peak_df <- fft_repeat_caller(fragments_repeats,
                        fragment_window = repeat_calling_algorithm_size_window_around_allele,
                        scan_peak_window = repeat_calling_algorithm_peak_assignment_scan_window)

      repeat_table_df <- data.frame(
        unique_id = fft_peak_df$unique_id,
        size = fft_peak_df$size,
        height = fft_peak_df$signal,
        calculated_repeats = (fft_peak_df$size - assay_size_without_repeat) / repeat_size,
        repeats = (fft_peak_df$size - assay_size_without_repeat) / repeat_size,
        off_scale = fft_peak_df$off_scale
      )

    }
    else if(repeat_calling_algorithm == "size_period"){

      #check to see that fragments repeats has trace data since that is required.
      if(is.null(fragments_repeats$trace_bp_df)){
        stop("size_period algorithim requires trace data. Use fsa samples rather than peak table is inputs into the pipeline.",
             call. = FALSE)
      }

      size_period_df <- size_period_repeat_caller(fragments_repeats,
                                       size_period = repeat_calling_algorithm_size_period,
                                       fragment_window = repeat_calling_algorithm_size_window_around_allele,
                                       scan_peak_window = repeat_calling_algorithm_peak_assignment_scan_window)

      repeat_table_df <- data.frame(
        unique_id = size_period_df$unique_id,
        size = size_period_df$size,
        height = size_period_df$signal,
        calculated_repeats = (size_period_df$size - assay_size_without_repeat) / repeat_size,
        repeats = (size_period_df$size - assay_size_without_repeat) / repeat_size,
        off_scale = size_period_df$off_scale
      )

    }
    else{
      stop(call. = FALSE,
           "Invalid repeat calling algorithim selected")
    }


    # Correct repeat length with positive controls
    if (correct_repeat_length == TRUE) {

      repeat_table_df$plate_id <- rep(fragments_repeats$plate_id, nrow(repeat_table_df))
      repeat_table_df$calculated_repeats <- stats::predict.lm(fragments_repeats$.__enclos_env__$private$correction_mod, repeat_table_df)
      repeat_table_df$repeats <- repeat_table_df$calculated_repeats
    }

    # Force the repeat units to be whole numbers
    if (force_whole_repeat_units == TRUE) {
      repeat_table_df$repeats <- np_repeat(
        size = repeat_table_df$size,
        main_peak_size = fragments_repeats$allele_1_size,
        main_peak_repeat = repeat_table_df$calculated_repeats[which(repeat_table_df$size == fragments_repeats$allele_1_size)],
        repeat_size = repeat_size
      )
    }

    # Finally save main peak repeat length and repeats data
    repeat_class$allele_1_repeat <- repeat_table_df$repeats[which(repeat_table_df$size == fragments_repeats$allele_1_size)]
    repeat_class$allele_2_repeat <- repeat_table_df$repeats[which(repeat_table_df$size == fragments_repeats$allele_2_size)]
    repeat_class$repeat_table_df <- repeat_table_df
  }


  # also calculate repeat length for the trace-level data if it exists
  if(!is.null(repeat_class$trace_bp_df)){
    if(correct_repeat_length == TRUE){

      fragments_repeats$trace_bp_df$plate_id <- rep(fragments_repeats$plate_id, nrow(fragments_repeats$trace_bp_df))
      repeat_class$trace_bp_df$calculated_repeats <- stats::predict.lm(fragments_repeats$.__enclos_env__$private$correction_mod, fragments_repeats$trace_bp_df)
    }
    else{
      repeat_class$trace_bp_df$calculated_repeats <- (fragments_repeats$trace_bp_df$size - assay_size_without_repeat) / repeat_size
    }
  }

  return(repeat_class)
}



# call_repeats ------------------------------------------------------------

#' Call Repeats for Fragments
#'
#' This function calls the repeat lengths for a list of fragments.
#'
#' @param fragments_list A list of fragments_repeats objects containing fragment data.
#' @param assay_size_without_repeat An integer specifying the assay size without repeat for repeat calling. Default is 87.
#' @param repeat_size An integer specifying the repeat size for repeat calling. Default is 3.
#' @param force_whole_repeat_units A logical value specifying if the peaks should be forced to be whole repeat units apart. Usually the peaks are slightly under the whole repeat unit if left unchanged.
#' @param repeat_length_correction A character specifying the repeat length correction method. Options: \code{"none"}, \code{"from_metadata"}, \code{"from_genemapper"}. Default is \code{"none"}.
#' @param repeat_calling_algorithm A character specifying the repeat calling algorithm. Options: \code{"simple"}, \code{"fft"}, or \code{"size_period"} (see details section for more information on these).
#' @param repeat_calling_algorithm_size_window_around_allele A numeric value for how big of a window around the tallest peak should be used to find the peak periodicity. Used for both \code{"fft"} and \code{"size_period"}. For \code{"fft"}, you want to make sure that this window is limited to where there are clear peaks. For \code{"size_period"}, it will not make a big difference.
#' @param repeat_calling_algorithm_peak_assignment_scan_window A numeric value for the scan window when assigning the peak. This is used for both \code{"fft"} and \code{"size_period"}. When the scan period is determined, the algorithm jumps to the predicted scan for the next peak. This value opens a window of the neighboring scans to pick the tallest in.
#' @param repeat_calling_algorithm_size_period A numeric value \code{"size_period"} algorithm to set the peak periodicity by bp size. This is the key variable to change for \code{"size_period"}. In fragment analysis, the peaks are usually slightly below the actual repeat unit size.
#'
#' @return A list of \code{"fragments_repeats"} objects with repeat data added.
#'
#' @details
#' The calculated repeat lengths are assigned to the corresponding peaks in the provided `fragments_repeats` object. The repeat lengths can be used for downstream instability analysis.
#'
#' The `simple` algorithm is just the repeat size calculated either directly, or when size standards are used to correct the repeat, it's the repeat length calculated from the model of bp vs repeat length.
#'
#' The `fft` or `size_period` algorithms both re-call the peaks based on empirically determined (`fft`) or specified (`size_period`) periodicity of the peaks. The main application of these algorithms is to solve the issue of contaminating peaks that are not expected in the expected regular pattern of peaks. The `fft` approach applies a fourier transform to the peak signal to determine the underlying periodicity of the signal. `size_period` is similar and simpler, where instead of automatically figuring out the periodicity we as users usually know the size distance between repeat units. We can use that known peroidicty to jump between peaks.
#'
#' The `force_whole_repeat_units` algorithm aims to correct for the systematic drift in fragment sizes that occurs. It calculates repeat lengths in a way that helps align peaks with the underlying repeat pattern, making the estimation of repeat lengths more reliable relative to the main peak. The calculated repeat lengths start from the main peak's repeat length and increases in increments of the specified `repeat_size`.
#'
#' @seealso [find_alleles()]
#'
#' @export
#'
#' @examples
#'
#' file_list <- instability::cell_line_fsa_list[c(90:92)]
#'
#' test_ladders <- find_ladders(file_list)
#'
#' fragments_list <- find_fragments(test_ladders,
#'   min_bp_size = 300
#' )
#'
#' test_alleles <- find_alleles(
#'   fragments_list = fragments_list
#' )
#'
#' # Simple conversion from bp size to repeat size
#' test_repeats <- call_repeats(
#'   fragments_list = test_alleles,
#'   repeat_calling_algorithm = "simple",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(test_repeats[1], xlim = c(120,170))
#'
#'
#' # use different algorithms to call the repeats to ensure only periodic peaks are called
#'
#' # fft to automatically find peak period
#' test_repeats_fft <- call_repeats(
#'   fragments_list = test_alleles,
#'   repeat_calling_algorithm = "fft",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(test_repeats_fft[1], xlim = c(120,170))
#'
#' # size_period to manually supply the peak period
#' test_repeats_size_period <- call_repeats(
#'   fragments_list = test_alleles,
#'   repeat_calling_algorithm = "size_period",
#'   repeat_calling_algorithm_size_period = 2.75,
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(test_repeats_size_period[1], xlim = c(120,170))
#'
#'
#' # Use force_whole_repeat_units algorithm to make sure called
#' # repeats are the exact number of bp apart
#'
#' test_repeats_whole_units <- call_repeats(
#'   fragments_list = test_alleles,
#'   force_whole_repeat_units = TRUE,
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#' plot_traces(test_repeats_whole_units[1], xlim = c(120,170))
#'
#' # correct repeat length from metadata
#'
#' test_alleles_metadata <- add_metadata(
#'   fragments_list = test_alleles,
#'   metadata_data.frame = instability::metadata
#' )
#'
#' test_repeats_corrected <- call_repeats(
#'   fragments_list = test_alleles_metadata,
#'   repeat_length_correction = "from_metadata"
#' )
#'
#'
#' plot_traces(test_repeats_corrected[1], xlim = c(120,170))
#'
#'
#'
#'
call_repeats <- function(
    fragments_list,
    assay_size_without_repeat = 87,
    repeat_size = 3,
    force_whole_repeat_units = FALSE,
    repeat_length_correction = "none",
    repeat_calling_algorithm = "simple",
    repeat_calling_algorithm_size_window_around_allele = repeat_size * 5,
    repeat_calling_algorithm_peak_assignment_scan_window = 3,
    repeat_calling_algorithm_size_period = repeat_size * 0.93) {
  # Check to see if repeats are to be corrected
  # if so, supply the model to each of the samples in the list
  if (repeat_length_correction %in% c("from_metadata", "from_genemapper")) {
    mod <- model_repeat_length(
      fragments_list = fragments_list,
      assay_size_without_repeat = assay_size_without_repeat,
      repeat_size = repeat_size,
      repeat_length_correction = repeat_length_correction
    )

    for (i in seq_along(fragments_list)) {
      fragments_list[[i]]$.__enclos_env__$private$correction_mod <- mod$correction_mods
      fragments_list[[i]]$.__enclos_env__$private$controls_repeats_df <- mod$controls_repeats_df
    }
  }

  # call repeats for each sample
  added_repeats <- lapply(
    fragments_list,
    function(x) {
      x <- add_repeats_helper(
        x,
        assay_size_without_repeat = assay_size_without_repeat,
        repeat_size = repeat_size,
        repeat_calling_algorithm = repeat_calling_algorithm,
        repeat_calling_algorithm_size_window_around_allele = repeat_calling_algorithm_size_window_around_allele,
        repeat_calling_algorithm_peak_assignment_scan_window = repeat_calling_algorithm_peak_assignment_scan_window,
        repeat_calling_algorithm_size_period = repeat_calling_algorithm_size_period,
        force_whole_repeat_units = force_whole_repeat_units,
        correct_repeat_length = ifelse(repeat_length_correction == "none", FALSE, TRUE)
      )

      return(x)
    }
  )

  return(added_repeats)
}





