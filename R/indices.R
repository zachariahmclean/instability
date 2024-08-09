######################## Helper functions ####################################
# instability index ---------------------------------------------------------
instability_index <- function(repeats,
                              heights,
                              index_peak_height,
                              index_peak_repeat,
                              peak_threshold,
                              abs_sum = FALSE) {
  # apply height threshold
  peak_over_threshold <- which(heights / index_peak_height > peak_threshold)
  repeats <- repeats[peak_over_threshold]
  heights <- heights[peak_over_threshold]

  # normalised peak height
  heights_normalised <- heights / sum(heights)

  # distance to index peak
  repeat_delta <- repeats - index_peak_repeat
  if (abs_sum == FALSE) {
    sum(heights_normalised * repeat_delta)
  } else if (abs_sum == TRUE) {
    sum(abs(heights_normalised * repeat_delta))
  }
}

# function for finding quantiles -----------------------------------------------

find_percentiles <- function(repeats,
                             heights,
                             index_peak_repeat,
                             type, # "percentile" or "repeat"
                             range,
                             col_preffix) {
  # if there are double peaks select the tallest of the peaks, otherwise approx interpolation doesn't work
  # also if there are no main peak caled filter out (for example samples used to call repeats but irrelevat for metrics)
  df_names <- paste(col_preffix, range, sep = "_")

  # Deal with case when there are no expansion peaks by returning 0s
  if (sum(repeats > index_peak_repeat) <= 1) {
    percentile_df <- as.data.frame(setNames(as.list(rep(0, length(range))), df_names))
  } else {
    unique_repeat_df <- aggregate(heights ~ repeats, FUN = max)
    cumsum_pct <- cumsum(unique_repeat_df$heights) / sum(unique_repeat_df$heights)
    repeat_delta <- unique_repeat_df$repeats - index_peak_repeat

    values <- vector("numeric", length(range))

    if (type == "percentile") {
      for (i in seq_along(range)) {
        values[[i]] <- approx(cumsum_pct,
          repeat_delta,
          xout = range[[i]],
          yleft = min(repeat_delta)
        )$y
      }
    } else if (type == "repeat") {
      for (i in seq_along(range)) {
        values[[i]] <- approx(repeat_delta,
          cumsum_pct,
          xout = range[[i]],
          yleft = min(cumsum_pct)
        )$y
      }
    }

    percentile_df <- as.data.frame(setNames(as.list(values), df_names))
  }

  return(percentile_df)
}


# skewness ------------------------------------------------------------------

fishers_skewness <- function(x, y) {
  mean_val <- sum(x * y)
  sd_val <- sqrt(sum(y * (x - mean_val)^2))

  skewness <- sum(y * (x - mean_val)^3) / sd_val^3

  return(skewness)
}


# kurtosis -----------------------------------------------------------------

fishers_kurtosis <- function(x, y) {
  mean_val <- sum(x * y)
  sd_val <- sqrt(sum(y * (x - mean_val)^2))

  kurtosis <- (sum(y * (x - mean_val)^4) / sd_val^4) - 3
  return(kurtosis)
}

# subsetting repeat table ---------------------------------------------------

repeat_table_subset <- function(repeat_table_df,
                                allele_1_height,
                                index_repeat,
                                peak_threshold,
                                window_around_main_peak) {
  # Filter to include only the peaks above the certain threshold
  # height threshold is set on the modal peak rather than the index peak
  repeat_table_df$peak_percent <- repeat_table_df$height / allele_1_height
  height_filtered_df <- repeat_table_df[which(repeat_table_df$peak_percent > peak_threshold), ]

  # Ensure window_around_main_peak is exactly length 2
  if (length(window_around_main_peak) != 2) {
    stop("window_around_main_peak must be a vector of length 2")
  }

  # Filter to include only peaks of a certain size
  lower_lim <- ifelse(is.na(window_around_main_peak[1]),
    min(height_filtered_df$repeats),
    index_repeat - abs(window_around_main_peak[1])
  )
  upper_lim <- ifelse(is.na(window_around_main_peak[1]),
    max(height_filtered_df$repeats),
    index_repeat + abs(window_around_main_peak[2])
  )
  size_filtered_df <- height_filtered_df[which(height_filtered_df$repeats >= lower_lim & height_filtered_df$repeats <= upper_lim), ]

  return(size_filtered_df)
}


###################### R6 Class Method Helpers ####################################

# Calculating metrics --------------------------------------------------------

compute_metrics <- function(fragments_repeats,
                            peak_threshold,
                            window_around_main_peak,
                            percentile_range,
                            repeat_range) {
  # filter dataset to user supplied thresholds
  size_filtered_df <- repeat_table_subset(
    repeat_table_df = fragments_repeats$repeat_table_df,
    allele_1_height = fragments_repeats$allele_1_height,
    index_repeat = fragments_repeats$index_repeat,
    peak_threshold = peak_threshold,
    window_around_main_peak = window_around_main_peak
  )

  if(!is.null(fragments_repeats$.__enclos_env__$private$index_samples)){

    control_weighted_mean_repeat <- sapply(fragments_repeats$.__enclos_env__$private$index_samples, function(x){
      control_filtered_df <- repeat_table_subset(
        repeat_table_df = x[[2]],
        allele_1_height = x[[2]][which(x[[2]]$repeats == x[[1]]), "height"],
        index_repeat = x[[1]],
        peak_threshold = peak_threshold,
        window_around_main_peak = window_around_main_peak
      )

      weighted.mean(control_filtered_df$repeats, control_filtered_df$height)
    })

    index_weighted_mean_repeat <- median(control_weighted_mean_repeat, na.rm = TRUE)
  } else{
    index_weighted_mean_repeat <- NA
  }




  # first subset to make some dataframe that are just for contractions or expansions
  size_filtered_df$repeat_delta_index_peak <- size_filtered_df$repeats - fragments_repeats$index_repeat
  expansion_filtered <- size_filtered_df[which(size_filtered_df$repeat_delta_index_peak >= 0), ]
  contraction_filtered <- size_filtered_df[which(size_filtered_df$repeat_delta_index_peak <= 0), ]

  # QCs
  QC_modal_peak_height <- if (fragments_repeats$allele_1_height > 500) {
    NA_character_
  } else if (fragments_repeats$allele_1_height > 100) {
    "Low"
  } else {
    "Extremly low"
  }

  QC_peak_number <- if (nrow(fragments_repeats$repeat_table_df) > 20) {
    NA_character_
  } else if (nrow(fragments_repeats$repeat_table_df) > 10) {
    "Low"
  } else {
    "Extremly low"
  }

  QC_off_scale <- if (any(fragments_repeats$repeat_table_df$off_scale)) {
    paste(
      "The following repeats were determined off scale (check ladder too, could be scans in any channel):",
      paste(round(fragments_repeats$repeat_table_df[which(fragments_repeats$repeat_table_df$off_scale), "repeats"]), collapse = ", ")
    )
  } else {
    NA_character_
  }

  # make a wide dataframe
  metrics <- data.frame(
    unique_id = fragments_repeats$unique_id,
    QC_comments = NA_character_,
    QC_modal_peak_height = QC_modal_peak_height,
    QC_peak_number = QC_peak_number,
    QC_off_scale = QC_off_scale,
    modal_peak_repeat = fragments_repeats$allele_1_repeat,
    modal_peak_height = fragments_repeats$allele_1_height,
    index_peak_repeat = fragments_repeats$index_repeat,
    index_peak_height = fragments_repeats$index_height,
    index_weighted_mean_repeat = index_weighted_mean_repeat,
    n_peaks_total = nrow(fragments_repeats$repeat_table_df),
    n_peaks_analysis_subset = nrow(size_filtered_df),
    n_peaks_analysis_subset_expansions = nrow(expansion_filtered),
    min_repeat = min(size_filtered_df$repeats),
    max_repeat = max(size_filtered_df$repeats),
    mean_repeat = mean(size_filtered_df$repeats),
    weighted_mean_repeat = weighted.mean(size_filtered_df$repeats, size_filtered_df$height),
    median_repeat = median(size_filtered_df$repeats),
    max_height = max(size_filtered_df$height),
    max_delta_neg = min(size_filtered_df$repeat_delta_index_peak),
    max_delta_pos = max(size_filtered_df$repeat_delta_index_peak),
    skewness = fishers_skewness(size_filtered_df$repeats, size_filtered_df$height),
    kurtosis = fishers_kurtosis(size_filtered_df$repeats, size_filtered_df$height),
    modal_repeat_delta = fragments_repeats$allele_1_repeat - fragments_repeats$index_repeat,
    average_repeat_gain = weighted.mean(size_filtered_df$repeats, size_filtered_df$height) - index_weighted_mean_repeat,
    instability_index = instability_index(
      repeats = size_filtered_df$repeats,
      heights = size_filtered_df$height,
      index_peak_height = fragments_repeats$allele_1_height,
      index_peak_repeat = fragments_repeats$index_repeat,
      peak_threshold = peak_threshold,
      abs_sum = FALSE
    ),
    instability_index_abs = instability_index(
      repeats = size_filtered_df$repeats,
      heights = size_filtered_df$height,
      index_peak_height = fragments_repeats$allele_1_height,
      index_peak_repeat = fragments_repeats$index_repeat,
      peak_threshold = peak_threshold,
      abs_sum = TRUE
    ),
    expansion_index = instability_index(
      repeats = expansion_filtered$repeats,
      heights = expansion_filtered$height,
      index_peak_height = fragments_repeats$allele_1_height,
      index_peak_repeat = fragments_repeats$index_repeat,
      peak_threshold = peak_threshold,
      abs_sum = FALSE
    ),
    contraction_index = instability_index(
      repeats = contraction_filtered$repeats,
      heights = contraction_filtered$height,
      index_peak_height = fragments_repeats$allele_1_height,
      index_peak_repeat = fragments_repeats$index_repeat,
      peak_threshold = peak_threshold,
      abs_sum = FALSE
    ),
    expansion_ratio = sum(expansion_filtered$peak_percent) - 1, # remove the main peak by subtracting 1
    contraction_ratio = sum(contraction_filtered$peak_percent) - 1
  )

  expansion_percentile <- find_percentiles(
    expansion_filtered$repeats,
    expansion_filtered$height,
    fragments_repeats$index_repeat,
    type = "percentile",
    range = percentile_range,
    col_preffix = "expansion_percentile"
  )

  expansion_repeat <- find_percentiles(
    expansion_filtered$repeats,
    expansion_filtered$height,
    fragments_repeats$index_repeat,
    type = "repeat",
    range = repeat_range,
    col_preffix = "expansion_percentile_for_repeat"
  )

  metrics <- cbind(metrics, expansion_percentile)
  metrics <- cbind(metrics, expansion_repeat)


  return(metrics)
}
