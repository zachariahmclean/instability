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
    index_weighted_mean_repeat = fragments_repeats$index_weighted_mean_repeat,
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
    modal_repeat_delta = fragments_repeats$allele_1_repeat - fragments_repeats$index_repeat,
    average_repeat_gain = weighted.mean(size_filtered_df$repeats, size_filtered_df$height) - fragments_repeats$index_weighted_mean_repeat,
    instabity_index_jml = instability_index(
      repeats = size_filtered_df$repeats,
      heights = size_filtered_df$height,
      index_peak_height = fragments_repeats$allele_1_height,
      index_peak_repeat = fragments_repeats$index_repeat,
      peak_threshold = peak_threshold,
      abs_sum = FALSE
    ),
    abs_index = instability_index(
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
    pps_index = sum(expansion_filtered$peak_percent) - 1, # remove the main peak from the pps by subracting 1
    som_mos_index = (sum(expansion_filtered$height) - fragments_repeats$allele_1_height) / fragments_repeats$allele_1_height,
    contration_index = instability_index(
      repeats = contraction_filtered$repeats,
      heights = contraction_filtered$height,
      index_peak_height = fragments_repeats$allele_1_height,
      index_peak_repeat = fragments_repeats$index_repeat,
      peak_threshold = peak_threshold,
      abs_sum = FALSE
    ),
    skewness = fishers_skewness(size_filtered_df$repeats, size_filtered_df$height),
    kurtosis = fishers_kurtosis(size_filtered_df$repeats, size_filtered_df$height)
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

  return(metrics)
}


############################# accessor-function helper ##################################


metrics_grouping_helper <- function(fragments_list,
                                    peak_threshold,
                                    window_around_main_peak) {
  # extract dataframe to identify controls
  sample_list <- lapply(fragments_list, function(x) {
    # we need to do this to find the weighted mean repeat for the index sample
    # therefore only do the computational work if it's actual the index samples
    if (x$metrics_baseline_control == TRUE) {
      filtered_df <- repeat_table_subset(
        repeat_table_df = x$repeat_table_df,
        allele_1_height = x$allele_1_height,
        index_repeat = x$allele_1_repeat,
        peak_threshold = peak_threshold,
        window_around_main_peak = window_around_main_peak
      )
    }

    data.frame(
      unique_id = x$unique_id,
      group_id = x$group_id,
      repeats = x$allele_1_repeat,
      metrics_baseline_control = x$metrics_baseline_control,
      weighted_mean_repeat = ifelse(x$metrics_baseline_control == TRUE,
        weighted.mean(filtered_df$repeats, filtered_df$height),
        NA_real_
      )
    )
  })
  groups_df <- do.call(rbind, sample_list)
  groups_list <- split(groups_df, groups_df$group_id)
  # do some quality control to make sure each group has an index sample
  for (i in seq_along(groups_list)) {
    controls <- groups_list[[i]][which(groups_list[[i]]$metrics_baseline_control == TRUE), ]
    n_controls <- nrow(controls)
    controls_missing_allele <- any(is.na(controls$repeats))
    if (n_controls == 0) {
      stop(paste0("Group '", names(groups_list)[[i]], "' has no 'metrics_baseline_control'"),
        call. = FALSE
      )
    } else if (controls_missing_allele == TRUE) {
      stop(paste0("Group '", names(groups_list)[[i]], "' control has no allele called. Grouped analysis won't work for these samples."),
        call. = FALSE
      )
    } else if (n_controls > 1) {
      message(paste0("Group '", names(groups_list)[[i]], "' has more than one 'metrics_baseline_control'. The median repeat of the assigned samples will be used to assign the index peak"))
    }

    ## pull out the index controls
    group_controls_df <- groups_list[[i]][which(groups_list[[i]]$metrics_baseline_control == TRUE), ]
    groups_list[[i]] <- data.frame(
      group_id = unique(group_controls_df$group_id),
      repeats = median(group_controls_df$repeats),
      weighted_mean_repeat = median(group_controls_df$weighted_mean_repeat),
      metrics_baseline_control = TRUE
    )
  }
  # merge the list back together
  controls_df <- do.call(rbind, groups_list)

  # add the index values to all samples
  fragments_indexed_list <- lapply(fragments_list, function(x) {
    # since the repeat size may not be an integer, need t find what the closest peak is to the control sample
    # delta between repeat of index sample and all repeats of sample
    index_delta <- x$repeat_table_df$repeats - controls_df[which(controls_df$group_id == x$group_id), "repeats"]
    closest_peak <- which(abs(index_delta) == min(abs(index_delta)))

    if (length(closest_peak) == 1) {
      x$index_repeat <- x$repeat_table_df$repeats[closest_peak]
      x$index_height <- x$repeat_table_df$height[closest_peak]
    } else {
      tallest_candidate <- closest_peak[which(x$repeat_table_df$height[closest_peak] == max(x$repeat_table_df$height[closest_peak]))]
      x$index_repeat <- x$repeat_table_df$repeats[tallest_candidate]
      x$index_height <- x$repeat_table_df$height[tallest_candidate]
    }

    x$index_weighted_mean_repeat <- controls_df[which(controls_df$group_id == x$group_id), "weighted_mean_repeat"]

    return(x)
  })

  return(fragments_indexed_list)
}

# metrics_override_helper ---------------------------------------------------------
metrics_override_helper <- function(fragments_list,
                                    index_override_dataframe) {
  lapply(fragments_list, function(x) {
    # if there is nothing to override, then just return the existing index values
    if (!any(index_override_dataframe$unique_id == x$unique_id)) {
      x$index_repeat <- x$index_repeat
      x$index_height <- x$index_height
    } else if (any(index_override_dataframe[, 1] == x$unique_id)) {
      index_delta <- x$repeat_table_df$repeats - index_override_dataframe[which(index_override_dataframe[, 1] == x$unique_id), 2]

      closest_peak <- which(abs(index_delta) == min(abs(index_delta)))
      if (length(closest_peak) == 1) {
        x$index_repeat <- x$repeat_table_df$repeats[closest_peak]
        x$index_height <- x$repeat_table_df$height[closest_peak]
      } else {
        tallest_candidate <- closest_peak[which(x$repeat_table_df$height[closest_peak] == max(x$repeat_table_df$height[closest_peak]))]
        x$index_repeat <- x$repeat_table_df$repeats[tallest_candidate]
        x$index_height <- x$repeat_table_df$height[tallest_candidate]
      }
    }
  })
}
