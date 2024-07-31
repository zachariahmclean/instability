

############################# accessor-function helper ##################################


metrics_grouping_helper <- function(fragments_list,
                                    peak_threshold,
                                    window_around_main_peak) {




  # need to put the dataframe and modal peak inside samples
















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


  index_override_dataframe <- as.data.frame(index_override_dataframe)

  if(!any(index_override_dataframe[, 1] %in% names(fragments_list))){
    missing_unique_ids <- which(!index_override_dataframe[, 1] %in% names(fragments_list))

    warning(call. = FALSE,
            paste0("The following unique ids from the index override data frame are not in the repeats list:",
                   paste0(index_override_dataframe[, 1], collapse = ", ")
            )
    )
  }

  lapply(fragments_list, function(x) {
    # if there is nothing to override, then just return the existing index values
    if (any(index_override_dataframe[, 1] == x$unique_id)) {
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
    return(x)
  })
}
