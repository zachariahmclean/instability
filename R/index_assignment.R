

############################# accessor-function helper ##################################


metrics_grouping_helper <- function(fragments_list,
                                    peak_threshold,
                                    window_around_main_peak) {

  # what we're doing here is pulling out the key data for all the samples that are metrics controls
  # each sample will then have the data for their appropriate control inserted inside
  # that can then be used in the calculation of instability metrics



  # make a list of dataframes and alleles
  group_ids <- sapply(fragments_list, function(x) x$group_id)
  unique_group_ids <- unique(group_ids)

  baseline_control_list <- vector("list", length(unique_group_ids))
  names(baseline_control_list) <- unique_group_ids

  for(i in seq_along(fragments_list)){
    if(fragments_list[[i]]$metrics_baseline_control == TRUE){

      #since there can be more than one control, make a list of them
      baseline_control_list[[fragments_list[[i]]$group_id]] <- c(
         baseline_control_list[[fragments_list[[i]]$group_id]],
         list(
            list(
              fragments_list[[i]]$allele_1_repeat,
              fragments_list[[i]]$repeat_table_df
             )
          )
       )
    }
  }

  # do some quality control

  for (i in seq_along(baseline_control_list)) {
    controls_missing_allele <- all(sapply(baseline_control_list[[fragments_list[[i]]$group_id]], function(x) is.na(x[[1]])))

    if (length(baseline_control_list[[i]]) == 0) {
      stop(paste0("Group '", names(baseline_control_list)[[i]], "' has no 'metrics_baseline_control'"),
           call. = FALSE
      )
    } else if (controls_missing_allele == TRUE) {
      stop(paste0("Group '", names(baseline_control_list)[[i]], "' control has no allele called. Grouped analysis won't work for these samples."),
           call. = FALSE
      )
    } else if (length(baseline_control_list[[i]]) > 1) {
      message(paste0("Group '", names(baseline_control_list)[[i]], "' has more than one 'metrics_baseline_control'. The median repeat of the assigned samples will be used to assign the index peak"))
    }
  }

  #loop over each sample and put data inside
  for(i in seq_along(fragments_list)){
    fragments_list[[i]]$.__enclos_env__$private$index_samples <- baseline_control_list[[fragments_list[[i]]$group_id]]

    control_index_median_repeat <- median(sapply(baseline_control_list[[fragments_list[[i]]$group_id]], function(x) x[[1]]))


    # since the repeat size may not be an integer, need t find what the closest peak is to the control sample
    # delta between repeat of index sample and all repeats of sample
    index_delta <- fragments_list[[i]]$repeat_table_df$repeats - control_index_median_repeat
    closest_peak <- which(abs(index_delta) == min(abs(index_delta)))

    if (length(closest_peak) == 1) {
      fragments_list[[i]]$index_repeat <- fragments_list[[i]]$repeat_table_df$repeats[closest_peak]
      fragments_list[[i]]$index_height <- fragments_list[[i]]$repeat_table_df$height[closest_peak]
    } else {
      tallest_candidate <- closest_peak[which(fragments_list[[i]]$repeat_table_df$height[closest_peak] == max(fragments_list[[i]]$repeat_table_df$height[closest_peak]))]
      fragments_list[[i]]$index_repeat <- fragments_list[[i]]$repeat_table_df$repeats[tallest_candidate]
      fragments_list[[i]]$index_height <- fragments_list[[i]]$repeat_table_df$height[tallest_candidate]
    }

  }

  return(fragments_list)
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
