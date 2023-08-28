## set class from peak table ----------------------------------------------------------------

#' Convert peak table into fragments objects
#'
#' A convenient function to turn a peak table (eg genemapper output) into a list, where each element is a sample of the chosen \code{fragment_class}.
#'
#' @param df data.frame
#' @param fragment_class Class that data.frame will be converted into. Choose between classes: HTT_fragments
#' @param data_format The format that the data.frame is in (for example a genemapper peak table). Choose between: genemapper5, generic
#' @param peak_size_col The column in the data.frame as a character string giving the peak size
#' @param peak_height_col The column in the data.frame as a character string giving the peak height
#' @param unique_id The column in the data.frame as a character string giving the unique sample id (often the file name)
#' @param dye_col The column in the data.frame as a character string indicating the dye channel
#' @param dye_channel A character string indicating the channel to extract data from. For example 6-FAM is often "B"
#' @param allele_col The column in the data.frame as a character string indicating the called alleles. This is often for when the peaks have been called in genemapper.
#' @param min_size_bp Numeric value indicating the minimum size of peak table to import
#' @param max_size_bp Numeric value indicating the maximum size of peak table to import
#'
#' @return A list, where each element is a sample of the chosen \code{fragment_class}
#' @export
#'
#' @examples
#'
peak_table_to_fragments <- function(df,
                                    fragment_class = "bp_fragments",
                                    data_format = NULL,
                                    peak_size_col = NULL,
                                    peak_height_col = NULL,
                                    unique_id = NULL,
                                    dye_col = NULL,
                                    dye_channel = NULL,
                                    allele_col = NULL,
                                    min_size_bp = 100,
                                    max_size_bp = 1000){
  #check to make sure that if the user supplies a column name, that it's actually in the dataframe
  if(any(!is.null(peak_size_col), !is.null(peak_height_col), !is.null(unique_id))){
    function_input_vector <- c(peak_size_col, peak_height_col, unique_id)
    function_input_name_vector <- c("peak_size_col", "peak_height_col", "unique_id")
    for (i in seq_along(function_input_vector)) {
      if(!any(names(df) == function_input_vector[[i]])){
        stop(paste0(function_input_name_vector[[i]], " input '", function_input_vector[[i]], "' was not detected as a column name in the supplied dataframe. Check column names and supply the right character string for the ",function_input_name_vector[[i]]," input"),
             call. = FALSE)
      }
    }
  }

  # chose the tidying function
  # Use the supplied user column names if given
  if (data_format == "genemapper5") {
    df2 <- clean_genemapper5(df,
                             peak_size_col = ifelse(length(peak_size_col) == 0, "Size", peak_size_col),
                             peak_height_col = ifelse(length(peak_height_col) == 0, "Height", peak_height_col),
                             unique_id = ifelse(length(unique_id) == 0, "Sample.File.Name", unique_id),
                             dye_col = ifelse(length(dye_col) == 0, "Dye.Sample.Peak", dye_col),
                             dye_channel = ifelse(length(dye_channel) == 0, "B", dye_channel),
                             allele_col = ifelse(length(allele_col) == 0, "Allele", allele_col))
  } else if(data_format == "generic"){
    df2 <- clean_generic(df,
                         peak_size_col = peak_size_col,
                         peak_height_col = peak_height_col,
                         unique_id = unique_id)
  } else {
    stop("Data format not recognised. Choose between: genemapper5, generic",
         call. = FALSE)
  }

  #filter size and split up into a list of fragments
  fragments_list <-
    lapply(split(df2, df2$unique_id),
           function(x) {
             #filter size
             df <- x[x$size > min_size_bp & x$size < max_size_bp & !is.na(x$size), , drop = FALSE]
             #check to see if all rows removed and give warning
             if(nrow(df) == 0){
               warning(paste0("Size filtering removed all rows for ", unique(x$unique_id)),
                       call. = FALSE)
             }

             bp_fragments$new(unique_id = unique(df$unique_id),
                              peak_data = df)
             })

  return(fragments_list)
}


## set class from repeats table ----------------------------------------------------------------

#' Convert repeat table into repeats objects
#'
#' A convenient function to turn a repeat table (eg summarized sequencing data) into a list, where each element is a sample of the chosen \code{repeat_class}.
#'
#' @param df data.frame
#' @param repeat_class Class that data.frame will be converted into. Choose between classes: HTT_repeats
#' @param repeat_col The column in the data.frame as a character string giving the repeat size
#' @param frequency_col The column in the data.frame as a character string giving the read count / height
#' @param unique_id The column in the data.frame as a character string giving the unique sample id
#'
#' @return A list, where each element is a sample of the chosen \code{fragment_class}
#' @export
#'
#' @examples
#'
repeat_table_to_repeats <- function(df,
                               repeat_class = "repeats_fragments",
                               repeat_col,
                               frequency_col,
                               unique_id){
  # validate inputs to give good errors to user
  ## check to make sure that if the user supplies a column name, that it's actually in the dataframe
  function_input_vector <- c(repeat_col, frequency_col, unique_id)
  function_input_name_vector <- c("repeat_col", "frequency_col", "unique_id")
  for (i in seq_along(function_input_vector)) {
    if(!any(names(df) == function_input_vector[[i]])){
      stop(paste0(function_input_name_vector[[i]], " input '", function_input_vector[[i]], "' was not detected as a column name in the supplied dataframe. Check column names and supply the right character string for the ",function_input_name_vector[[i]]," input"),
           call. = FALSE)
    }
  }

  names(df)[names(df) == repeat_col] <- 'repeats'
  names(df)[names(df) == frequency_col] <- 'height'
  names(df)[names(df) == unique_id] <- 'unique_id'


  if (repeat_class == "repeats_fragments") {
    repeats_list <-
      lapply(split(df, df$unique_id),
             function(x) {
               repeats_fragments$new(unique_id = unique(x$unique_id),
                                     repeat_data = x)
             })
  }

  return(repeats_list)
}

# add metadata ------------------------------------------------------------

#' Add metadata to a list of samples
#'
#' @param fragments_list a list of fragments class (can be either a fragments or repeats class)
#' @param metadata_data.frame A data.frame containing a single row per sample with a unique id matching the samples in the fragments_list
#' @param unique_id The column in the data.frame as a character string giving the unique sample id
#' @param plate_id The column in the data.frame as a character string giving the plate or batch samples were run in
#' @param sample_group_id The column in the data.frame as a character string indicating the sample grouping. This is used in the \code{calculate_instability_metrics()} function to set a common index peak among a group of samples. For example this could be a single mouse.
#' @param metrics_baseline_control The column in the data.frame as a character string indicating if the sample should be used to set the index peak. For example, the tail of the mouse. The values within the column should be logical.
#' @param repeat_positive_control_TF The column in the data.frame as a character string indicating if samples are positive controls that are to be used to correct repeat length. The values within the column should be logical.
#' @param repeat_positive_control_length The column in the data.frame as a character string indicating the repeat size of the expanded allele for the positive control samples. The values within the column should be numeric.
#'
#' @return A list, where each element is a fragments or repeats class
#' @export
#'
#' @examples
#'
add_metadata <- function(fragments_list,
                         metadata_data.frame,
                         unique_id = "sample_file_name",
                         plate_id = "plate_id",
                         sample_group_id = "sample_group_id",
                         metrics_baseline_control = "metrics_baseline_control_TF",
                         repeat_positive_control_TF = "repeat_positive_control_TF",
                         repeat_positive_control_length = "repeat_positive_control_length") {

  # validate inputs to give good errors to user
  ## check to make sure that if the user supplies a column name, that it's actually in the dataframe they supplied
  function_input_vector <- c(plate_id, sample_group_id, unique_id, repeat_positive_control_TF,
                             repeat_positive_control_length, metrics_baseline_control)
  function_input_name_vector <- c("plate_id", "sample_group_id", "unique_id", "repeat_positive_control_TF",
                                  "repeat_positive_control_length", "metrics_baseline_control")
  for (i in seq_along(function_input_vector)) {
    if(!any(names(metadata_data.frame) == function_input_vector[[i]])){
      stop(paste0(function_input_name_vector[[i]], " input '", function_input_vector[[i]], "' was not detected as a column name in the supplied dataframe. Check column names and supply the right character string for the ",function_input_name_vector[[i]]," input"),
           call. = FALSE)
    }
  }
  ## check if user has any duplicated unique ids
  supplied_ids <- metadata_data.frame[,unique_id, drop = TRUE]
  if(anyDuplicated(supplied_ids) != 0){
    stop(paste0(unique_id, " does not contain unique sample ids. The metadata must have one row per unique sample id."),
         call. = FALSE)
  }
  ## Give warning if samples don't have metadata
  not_in_metadata <- which(!names(fragments_list) %in% supplied_ids)
  if(length(not_in_metadata) > 0){
    warning(paste0("The following samples do not have a corresponding unique id in the metadata: ",
                   paste0(names(fragments_list)[not_in_metadata], collapse = ", ")),
            call. = FALSE)
  }

  ## Give warning if user tries to give a of metadata but it's not in sample list
  not_in_samples <- which(!supplied_ids %in% names(fragments_list))
  if(length(not_in_samples) > 0){
    warning(paste0("The following unique ids in the metadata file do not have a corresponding sample: ",
                   paste0(supplied_ids[not_in_samples], collapse = ", ")),
            call. = FALSE)
  }

  metadata_added <- lapply(fragments_list,
                           function(x) {
                             x$add_metadata(
                               metadata_data.frame = metadata_data.frame,
                               unique_id = unique_id,
                               plate_id = plate_id,
                               sample_group_id = sample_group_id,
                               repeat_positive_control_TF = repeat_positive_control_TF,
                               repeat_positive_control_length = repeat_positive_control_length,
                               metrics_baseline_control = metrics_baseline_control
                             )
                           })
}




# find alleles ------------------------------------------------------------


#' Find modal peaks
#'
#' Find the modal peaks across the list of samples
#' @param fragments_list a list of fragments class (can be either a fragments or repeats class)
#' @param number_of_peaks_to_return Numeric value of 1 or 2 indicating if one (mouse) or two (human) modal peaks should be identified, respectively.
#' @param peak_selection_method A heuristic to help chose the peak region that should be used to find the expanded allele. "tallest" will prioritize the peak region that has the tallest modal peak, while "largest" will pick the peak region that has the largest size. "largest" is often useful for when using in vitro cell models that have extremly long expanded alleles.
#' @param peak_region_size_gap_threshold Numeric value to indicate how large the gap should be between peaks before they are no longer included in the same peak regions
#' @param peak_region_height_threshold_multiplier Numeric value to set the peak height required to be in the same peak region. The threshold is the mean of all of the peaks and this is a multiplier that adjusts that threshold. For example, 0.5 will set the height threshold to be half of the mean peak height.
#'
#' @return A list, where each element is a fragments or repeats class
#' @export
#'
#' @examples
find_alleles <- function(fragments_list,
                         number_of_peaks_to_return = 2,
                         peak_region_size_gap_threshold = 6,
                         peak_region_height_threshold_multiplier = 1){

  main_peaks <- lapply(fragments_list, function(x){

    x$find_main_peaks(number_of_peaks_to_return = number_of_peaks_to_return,
                      peak_region_size_gap_threshold = peak_region_size_gap_threshold,
                      peak_region_height_threshold_multiplier = peak_region_height_threshold_multiplier)

  })
}


# call_repeats ------------------------------------------------------------

#' Call repeats
#'
#' Use base-pair level data and the identified modal peaks to call repeat size
#'
#' @param fragments_list  a list of fragments class (can be either a fragments or repeats class)
#' @param repeat_algorithm An optional algorithm that can be used to call the repeat size. Either "nearest_peak" (see description for more details) or "none" to just use the base-pair size and arithmetic to calculate the repeat.
#' @param assay_size_without_repeat A numeric value indicating how many flanking bases there are around the repeat. Required when repeat size is to be directly calculted.
#' @param repeat_size The number of nucleotides in each repeat unit. For example a triplet repeat would be 3
#' @param repeat_length_correction Indicating whether sample metadata should be used to correct the repeat length across the samples. This requires metadata added via the \code{add_metadata()} function. Options are either "none", "from_metadata" (added via metadata on selected samples) or "from_genemapper" (when peak table exported from genemapper has expanded allele manually called)
#'
#' @return A list, where each element is a repeats class
#' @export
#'
#' @examples
#'
call_repeats <- function(fragments_list,
                         repeat_algorithm = "nearest_peak",
                         assay_size_without_repeat = 87,
                         repeat_size = 3,
                         repeat_length_correction = "none" # or "from_metadata" or "from_genemapper"
) {
  # Check to see if repeats are to be corrected
  # if so, supply the model to each of the samples in the list
  if (repeat_length_correction %in% c("from_metadata", "from_genemapper")) {
    mod <- model_repeat_length(fragments_list = fragments_list,
                               assay_size_without_repeat = assay_size_without_repeat,
                               repeat_size = repeat_size,
                               repeat_length_correction = repeat_length_correction)

    for (i in seq_along(fragments_list)) {
      fragments_list[[i]]$.__enclos_env__$private$correction_mod <- mod$correction_mods
      fragments_list[[i]]$.__enclos_env__$private$controls_repeats_df <- mod$controls_repeats_df
    }
  }

  # call repeats for each sample
    added_repeats <- lapply(fragments_list,
                            function(x) {
                              x$add_repeats(
                                repeat_algorithm = repeat_algorithm,
                                assay_size_without_repeat = assay_size_without_repeat,
                                repeat_size = repeat_size,
                                correct_repeat_length = ifelse(repeat_length_correction == "none", FALSE, TRUE)
                              )
                            })

    return(added_repeats)

}



# Calculate metrics -------------------------------------------------------

#' calculate instability metrics
#'
#' Calculate repeat metrics from a list of repeats objects
#' @param fragments_list a list of repeats class
#' @param grouped A logical value for whether the samples should be grouped to have a common index peak. For example, using the mouse tail to send the repeat size of the inherited allele. This requires metadata added via the \code{add_metadata()} function.
#' @param peak_threshold A threshold for the height of the peaks that should be includes for the calculations. Always proportional to the modal peak height of the expanded allele.
#' @param window_around_main_peak A numeric vector with a length of two to limit set a threshold of size around the index peak. The first number being the the number of repeats before the index peak and the second being the number of repeats after. For example, c(-5, 40) when the index peak is 100 would limit the analysis to 95 to 140 repeats.
#' @param percentile_range A numeric vector of any length to indicate which percentiles should be calculates.
#' @param index_override_dataframe A data.frame which can be used to manually set the index peak. First column being the unique sample id and second column the desired index peak. The closest peak in the sample to the desired index peak will then be selected.
#'
#' @return data.frame
#' @export
#'
#' @examples
calculate_instability_metrics <- function(fragments_list,
                              grouped = FALSE,
                              peak_threshold = 0.05,
                              # note the lower lim should be a negative value
                              window_around_main_peak = c(NA, NA),
                              percentile_range = c(0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99),
                              repeat_range = c( 1,2,3,4,seq(6,20,2)),
                              index_override_dataframe = NULL){
 # is it grouped and the index peak needs to be determined from another sample?
 if(grouped == TRUE){
   fragments_list <- metrics_grouping_helper(fragments_list = fragments_list,
                                             peak_threshold = peak_threshold,
                                             window_around_main_peak = window_around_main_peak)
 } else if(is.null(index_override_dataframe)){
    #this is to make sure that we use the modal peak as the index peak
   #fixes cases where index peak has been assigned previously and we need to make sure it's the modal peak
   fragments_list <- lapply(fragments_list, function(x){
     x$index_repeat <- x$allele_1_repeat
     x$index_height <- x$allele_1_height
     x$index_weighted_mean_repeat <- NA_real_
     return(x)
     })

  }

  # override index peak with manually supplied values
  if(!is.null(index_override_dataframe)){
    fragments_list <- metrics_override_helper(
      fragments_list = fragments_list,
      index_override_dataframe = index_override_dataframe)
  }

  # calculate metrics
  metrics_list <- lapply(fragments_list, function(x){
    x$instability_metrics(peak_threshold = peak_threshold,
                          window_around_main_peak = window_around_main_peak,
                          percentile_range = percentile_range,
                          repeat_range = repeat_range)
  })
  metrics <- do.call(rbind, metrics_list)

  #add back in any samples that failed to calculate metrics (they are returned as NULL and therefore not in the dataframe)
  missing_samples <- names(fragments_list)[!names(fragments_list) %in% metrics$unique_id]
  if(length(missing_samples) > 0){
    metrics[nrow(metrics) + seq_along(missing_samples), "unique_id"] <- missing_samples
    rownames(metrics) <- metrics$unique_id
    metrics$QC_comments <- ifelse(metrics$unique_id %in% missing_samples, "metrics could not be calculated", NA_character_)
  }

  return(metrics)

}


# Extract alleles -------------------------------------------------------

#' Extract modal peaks
#'
#' Extract the modal peak information from each sample in a list of fragments
#' @param fragments_list  a list of fragments class (can be either a fragments or repeats class)
#'
#' @return data.frame
#' @export
#'
#' @examples
extract_alleles <- function(fragments_list) {

  suppressWarnings(

    if(grepl("bp", class(fragments_list[[1]])[1])){

      extracted <- lapply(fragments_list, function(x){

        data.frame(unique_id = rep(x$unique_id,2),
                   size = c(x$allele_2_size, x$allele_1_size),
                   height = c(x$allele_2_height, x$allele_1_height),
                   peak_allele = c(1, 2))
      })
      extracted_df <- do.call(rbind, extracted)

    }else if(grepl("repeats", class(fragments_list[[1]])[1])){
      extracted <- lapply(fragments_list, function(x){

        data.frame(unique_id = rep(x$unique_id,2),
                   repeats = c(x$allele_2_repeat, x$allele_1_repeat),
                   height = c(x$allele_2_height, x$allele_1_height),
                   peak_allele = c(1, 2))
      })
      extracted_df <- do.call(rbind, extracted)

    }

  )

  return(extracted_df)

}

# Extract fragments -------------------------------------------------------

#' Extract all fragments
#'
#' Extract the peak data from each sample in a list of fragments
#' @param fragments_list a list of fragments class (can be either a fragments or repeats class)
#'
#' @return data.frame
#' @export
#'
#' @examples
extract_fragments <- function(fragments_list) {

  suppressWarnings(

    if(grepl("bp", class(fragments_list[[1]])[1])){
      extracted <- lapply(fragments_list, function(x){
        if(is.null(x$peak_data)){
          return(NULL)
        } else{
          df_length <- nrow(x$peak_data)
          data.frame(unique_id = rep(x$unique_id,df_length),
                     main_peak_size = rep(x$allele_1_size,df_length),
                     main_peak_height = rep(x$allele_1_height, df_length),
                     height = x$peak_data$height,
                     size = x$peak_data$size,
                     peak_region = x$.__enclos_env__$private$peak_regions)
        }
      })
    }else if(grepl("repeats", class(fragments_list[[1]])[1])){
      extracted <- lapply(fragments_list, function(x){
        if(is.null(x$repeat_data)){
          return(NULL)
        } else{
          df_length <- nrow(x$repeat_data)
          data.frame(unique_id = rep(x$unique_id,df_length),
                     main_peak_repeat = rep(x$allele_1_repeat,df_length),
                     main_peak_height = rep(x$allele_1_height, df_length),
                     height = x$repeat_data$height,
                     repeats = x$repeat_data$repeats,
                     peak_region = x$.__enclos_env__$private$peak_regions)
        }
      })
    }

  )

  extracted_df <- do.call(rbind, extracted)


  return(extracted_df)

}

# remove fragments -------------------------------------------------------

#' Remove samples from list
#'
#' A convenient function to remove samples from a list.
#' @param fragments_list a list of fragments class (can be either a fragments or repeats class)
#'
#' @return data.frame
#' @export
#'
#' @examples
remove_fragments <- function(fragments_list,
                             samples_to_remove){

  unique_ids <- vector("numeric", length = length(fragments_list))
  for (i in seq_along(fragments_list)) {
    unique_ids[[i]] <- fragments_list[[i]]$unique_id
  }
  samples_removed <- fragments_list
  suppressWarnings(
    samples_removed[which(unique_ids %in% samples_to_remove)] <- NULL
  )
  return(samples_removed)
}




# plot fragment data -------------------------------------------------------

#' Plot peak data
#'
#' Plots data from a list of fragments
#' @param fragments_list  A list of fragments class (can be either a fragments or repeats class)
#' @param sample_subset A character vector of unique sample ids of samples to pick out and plot
#' @param zoom_in_on_peak A numeric value for which allele should be zoomed in on. Must be 1 or 2. NA will show the whole trace.
#' @param zoom_range A numeric vector with a length of 2 to indicate the size around zoomed in peak to include. For example c(20, 20).
#' @param show_peak_regions A logical value to indicate whether the peak regions should be included in the plot
#' @param col_width A numeric value to for the column widths
#' @param facet_nrow A numeric value for how many faceting rows. Is passed to 'nrow' in \code{facet_wrap()}.
#' @param facet_ncol A numeric value for how many faceting columns. Is passed to 'ncol' in \code{facet_wrap()}.
#' @param facet_scales A character string for setting the axis scales. Is passed to 'scales' in \code{facet_wrap()}.
#'
#' @return data.frame
#' @export
#'
#' @examples
plot_fragments <- function(fragments_list,
                           sample_subset = NULL,
                           zoom_in_on_peak = NA, # 1 or 2
                           zoom_range = NULL, #must be numeric vector with 2 values, first one negative
                           show_peak_regions = FALSE,
                           col_width = NULL,
                           facet_nrow = NULL,
                           facet_ncol = NULL,
                           facet_scales = "fixed"
){
  suppressWarnings(
    #figure out which column to use for plotting and use default settings for different class
    if(grepl("bp", class(fragments_list[[1]])[1])){
      x_axis_col <- "size"
      defined_col_width <- ifelse(is.null(col_width), 1.5, col_width)
      defined_zoom_range <- if(is.null(zoom_range)) c(-60, 60) else zoom_range
    } else if(grepl("repeats", class(fragments_list[[1]])[1])){
      x_axis_col <- "repeats"
      defined_col_width <- ifelse(is.null(col_width), 0.6, col_width)
      defined_zoom_range <- if(is.null(zoom_range)) c(-20, 20) else zoom_range
    } else {
      stop("Invalid input. The 'fragments_list' must contain a list of 'fragments' objects",
           call. = FALSE)
    }
  )

  #subset data
  if(!is.null(sample_subset)){
    fragments_list <- fragments_list[which(names(fragments_list) %in% sample_subset)]
  }
  modal_peaks <- extract_alleles(fragments_list)
  all_peaks <- extract_fragments(fragments_list)
  #zoom in on data if desired
  if(zoom_in_on_peak %in% c(1,2)){
    modal_peaks <- modal_peaks[which(modal_peaks$peak_allele == zoom_in_on_peak), ,drop = FALSE]
    all_peaks_by <- lapply(split(all_peaks, all_peaks$unique_id),
                           function(df){
                             upper_lim <- abs(defined_zoom_range[2]) + modal_peaks[which(modal_peaks$unique_id == df$unique_id[[1]]), x_axis_col]
                             lower_lim <- modal_peaks[which(modal_peaks$unique_id == df$unique_id[[1]]), x_axis_col] - abs(defined_zoom_range[1])
                             df[which(df[[x_axis_col]] < upper_lim & df[[x_axis_col]] > lower_lim), ,drop = FALSE]
                           })
    all_peaks <- do.call(rbind, all_peaks_by)
  }
  else if(!is.na(zoom_in_on_peak)){
    stop("zoom_in_on_peak must be either '1' or '2'",
         call. = FALSE)
  }

  #extract peak region data if required
  if(show_peak_regions == TRUE){
    pr_df <- all_peaks[which(!is.na(all_peaks$peak_region)), ]
    # make new column name to split by
    max_heights <- aggregate(data = pr_df, height ~ unique_id, FUN = max)
    pr_df$unique_id_peak_region <- paste(pr_df$unique_id, pr_df$peak_region, sep = "_")
    pr_df_by <- lapply(split(pr_df, pr_df$unique_id_peak_region),
                       function(df){
                         data.frame(unique_id = df$unique_id[[1]],
                                    peak_region = df$peak_region[[1]],
                                    xmin = min(df[x_axis_col]),
                                    xmax = max(df[x_axis_col]),
                                    ymax = max_heights[which(max_heights$unique_id == unique(df$unique_id)), "height"])
                       })
    pr_summarised_df <- do.call(rbind, pr_df_by)
  }

  # compose plot
  gg <- ggplot2::ggplot(all_peaks)

  if(show_peak_regions == TRUE){
    gg <- gg +
      ggplot2::geom_rect(data = pr_summarised_df,
                         ggplot2::aes(xmin = xmin,
                    xmax = xmax,
                    ymin = 0,
                    ymax = ymax,
                    fill = as.factor(peak_region)),
                alpha = 0.5) +
      ggplot2::scale_fill_viridis_d(begin = 0.3)
  }

  gg <- gg +
    ggplot2::geom_col(ggplot2::aes(x = eval(parse(text = x_axis_col)),
                 y = height,
                 text = paste('Unique id: ', unique_id,
                              '<br>Repeat:', round(eval(parse(text = x_axis_col)), digits = 2),
                              '<br>Height:', height)),
             width = defined_col_width) +
    ggplot2::geom_point(data = modal_peaks,
                        ggplot2::aes(x = eval(parse(text = x_axis_col)),
                   y = height,
                   colour = fct_rev(as.character(peak_allele)),
                   text = paste('Unique id: ', unique_id,
                                '<br>Repeat:', round(eval(parse(text = x_axis_col)), digits = 2),
                                '<br>Height:', height)),
               size = 1) +
    ggplot2::facet_wrap(ggplot2::vars(unique_id),
               nrow = facet_nrow,
               ncol = facet_ncol,
               scales = facet_scales) +
    ggplot2::labs(x = x_axis_col,
         fill = "Peak Region",
         colour = NULL)
  return(gg)
}


# plot repeat correction model --------------------------------------------




plot_repeat_correction_model <- function(
    fragments_list
){
  ###may want to check to see if all models in this list are the same
  first_model_df <- fragments_list[[1]]$.__enclos_env__$private$correction_mod$model
  identical_model_test <- vector("logical", length(fragments_list))
  for (i in seq_along(fragments_list)) {
    identical_model_test[i] <- identical(first_model_df, fragments_list[[i]]$.__enclos_env__$private$correction_mod$model)
  }
  if(!all(identical_model_test)){
    stop("The supplied fragments list must come from the same 'call_repeats' function output",
         call. = FALSE)
  }

  controls_repeats_df <- fragments_list[[1]]$.__enclos_env__$private$controls_repeats_df

  ggplot2::ggplot(controls_repeats_df) +
    ggplot2::geom_point(ggplot2::aes(size, validated_repeats, colour = unique_id),
               size = 2, shape = 21) +
    ggplot2::geom_line(ggplot2::aes(size,predicted_repeat),
              alpha = 0.5, colour = "blue") +
    ggplot2::facet_wrap(ggplot2::vars(plate_id)) +
    ggplot2::labs(y = "User supplied repeat length")



}






