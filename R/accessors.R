

# read_fsa ----------------------------------------------------------------

#' Read fsa file
#'
#' This function is just a wrapper of seqinr::read.abif where each samples is an
#' element of a list
#'
#' @param files folder path containing the fsa files
#'
#' @return
#' @export
#' @importFrom seqinr read.abif
#'
#' @examples
read_fsa <- function(files){
  #make sure file extension is fsa
  unique_file_ext = unique(tools::file_ext(files))
  if(length(unique_file_ext) > 1){
    stop("Files must be only be .fsa")
  }
  if(unique_file_ext != "fsa"){
    stop("Files must be .fsa")
  }

  #read in samples
  file_list <- list()
  for(i in seq_along(files)){
    file_list[[i]] <- seqinr::read.abif(files[i])
  }

  names(file_list) <- basename(files)

  return(file_list)
}


# ladder ------------------------------------------------------------------


#' Find the ladder peaks in the fsa file and use that to call bp size
#'
#' @param fsa_list list from 'read_fsa' function
#' @param ladder_channel string: which channel in the fsa file contains the
#'        ladder signal
#' @param signal_channel string: which channel in the fsa file contains the data
#'        signal
#' @param ladder_sizes numeric vector: bp sizes of ladder used in fragment analysis.
#'        defaults to GeneScan™ 500 LIZ™
#' @param hq_ladder logical: If TRUE, c(35, 250, 340) will be dropped from ladder
#' @param spike_location numeric: indicate the scan number of the large spike peak
#' @param smoothing_window numeric: scan window size for smoothing ladder signal
#' @param max_combinations numeric: what is the maximum number of ladder
#'        combinations that should be tested
#' @param ladder_selection_window numeric: in the ladder assigning algorithm,
#'        the we iterate through the scans in blocks and test their linear fit.
#'        In fragment analysis, we can assume that the ladder is linear over a
#'        short distance. This value defines how large that block of peaks should be.
#'
#' @return list of fragments_trace objects
#' @export
#'
#' @examples
find_ladders <- function(fsa_list,
                             ladder_channel = "DATA.105",
                             signal_channel = "DATA.1",
                             #sample_id_channel = 'SpNm.1',
                             ladder_sizes = NULL,
                             hq_ladder=TRUE,
                             spike_location = NULL,
                             smoothing_window = 5,
                             max_combinations = 2500000,
                             ladder_selection_window = 5){


  ladder_list <- vector("list", length(fsa_list))
  for(i in seq_along(fsa_list)){
    ladder_list[[i]] <- fragments_trace$new(unique_id = names(fsa_list[i]))
    ladder_list[[i]]$raw_ladder <- fsa_list[[i]]$Data[[ladder_channel]]
    ladder_list[[i]]$raw_data <- fsa_list[[i]]$Data[[signal_channel]]
    ladder_list[[i]]$scan <- 0:(length(fsa_list[[i]]$Data[[signal_channel]])- 1)
  }

  pb <- txtProgressBar(min = 0, max = length(ladder_list), style = 3)

  #find ladder for each sample and use that to predict bp size
  for (i in seq_along(ladder_list)) {

    #ladder

    ladder_fit <- fit_ladder(
      ladder = ladder_list[[i]]$raw_ladder,
      scans = ladder_list[[i]]$scan,
      ladder_sizes = ladder_sizes,
      hq_ladder = hq_ladder,
      spike_location = spike_location,
      smoothing_window = smoothing_window,
      max_combinations = max_combinations,
      ladder_selection_window = ladder_selection_window)

    ladder_list[[i]]$ladder_df <- ladder_fit

    # predict bp size
    ladder_list[[i]]$generate_mod_parameters()
    data_bp <- ladder_list[[i]]$predict_size()

    ladder_list[[i]]$trace_bp_df <- data.frame(
      unique_id = rep(ladder_list[[i]]$unique_id, length(ladder_list[[i]]$scan)),
      scan = ladder_list[[i]]$scan,
      size = ladder_list[[i]]$predict_size(),
      signal = ladder_list[[i]]$raw_data,
      ladder_signal = ladder_list[[i]]$raw_ladder
    )

    #make a warning if one of the ladder modes is bad

    rsq <- sapply(ladder_list[[i]]$mod_parameters, function(x) suppressWarnings(summary(x$mod)$r.squared))
    if(any(rsq < 0.998)){
      size_ranges <- sapply(ladder_list[[i]]$mod_parameters, function(x) x$mod$model$yi)
      size_ranges <- size_ranges[ , which(rsq < 0.998), drop = FALSE]
      size_ranges_vector <- vector('numeric', ncol(size_ranges))
      for (j in seq_along(size_ranges_vector)) {
        size_ranges_vector[j] <- paste0(size_ranges[1,j],"-",size_ranges[3,j])

      }
      warning(call. = FALSE,
              paste("sample", ladder_list[[i]]$unique_id, "has badly fitting ladder for bp sizes:",
                    paste0(size_ranges_vector, collapse = ", ")))
    }


    setTxtProgressBar(pb, i)
  }

  names(ladder_list) <- names(fsa_list)

  return(ladder_list)

}



#' Title
#'
#' @param fragments_trace_list
#' @param unique_ids
#' @param size_threshold
#' @param size_tolerance
#' @param rsq_threshold
#'
#' @return
#' @export
#'
#' @examples
fix_ladders <- function(fragments_trace_list,
                        unique_ids,
                        size_threshold = 60,
                        size_tolerance = 2.5,
                        rsq_threshold = 0.9985){
  fragments_trace_list_2 <- vector("list", length(fragments_trace_list))
  for (i in seq_along(fragments_trace_list)) {
    if(fragments_trace_list[[i]]$unique_id %in% unique_ids){
      fragments_trace_list_2[[i]] <- fragments_trace_list[[i]]$ladder_correction_auto(size_threshold = size_threshold,
                                                                  size_tolerance = size_tolerance,
                                                                  rsq_threshold = rsq_threshold)
    }
    else{
      fragments_trace_list_2[[i]] <- fragments_trace_list[[i]]
      }
  }

  names(fragments_trace_list_2) <- names(fragments_trace_list)

  return(fragments_trace_list_2)

}




#' Title
#'
#' @param fragments_trace_list
#'
#' @return
#' @export
#'
#' @examples
extract_trace_table <- function(fragments_trace_list){

  #turn the output into a dataframe
  plate_list <- lapply(fragments_trace_list, function(x){
    x$trace_bp_df
  })

  plate_combined_df <- do.call(rbind, plate_list)

  return(plate_combined_df)

}




## set class from peak table ----------------------------------------------------------------

#' Convert Peak Table to Fragments_repeats class
#'
#' This function converts a peak table data frame into a list of fragments_repeats objects.
#'
#' @param df A data frame containing the peak data.
#' @param data_format The format that the data frame is in (for example, a genemapper peak table). Choose between: genemapper5, generic.
#' @param unique_id A character string specifying column name giving the unique sample id (often the file name).
#' @param peak_size_col A character string specifying column name giving the peak size.
#' @param peak_height_col A character string specifying column name giving the peak height.
#' @param min_size_bp Numeric value indicating the minimum size of the peak table to import.
#' @param max_size_bp Numeric value indicating the maximum size of the peak table to import.
#' @param dye_col Genemapper specific. A character string specifying column name indicating the dye channel.
#' @param dye_channel Genemapper specific. A character string indicating the channel to extract data from. For example, 6-FAM is often "B".
#' @param allele_col Genemapper specific. A character string specifying column name indicating the called alleles. This is often used when the peaks have been called in genemapper.
#'
#' @return A list of fragments_repeats. objects.
#'
#' @details This function takes a peak table data frame (eg. Genemapper output) and converts it into a list of fragment objects.
#' The function supports different data formats and allows specifying column names for various attributes.
#'
#' @seealso \code{\link{repeat_table_to_repeats}}
#'
#' @examples
#'
#'   gm_raw <- instability::example_data
#'
#'   test_fragments <- peak_table_to_fragments(
#'     gm_raw,
#'     data_format = "genemapper5",
#'     dye_channel = "B",
#'     min_size_bp = 400)
#'
#' @export
peak_table_to_fragments <- function(df,
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

             new_fragments_repeats <- fragments_repeats$new(unique_id = unique(x$unique_id))
             new_fragments_repeats$peak_table_df <- df

             return(new_fragments_repeats)
             })

  return(fragments_list)
}


## set class from repeats table ----------------------------------------------------------------

#' Convert Repeat Table to Repeats Fragments
#'
#' This function converts a repeat table data frame into a list of fragments_repeats. class.
#'
#' @param df A data frame containing the repeat data.
#' @param unique_id A character string indicating the column name for unique identifiers.
#' @param repeat_col A character string indicating the column name for the repeats.
#' @param frequency_col A character string indicating the column name for the repeat frequencies.
#'
#' @return A list of fragments_repeats.
#'
#' @details This function takes a repeat table data frame and converts it into a list of repeats fragments.
#' The function allows specifying column names for repeats, frequencies, and unique identifiers.
#' @export
#'
#' @examples
#'   repeat_table <- instability::example_data_repeat_table
#'   test_fragments <- repeat_table_to_repeats(
#'     repeat_table,
#'     repeat_col = "repeats",
#'     frequency_col = "height",
#'     unique_id = "unique_id"
#'     )

repeat_table_to_repeats <- function(df,
                                    unique_id,
                                    repeat_col,
                                    frequency_col
                               ){
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

  repeats_list <- lapply(split(df, df$unique_id),
             function(x) {
               new_fragments_repeats <- fragments_repeats$new(unique_id = unique(x$unique_id))
               new_fragments_repeats$repeat_table_df <- x
               return(new_fragments_repeats)
             })


  return(repeats_list)
}

# add metadata ------------------------------------------------------------

#' Add Metadata to Fragments List
#'
#' This function adds metadata information to a list of fragments.
#'
#' @param fragments_list A list of fragment objects to which metadata will be added.
#' @param metadata_data.frame A data frame containing the metadata information.
#' @param unique_id A character string indicating the column name for unique sample identifiers in the metadata.
#' @param plate_id A character string indicating the column name for plate identifiers in the metadata.
#' @param group_id A character string indicating the column name for sample group identifiers in the metadata.
#' @param metrics_baseline_control A character string indicating the column name for baseline control indicators in the metadata.
#' @param size_standard A character string indicating the column name for repeat positive control indicators in the metadata.
#' @param size_standard_repeat_length A character string indicating the column name for repeat positive control lengths in the metadata.
#'
#' @return A modified list of fragment objects with added metadata.
#'
#' @details This function adds specified metadata attributes to each fragment in the list.
#' It matches the unique sample identifiers from the fragments list with those in the metadata data frame.
#'
#' @export
#'
#' @examples
#'
#' gm_raw <- instability::example_data
#' metadata <- instability::metadata
#'
#' test_fragments <- peak_table_to_fragments(gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B")
#'
#' test_metadata <- add_metadata(
#'   fragments_list = test_fragments,
#'   metadata_data.frame = metadata,
#'   unique_id = "unique_id",
#'   plate_id = "plate_id",
#'   group_id = "cell_line",
#'   metrics_baseline_control = "metrics_baseline_control_TF",
#'   size_standard = "size_standard",
#'   size_standard_repeat_length = "size_standard_repeat_length")

add_metadata <- function(fragments_list,
                         metadata_data.frame,
                         unique_id = "sample_file_name",
                         plate_id = "plate_id",
                         group_id = "group_id",
                         metrics_baseline_control = "metrics_baseline_control_TF",
                         size_standard = "size_standard",
                         size_standard_repeat_length = "size_standard_repeat_length") {

  # validate inputs to give good errors to user
  ## check to make sure that if the user supplies a column name, that it's actually in the dataframe they supplied
  function_input_vector <- c(plate_id, group_id, unique_id, size_standard,
                             size_standard_repeat_length, metrics_baseline_control)
  function_input_name_vector <- c("plate_id", "group_id", "unique_id", "size_standard",
                                  "size_standard_repeat_length", "metrics_baseline_control")
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
                               group_id = group_id,
                               size_standard = size_standard,
                               size_standard_repeat_length = size_standard_repeat_length,
                               metrics_baseline_control = metrics_baseline_control
                             )
                           })
}



# find fragments ------------------------------------------------------------



#' Title
#'
#' @param fragments_trace_list
#' @param smoothing_window
#' @param minumum_peak_signal
#' @param min_bp_size
#' @param max_bp_size
#'
#' @return
#' @export
#'
#' @examples
find_fragments <- function(fragments_trace_list,
                           smoothing_window = 5,
                           minumum_peak_signal = 20,
                           min_bp_size = 100,
                           max_bp_size = 1000){

  fragments_list <- lapply(fragments_trace_list, function(x){

    x$call_peaks(
      smoothing_window = smoothing_window,
      minumum_peak_signal = minumum_peak_signal,
      min_bp_size = min_bp_size,
      max_bp_size = max_bp_size
    )

    new_fragments_repeats <- fragments_repeats$new(unique_id = x$unique_id)
    new_fragments_repeats$trace_bp_df <- x$trace_bp_df
    new_fragments_repeats$peak_table_df <- x$peak_table_df
    new_fragments_repeats <- transfer_metadata_helper(x, new_fragments_repeats)

    #placeholder for function that transfers metadata

    return(new_fragments_repeats)

  })
}

# find alleles ------------------------------------------------------------


#' Find Alleles in Fragments List
#'
#' This function identifies main alleles within each fragment in a list of fragments.
#'
#' @param fragments_list A list of fragment objects containing peak data.
#' @param number_of_peaks_to_return Number of main peaks to be returned for each fragment. Default is 2.
#' @param peak_region_size_gap_threshold Gap threshold for identifying peak regions. The peak_region_size_gap_threshold is a parameter used to determine the maximum allowed gap between peak sizes within a peak region. Adjusting this parameter affects the size range of peaks that can be grouped together in a region. A smaller value makes it more stringent, while a larger value groups peaks with greater size differences, leading to broader peak regions that may encompass wider size ranges.
#' @param peak_region_height_threshold_multiplier Multiplier for the peak height threshold. The peak_region_height_threshold_multiplier parameter allows adjusting the threshold for identifying peak regions based on peak heights. Increasing this multiplier value will result in higher thresholds, making it more stringent to consider peaks as part of a peak region. Conversely, reducing the multiplier value will make the criteria less strict, potentially leading to more peaks being grouped into peak regions. It's important to note that this parameter's optimal value depends on the characteristics of the data and the specific analysis goals. Choosing an appropriate value for this parameter can help in accurately identifying meaningful peak regions in the data.
#'
#' @return A list of fragments with identified main alleles.
#'
#' @details This function finds the main alleles for each fragment in the list by identifying clusters of peaks ("peak regions")
#' with the highest signal intensities. This is based on the idea that PCR amplicons of repeats have broad peaks and PCR artififacts that help identifying the alleles.
#'  The number of peaks to be returned, and the parameters for identifying peak regions can be customized.
#'  It's important to note that both peak_region_height_threshold_multiplier and peak_region_size_gap_threshold influence the criteria for identifying peak regions, and finding the right balance between them is crucial.
#' @export
#'
#' @examples
#'  gm_raw <- instability::example_data
#'
#' test_fragments <- peak_table_to_fragments(
#'   gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B")
#'
#' test_alleles <- find_alleles(
#'   fragments_list = test_fragments,
#'   number_of_peaks_to_return = 2,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_height_threshold_multiplier = 1)

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

#' Call Repeats for Fragments
#'
#' This function calls the repeat lengths for a list of fragments.
#'
#' @param fragments_list A list of fragments_repeats objects containing fragment data.
#' @param repeat_algorithm A character specifying the repeat calling algorithm. Options: \code{"simple"} or \code{"nearest_peak"}.
#' @param assay_size_without_repeat An integer specifying the assay size without repeat for repeat calling. Default is 87.
#' @param repeat_size An integer specifying the repeat size for repeat calling. Default is 3.
#' @param repeat_length_correction A character specifying the repeat length correction method. Options: \code{"none"}, \code{"from_metadata"}, \code{"from_genemapper"}. Default is \code{"none"}.
#'
#' @return A list of \code{\link{fragments_repeats}} objects with repeat data added.
#'
#' #' @details
#' The calculated repeat lengths are assigned to the corresponding peaks in the provided `fragments_repeats` object. The repeat lengths can be used for downstream instability analysis.
#'
#' The `simple` algorithm is just the repeat size calculated either directly, or when size standards are used to correct the repeat, it's the repeat length calculated from the model of bp vs repeat length.
#'
#' The `nearest_peak` algorithm aims to model and correct for the systematic variation in fragment sizes that occurs over base pairs. It calculates repeat lengths in a way that helps align peaks with the underlying repeat pattern, making the estimation of repeat lengths more reliable relative to the main peak. The calculated repeat lengths start from the main peak's repeat length (repeat length of the main peak calculated with simple algorithm described above) and increase in increments of the specified `repeat_size`. This approach is particularly useful for mitigating the impact of size measurement underestimate, often referred to as "drift," that can occur over base pairs. By ensuring that the calculated repeat lengths are evenly spaced apart by a fixed amount (`repeat_size`), this algorithm helps stabilize the estimation of repeat lengths across peaks, leading to more consistent results.
#'
#'
#' @seealso [instability::find_main_peaks()]
#'
#' @export
#'
#' @examples
#'
#' gm_raw <- instability::example_data
#' metadata <- instability::metadata
#'
#' test_fragments <- peak_table_to_fragments(
#'   gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B")
#'
#' test_alleles <- find_alleles(
#'   fragments_list = test_fragments,
#'   number_of_peaks_to_return = 2,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_height_threshold_multiplier = 1)
#'
#' # Simple conversion from bp size to repeat size
#' test_repeats <- call_repeats(
#'   fragments_list = test_alleles,
#'   repeat_algorithm = "simple",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#'
#' # Use nearest peak algorithm to make sure called repeats are the exact number
#' # of bp apart
#'
#' test_repeats_np <- call_repeats(
#'   fragments_list = test_alleles,
#'   repeat_algorithm = "nearest_peak",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3
#' )
#'
#'
#' # correct repeat length from metadata
#'
#' test_alleles_metadata <- add_metadata(
#'   test_alleles, metadata,
#'   group_id = "cell_line",
#'   unique_id = "unique_id")
#'
#' test_repeats_corrected <- call_repeats(
#'   fragments_list = test_alleles_metadata,
#'   repeat_algorithm = "simple",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3,
#'   repeat_length_correction = "from_metadata"
#' )

call_repeats <- function(fragments_list,
                         repeat_algorithm = "simple",
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

#' Calculate Repeat Instability Metrics
#'
#' This function computes instability metrics from a list of fragments_repeats data objects.
#'
#' @param fragments_list A list of "fragments_repeats" class objects representing fragment data.
#' @param grouped Logical value indicating whether samples should be grouped to share a common index peak. Useful for cases like inferring repeat size of inherited alleles from mouse tail data. Requires metadata via \code{add_metadata()}.
#' @param peak_threshold The threshold for peak heights to be considered in the calculations, relative to the modal peak height of the expanded allele.
#' @param window_around_main_peak A numeric vector (length = 2) defining the range around the index peak. First number specifies repeats before the index peak, second after. For example, \code{c(-5, 40)} around an index peak of 100 would analyze repeats 95 to 140.
#' @param percentile_range A numeric vector of percentiles to compute (e.g., c(0.5, 0.75, 0.9, 0.95)).
#' @param repeat_range A numeric vector specifying ranges of repeats for the inverse quantile computation.
#' @param index_override_dataframe A data.frame to manually set index peaks. Column 1: unique sample IDs, Column 2: desired index peaks. Closest peak in each sample is selected.
#'
#' @return A data.frame with calculated instability metrics for each sample.
#'
#' @export
#'
#' @examples
#'  gm_raw <- instability::example_data
#'  metadata <- instability::metadata
#'
#'  test_fragments <- peak_table_to_fragments(gm_raw,
#'    data_format = "genemapper5",
#'    dye_channel = "B",
#'    min_size_bp = 400)
#'
#'  test_metadata <- add_metadata(
#'    fragments_list = test_fragments,
#'    metadata_data.frame = metadata,
#'    unique_id = "unique_id",
#'    plate_id = "plate_id",
#'    group_id = "cell_line",
#'    metrics_baseline_control = "metrics_baseline_control_TF",
#'    size_standard = "size_standard",
#'    size_standard_repeat_length = "size_standard_repeat_length")
#'
#'  test_alleles <- find_alleles(
#'    fragments_list = test_metadata,
#'    number_of_peaks_to_return = 1,
#'    peak_region_size_gap_threshold = 6,
#'    peak_region_height_threshold_multiplier = 1)
#'
#'
#'  test_repeats <- call_repeats(
#'    fragments_list = test_alleles,
#'    repeat_algorithm = "simple",
#'    assay_size_without_repeat = 87,
#'    repeat_size = 3,
#'    repeat_length_correction = "none"
#'   )
#'
#'   # grouped metrics
#'   # uses t=0 samples as indicated in metadata
#'  test_metrics_grouped <- calculate_instability_metrics(
#'    fragments_list = test_repeats,
#'    grouped = TRUE,
#'    peak_threshold = 0.05,
#'    window_around_main_peak = c(-40, 40),
#'    percentile_range = c(0.5, 0.75, 0.9, 0.95),
#'    repeat_range = c(2, 5, 10, 20))

calculate_instability_metrics <- function(fragments_list,
                              grouped = FALSE,
                              peak_threshold = 0.05,
                              # note the lower lim should be a negative value
                              window_around_main_peak = c(NA, NA),
                              percentile_range = c(0.5, 0.75, 0.9, 0.95),
                              repeat_range = c(2, 5, 10, 20),
                              index_override_dataframe = NULL){
 # is it grouped and the index peak needs to be determined from another sample?
 if(grouped == TRUE){
   fragments_list <- metrics_grouping_helper(
     fragments_list = fragments_list,
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

  #add back in any samples that were removed earlier or failed to calculate metrics (they are returned as NULL and therefore not in the dataframe)
  missing_samples <- names(fragments_list)[!names(fragments_list) %in% metrics$unique_id]
  if(length(missing_samples) > 0){
    metrics[nrow(metrics) + seq_along(missing_samples), "unique_id"] <- missing_samples
    rownames(metrics) <- metrics$unique_id
    metrics$QC_comments <- ifelse(metrics$unique_id %in% missing_samples, "metrics could not be calculated", NA_character_)
  }

  return(metrics)

}


# Extract alleles -------------------------------------------------------

#' Extract Modal Peaks
#'
#' Extracts modal peak information from each sample in a list of fragments.
#'
#' @param fragments_list A list of "fragments" class objects (can be either "fragments_repeats" or "fragments_repeats" class).
#'
#' @return A data.frame containing modal peak information for each sample.
#' @export
#'
#' @examples
#'  gm_raw <- instability::example_data
#'
#'  test_fragments <- peak_table_to_fragments(gm_raw,
#'    data_format = "genemapper5",
#'    dye_channel = "B",
#'    min_size_bp = 400)
#'
#'  test_alleles <- find_alleles(
#'    fragments_list = test_fragments,
#'    number_of_peaks_to_return = 1,
#'    peak_region_size_gap_threshold = 6,
#'    peak_region_height_threshold_multiplier = 1)
#'
#'   extract_alleles(test_alleles)
#'
extract_alleles <- function(fragments_list) {


  extracted <- lapply(fragments_list, function(x){

    data.frame(unique_id = rep(x$unique_id,2),
               size = c(x$allele_1_size, x$allele_2_size),
               repeats = c(x$allele_1_repeat, x$allele_2_repeat),
               height = c(x$allele_1_height, x$allele_2_height),
               peak_allele = c(1, 2))
  })
  extracted_df <- do.call(rbind, extracted)

  return(extracted_df)

}

# Extract fragments -------------------------------------------------------

#' Extract All Fragments
#'
#' Extracts peak data from each sample in a list of fragments.
#'
#' @param fragments_list A list of "fragments" class objects (can be either "fragments_repeats" or "fragments_repeats" class).
#'
#' @return A data.frame containing peak data for each sample.
#' @export
#'
#' @examples
#'  gm_raw <- instability::example_data
#'  metadata <- instability::metadata
#'
#'  test_fragments <- peak_table_to_fragments(gm_raw,
#'    data_format = "genemapper5",
#'    dye_channel = "B",
#'    min_size_bp = 400)
#'
#'  test_metadata <- add_metadata(
#'    fragments_list = test_fragments,
#'    metadata_data.frame = metadata,
#'    unique_id = "unique_id",
#'    plate_id = "plate_id",
#'    group_id = "cell_line",
#'    metrics_baseline_control = "metrics_baseline_control_TF",
#'    size_standard = "size_standard",
#'    size_standard_repeat_length = "size_standard_repeat_length")
#'
#'  test_alleles <- find_alleles(
#'    fragments_list = test_metadata,
#'    number_of_peaks_to_return = 1,
#'    peak_region_size_gap_threshold = 6,
#'    peak_region_height_threshold_multiplier = 1)
#'
#'  test_repeats <- call_repeats(
#'    fragments_list = test_alleles,
#'    repeat_algorithm = "simple",
#'    assay_size_without_repeat = 87,
#'    repeat_size = 3,
#'    repeat_length_correction = "none"
#'   )
#'
#'   extract_alleles(test_repeats)
#'
extract_fragments <- function(fragments_list) {

  suppressWarnings(
    extracted <- lapply(fragments_list, function(x){
      if(is.null(x$peak_table_df) & is.null(x$repeat_table_df)){
        return(NULL)
      }
      else if(!is.null(x$peak_table_df) & is.null(x$repeat_table_df)){
        df_length <- nrow(x$peak_table_df)
        data.frame(unique_id = rep(x$unique_id,df_length),
                   main_peak_size = rep(x$allele_1_size,df_length),
                   main_peak_height = rep(x$allele_1_height, df_length),
                   height = x$peak_table_df$height,
                   size = x$peak_table_df$size,
                   peak_region = x$.__enclos_env__$private$peak_regions)
      }
      else if(!is.null(x$repeat_table_df)){
        df_length <- nrow(x$repeat_table_df)
        data.frame(unique_id = rep(x$unique_id,df_length),
                   main_peak_repeat = rep(x$allele_1_repeat,df_length),
                   main_peak_height = rep(x$allele_1_height, df_length),
                   height = x$repeat_table_df$height,
                   repeats = x$repeat_table_df$repeats,
                   peak_region = x$.__enclos_env__$private$peak_regions)
      }
    })
  )
  extracted_df <- do.call(rbind, extracted)


  return(extracted_df)

}

# remove fragments -------------------------------------------------------

#' Remove Samples from List
#'
#' A convenient function to remove specific samples from a list of fragments.
#'
#' @param fragments_list A list of "fragments" class objects (can be either "fragments_repeats" or "fragments_repeats" class).
#' @param samples_to_remove A character vector containing the unique IDs of the samples to be removed.
#'
#' @return A modified list of fragments with the specified samples removed.
#' @export
#'
#' @examples
#' gm_raw <- instability::example_data
#' metadata <- instability::metadata
#'
#' test_fragments <- peak_table_to_fragments(
#'   gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B")
#'
#' all_fragment_names <- names(test_fragments)
#'
#' #pull out unique ids of samples to remove
#' samples_to_remove <- all_fragment_names[c(1,5,10)]
#'
#' samples_removed <- remove_fragments(test_fragments, samples_to_remove)
#'
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


# plot ladder -------------------------------------------------------------

#' Title
#'
#' @param fragments_trace_list
#' @param n_facet_col
#' @param xlim
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
plot_ladders <- function(fragments_trace_list,
                         n_facet_col = 2,
                         xlim = NULL,
                         ylim = NULL) {

  par(mfrow=c(ceiling(length(fragments_trace_list)/n_facet_col), n_facet_col)) # Adjust layout as needed
  for (i in seq_along(fragments_trace_list)) {
    fragments_trace_list[[i]]$plot_ladder(
      xlim = xlim,
      ylim = ylim)
  }
  par(mfrow=c(1, 1)) # Reset the layout
}

# plot traces -------------------------------------------------------------

#' Title
#'
#' @param fragments_list
#' @param show_peaks
#' @param n_facet_col
#' @param xlim
#' @param ylim
#'
#' @return
#' @export
#'
#' @examples
plot_traces <- function(fragments_list,
                        show_peaks = TRUE,
                         n_facet_col = 2,
                         xlim = NULL,
                         ylim = NULL) {

  par(mfrow=c(ceiling(length(fragments_list)/n_facet_col), n_facet_col)) # Adjust layout as needed
  for (i in seq_along(fragments_list)) {
    fragments_list[[i]]$plot_trace(
      show_peaks = show_peaks,
      xlim = xlim,
      ylim = ylim)
  }
  par(mfrow=c(1, 1)) # Reset the layout
}

# plot fragment data -------------------------------------------------------

#' Plot Peak Data
#'
#' Plots peak data from a list of fragments.
#'
#' @param fragments_list A list of "fragments" class objects (can be either "fragments_repeats" or "fragments_repeats" class).
#' @param sample_subset A character vector of unique sample IDs of samples to pick out and plot.
#' @param zoom_in_on_peak A numeric value (1 or 2) for which allele's peak should be zoomed in on. Use NA to show the whole trace.
#' @param zoom_range A numeric vector with a length of 2 to indicate the size range around the zoomed-in peak. For example, c(-20, 20).
#' @param show_peak_regions A logical value indicating whether peak regions should be included in the plot.
#' @param col_width A numeric value for the column widths of the plotted bars.
#' @param facet_nrow A numeric value indicating the number of rows for faceting in the plot (passed to 'nrow' in \code{facet_wrap()}).
#' @param facet_ncol A numeric value indicating the number of columns for faceting in the plot (passed to 'ncol' in \code{facet_wrap()}).
#' @param facet_scales A character string for setting the axis scales for faceting (passed to 'scales' in \code{facet_wrap()}).
#'
#' @return A base R plot object displaying the peak data.
#' @export
#'
#' @examples
#' gm_raw <- instability::example_data
#'
#' test_fragments <- peak_table_to_fragments(gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B")
#'
#' test_alleles <- find_alleles(
#'   fragments_list = test_fragments,
#'   number_of_peaks_to_return = 2,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_height_threshold_multiplier = 1)
#'
#' plot_fragments(test_alleles,
#'   names(test_alleles)[1:9])

plot_fragments <- function(fragments_list,
                           n_facet_col = 2,
                           xlim = NULL,
                           ylim = NULL) {

  par(mfrow=c(ceiling(length(fragments_list)/n_facet_col), n_facet_col)) # Adjust layout as needed
  for (i in seq_along(fragments_list)) {
    fragments_list[[i]]$plot_fragments(
      xlim = xlim,
      ylim = ylim)
  }
  par(mfrow=c(1, 1)) # Reset the layout
}


# plot repeat correction model --------------------------------------------



#' Plot Repeat Correction Model
#'
#' Plots the results of the repeat correction model for a list of fragments.
#'
#' @param fragments_list A list of fragments_repeats class objects obtained from the 'call_repeats' function when the 'repeat_length_correction' was either 'from_metadata' or 'from_genemapper'.
#'
#' @return A base R graphic object displaying the repeat correction model results.
#' @export
#'
#' @examples
#' gm_raw <- instability::example_data
#' metadata <- instability::metadata
#'
#' test_fragments <- peak_table_to_fragments(
#'   gm_raw,
#'   data_format = "genemapper5",
#'   dye_channel = "B")
#'
#' test_alleles <- find_alleles(
#'   fragments_list = test_fragments,
#'   number_of_peaks_to_return = 2,
#'   peak_region_size_gap_threshold = 6,
#'   peak_region_height_threshold_multiplier = 1)
#'
#' test_alleles_metadata <- add_metadata(
#'   test_alleles, metadata,
#'   group_id = "cell_line",
#'   unique_id = "unique_id")
#'
#' test_repeats_corrected <- call_repeats(
#'   fragments_list = test_alleles_metadata,
#'   repeat_algorithm = "simple",
#'   assay_size_without_repeat = 87,
#'   repeat_size = 3,
#'   repeat_length_correction = "from_metadata"
#' )
#'
#' plot_repeat_correction_model(test_repeats_corrected)
#'
#'


plot_repeat_correction_model <- function(fragments_list) {
  # Check if all models in the list are the same
  first_model_df <- fragments_list[[1]]$.__enclos_env__$private$correction_mod$model
  identical_model_test <- logical(length(fragments_list))

  for (i in seq_along(fragments_list)) {
    identical_model_test[i] <- identical(first_model_df, fragments_list[[i]]$.__enclos_env__$private$correction_mod$model)
  }

  if (!all(identical_model_test)) {
    stop("The supplied fragments list must come from the same 'call_repeats' function output", call. = FALSE)
  }

  controls_repeats_df <- fragments_list[[1]]$.__enclos_env__$private$controls_repeats_df

  # Plotting
  par(mfrow = c(1, length(unique(controls_repeats_df$plate_id))))
  for (plate_id in unique(controls_repeats_df$plate_id)) {
    plate_data <- subset(controls_repeats_df, plate_id == plate_id)

    # Generating unique colors for each unique value in unique_id
    unique_ids <- unique(plate_data$unique_id)
    colors <- rainbow(length(unique_ids))
    id_color_map <- setNames(colors, unique_ids)

    plot(plate_data$size, plate_data$validated_repeats, pch = 21, col = id_color_map[plate_data$unique_id],
         cex = 2, main = paste("Plate ID:", plate_id), xlab = "Size", ylab = "User supplied repeat length")

    lm_model <- lm(validated_repeats ~ size, data = plate_data)
    abline(lm_model, col = "blue")
  }
}


