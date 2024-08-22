transfer_metadata_helper <- function(old_fragment,
                                     new_fragment) {
  metadata_names <- c(
    "unique_id",
    "plate_id",
    "group_id",
    "size_standard",
    "size_standard_sample_id",
    "size_standard_repeat_length",
    "metrics_baseline_control"
  )


  for (name in metadata_names) {
    eval(parse(
      text = paste0(
        "new_fragment$",
        name,
        "<- old_fragment$",
        name
      )
    ))
  }
  return(new_fragment)
}

# add metadata ------------------------------------------------------------

#' Add Metadata to Fragments List
#'
#' This function adds metadata information to a list of fragments.
#'
#' @param fragments_list A list of fragment objects to which metadata will be added.
#' @param metadata_data.frame A data frame containing the metadata information.
#' @param unique_id (required) A character string indicating the column name for unique sample identifiers in the metadata.
#' @param plate_id (optional) A character string indicating the column name for plate identifiers in the metadata.  To skip, provide NA.
#' @param group_id (optional) A character string indicating the column name for sample group identifiers in the metadata. To skip, provide NA.
#' @param metrics_baseline_control (optional) A character string indicating the column name for baseline control indicators in the metadata. To skip, provide NA.
#' @param size_standard (optional) A character string indicating the column name for repeat positive control indicators in the metadata. To skip, provide NA.
#' @param size_standard_sample_id (optional) A character string indicating the column name for an id of the size standard. For example, a sample code. This is so that checks can be done that the sample has the same modal peak across different fragment analysis runs. To skip, provide NA.
#' @param size_standard_repeat_length (optional) A character string indicating the column name for repeat positive control lengths in the metadata. To skip, provide NA.
#'
#' @return A modified list of fragment objects with added metadata.
#'
#' @details This function adds specified metadata attributes to each fragment in the list.
#' It matches the unique sample identifiers from the fragments list with those in the metadata data frame.
#' To skip any of the optional columns, make parameter NA.
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
#'   dye_channel = "B"
#' )
#'
#' test_metadata <- add_metadata(
#'   fragments_list = test_fragments,
#'   metadata_data.frame = metadata,
#'   unique_id = "unique_id",
#'   plate_id = "plate_id",
#'   group_id = "group_id",
#'   metrics_baseline_control = "metrics_baseline_control",
#'   size_standard = "size_standard",
#'   size_standard_repeat_length = "size_standard_repeat_length"
#' )
#'
#' # skip unwanted metadata by using NA
#'
#' test_metadata_skipped <- add_metadata(
#'   fragments_list = test_fragments,
#'   metadata_data.frame = metadata,
#'   unique_id = "unique_id",
#'   plate_id = "plate_id",
#'   group_id = "group_id",
#'   metrics_baseline_control = "metrics_baseline_control",
#'   size_standard = NA,
#'   size_standard_repeat_length = NA
#' )
#'
#'
#'
#'
#'
add_metadata <- function(
    fragments_list,
    metadata_data.frame,
    unique_id = "unique_id",
    plate_id = "plate_id",
    group_id = "group_id",
    metrics_baseline_control = "metrics_baseline_control",
    size_standard = "size_standard",
    size_standard_sample_id = "size_standard_sample_id",
    size_standard_repeat_length = "size_standard_repeat_length") {
  # validate inputs to give good errors to user
  ## check to make sure that if the user supplies a column name, that it's actually in the dataframe they supplied
  function_input_vector <- c(
    plate_id, group_id, unique_id, size_standard,
    size_standard_repeat_length, metrics_baseline_control
  )
  function_input_name_vector <- c(
    "plate_id", "group_id", "unique_id", "size_standard",
    "size_standard_repeat_length", "metrics_baseline_control"
  )
  for (i in seq_along(function_input_vector)) {
    if (!any(names(metadata_data.frame) == function_input_vector[[i]]) & !is.na(function_input_vector[[i]])) {
      stop(paste0(function_input_name_vector[[i]], " input '", function_input_vector[[i]], "' was not detected as a column name in the supplied dataframe. Check column names and supply the right character string for the ", function_input_name_vector[[i]], " input. If you don't want to add this metadata category, set '", function_input_name_vector[[i]], " = NA'"),
        call. = FALSE
      )
    }
  }
  ## check if user has any duplicated unique ids
  supplied_ids <- metadata_data.frame[, unique_id, drop = TRUE]
  if (anyDuplicated(supplied_ids) != 0) {
    stop(paste0(unique_id, " does not contain unique sample ids. The metadata must have one row per unique sample id."),
      call. = FALSE
    )
  }
  ## Give warning if samples don't have metadata
  not_in_metadata <- which(!names(fragments_list) %in% supplied_ids)
  if (length(not_in_metadata) > 0) {
    warning(
      paste0(
        "The following samples do not have a corresponding unique id in the metadata: ",
        paste0(names(fragments_list)[not_in_metadata], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ## Give warning if user tries to give a of metadata but it's not in sample list
  not_in_samples <- which(!supplied_ids %in% names(fragments_list))
  if (length(not_in_samples) > 0) {
    warning(
      paste0(
        "The following unique ids in the metadata file do not have a corresponding sample: ",
        paste0(supplied_ids[not_in_samples], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  metadata_added <- lapply(
    fragments_list,
    function(fragments) {

      # make sure dataframe, not tibble
      metadata_data.frame <- as.data.frame(metadata_data.frame)
    
      # filter for row of sample
      sample_metadata <- metadata_data.frame[which(metadata_data.frame[unique_id] == fragments$unique_id), , drop = FALSE]
    
      # add metadata to slots, checking if parameters are NA
      fragments$plate_id <- if (!is.na(plate_id)) as.character(sample_metadata[[plate_id]]) else NA_character_
      fragments$group_id <- if (!is.na(group_id)) as.character(sample_metadata[[group_id]]) else NA_character_
      fragments$size_standard <- if (!is.na(size_standard)) {
        ifelse(is.na(sample_metadata[[size_standard]]) || !as.logical(sample_metadata[[size_standard]]), FALSE, TRUE)
      } else {
        FALSE
      }
      fragments$size_standard_sample_id <- if (!is.na(size_standard_sample_id)) as.character(sample_metadata[[size_standard_sample_id]]) else NA_character_
      fragments$size_standard_repeat_length <- if (!is.na(size_standard_repeat_length)) as.double(sample_metadata[[size_standard_repeat_length]]) else NA_real_
      fragments$metrics_baseline_control <- if (!is.na(metrics_baseline_control)) {
        ifelse(is.na(sample_metadata[[metrics_baseline_control]]) || !as.logical(sample_metadata[[metrics_baseline_control]]), FALSE, TRUE)
      } else {
        FALSE
      }
    
      return(fragments)
    }
  )

  return(metadata_added)
}

