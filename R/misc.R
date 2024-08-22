print_helper <- function(fragment,
                         exclude = NULL,
                         sample_attrs) {
  class_name <- class(fragment)[1]
  unique_id <- fragment$unique_id

  cat(paste0("\033[1;34m< ", class_name, " object >\033[0m\n"))
  cat("\033[1;36m-----------------------------\033[0m\n")

  # Section: Sample Attributes
  for (attr in sample_attrs) {
    if (attr %in% names(fragment)) {
      value <- fragment[[attr]]
      cat(sprintf("\033[1;33m%-30s\033[0m", attr))

      if (is.null(value)) {
        cat("NULL\n")
      } else if (is.logical(value)) {
        cat(ifelse(value, "\033[1;32mTRUE\033[0m", "\033[1;31mFALSE\033[0m"), "\n")
      } else if (is.numeric(value)) {
        if (length(value) == 1) {
          cat(format(value), "\n")
        } else {
          cat(format(paste("numeric vector length", length(value))), "\n")
        }
      } else if (is.character(value)) {
        if (length(value) == 1) {
          cat(format(value), "\n")
        } else {
          cat(format(paste("character vector length", length(value))), "\n")
        }
      } else {
        cat(class(value), "\n")
      }
    }
  }

  cat("\033[1;36m-----------------------------\033[0m\n")

  slot_names <- ls(fragment, all.names = TRUE)
  all_exclusions <- c(
    exclude,
    sample_attrs,
    slot_names[which(sapply(slot_names, function(name) class(fragment[[name]]) == "function"))],
    ".__active__",
    ".__enclos_env__"
  )
  slot_names <- setdiff(slot_names, all_exclusions)

  for (name in slot_names) {
    value <- fragment[[name]]
    class_value <- class(value)

    cat(sprintf("\033[1m%-30s\033[0m", name))

    if (is.null(value)) {
      cat("NULL\n")
    } else if (is.logical(value)) {
      cat(ifelse(value, "\033[1;32mTRUE\033[0m", "\033[1;31mFALSE\033[0m"), "\n")
    } else if (is.numeric(value)) {
      if (length(value) == 1) {
        cat(format(value), "\n")
      } else {
        cat(format(paste("numeric vector length", length(value))), "\n")
      }
    } else if (is.character(value)) {
      if (length(value) == 1) {
        cat(format(value), "\n")
      } else {
        cat(format(paste("character vector length", length(value))), "\n")
      }
    } else if (is.data.frame(value)) {
      cat(sprintf("data.frame: %d rows × %d cols\n", nrow(value), ncol(value)))
    } else {
      cat(class_value, "\n")
    }
  }

  invisible(fragment)
}


# remove fragments -------------------------------------------------------

#' Remove Samples from List
#'
#' A convenient function to remove specific samples from a list of fragments.
#'
#' @param fragments_list A list of fragments_repeats objects containing fragment data.
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
#'   dye_channel = "B"
#' )
#'
#' all_fragment_names <- names(test_fragments)
#'
#' # pull out unique ids of samples to remove
#' samples_to_remove <- all_fragment_names[c(1, 5, 10)]
#'
#' samples_removed <- remove_fragments(test_fragments, samples_to_remove)
#'
remove_fragments <- function(
    fragments_list,
    samples_to_remove) {
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



# qmd templates -----------------------------------------------------------


#' Generate a Quarto file that has the instability pipeline preset
#'
#' @param file_name Name of file to create
#' @param size_standards Indicates if the functionality for correcting repeat size using size standards be included in the pipeline. See \code{\link{add_metadata}} & \code{\link{call_repeats}} for more info.
#' @param samples_grouped Indicates if the functionality for grouping samples for metrics calculations should be included in the pipeline. See \code{\link{add_metadata}} & \code{\link{assign_index_peaks}} for more info.
#'
#' @return A Quarto template file
#' @export
#'
#' @importFrom utils file.edit
#'
#' @examples
#'
#' if (interactive()) {
#'   generate_instability_template("test")
#' }
#'
#'
generate_instability_template <- function(
  file_name = NULL,
  size_standards = TRUE,
  samples_grouped = TRUE) {

comment_out_lines <- function(content, start_pattern, end_pattern, comment_message = NULL) {
  start_idx <- grep(start_pattern, content)
  end_idx <- grep(end_pattern, content)[grep(end_pattern, content) > start_idx][1]

  if (length(start_idx) > 0 && length(end_idx) > 0) {
    content[(start_idx + 1):(end_idx - 1)] <- paste("#", content[(start_idx + 1):(end_idx - 1)])
    if (!is.null(comment_message)) {
      content <- append(content, comment_message, after = start_idx - 1)
    }
  }

  return(content)
}

if (is.null(file_name)) {
  stop("You must provide a valid file_name")
}

source_file <- system.file("extdata/_extensions/template.qmd", package = "instability")

if (!file.exists(source_file)) {
  stop(paste("Source file does not exist:", source_file))
}

template_content <- readLines(source_file)

if (!size_standards) {
  template_content <- gsub('metadata\\$plate_id <- metadata\\$plate_id', '# metadata$plate_id <- metadata$plate_id', template_content)
  template_content <- gsub('metadata\\$size_standard <- metadata\\$size_standard', '# metadata$size_standard <- metadata$size_standard', template_content)
  template_content <- gsub('metadata\\$size_standard_repeat_length <', '# metadata\\$size_standard_repeat_length <', template_content)
  template_content <- gsub('repeat_length_correction = "from_metadata"', 'repeat_length_correction = "none"', template_content)
  template_content <- gsub('plate_id = "plate_id"', 'plate_id = NA', template_content)
  template_content <- gsub('size_standard = "size_standard"', 'size_standard = NA', template_content)
  template_content <- gsub('size_standard_repeat_length = "size_standard_repeat_length"', 'size_standard_repeat_length = NA', template_content)
}

if (!samples_grouped) {
  template_content <- gsub('metadata\\$group_id <- metadata\\$group_id', '# metadata$group_id <- metadata$group_id', template_content)
  template_content <- gsub('metadata\\$metrics_baseline_control <- metadata\\$metrics_baseline_control', '# metadata$metrics_baseline_control <- metadata$metrics_baseline_control', template_content)
  template_content <- gsub('grouped = TRUE', 'grouped = FALSE', template_content)
  template_content <- gsub('group_id = "group_id"', 'group_id = NA', template_content)
  template_content <- gsub('metrics_baseline_control = "metrics_baseline_control"', 'metrics_baseline_control = NA', template_content)
}

if (!size_standards & !samples_grouped) {
  template_content <- gsub('metadata <- read.csv\\("")', '# metadata <- read.csv("")', template_content)
  template_content <- gsub('fragments_list = metadata_added_list', 'fragments_list = peak_list', template_content)

  # Comment out block of input section
  template_content <- comment_out_lines(template_content, '#Provide the appropriate metadata below by replacing the placeholders', "\\`\\`\\`")

  # Comment out the Add metadata section
  template_content <- comment_out_lines(template_content, "\\`\\`\\`\\{r Add metadata\\}", "\\`\\`\\`", "metadata not used")
}

writeLines(template_content, paste0(file_name, ".qmd"))
file.edit(paste0(file_name, ".qmd"))
}

