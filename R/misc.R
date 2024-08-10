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


# metadata ----------------------------------------------------------------

add_metadata_helper <- function(
    fragments,
    metadata_data.frame,
    unique_id,
    plate_id,
    group_id,
    size_standard,
    size_standard_repeat_length,
    metrics_baseline_control) {
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
  fragments$size_standard_repeat_length <- if (!is.na(size_standard_repeat_length)) as.double(sample_metadata[[size_standard_repeat_length]]) else NA_real_
  fragments$metrics_baseline_control <- if (!is.na(metrics_baseline_control)) {
    ifelse(is.na(sample_metadata[[metrics_baseline_control]]) || !as.logical(sample_metadata[[metrics_baseline_control]]), FALSE, TRUE)
  } else {
    FALSE
  }

  return(fragments)
}

transfer_metadata_helper <- function(old_fragment,
                                     new_fragment) {
  metadata_names <- c(
    "unique_id",
    "plate_id",
    "group_id",
    "size_standard",
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

