# Fragments class ---------------------------------------------------------

fragments <- R6::R6Class("fragments",
  public = list(
    unique_id = NA_character_,
    plate_id = NA_character_,
    group_id = NA_character_,
    metrics_baseline_control = FALSE,
    size_standard = FALSE,
    size_standard_repeat_length = NA_real_,
    initialize = function(unique_id) {
      if (length(unique_id) != 1) stop("Fragments must have a single unique id", call. = FALSE)
      self$unique_id <- unique_id
    },
    add_metadata = function(metadata_data.frame,
                            unique_id,
                            plate_id,
                            group_id,
                            size_standard,
                            size_standard_repeat_length,
                            metrics_baseline_control) {
      # clone the class so that it doesn't modify in place
      self2 <- self$clone()

      self2 <- add_metadata_helper(
        fragments = self2,
        metadata_data.frame = metadata_data.frame,
        unique_id = unique_id,
        plate_id = plate_id,
        group_id = group_id,
        size_standard = size_standard,
        size_standard_repeat_length = size_standard_repeat_length,
        metrics_baseline_control = metrics_baseline_control
      )

      return(self2)
    },
    print = function(...) {
      print_helper(self,
        sample_attrs = c("unique_id", "plate_id", "group_id", "metrics_baseline_control", "size_standard", "size_standard_repeat_length")
      )
    },
    plot_trace = function(show_peaks = TRUE,
                          x_axis = NULL,
                          ylim = NULL,
                          xlim = NULL,
                          height_color_threshold = 0.05,
                          plot_title = NULL) {
      if (is.null(self$trace_bp_df)) {
        stop(
          call. = FALSE,
          paste(self$unique_id, "This sample does not have trace data. Use fsa files as inputs to pipeline to plot trace.")
        )
      }

      #there must be a simpler way of the following if else below
      if (is.null(x_axis) && is.null(self$repeat_table_df)) {
        data <- self$trace_bp_df
        data$x <- data$size
        x_axis_label <- "Size"
      } else if (is.null(x_axis) && !is.null(self$repeat_table_df)) {
        data <- self$trace_bp_df
        data$x <- data$calculated_repeats
        x_axis_label <- "Repeats"
      } else if (x_axis == "size") {
        data <- self$trace_bp_df
        data$x <- data$size
        x_axis_label <- "Size"
      } else {
        data <- self$trace_bp_df
        data$x <- data$calculated_repeats
        x_axis_label <- "Repeats"
      }

      if (!is.null(xlim)) {
        data <- data[which(data$x < xlim[2] & data$x > xlim[1]), ]
      }

      plot(data$x,
        data$signal,
        main = ifelse(is.null(plot_title), self$unique_id, plot_title),
        type = "l",
        xlab = x_axis_label,
        ylab = "Signal",
        ylim = ylim
      )


      if (any(data$off_scale)) {
        abline(v = data[which(data$off_scale), "x"], col = adjustcolor("red", alpha = 0.3), lwd = 2.5)
      }

      # add points onto plot showing peaks
      if (!is.null(self$peak_table_df) && show_peaks) {
        if (is.null(x_axis) && is.null(self$repeat_table_df)) {
          peak_table <- self$peak_table_df
          peak_table$x <- peak_table$size
        } else if (is.null(x_axis) && !is.null(self$repeat_table_df)) {
          peak_table <- self$repeat_table_df
          peak_table$x <- peak_table$repeats
        } else if (x_axis == "size") {
          peak_table <- self$peak_table_df
          peak_table$x <- peak_table$size
        } else {
          peak_table <- self$repeat_table_df
          peak_table$x <- peak_table$repeats
        }

        #exit early if the peak table is empty
        if(nrow(peak_table) == 0){
          return()
        }

        if (!is.null(xlim)) {
          peak_table <- peak_table[which(peak_table$x < xlim[2] & peak_table$x > xlim[1]), ]
        }

        tallest_peak_height <- peak_table[which(peak_table$height == max(peak_table$height)), "height"]
        tallest_peak_x <- peak_table[which(peak_table$height == tallest_peak_height), "x"]
        if (!is.null(self$allele_1_height) && !is.na(self$allele_1_height)) {
          tallest_peak_height <- self$allele_1_height
          #find the tallest peak x axis position
          if (is.null(x_axis) && is.na(self$allele_1_repeat)) {
            tallest_peak_x <- self$allele_1_size
          } else if (is.null(x_axis) && !is.na(self$allele_1_repeat)) {
            tallest_peak_x <- self$allele_1_repeat
          } else if (x_axis == "size") {
            tallest_peak_x <- self$allele_1_size
          } else {
            tallest_peak_x <- self$allele_1_repeat
          }
        }

        peaks_above <- peak_table[which(peak_table$height > tallest_peak_height * height_color_threshold), ]
        peaks_below <- peak_table[which(peak_table$height < tallest_peak_height * height_color_threshold), ]

        # Adding peaks
        points(peaks_above$x,
          peaks_above$height,
          col = "blue"
        )
        points(peaks_below$x,
          peaks_below$height,
          col = "purple"
        )
        points(tallest_peak_x,
               tallest_peak_height,
               col = "green")

        # Draw horizontal dotted lines to connect repeats to their actual place on the plot
        if (!is.null(peak_table$repeats) && !is.null(peak_table$calculated_repeats)) {
          for (i in 1:nrow(peak_table)) {
            segments(x0 = peak_table$repeats[i],
                     y0 = peak_table$height[i],
                     x1 = peak_table$calculated_repeats[i],
                     y1 = peak_table$height[i],
                     lty = 2)
          }
        }

      }


      if(!is.null(self$index_repeat) && !is.na(self$index_repeat)){
        abline(v = self$index_repeat, col = "black", lwd = 2, lty = 3)
      }

    }
  ),
  private = list(
    min_bp_size = NULL,
    max_bp_size = NULL,
    find_main_peaks_used = FALSE,
    peak_regions = NA_real_,
    repeats_not_called_reason = NA_character_,
    validated_peaks_df = NULL,
    correction_mod = NULL,
    controls_repeats_df = NULL,
    assigned_index_peak_used = FALSE,
    index_samples = NULL
  )
)



# fragments_trace class ------------------------------------------------------------


fragments_trace <- R6::R6Class(
  "fragments_trace",
  inherit = fragments,
  public = list(
    unique_id = NULL,
    raw_ladder = NULL,
    raw_data = NULL,
    scan = NULL,
    off_scale_scans = NULL,
    ladder_df = NULL,
    trace_bp_df = NULL,
    peak_table_df = NULL,
    mod_parameters = NULL,
    find_ladder = function(fsa,
                           ladder_channel = "DATA.105",
                           signal_channel = "DATA.1",
                           ladder_sizes = c(50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                           spike_location = NULL,
                           zero_floor = FALSE,
                           scan_subset = NULL,
                           smoothing_window = 21,
                           minimum_peak_signal = NULL,
                           max_combinations = 2500000,
                           ladder_selection_window = 5,
                           show_progress_bar = TRUE) {

      if(any(sapply(list(self$raw_ladder, self$raw_data, self$scan, self$off_scale_scans), is.null))){
        self$raw_ladder <- fsa$Data[[ladder_channel]]
        self$raw_data <- fsa$Data[[signal_channel]]
        self$scan <- 0:(length(fsa$Data[[signal_channel]]) - 1)
        self$off_scale_scans <- fsa$Data$OfSc.1
      }

      self2 <- find_ladder_helper(
        fragments_trace = self,
        ladder_channel = ladder_channel,
        signal_channel = signal_channel,
        ladder_sizes = ladder_sizes,
        spike_location = spike_location,
        zero_floor = zero_floor,
        scan_subset = scan_subset,
        smoothing_window = smoothing_window,
        minimum_peak_signal = minimum_peak_signal,
        max_combinations = max_combinations,
        ladder_selection_window = ladder_selection_window
      )

      return(self2)
    },
    call_peaks = function(smoothing_window = 4,
                          minimum_peak_signal = 20,
                          min_bp_size = 100,
                          max_bp_size = 1000,
                          ...) {
      df <- find_fragment_peaks(self$trace_bp_df,
        smoothing_window = smoothing_window,
        minimum_peak_signal = minimum_peak_signal,
        ...
      )

      df$unique_id <- rep(self$unique_id, nrow(df))
      self$peak_table_df <- df[which(df$size > min_bp_size & df$size < max_bp_size), ]

      invisible(self)
    },
    ladder_correction_auto = function(size_threshold = 60,
                                      size_tolerance = 2.5,
                                      rsq_threshold = 0.9985) {
      self2 <- ladder_self_mod_predict(self,
        size_threshold = size_threshold,
        size_tolerance = size_tolerance,
        rsq_threshold = rsq_threshold
      )


      return(self2)
    },
    ladder_correction_manual = function(replacement_ladder_df) {
      fixed_fragments_trace <- ladder_fix_helper(
        self,
        replacement_ladder_df = replacement_ladder_df
      )

      return(fixed_fragments_trace)
    },
    plot_ladder = function(xlim = NULL, ylim = NULL,
                           plot_title = NULL) {
      # Scatter plot
      plot(self$trace_bp_df$scan, self$trace_bp_df$ladder_signal,
        xlab = "Scan", ylab = "Ladder Signal",
        main = ifelse(is.null(plot_title), self$unique_id, plot_title),
        type = "l",
        xlim = xlim,
        ylim = ylim
      )

      # Adding text
      text(self$ladder_df$scan, rep(max(self$trace_bp_df$ladder_signal) / 3, nrow(self$ladder_df)),
        labels = self$ladder_df$size,
        adj = 0.5, cex = 0.7, srt = 90
      )

      # Adding vertical lines with transparency
      for (i in 1:nrow(self$ladder_df)) {
        abline(
          v = self$ladder_df$scan[i],
          lty = 3,
          col = rgb(1, 0, 0, alpha = 0.3)
        )
      }
    }
  )
)


# repeats -------------------------------
# responsibility if this class is to calculate the instability metrics from non-continuous data



fragments_repeats <- R6::R6Class(
  "fragments_repeats",
  inherit = fragments,
  public = list(
    allele_1_size = NA_real_,
    allele_1_repeat = NA_real_,
    allele_1_height = NA_real_,
    allele_2_size = NA_real_,
    allele_2_repeat = NA_real_,
    allele_2_height = NA_real_,
    trace_bp_df = NULL,
    peak_table_df = NULL,
    repeat_table_df = NULL,
    index_repeat = NA_real_,
    index_height = NA_real_,
    find_main_peaks = function(number_of_peaks_to_return = 1,
                               peak_region_size_gap_threshold = 6,
                               peak_region_height_threshold_multiplier = 1) {
      # clone the class so that it doesn't modify in place
      self2 <- self$clone()

      self2 <- find_main_peaks_helper(
        fragments_repeats_class = self2,
        number_of_peaks_to_return = number_of_peaks_to_return,
        peak_region_size_gap_threshold = peak_region_size_gap_threshold,
        peak_region_height_threshold_multiplier = peak_region_height_threshold_multiplier
      )


      # finally, indicate in the private part of the class that this function has been used since that is required for next steps
      self2$.__enclos_env__$private$find_main_peaks_used <- TRUE


      return(self2)
    },
    add_repeats = function(assay_size_without_repeat = 87,
                           repeat_size = 3,
                           repeat_calling_algorithm  = "simple",
                           repeat_calling_algorithm_size_window_around_allele =  repeat_size * 5,
                           repeat_calling_algorithm_peak_assignment_scan_window = 3,
                           repeat_calling_algorithm_size_period = repeat_size * 0.93 ,
                           force_whole_repeat_units = FALSE,
                           correct_repeat_length = FALSE){
      repeat_class <- add_repeats_helper(
        self,
        assay_size_without_repeat = assay_size_without_repeat,
        repeat_size = repeat_size,
        repeat_calling_algorithm = repeat_calling_algorithm,
        repeat_calling_algorithm_size_window_around_allele = repeat_calling_algorithm_size_window_around_allele,
        repeat_calling_algorithm_peak_assignment_scan_window = repeat_calling_algorithm_peak_assignment_scan_window,
        repeat_calling_algorithm_size_period = repeat_calling_algorithm_size_period,
        force_whole_repeat_units = force_whole_repeat_units,
        correct_repeat_length = correct_repeat_length
      )

      return(repeat_class)
    },
    instability_metrics = function(peak_threshold = 0.05,
                                   window_around_index_peak = c(NA, NA), # note the lower lim should be a negative value
                                   percentile_range = c(
                                     0.01, 0.05, 0.1, 0.2, 0.3,
                                     0.4, 0.5, 0.6, 0.7, 0.8,
                                     0.9, 0.95, 0.99
                                   ),
                                   repeat_range = c(
                                     1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20
                                   )) {
      # check to make sure all the required inputs for the function have been given
      if(private$assigned_index_peak_used == FALSE){
        stop(paste0(self$unique_id, " requires an index peak to calculate repeat instability metrics. Use 'assign_index_peaks' to set the index peaks."),
             call. = FALSE
        )
      }else if (is.na(self$allele_1_repeat)) {
        message(paste0(self$unique_id, ": metrics not calculated (no main peaks in sample)"))
        return(NULL)
      }

      # set index repeat if it hasn't already been set
      if (is.na(self$index_repeat)) {
        self$index_repeat <- self$allele_1_repeat
        self$index_height <- self$allele_1_height
      }

      # compute metrics
      metrics <- compute_metrics(
        self,
        peak_threshold = peak_threshold,
        window_around_index_peak = window_around_index_peak,
        percentile_range = percentile_range,
        repeat_range = repeat_range
      )

      return(metrics)
    },
    plot_fragments = function( # show_peak_regions = FALSE,
                              ylim = NULL,
                              xlim = NULL,
                              plot_title = NULL) {
      if (is.null(self$repeat_table_df)) {
        data <- self$peak_table_df
        data$x <- data$size
      } else {
        data <- self$repeat_table_df
        data$x <- data$repeats
      }

      if (nrow(data) == 0) {
        plot.new()
        title(main = self$unique_id)
        return()
      }

      if (!is.null(xlim)) {
        if (length(xlim == 2) & is(xlim, "numeric")) {
          data <- data[which(data$x < xlim[2] & data$x > xlim[1]), ]
        } else {
          stop(
            call. = FALSE,
            "xlim must be a numeric vector with length of 2"
          )
        }
      }


      allele_1_mode <- ifelse(is.null(self$repeat_table_df), round(self$allele_1_size), round(self$allele_1_repeat))
      allele_2_mode <- ifelse(is.null(self$repeat_table_df), round(self$allele_2_size), round(self$allele_2_repeat))

      # Fill missing y values with zeros
      rounded_x <- round(data$x)
      all_x_values <- seq(min(rounded_x), max(rounded_x))
      y_values <- rep(0, length(all_x_values))
      for (i in seq_along(rounded_x)) {
        y_values[which(all_x_values == rounded_x[i])] <- data[which(data$x == data$x[i]), "height"]
      }



      barplot(
        names.arg = all_x_values,
        height = y_values,
        main = ifelse(is.null(plot_title), self$unique_id, plot_title),
        xlab = ifelse(is.null(self$repeat_table_df), "Size", "Repeat"),
        ylab = "Signal",
        ylim = ylim,
        beside = TRUE,
        col = sapply(all_x_values, function(x) if (!is.na(allele_1_mode) && x == allele_1_mode) "red" else if (!is.na(allele_2_mode) && x == allele_2_mode) "blue" else "gray")
      )

      # why doesn't the following code work? can't even get rectangle to show up on its own
      # if(show_peak_regions == TRUE && !is.na(allele_1_mode)) {
      #   data$peak_region <- private$peak_regions
      #   unique_peak_regions <- unique(na.omit(private$peak_regions))
      #   for (i in seq_along(unique_peak_regions)) {
      #     tmp_df <- data[which(data$peak_region == unique_peak_regions[i]), ]
      #     rect(min(tmp_df$x), 0, max(tmp_df$x), max(tmp_df$height), col = palette.colors(palette  = "ggplot2")[2:8][i])
      #   }
      #
      # }
    }
  )
)
