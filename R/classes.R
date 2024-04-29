# Fragments class ---------------------------------------------------------

fragments <- R6::R6Class("fragments", public = list(
  unique_id = NA_character_,
  plate_id = NA_character_,
  group_id = NA_character_,
  metrics_baseline_control = FALSE,
  size_standard = FALSE,
  size_standard_repeat_length = NA_real_,

  initialize = function(unique_id) {
    if(length(unique_id) != 1) stop("Fragments must have a single unique id", call. = FALSE)
    self$unique_id <- unique_id
  },

  add_metadata = function(metadata_data.frame,
                          unique_id,
                          plate_id,
                          group_id,
                          size_standard,
                          size_standard_repeat_length,
                          metrics_baseline_control){

    # clone the class so that it doesn't modify in place
    self2 <- self$clone()

    # filter for row of sample
    sample_metadata <- metadata_data.frame[which(metadata_data.frame[unique_id] == self2$unique_id), ,drop = FALSE]

    # add metadata to slots
    self2$plate_id <- ifelse(!is.na(plate_id), as.character(sample_metadata[plate_id]), NA_character_)
    self2$group_id <- ifelse(!is.na(group_id), as.character(sample_metadata[group_id]), NA_character_)
    self2$size_standard <- ifelse(!is.na(size_standard),
                                  ifelse(is.na(sample_metadata[size_standard]) || !as.logical(sample_metadata[size_standard]), FALSE, TRUE),
                                  FALSE)
    self2$size_standard_repeat_length <- ifelse(!is.na(size_standard_repeat_length), as.double(sample_metadata[size_standard_repeat_length]), NA_real_)
    self2$metrics_baseline_control <- ifelse(!is.na(metrics_baseline_control),
                                             ifelse(is.na(sample_metadata[metrics_baseline_control]) || !as.logical(sample_metadata[metrics_baseline_control]), FALSE, TRUE),
                                             FALSE)

    return(self2)
  },
  print = function(...) {
    print_helper(self,
                 sample_attrs = c("unique_id", "plate_id", "group_id", "metrics_baseline_control", "size_standard", "size_standard_repeat_length")
    )
  },
  plot_trace = function(show_peaks = TRUE,
                         ylim = NULL,
                         xlim = NULL){
    if(is.null(self$trace_bp_df)){
      stop(call. = FALSE,
           paste(self$unique_id, "This sample does not have trace data, so cannot make plot"))
    }

    plot(self$trace_bp_df$size,
         self$trace_bp_df$signal,
         main = self$unique_id,
         type = "l",
         xlab = "Size",
         ylab = "Signal",
         ylim = ylim,
         xlim = xlim)

    if(!is.null(self$peak_table_df) && show_peaks){
      # Adding peaks
      points(self$peak_table_df$size,
             self$peak_table_df$height,
             col = "blue")
    }
  }


),
private = list(
  find_main_peaks_used = FALSE,
  repeats_not_called_reason = NA_character_,
  validated_peaks_df = NULL,
  correction_mod = NULL,
  controls_repeats_df = NULL,
  peak_regions = NA_real_
))



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

    #model related
    mod_parameters = NULL,
    generate_mod_parameters = function() {
      # Perform any necessary calculations to fit the model and save the parameters
      ladder_df <- self$ladder_df[which(!is.na(self$ladder_df$size)), ]
      ladder_df <- ladder_df[which(!is.na(ladder_df$scan)), ]
      self$mod_parameters <- local_southern_fit(ladder_df$scan, ladder_df$size)
      invisible(self)
    },
    predict_size = function() {
      # Predict fragment sizes for new data points
      predicted_sizes <- local_southern_predict(local_southern_fit = self$mod_parameters , scans = self$scan)

      return(predicted_sizes)
    },
    call_peaks = function(smoothing_window = 4,
                          minimum_peak_signal = 20,
                          min_bp_size = 100,
                          max_bp_size = 1000,
                          ...){
      df <- find_fragment_peaks(self$trace_bp_df,
                                smoothing_window = smoothing_window,
                                minimum_peak_signal = minimum_peak_signal,
                                ...)

      df$unique_id <- rep(self$unique_id, nrow(df))
      self$peak_table_df <- df[which(df$size > min_bp_size & df$size < max_bp_size), ]

      invisible(self)
    },
    ladder_correction_auto = function(size_threshold = 60,
                                      size_tolerance = 2.5,
                                      rsq_threshold = 0.9985){

      self2 <- ladder_self_mod_predict(self,
                              size_threshold = size_threshold,
                              size_tolerance = size_tolerance,
                              rsq_threshold = rsq_threshold)


      return(self2)
    },
    ladder_correction_manual = function(replacement_ladder_df){
      fixed_fragments_trace <- ladder_fix_helper(
        self,
        replacement_ladder_df = replacement_ladder_df
      )

      return(fixed_fragments_trace)
    },
    plot_ladder = function(xlim = NULL, ylim = NULL){
        # Scatter plot
        plot(self$trace_bp_df$scan, self$trace_bp_df$ladder_signal,
             xlab = "Scan", ylab = "Ladder Signal",
             main = self$unique_id,
             pch = 16,
             xlim = xlim,
             ylim = ylim)

        # Adding text
        text(self$ladder_df$scan, rep(max(self$trace_bp_df$ladder_signal) / 3, nrow(self$ladder_df)),
             labels = self$ladder_df$size,
             adj = 0.5, cex = 0.7, srt = 90)

        # Adding vertical lines with transparency
        for (i in 1:nrow(self$ladder_df)) {
          abline(v = self$ladder_df$scan[i],
                 lty = 3,
                 col = rgb(1, 0, 0, alpha = 0.3))
        }



    }
  ))


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
    index_weighted_mean_repeat = NA_real_,

    find_main_peaks = function(number_of_peaks_to_return = 2,
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


      #finally, indicate in the private part of the class that this function has been used since that is required for next steps
      self2$.__enclos_env__$private$find_main_peaks_used <- TRUE


      return(self2)

    },
    add_repeats = function(repeat_algorithm = "simple", # "simple" or "nearest_peak"
                           assay_size_without_repeat = 87,
                           repeat_size = 3,
                           correct_repeat_length = FALSE){

      repeat_class <- add_repeats_helper(self,
                                         repeat_algorithm = repeat_algorithm,
                                         assay_size_without_repeat = assay_size_without_repeat,
                                         repeat_size = repeat_size,
                                         correct_repeat_length = correct_repeat_length)

      return(repeat_class)
    },

    instability_metrics = function(peak_threshold = 0.05,
                                   window_around_main_peak = c(NA, NA), # note the lower lim should be a negative value
                                   percentile_range = c(
                                     0.01, 0.05, 0.1, 0.2, 0.3,
                                     0.4, 0.5, 0.6, 0.7, 0.8,
                                     0.9, 0.95, 0.99),
                                   repeat_range = c(
                                     1,2,3,4,6,8,10,12,14,16,18,20
                                   )
    ) {

      # check to make sure all the required inputs for the function have been given
      if(private$find_main_peaks_used == FALSE){
        stop(paste0(self$unique_id, " requires main alleles to be identified before repeats can be called. Find alleles using 'find_main_peaks()' whitin the class, or use the 'find_alleles()' accesesor to find the main peaks across a list of 'HTT_fragments' objects"),
             call. = FALSE)
      } else if(is.na(self$allele_1_repeat)){
        message(paste0(self$unique_id, ": metrics not calculated (no main peaks in sample)"))
        return(NULL)
      }

      #set index repeat if it hasn't already been set
      if(is.na(self$index_repeat)){
        self$index_repeat <- self$allele_1_repeat
        self$index_height <- self$allele_1_height
      }

      # compute metrics
      metrics <- compute_metrics(
        self,
        peak_threshold = peak_threshold,
        window_around_main_peak = window_around_main_peak,
        percentile_range = percentile_range,
        repeat_range = repeat_range)

      return(metrics)

    },
    plot_fragments = function(
    #show_peak_regions = FALSE,
                              ylim = NULL,
                              xlim = NULL){

      if(is.null(self$repeat_table_df)){
        data <- self$peak_table_df
        data$x <- data$size
      }
      else{
        data <- self$repeat_table_df
        data$x <- data$repeats
      }

      if(nrow(data) == 0){
        plot.new()
        title(main = self$unique_id)
        return()
      }

      if(!is.null(xlim)){
        if(length(xlim == 2) & class(xlim) == "numeric"){
          data <- data[which(data$x < xlim[2] & data$x > xlim[1]), ]
        }
        else{
          stop(call. = FALSE,
               "xlim must be a numeric vector with length of 2")
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
        main = self$unique_id,
        xlab = ifelse(is.null(self$repeat_table_df), "Size", "Repeat"),
        ylab = "Signal",
        ylim = ylim,
        beside = TRUE,
        col = sapply(all_x_values, function(x) if(!is.na(allele_1_mode) && x == allele_1_mode) "red" else if(!is.na(allele_2_mode) && x == allele_2_mode) "blue" else "gray")
      )

      #why doesn't the following code work? can't even get rectangle to show up on its own
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

