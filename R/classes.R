# Fragments class ---------------------------------------------------------

fragments <- R6::R6Class("fragments",
  public = list(
    unique_id = NA_character_,
    plate_id = NA_character_,
    group_id = NA_character_,
    metrics_baseline_control = FALSE,
    size_standard = FALSE,
    size_standard_sample_id = NA_character_,
    size_standard_repeat_length = NA_real_,
    initialize = function(unique_id) {
      if (length(unique_id) != 1) stop("Fragments must have a single unique id", call. = FALSE)
      self$unique_id <- unique_id
    },
    print = function(...) {
      print_helper(self,
        sample_attrs = c("unique_id", "plate_id", "group_id", "metrics_baseline_control", "size_standard", "size_standard_sample_id", "size_standard_repeat_length")
      )
    },
    plot_trace = function(show_peaks = TRUE,
                          x_axis = NULL,
                          ylim = NULL,
                          xlim = NULL,
                          height_color_threshold = 0.05,
                          plot_title = NULL) {
      plot_trace_helper(
        fragments = self,
        show_peaks = show_peaks,
        x_axis = x_axis,
        ylim = ylim,
        xlim = xlim,
        height_color_threshold = height_color_threshold,
        plot_title = plot_title)

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
# responsibility if this class is to process continuous scan level data


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
    plot_ladder = function(xlim = NULL, ylim = NULL,
                           plot_title = NULL) {
      plot_ladder_helper(self,
                         xlim = xlim, ylim = ylim,
                         plot_title = plot_title)
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
    plot_fragments = function(ylim = NULL,
                              xlim = NULL,
                              plot_title = NULL) {
      plot_fragments_helper(self,
                            ylim = ylim,
                            xlim = xlim,
                            plot_title = plot_title)

    }
  )
)
