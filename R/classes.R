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
    repeat_correction_mod = NULL,
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
    fsa = NULL,
    raw_ladder = NULL,
    raw_data = NULL,
    scan = NULL,
    off_scale_scans = NULL,
    ladder_df = NULL,
    trace_bp_df = NULL,
    peak_table_df = NULL,
    local_southern_mod = NULL,
    initialize = function(
      unique_id, 
      fsa_file,
      ladder_channel,
      signal_channel) {
        if (length(unique_id) != 1) stop("Fragments must have a single unique id", call. = FALSE)
        self$unique_id <- unique_id
        self$fsa <- fsa_file
        self$raw_ladder <- self$fsa$Data[[ladder_channel]]
        self$raw_data <- self$fsa$Data[[signal_channel]]
        self$scan <- 0:(length(self$fsa$Data[[signal_channel]]) - 1)
        self$off_scale_scans <- self$fsa$Data$OfSc.1
    },
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
    trace_bp_df = NULL,
    peak_table_df = NULL,
    repeat_table_df = NULL,
    get_alleles = function(){
      alleles <- list(
        allele_1_size = private$allele_1_size,
        allele_1_height = private$allele_1_height,
        allele_1_repeat = private$allele_1_repeat,
        allele_2_size = private$allele_2_size,
        allele_2_height = private$allele_2_height,
        allele_2_repeat = private$allele_2_repeat
      )
      return(alleles)
    },
    set_allele = function(allele, unit, value){

      if(!is.na(value)){
        if(is.null(self$repeat_table_df)){
          if(unit != "size") stop("Only size can be used to set alleles if repeats have not been called", call. = FALSE )
          df <- self$peak_table_df
        } else{
          df <- self$repeat_table_df
        }

        size_diff <- df[[unit]]- value
        allele_df <- df[which.min(abs(size_diff)), , drop = FALSE]

        if(nrow(allele_df) > 1){
          stop("More than one peak was selected with the value provided", call. = FALSE)
        }

        # Ensure the allele is either 1 or 2
        if (!(allele %in% c(1, 2))) {
          stop("Invalid 'allele' input. Please select between 1 or 2", call. = FALSE)
        }
      }
      # Dynamically construct the variable names and assign values
      private[[paste0("allele_", allele, "_size")]] <- ifelse(!is.na(value), allele_df$size, NA_real_)
      private[[paste0("allele_", allele, "_height")]] <- ifelse(!is.na(value), allele_df$height, NA_real_)
      private[[paste0("allele_", allele, "_repeat")]] <- ifelse(!is.null(self$repeat_table_df) && !is.na(value), allele_df$repeats, NA_real_)      
      private$find_main_peaks_used <- TRUE

      invisible(self)
    },
    get_index_peak = function(){
      index <- list(
        index_repeat = private$index_repeat,
        index_height = private$index_height
      )
      return(index)
    },
    set_index_peak = function(value){
      if(is.null(self$repeat_table_df)){
        stop("Index assignment requires repeats to be called", call. = FALSE )
      }

      if(!is.na(value)){
        size_diff <- self$repeat_table_df$repeats- value
        index_df <- self$repeat_table_df[which.min(abs(size_diff)), , drop = FALSE]

        if(nrow(index_df) > 1){
          stop("More than one peak was selected with the value provided", call. = FALSE)
        }
      }
      private$index_repeat <- ifelse(!is.na(value), index_df$repeats, NA_real_)
      private$index_height <- ifelse(!is.na(value), index_df$height, NA_real_)
      private$assigned_index_peak_used <- TRUE
      
      invisible(self)
    },
    plot_fragments = function(ylim = NULL,
                              xlim = NULL,
                              plot_title = NULL) {
      plot_fragments_helper(self,
                            ylim = ylim,
                            xlim = xlim,
                            plot_title = plot_title)

    }
  ),
  private = list(
    allele_1_size = NA_real_,
    allele_1_repeat = NA_real_,
    allele_1_height = NA_real_,
    allele_2_size = NA_real_,
    allele_2_repeat = NA_real_,
    allele_2_height = NA_real_,
    index_repeat = NA_real_,
    index_height = NA_real_
  )
)
