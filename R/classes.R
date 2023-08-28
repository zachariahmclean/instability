

# Fragments class ---------------------------------------------------------

fragments <- R6::R6Class("fragments", public = list(
  unique_id = NA_character_,
  plate_id = NA_character_,
  group_id = NA_character_,
  size_standard = FALSE,
  size_standard_repeat_length = NA_real_,
  metrics_baseline_control = FALSE,
  peak_data = NULL,
  repeat_data = NULL,

  initialize = function(unique_id,
                        peak_data = NULL,
                        repeat_data = NULL,
                        plate_id = NA_character_,
                        group_id = NA_character_,
                        size_standard = FALSE,
                        size_standard_repeat_length = NA_real_,
                        metrics_baseline_control = FALSE) {

    stopifnot(length(unique_id) == 1)

    self$unique_id <- unique_id
    self$plate_id <- plate_id
    self$group_id <- group_id
    self$size_standard <- size_standard
    self$size_standard_repeat_length <- size_standard_repeat_length
    self$metrics_baseline_control <- metrics_baseline_control
    self$peak_data <- peak_data
    self$repeat_data <- repeat_data
  },

  add_metadata = function(metadata_data.frame,
                          unique_id = "sample_file_name",
                          plate_id = "plate_id",
                          sample_group_id = "sample_group_id",
                          repeat_positive_control_TF = "repeat_positive_control_TF",
                          repeat_positive_control_length = "repeat_positive_control_length",
                          metrics_baseline_control = "metrics_baseline_control_TF"){
    # validate inputs to give good errors to user
    function_input_vector <- c(unique_id, plate_id, sample_group_id, repeat_positive_control_TF, repeat_positive_control_length, metrics_baseline_control)
    function_input_name_vector <- c("unique_id", "plate_id", "sample_group_id", "repeat_positive_control_TF", "repeat_positive_control_length", "metrics_baseline_control")
    for (i in seq_along(function_input_vector)) {
      if(!any(names(metadata_data.frame) == function_input_vector[[i]])){
        stop(paste0(function_input_name_vector[[i]], " input '", function_input_vector[[i]], "' was not detected as a column name in the 'metadata_data.frame'. Check column names and supply the right character string for the ",function_input_name_vector[[i]]," input"),
             call. = FALSE)
      }
    }

    # clone the class so that it doesn't modify in place
    self2 <- self$clone()

    # filter for row of sample
    sample_metadata <- metadata_data.frame[which(metadata_data.frame[unique_id] == self2$unique_id), ,drop = FALSE]

    # add metadata to slots
    self2$plate_id <- as.character(sample_metadata[plate_id])
    self2$group_id <- as.character(sample_metadata[sample_group_id])
    self2$size_standard <- as.logical(sample_metadata[repeat_positive_control_TF]) #give a better error if this coercion isn't possible
    self2$size_standard_repeat_length <- as.double(sample_metadata[repeat_positive_control_length])
    self2$metrics_baseline_control <- as.character(sample_metadata[metrics_baseline_control])

    return(self2)
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

# HTT_fragments -------------------------------
#responsibility of this class is to turn the bp level peak table into called repeats

bp_fragments <- R6::R6Class("bp_fragments",
                         inherit = fragments,
                         public = list(
  allele_1_size = NA_real_,
  allele_1_height = NA_real_,
  allele_2_size = NA_real_,
  allele_2_height = NA_real_,

  find_main_peaks = function(number_of_peaks_to_return = 2,
                             peak_region_size_gap_threshold = 6,
                             peak_region_height_threshold_multiplier = 1){
    # clone the class so that it doesn't modify in place
    self2 <- self$clone()

    self2 <- find_main_peaks_helper(
      fragments_class = self2,
      fragment_sizes = self2$peak_data$size,
      fragment_heights = self2$peak_data$height,
      data_type = "bp_size",
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
  }
  )
)



# repeats -------------------------------
# responsibility if this class is to calculate the instability metrics from a repeat table



repeats_fragments <- R6::R6Class(
  "repeats_fragments",
  inherit = fragments,
  public = list(
    allele_1_repeat = NA_real_,
    allele_1_height = NA_real_,
    allele_2_repeat = NA_real_,
    allele_2_height = NA_real_,
    index_repeat = NA_real_,
    index_height = NA_real_,
    index_weighted_mean_repeat = NA_real_,

    find_main_peaks = function(number_of_peaks_to_return = 2,
                               peak_region_size_gap_threshold = 2, #this is the only thing different with the find main peaks above
                               peak_region_height_threshold_multiplier = 1) {
      # clone the class so that it doesn't modify in place
      self2 <- self$clone()

      self2 <- find_main_peaks_helper(
        fragments_class = self2,
        fragment_sizes = self2$repeat_data$repeats,
        fragment_heights = self2$peak_data$height,
        data_type = "repeat",
        number_of_peaks_to_return = number_of_peaks_to_return,
        peak_region_size_gap_threshold = peak_region_size_gap_threshold,
        peak_region_height_threshold_multiplier = peak_region_height_threshold_multiplier
      )


      #finally, indicate in the private part of the class that this function has been used since that is required for next steps
      self2$.__enclos_env__$private$find_main_peaks_used <- TRUE


      return(self2)

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

      # filter dataset to user supplied thresholds
      filtered_df <- repeat_table_subset(repeat_data = self$repeat_data,
                                         allele_1_height = self$allele_1_height,
                                         index_repeat =  self$index_repeat,
                                         peak_threshold = peak_threshold,
                                         window_around_main_peak = window_around_main_peak)

      # compute metrics
      metrics <- compute_metrics(
        self,
        filtered_df,
        peak_threshold = peak_threshold,
        window_around_main_peak = window_around_main_peak,
        percentile_range = percentile_range,
        repeat_range = repeat_range)

      return(metrics)

    }
  )
)

