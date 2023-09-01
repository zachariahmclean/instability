######################## Helper functions ####################################
# np_repeat algorithim -----------------------------------------------------------
np_repeat <- function(size,
                      main_peak_size,
                      main_peak_repeat,
                      repeat_size) {

  # For loop to get values
  # first find the distance to main peak and the peak before
  size_delta_from_main_peak <- size - main_peak_size
  # create a numeric vector of length of peak to store values
  np_repeat <- vector(mode = 'numeric', length = length(size))
  # note that the loop has to start in the middle at the main peak since it's looking for the previous value in the new vector
  for (i in seq_along(size)) {
    #set main peak
    if (i == which(size_delta_from_main_peak == 0)) {
      np_repeat[[i]] <- main_peak_repeat
    }
    #calculate for peaks greater than main peak
    if (i > which(size_delta_from_main_peak == 0)) {
      # calculate size distance to nearest peak in cag length,
      # add that on to previous cag, then round to whole cag
      np_repeat[[i]] <-
        np_repeat[[i - 1]] + round((size[[i]] - size[[i-1]]) / repeat_size)
    }
  }
  # --- reverse order and find smaller repeats ---
  size_rev <- rev(size)
  np_repeat_rev <-rev(np_repeat)
  size_delta_from_main_peak_rev <- rev(size_delta_from_main_peak)
  #second loop that goes along the vector in the from larger to smaller bp peaks
  for (i in seq_along(size_rev)) {
    #calculate peaks greater than main peak (smaller bp peaks since reversed)
    if (i > which(size_delta_from_main_peak_rev == 0)) {
      # note this is a subtraction since the size_delta_peak_after is a negative value
      np_repeat_rev[[i]] <-
        np_repeat_rev[[i - 1]] + round((size_rev[[i]] - size_rev[[i-1]]) / repeat_size)
    }
  }
  return(rev(np_repeat_rev))
}

# repeat length correction -------------------------------------------------------

model_repeat_length <- function(fragments_list,
                                repeat_size,
                                assay_size_without_repeat,
                                repeat_length_correction) {

  calling_close_neighbouring_repeats <- function(controls_fragments){
    #use np_repeat to accurately call the repeat length of the neighboring peaks
    #extract a dataframe of the called repeats that can then be used to make a model
    controls_fragments_df_list <- lapply(controls_fragments, function(x){
      df_length <- nrow(x$peak_data)
      #identify peaks close to modal peak and at least 20% as high
      main_peak_delta <- x$peak_data$size - x$allele_1_size
      height_prop <- x$peak_data$height / x$allele_1_height
      peak_cluster <- vector("logical", length = nrow(x$peak_data))
      for (i in seq_along(main_peak_delta)) {
        if(abs(main_peak_delta[[i]]) < 30 & height_prop[[i]] > 0.2){
          peak_cluster[[i]] <- TRUE
        } else{
          peak_cluster[[i]] <- FALSE
        }
      }
      cluster_df <- x$peak_data[peak_cluster,]
      cluster_df_length <- nrow(cluster_df)
      # use np_repeat method to accurately call the neighboring repeats
      data.frame(unique_id = rep(x$unique_id, cluster_df_length),
                 size = cluster_df$size,
                 validated_repeats = np_repeat(
                   size = cluster_df$size,
                   main_peak_size = x$allele_1_size,
                   main_peak_repeat = x$size_standard_repeat_length,
                   repeat_size = repeat_size),
                 plate_id = rep(x$plate_id, cluster_df_length)
      )
    })

    controls_repeats_df <- do.call(rbind, controls_fragments_df_list)
  }

  # correct repeats
  if(repeat_length_correction == "from_metadata"){
    ## first pull out a dataframe for all samples with a column that indicates if it's a positive control or not
    extracted <- lapply(fragments_list, function(x) {
      data.frame(
        unique_id = x$unique_id,
        size_standard = x$size_standard,
        allele_1_size = x$allele_1_size,
        plate_id = x$plate_id
      )
    })
    extracted_df <- do.call(rbind, extracted)

    # Check to see if there are controls, if there are none, give error
    if(!any(extracted_df$size_standard == TRUE)){
      stop("No repeat-length control samples were detected. Ensure that the metadata has been added to the samples with 'add_metadata()' and check your metadata to make sure 'TRUE' is indicated in the appropriate column to indicate samples that are to be used for predicting the repeat length",
           call. = FALSE)
    }
    # pull out the controls
    controls_df <- extracted_df[which(extracted_df$size_standard == TRUE), , drop = FALSE]
    controls_fragments <- fragments_list[which(names(fragments_list) %in% controls_df$unique_id)]
    controls_repeats_df <- calling_close_neighbouring_repeats(controls_fragments)

  } else if(repeat_length_correction == "from_genemapper"){
    # Do some checks by identify samples that have genemapper called alleles
    controls_samples_list <- lapply(fragments_list, function(x) x$peak_data[which(!is.na(x$peak_data$allele)), ] )
    controls_samples_df <- do.call(rbind, controls_samples_list)
    if(nrow(controls_samples_df) == 0){
      stop(paste("Correction could not go ahead because no genemapper alleles could be indetified"),
           call. = FALSE)
    }

    # pick the closest peak to the main peak size and temporarily make that allele_1
    controls_fragments <- lapply(fragments_list[which( names(fragments_list) %in% unique(controls_samples_df$unique_id))],
                                 function(x){
                                   genemapper_alleles <- controls_samples_df[which(controls_samples_df$unique_id == x$unique_id), ]
                                   allele_1_delta_abs <- abs(genemapper_alleles$size - x$allele_1_size)
                                   closest_to_allele_1 <- which(allele_1_delta_abs == min(allele_1_delta_abs))
                                   selected_genemapper_allele <- genemapper_alleles[closest_to_allele_1[1],]

                                   #make sure it doesn't modify in place and mess up the selection of the real main peak
                                   y <- x$clone()
                                   y$allele_1_size <- selected_genemapper_allele$size
                                   y$size_standard <- TRUE
                                   y$size_standard_repeat_length <- selected_genemapper_allele$allele
                                   return(y)
                                 })

    controls_repeats_df <- calling_close_neighbouring_repeats(controls_fragments)

  }

  # Check to see if there are controls for each plate, if there are no controls for a plate, give error
  all_plate_ids <- lapply(fragments_list, function(x) x$plate_id )
  control_plate_ids <- unique(controls_repeats_df$plate_id)
  if(length(unique(control_plate_ids)) != length(unique(all_plate_ids))){
    plates_missing_controls <- paste0(all_plate_ids[which(!all_plate_ids %in% control_plate_ids)], collapse = ", ")
    stop(paste("Plate(s)", plates_missing_controls, "have no repeat-length control samples"),
         call. = FALSE)
  }

  message(paste0("Repeat correction model: ", length(unique(controls_repeats_df$unique_id)), " samples used to build model"))

  # Can now make a model based on the bp size and the known repeat size
  if(length(unique(controls_repeats_df$plate_id)) == 1){
    # when there's only one plate just set up simple lm
    correction_mods <- stats::lm(validated_repeats ~ size, data = controls_repeats_df)
    repeat_bp_size <- round(1/correction_mods$coefficients[2], 2)
    message(paste0("Repeat correction model: ", repeat_bp_size, " bp increase per repeat"))

  } else{
    # when there are multiple samples a linear model can be made using the modal peak and the known repeat length of the modal peak
    correction_mods <- lm(validated_repeats ~ size*plate_id, data = controls_repeats_df)
  }

  # check to see if any samples look off

  controls_repeats_df$predicted_repeat <- stats::predict.lm(correction_mods, controls_repeats_df)
  controls_repeats_df$residuals <- correction_mods$residuals
  message(paste0("Repeat correction model: Average repeat residual ", round(mean(controls_repeats_df$residuals), 10)))

  if(any(abs(controls_repeats_df$residuals) > 0.5)){

    message("Repeat correction model: Warning! The following samples may be off and need investigaion")

    samples_all_controls <- unique(controls_repeats_df$unique_id)
    samples_high_diff <- unique(controls_repeats_df[which(abs(controls_repeats_df$residuals ) > 0.5), "unique_id"])
    for (i in seq_along(samples_high_diff)) {
      sample_id <- samples_high_diff[i]
      sample_control_df <- controls_repeats_df[which(controls_repeats_df$unique_id == sample_id), ]
      sample_control_peaks_n <- nrow(sample_control_df)
      sample_control_peaks_off_df <- sample_control_df[which(abs(sample_control_df$residuals ) > 0.5),]
      sample_control_peaks_off_df_n <- nrow(sample_control_peaks_off_df)

      message(paste0(sample_id," has ", sample_control_peaks_off_df_n, "/", sample_control_peaks_n, " peaks used for making model with high residual repeat size (average residual ", round(mean(sample_control_df$residuals), 2), " repeats)"))

    }
  }
  return(list(correction_mods = correction_mods, controls_repeats_df = controls_repeats_df))
}


##################### R6 Class Method Helpers ##################################

add_repeats_helper <- function(bp_fragments,
                               repeat_algorithm,
                               assay_size_without_repeat,
                               repeat_size,
                               correct_repeat_length){
  ##
  ###in this function, need to do the following in the correct order
  #### 1) calculate repeats by assay_size_without_repeat and repeat size
  #### 2) correct repeat length with positive control
  #### 3) use a method to calculate repeats

  # check to make sure all the required inputs for the function have been given
  if(bp_fragments$.__enclos_env__$private$find_main_peaks_used == FALSE){
    stop(paste0(bp_fragments$unique_id, " requires main alleles to be identified before repeats can be called. Find alleles using 'find_main_peaks()' whitin the class, or use the 'find_alleles()' accesesor to find the main peaks across a list of 'HTT_fragments' objects"),
         call. = FALSE)
  } else if(correct_repeat_length == TRUE & is.null(bp_fragments$.__enclos_env__$private$correction_mod)){
    stop("Correcting the repeat length requires a model based on positive controls, so 'correct_repeat_length' & 'correction_mod' inputs are not meant for users to directly use. To correct the repeat length, you need to work on the 'HTT_fragments' objects in a list format and use accessor functions. On a list of 'HTT_fragments' objects, i) use 'add_metadata()' to indicate which samples are positive controls, and ii) use 'find_alleles()' accesesor function to call and correct repeat lengths across all samples",
         call. = FALSE)
  }

  #only continue from here if main peaks were successfully found, otherwise, don't return repeat data (ie it can stay NULL)
  if(is.na(bp_fragments$allele_1_size) | is.na(bp_fragments$allele_1_height)){
    bp_fragments$.__enclos_env__$private$repeats_not_called_reason <- "No main peaks"
    warning(paste0(bp_fragments$unique_id, ": repeats were not called (no main peaks in sample)"),
            call. = FALSE)
    repeat_class <- repeats_fragments$new(
      unique_id = bp_fragments$unique_id)
    #populate with empty dataframe to help the rest of the pipeline
    repeat_class$repeat_data <- data.frame(
      unique_id = character(),
      size = numeric(),
      height = numeric(),
      repeats = numeric()
      )

  }
  else{
    #need to save info that isn't cloned over
    repeat_data <-  data.frame(unique_id = bp_fragments$peak_data$unique_id,
                               size = bp_fragments$peak_data$size,
                               height = bp_fragments$peak_data$height,
                               repeats = (bp_fragments$peak_data$size - assay_size_without_repeat) / repeat_size)

    # Correct repeat length with positive controls
    if(correct_repeat_length == TRUE){

      repeat_data$plate_id <- rep(bp_fragments$plate_id, nrow(repeat_data))
      repeat_data$repeats <- stats::predict.lm(bp_fragments$.__enclos_env__$private$correction_mod, repeat_data)

    }

    # use different repeat calling algorithm
    if(repeat_algorithm == "simple"){
      # don't need to do anything

    } else if(repeat_algorithm == "nearest_peak"){
      repeat_data$repeats <- np_repeat(
        size = bp_fragments$peak_data$size,
        main_peak_size = bp_fragments$allele_1_size,
        main_peak_repeat = repeat_data$repeats[which(bp_fragments$peak_data$size == bp_fragments$allele_1_size)],
        repeat_size = repeat_size)
    }

    ## initialize new class
    repeat_class <- repeats_fragments$new(
      unique_id = bp_fragments$unique_id)

    # Finally save main peak repeat length and repeats data
    repeat_class$allele_1_repeat <- repeat_data$repeats[which(bp_fragments$peak_data$size == bp_fragments$allele_1_size)]
    repeat_class$allele_2_repeat <- repeat_data$repeats[which(bp_fragments$peak_data$size == bp_fragments$allele_2_size)]
    repeat_class$repeat_data <- repeat_data
  }

  # transfer over public items
  public_items <- names(bp_fragments)
  public_items_classes <- as.vector(sapply(bp_fragments, class))
  public_items_to_transfer <- public_items[which(public_items_classes != "function" & public_items != ".__enclos_env__"& public_items != "repeat_data" )]
  public_items_repeats <- names(repeat_class)
  public_items_to_transfer <- public_items_repeats[which(public_items_repeats %in% public_items_to_transfer)]

  for (i in seq_along(public_items_to_transfer)) {
    repeat_class[[ public_items_to_transfer[i] ]] <- bp_fragments[[ public_items_to_transfer[i] ]]
  }

  # transfer over private items
  private_items <- names(bp_fragments$.__enclos_env__$private)
  private_items_classes <- as.vector(sapply(bp_fragments$.__enclos_env__$private, class))
  private_items_to_transfer <- private_items[which(private_items_classes != "function")]
  private_items_repeats <- names(repeat_class$.__enclos_env__$private)
  private_items_to_transfer <- private_items_repeats[which(private_items_repeats %in% private_items_to_transfer)]

  for (i in seq_along(private_items_to_transfer)) {
    repeat_class$.__enclos_env__$private[[ private_items_to_transfer[i] ]] <- bp_fragments$.__enclos_env__$private[[ private_items_to_transfer[i] ]]
  }

  return(repeat_class)
}
