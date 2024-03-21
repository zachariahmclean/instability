
# peak_table_to_fragments -------------------------------------------------

testthat::test_that("peak_table_to_fragments",{

  gm_raw <- instability::example_data

  test_fragments <- peak_table_to_fragments(
    gm_raw,
    data_format = "genemapper5",
    peak_size_col = "Size",
    peak_height_col = "Height",
    unique_id = "Sample.File.Name",
    dye_col = "Dye.Sample.Peak",
    dye_channel = "B",
    allele_col = "Allele",
    min_size_bp = 100,
    max_size_bp = 1000
  )

  test_sample <- test_fragments[[1]]

  test_fragments_classes <- vector("character", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_classes[i] <- class(test_fragments[[i]])[1]
  }

  test_fragments_ids <- vector("character", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_ids[i] <- test_fragments[[i]]$unique_id
  }

  # it's a list
  testthat::expect_true(class(test_fragments) == "list")
  # we only expect one class
  testthat::expect_true(length(unique(test_fragments_classes)) == 1)
  # everything in the list is the class we expect
  testthat::expect_true(unique(test_fragments_classes) == "fragments_repeats")
  # no missing unique ids
  testthat::expect_false(any(is.na(test_fragments_ids)))

})

# repeat_table_to_fragments -------------------------------------------------

testthat::test_that("repeat_table_to_fragments",{
  repeat_table <- instability::example_data_repeat_table

  test_fragments <- repeat_table_to_repeats(
    repeat_table,
    repeat_col = "repeats",
    frequency_col = "height",
    unique_id = "unique_id"
  )

  test_sample <- test_fragments[[1]]

  test_fragments_classes <- vector("character", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_classes[i] <- class(test_fragments[[i]])[1]
  }

  test_fragments_ids <- vector("character", length(test_fragments))
  for (i in seq_along(test_fragments)) {
    test_fragments_ids[i] <- test_fragments[[i]]$unique_id
  }

  # it's a list
  testthat::expect_true(class(test_fragments) == "list")
  # we only expect one class
  testthat::expect_true(length(unique(test_fragments_classes)) == 1)
  # everything in the list is the class we expect
  testthat::expect_true(unique(test_fragments_classes) == "fragments_repeats")
  # no missing unique ids
  testthat::expect_false(any(is.na(test_fragments_ids)))


})


# add metadata -------------------------------------------------

testthat::test_that("add_metadata", {

  gm_raw <- instability::example_data
  metadata <- instability::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(gm_raw,
                                            data_format = "genemapper5",
                                            # peak_size_col = "size",
                                            # peak_height_col = "signal",
                                            dye_channel = "B")

  test_metadata <- add_metadata(
    fragments_list = test_fragments,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    plate_id = "plate_id",
    group_id = "cell_line",
    metrics_baseline_control = "metrics_baseline_control_TF",
    size_standard = "repeat_positive_control_TF",
    size_standard_repeat_length = "repeat_positive_control_length")


  # plate id assigned
  test_fragments_plate_id <- vector("character", length(test_metadata))
  for (i in seq_along(test_metadata)) {
    test_fragments_plate_id[i] <- test_metadata[[i]]$plate_id
  }

  testthat::expect_true(all(test_fragments_plate_id == 20230414))

  # sample group id assigned
  test_fragments_group_id <- vector("character", length(test_metadata))
  for (i in seq_along(test_metadata)) {
    test_fragments_group_id[i] <- test_metadata[[i]]$group_id
  }

  testthat::expect_true(all(!is.na(test_fragments_group_id)))

  #index samples assigned
  test_fragments_index <- vector("character", length(test_metadata))
  for (i in seq_along(test_metadata)) {
    test_fragments_index[i] <- test_metadata[[i]]$metrics_baseline_control
  }

  index_samples <- which(metadata$metrics_baseline_control_TF == TRUE)

  testthat::expect_true(all(as.logical(test_fragments_index[index_samples])))
  testthat::expect_true(unique(test_fragments_index[which(!seq_along(test_fragments_index) %in% index_samples)]) == "NA")

  #  sizing controls assigned
  test_fragments_repeat_sizing <- vector("character", length(test_metadata))
  for (i in seq_along(test_metadata)) {
    test_fragments_repeat_sizing[i] <- test_metadata[[i]]$size_standard
  }

  repeat_sizing_samples <- which(metadata$size_standard == TRUE)

  testthat::expect_true(all(as.logical(test_fragments_repeat_sizing[repeat_sizing_samples])))
  testthat::expect_true(all(is.na(test_fragments_repeat_sizing[which(!seq_along(test_fragments_repeat_sizing) %in% repeat_sizing_samples)])))

  # size values assigned
  test_fragments_repeat_sizing_value <- vector("character", length(test_metadata))
  for (i in seq_along(test_metadata)) {
    test_fragments_repeat_sizing_value[i] <- test_metadata[[i]]$size_standard_repeat_length
  }

  testthat::expect_true(test_fragments_repeat_sizing_value[repeat_sizing_samples[1]] == 113 & test_fragments_repeat_sizing_value[repeat_sizing_samples[2]] == 115 )
  testthat::expect_true(all(is.na(test_fragments_repeat_sizing_value[which(!seq_along(test_fragments_repeat_sizing_value) %in% repeat_sizing_samples)])))


})


# find alleles ---------------------------------

testthat::test_that("find_alleles", {
  gm_raw <- instability::example_data
  metadata <- instability::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(gm_raw,
                                            data_format = "genemapper5",
                                            dye_channel = "B")

  test_alleles <- find_alleles(
    fragments_list = test_fragments,
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1)

  allele_1_size <- vector("numeric", length(test_alleles))
  for (i in seq_along(test_alleles)) {
    allele_1_size[i] <- test_alleles[[i]]$allele_1_size
  }

  testthat::expect_true(all(!is.na(allele_1_size)))


})

# call repeats ---------------------------------------

testthat::test_that("call_repeats", {
  gm_raw <- instability::example_data
  metadata <- instability::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(gm_raw,
                                            data_format = "genemapper5",
                                            dye_channel = "B")

  test_alleles <- find_alleles(
    fragments_list = test_fragments,
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1)


  test_repeats <- call_repeats(
    fragments_list = test_alleles,
    repeat_algorithm = "simple",
    assay_size_without_repeat = 87,
    repeat_size = 3,
    repeat_length_correction = "none"
  )


  test_repeats_class <- vector("numeric", length(test_repeats))
  for (i in seq_along(test_repeats)) {
    test_repeats_class[i] <- class(test_repeats[[i]])[1]
  }


  testthat::expect_true(all(unique(test_repeats_class) == "fragments_repeats"))


  # nearest peak algo

  test_repeats_np <- call_repeats(
    fragments_list = test_alleles,
    repeat_algorithm = "nearest_peak",
    assay_size_without_repeat = 87,
    repeat_size = 3,
    repeat_length_correction = "none"
  )


  test_repeats_np_dif <- vector("list", length(test_repeats_np))
  for (i in seq_along(test_repeats_np)) {

    repeat_sizes <- test_repeats_np[[i]]$repeat_table_df$repeats

    lag <- vector("numeric", length(repeat_sizes))
    for (j in 2:length(repeat_sizes)) {
      lag[j] <- repeat_sizes[j] - repeat_sizes[j-1]
    }

    test_repeats_np_dif[[i]] <- lag
  }

  all_integers <- sapply(test_repeats_np_dif, function(x) all(round(x, 10) %in% 0:200))

  testthat::expect_true(all(all_integers))

  # correct repeat length

  test_alleles_metadata <- add_metadata(test_alleles, metadata,
                                        group_id = "cell_line",
                                        unique_id = "unique_id")

  test_repeats_corrected <- call_repeats(
    fragments_list = test_alleles_metadata,
    repeat_algorithm = "simple",
    assay_size_without_repeat = 87,
    repeat_size = 3,
    repeat_length_correction = "from_metadata"
  )

  mod_coefficients <- test_repeats_corrected[[1]]$.__enclos_env__$private$correction_mod$coefficients

  testthat::expect_true(round(mod_coefficients[2],3) == 0.337)

})


testthat::test_that("call_repeats with correction from genemapper alleles", {

  gm_raw <- instability::example_data_genemapper_alleles
  metadata <- instability::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(gm_raw,
                                            data_format = "genemapper5",
                                            dye_channel = "B")

  test_alleles <- find_alleles(
    fragments_list = test_fragments,
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1)

  test_repeats_corrected <- call_repeats(
    fragments_list = test_alleles,
    repeat_algorithm = "simple",
    assay_size_without_repeat = 87,
    repeat_size = 3,
    repeat_length_correction = "from_genemapper"
  )

  mod_coefficients <- test_repeats_corrected[[1]]$.__enclos_env__$private$correction_mod$coefficients

  testthat::expect_true(round(mod_coefficients[2],3) == 0.337)

})


# metrics ---------------------------------------

testthat::test_that("calculate metrics", {
  gm_raw <- instability::example_data
  metadata <- instability::metadata
  # Save raw data as a fragment class
  suppressWarnings({

    test_fragments <- peak_table_to_fragments(gm_raw,
                                              data_format = "genemapper5",
                                              dye_channel = "B",
                                              min_size_bp = 400)
  })


  suppressWarnings({
    test_metadata <- add_metadata(
      fragments_list = test_fragments,
      metadata_data.frame = metadata,
      unique_id = "unique_id",
      plate_id = "plate_id",
      group_id = "cell_line",
      metrics_baseline_control = "metrics_baseline_control_TF",
      size_standard = "repeat_positive_control_TF",
      size_standard_repeat_length = "repeat_positive_control_length")

    test_alleles <- find_alleles(
      fragments_list = test_metadata,
      number_of_peaks_to_return = 1,
      peak_region_size_gap_threshold = 6,
      peak_region_height_threshold_multiplier = 1)


    test_repeats <- call_repeats(
      fragments_list = test_alleles,
      repeat_algorithm = "simple",
      assay_size_without_repeat = 87,
      repeat_size = 3,
      repeat_length_correction = "none"
    )

    # grouped

    test_metrics_grouped <- calculate_instability_metrics(
      fragments_list = test_repeats,
      grouped = TRUE,
      peak_threshold = 0.05,
      # note the lower lim should be a negative value
      window_around_main_peak = c(-40, 40),
      percentile_range = c(0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99),
      repeat_range = c( 1,2,3,4,seq(6,20,2)),
      index_override_dataframe = NULL)
  })



  testthat::expect_true(round(mean(test_metrics_grouped$expansion_index, na.rm = TRUE), 3) == 6.673)
  testthat::expect_true(round(mean(test_metrics_grouped$average_repeat_gain, na.rm = TRUE), 3) == 4.238)
  testthat::expect_true(round(mean(test_metrics_grouped$skewness, na.rm = TRUE), 5) == -0.01007)

  # ungrouped


  suppressWarnings({
    test_repeats <- call_repeats(
      fragments_list = test_alleles,
      repeat_algorithm = "simple",
      assay_size_without_repeat = 87,
      repeat_size = 3,
      repeat_length_correction = "none"
    )

    test_metrics_ungrouped <- calculate_instability_metrics(
      fragments_list = test_repeats,
      grouped = FALSE,
      peak_threshold = 0.05,
      # note the lower lim should be a negative value
      window_around_main_peak = c(-40, 40),
      percentile_range = c(0.01, 0.05, seq(0.1, 0.9, 0.1), 0.95, 0.99),
      repeat_range = c( 1,2,3,4,seq(6,20,2)),
      index_override_dataframe = NULL)
    })



  testthat::expect_true(round(mean(test_metrics_ungrouped$expansion_index, na.rm = TRUE), 3) == 4.899)
  testthat::expect_true(all(is.na(test_metrics_ungrouped$average_repeat_gain)))
  testthat::expect_true(round(mean(test_metrics_ungrouped$skewness, na.rm = TRUE), 5) == -0.01007)



})


# remove fragments --------------------------------------------------------

testthat::test_that("remove fragments",{

  gm_raw <- instability::example_data
  metadata <- instability::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(
    gm_raw,
    data_format = "genemapper5",
    dye_channel = "B")

  all_fragment_names <- names(test_fragments)
  samples_to_remove <- all_fragment_names[c(1,length(all_fragment_names))]

  samples_removed <- remove_fragments(test_fragments, samples_to_remove)

  testthat::expect_true(all(!names(samples_removed) %in% samples_to_remove))

})


# plot fragments ----------------------------------------------------------

testthat::test_that("plot fragments",{

  # gm_raw <- instability::example_data
  # metadata <- instability::metadata
  # # Save raw data as a fragment class
  #
  # test_fragments <- peak_table_to_fragments(gm_raw,
  #                                           data_format = "genemapper5",
  #                                           dye_channel = "B")
  #
  # test_alleles <- find_alleles(
  #   fragments_list = test_fragments,
  #   number_of_peaks_to_return = 2,
  #   peak_region_size_gap_threshold = 6,
  #   peak_region_height_threshold_multiplier = 1)
  #
  # gg <- plot_fragments(test_alleles,
  #                names(test_alleles)[1:9])
  #

  ## come up with tests


})
