
# init --------------------------------------------------------------------

testthat::test_that("fragments class initialization", {

 test_fragments <- fragments$new(
    unique_id = "test"
  )

  testthat::expect_equal(
    class(test_fragments)[1], "fragments"
  )
})


# add_metadata ------------------------------------------------------------


testthat::test_that("add_metadata function", {

  gm_raw <- instability::example_data
  metadata <- instability::metadata

  test_df <- gm_raw[which(gm_raw$Sample.File.Name == metadata$unique_id[1]), ]
  metadata <- metadata[which(metadata$unique_id == metadata$unique_id[1]), ]

  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(test_df,
                                       data_format = "genemapper5",
                                       dye_channel = "B")

  test_metadata <- test_fragments[[1]]$add_metadata(
    metadata,
    unique_id = "unique_id",
    plate_id = "plate_id",
    sample_group_id = "cell_line",
    repeat_positive_control_TF = "repeat_positive_control_TF",
    repeat_positive_control_length = "repeat_positive_control_length",
    metrics_baseline_control = "metrics_baseline_control_TF"
  )

  testthat::expect_true(!is.na(test_metadata$unique_id))
  testthat::expect_true(!is.na(test_metadata$plate_id))
  testthat::expect_true(!is.na(test_metadata$group_id))

})


# find_main_peaks ---------------------------------------------------------

testthat::test_that("main peaks bp_fragments",{

  gm_raw <- instability::example_data
  test_sample <- "20230413_A01.fsa"

  test_df <- gm_raw[which(gm_raw$Sample.File.Name == test_sample), ]

  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(test_df,
                                       data_format = "genemapper5",
                                       # peak_size_col = "size",
                                       # peak_height_col = "signal",
                                       dye_channel = "B")

  test_main_peaks <- test_fragments[[1]]$find_main_peaks(
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1
  )

  testthat::expect_true(!is.na(test_main_peaks$allele_2_size))
  testthat::expect_true(!is.na(test_main_peaks$allele_1_size))
  testthat::expect_true(!is.na(test_main_peaks$allele_2_height))
  testthat::expect_true(!is.na(test_main_peaks$allele_1_height))
  testthat::expect_true(test_main_peaks$allele_2_size == 131.94,
                        test_main_peaks$allele_1_size == 465.82)
})


# add_repeats ---------------------------------------------------------

testthat::test_that("add_repeats",{

  gm_raw <- instability::example_data
  test_sample <-  "20230413_A01.fsa"

  test_df <- gm_raw[which(gm_raw$Sample.File.Name == test_sample), ]

  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(test_df,
                                            data_format = "genemapper5",
                                            # peak_size_col = "size",
                                            # peak_height_col = "signal",
                                            dye_channel = "B")

  test_main_peaks <- test_fragments[[1]]$find_main_peaks(
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1
  )

  test_repeats_fragments_simple <- test_main_peaks$add_repeats(
    repeat_algorithm = "simple",
    assay_size_without_repeat = 87,
    repeat_size = 3,
    correct_repeat_length = FALSE)


  test_repeats_fragments_np <- test_main_peaks$add_repeats(
    repeat_algorithm = "nearest_peak",
    assay_size_without_repeat = 87,
    repeat_size = 3,
    correct_repeat_length = FALSE)

  testthat::expect_true(nrow(test_repeats_fragments_simple$repeat_data) == 40)
  testthat::expect_true(round(test_repeats_fragments_simple$repeat_data[40,4],2) == 142.43)
})

# find_main_peaks for repeats class ------------------------------------------



testthat::test_that("add_repeats",{

  gm_raw <- instability::example_data
  test_sample <- "20230413_A01.fsa"

  test_df <- gm_raw[which(gm_raw$Sample.File.Name == test_sample), ]

  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(test_df,
                                            data_format = "genemapper5",
                                            # peak_size_col = "size",
                                            # peak_height_col = "signal",
                                            dye_channel = "B")

  test_main_peaks <- test_fragments[[1]]$find_main_peaks(
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1
  )

  test_repeats_fragments_simple <- test_main_peaks$add_repeats(
    repeat_algorithm = "simple",
    assay_size_without_repeat = 87,
    repeat_size = 3,
    correct_repeat_length = FALSE)

  test_repeats_main_peaks <- test_repeats_fragments_simple$find_main_peaks(
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 2, #this is the only thing different with the find main peaks above
    peak_region_height_threshold_multiplier = 1
  )

  testthat::expect_true(!is.na(test_repeats_main_peaks$allele_2_repeat))
  testthat::expect_true(!is.na(test_repeats_main_peaks$allele_1_repeat))
  testthat::expect_true(!is.na(test_repeats_main_peaks$allele_2_height))
  testthat::expect_true(!is.na(test_repeats_main_peaks$allele_1_height))
  testthat::expect_true(round(test_repeats_main_peaks$allele_2_repeat, 2) == 14.98,
                        round(test_repeats_main_peaks$allele_1_repeat, 2) == 126.27)
})

# instability metrics -------------------------------------------------

testthat::test_that("instability metrics",{


  gm_raw <- instability::example_data
  test_sample <- "20230413_A01.fsa"

  test_df <- gm_raw[which(gm_raw$Sample.File.Name == test_sample), ]

  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(test_df,
                                            data_format = "genemapper5",
                                            # peak_size_col = "size",
                                            # peak_height_col = "signal",
                                            dye_channel = "B")

  test_main_peaks <- test_fragments[[1]]$find_main_peaks(
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1
  )

  test_repeats_fragments_simple <- test_main_peaks$add_repeats(
    repeat_algorithm = "simple",
    assay_size_without_repeat = 87,
    repeat_size = 3,
    correct_repeat_length = FALSE)

  test_metrics <- test_repeats_fragments_simple$instability_metrics(
    peak_threshold = 0.05,
    window_around_main_peak = c(-40, 40), # note the lower lim should be a negative value
    percentile_range = c(
      0.01, 0.05, 0.1, 0.2, 0.3,
      0.4, 0.5, 0.6, 0.7, 0.8,
      0.9, 0.95, 0.99)
  )

  testthat::expect_true(round(test_metrics$weighted_mean_repeat, 2) == 126.79)
  testthat::expect_true(round(test_metrics$instabity_index_jml, 2) == 0.51)
  testthat::expect_true(round(test_metrics$abs_index , 2) == 4.25)
  testthat::expect_true(round(test_metrics$expansion_index, 2) == 4.38)
  testthat::expect_true(round(test_metrics$pps_index, 2) == 6.5)
  testthat::expect_true(round(test_metrics$contration_index, 2) == -3.54)
  testthat::expect_true(round(test_metrics$expansion_percentile_0.9, 2) == 9.53)

})

