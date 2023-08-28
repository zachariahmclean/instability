# percentiles

testthat::test_that("percentiles",{

  gm_raw <- read.csv("data/example_data.txt", sep = "\t")
  test_sample <- "20230413_A01.fsa"

  test_df <- gm_raw[which(gm_raw$Sample.File.Name == test_sample), ]

  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(test_df,
                                            fragment_class = "bp_fragments",
                                            data_format = "genemapper5",
                                            # peak_size_col = "size",
                                            # peak_height_col = "signal",
                                            dye_channel = "B",
                                            min_size_bp = 350)

  test_main_peaks <- test_fragments[[1]]$find_main_peaks(
    number_of_peaks_to_return = 1,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1
  )

  test_repeats_fragments_simple <- test_main_peaks$add_repeats(
    repeat_algorithm = "simple",
    assay_size_without_repeat = 87,
    repeat_size = 3,
    correct_repeat_length = FALSE)


  test_distribution_df <- test_repeats_fragments_simple$repeat_data
  test_distribution_df <- test_distribution_df[which(test_distribution_df$repeats > 113), ]


  percentiles <- find_percentiles(
    repeats = test_distribution_df$repeats,
    heights = test_distribution_df$height,
    index_peak_repeat = test_repeats_fragments_simple$allele_1_repeat,
    type = "percentile", #"percentile" or "repeat"
    range = seq(0.1,0.99,.10),
    col_preffix = "percentile"
  )

  repeat_test <- find_percentiles(
    repeats = test_distribution_df$repeats,
    heights = test_distribution_df$height,
    index_peak_repeat = test_repeats_fragments_simple$allele_1_repeat,
    type = "repeat", #"percentile" or "repeat"
    range = percentiles$percentile_0.2,
    col_preffix = "repeat"
  )

  #the values go full cirle
  testthat::expect_true(repeat_test[[1]] == 0.2)

})
