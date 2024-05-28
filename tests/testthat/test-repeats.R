

testthat::test_that("call_repeats", {
  gm_raw <- instability::example_data
  metadata <- instability::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(gm_raw,
                                            data_format = "genemapper5",
                                            dye_channel = "B"
  )

  test_alleles <- find_alleles(
    fragments_list = test_fragments,
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1
  )

  suppressMessages(
    test_repeats <- call_repeats(
      fragments_list = test_alleles,
      repeat_calling_algorithm = "simple",
      assay_size_without_repeat = 87,
      repeat_size = 3,
      repeat_length_correction = "none"
    )
  )


  test_repeats_class <- vector("numeric", length(test_repeats))
  for (i in seq_along(test_repeats)) {
    test_repeats_class[i] <- class(test_repeats[[i]])[1]
  }


  testthat::expect_true(all(unique(test_repeats_class) == "fragments_repeats"))


  # nearest peak algo

  suppressMessages(
    test_repeats_np <- call_repeats(
      fragments_list = test_alleles,
      force_whole_repeat_units = TRUE,
      assay_size_without_repeat = 87,
      repeat_size = 3,
      repeat_length_correction = "none"
    )
  )


  test_repeats_np_dif <- vector("list", length(test_repeats_np))
  for (i in seq_along(test_repeats_np)) {
    repeat_sizes <- test_repeats_np[[i]]$repeat_table_df$repeats

    lag <- vector("numeric", length(repeat_sizes))
    for (j in 2:length(repeat_sizes)) {
      lag[j] <- repeat_sizes[j] - repeat_sizes[j - 1]
    }

    test_repeats_np_dif[[i]] <- lag
  }

  all_integers <- sapply(test_repeats_np_dif, function(x) all(round(x, 10) %in% 0:200))

  testthat::expect_true(all(all_integers))

  # correct repeat length

  test_alleles_metadata <- add_metadata(test_alleles, metadata,
                                        group_id = "cell_line",
                                        unique_id = "unique_id",
                                        size_standard = "repeat_positive_control_TF",
                                        size_standard_repeat_length = "repeat_positive_control_length"
  )

  suppressMessages(
    suppressMessages(
      test_repeats_corrected <- call_repeats(
        fragments_list = test_alleles_metadata,
        repeat_calling_algorithm = "simple",
        assay_size_without_repeat = 87,
        repeat_size = 3,
        repeat_length_correction = "from_metadata"
      )
    )
  )

  mod_coefficients <- test_repeats_corrected[[1]]$.__enclos_env__$private$correction_mod$coefficients

  testthat::expect_true(round(mod_coefficients[2], 3) == 0.337)
})


testthat::test_that("call_repeats with correction from genemapper alleles", {
  gm_raw <- instability::example_data_genemapper_alleles
  metadata <- instability::metadata
  # Save raw data as a fragment class

  test_fragments <- peak_table_to_fragments(gm_raw,
                                            data_format = "genemapper5",
                                            dye_channel = "B"
  )

  test_alleles <- find_alleles(
    fragments_list = test_fragments,
    number_of_peaks_to_return = 2,
    peak_region_size_gap_threshold = 6,
    peak_region_height_threshold_multiplier = 1
  )

  suppressMessages(
    test_repeats_corrected <- call_repeats(
      fragments_list = test_alleles,
      repeat_calling_algorithm = "simple",
      assay_size_without_repeat = 87,
      repeat_size = 3,
      repeat_length_correction = "from_genemapper"
    )
  )

  mod_coefficients <- test_repeats_corrected[[1]]$.__enclos_env__$private$correction_mod$coefficients

  testthat::expect_true(round(mod_coefficients[2], 3) == 0.337)
})


test_that("fft", {

  suppressWarnings(
    test_ladders <- find_ladders(cell_line_fsa_list[1],
                                 ladder_channel = "DATA.105",
                                 signal_channel = "DATA.1",
                                 ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                                 max_combinations = 2500000,
                                 ladder_selection_window = 5,
                                 show_progress_bar = FALSE
    )
  )

  # plot_ladders(test_ladders[1:9], n_facet_col = 3,
  #              xlim = c(1000, 4800),
  #              ylim = c(0, 15000))


  # # Start a PDF device
  # pdf(file = "C:/Users/zlm2/Downloads/ladder.pdf", width = 12, height = 6) # Set width and height as desired
  #
  # # Loop through the list of plots
  # for (i in seq_along(test_ladders)) {
  #   test_ladders[[i]]$plot_ladder(xlim = c(1400, 4500))
  # }
  #
  # # Close the PDF device
  # dev.off()


  peak_list <- find_fragments(test_ladders,
                              minimum_peak_signal = 20,
                              min_bp_size = 300
  )


  fragment_alleles <- find_alleles(
    fragments_list = peak_list,
    number_of_peaks_to_return = 1
  )


  suppressMessages(
    suppressWarnings(
      test_repeats_size_fft <- call_repeats(
        fragments_list = fragment_alleles,
        repeat_calling_algorithm = "fft"
      )
    )
  )

  suppressMessages(
    suppressWarnings(
      test_repeats_size_simple <- call_repeats(
        fragments_list = fragment_alleles,
        repeat_calling_algorithm = "simple"
      )
    )
  )


  testthat::expect_true(
    identical( test_repeats_size_fft[[1]]$repeat_table_df[which(test_repeats_size_fft[[1]]$repeat_table_df$repeats > 117 & test_repeats_size_fft[[1]]$repeat_table_df$repeats <140),"repeats"],
               test_repeats_size_simple[[1]]$repeat_table_df[which(test_repeats_size_simple[[1]]$repeat_table_df$repeats > 117 & test_repeats_size_simple[[1]]$repeat_table_df$repeats <140),"repeats"]
    )
  )

})



test_that("repeat period", {

  suppressWarnings(
    test_ladders <- find_ladders(cell_line_fsa_list[1],
                                 ladder_channel = "DATA.105",
                                 signal_channel = "DATA.1",
                                 ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                                 max_combinations = 2500000,
                                 ladder_selection_window = 5,
                                 show_progress_bar = FALSE
    )
  )

  # plot_ladders(test_ladders[1:9], n_facet_col = 3,
  #              xlim = c(1000, 4800),
  #              ylim = c(0, 15000))


  # # Start a PDF device
  # pdf(file = "C:/Users/zlm2/Downloads/ladder.pdf", width = 12, height = 6) # Set width and height as desired
  #
  # # Loop through the list of plots
  # for (i in seq_along(test_ladders)) {
  #   test_ladders[[i]]$plot_ladder(xlim = c(1400, 4500))
  # }
  #
  # # Close the PDF device
  # dev.off()


  peak_list <- find_fragments(test_ladders,
                              minimum_peak_signal = 20,
                              min_bp_size = 300
  )


  fragment_alleles <- find_alleles(
    fragments_list = peak_list,
    number_of_peaks_to_return = 1
  )

  suppressMessages(
    suppressWarnings(
      test_repeats_size_period <- call_repeats(
        fragments_list = fragment_alleles,
        repeat_calling_algorithm = "size_period"
      )
    )
  )

  suppressMessages(
    suppressWarnings(
      test_repeats_size_simple <- call_repeats(
        fragments_list = fragment_alleles,
        repeat_calling_algorithm = "simple"
      )
    )
  )


  testthat::expect_true(
    identical( test_repeats_size_period[[1]]$repeat_table_df[which(test_repeats_size_period[[1]]$repeat_table_df$repeats > 117 & test_repeats_size_period[[1]]$repeat_table_df$repeats <140),"repeats"],
               test_repeats_size_simple[[1]]$repeat_table_df[which(test_repeats_size_simple[[1]]$repeat_table_df$repeats > 117 & test_repeats_size_simple[[1]]$repeat_table_df$repeats <140),"repeats"]
               )
  )


})
