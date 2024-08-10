

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

  test_alleles_metadata <- add_metadata(test_alleles, metadata
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



testthat::test_that("full pipline repeat size algo", {
  suppressWarnings(
    test_ladders <- find_ladders(cell_line_fsa_list,
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

  fragment_metadata <- add_metadata(
    fragments_list = peak_list,
    metadata_data.frame = metadata
  )

  fragment_alleles <- find_alleles(
    fragments_list = fragment_metadata,
    number_of_peaks_to_return = 1
  )

  suppressMessages(
    suppressWarnings(
      test_repeats <- call_repeats(
        fragments_list = fragment_alleles,
        repeat_length_correction = "from_metadata",
        repeat_calling_algorithm = "size_period"
      )
    )
  )

  # plot_traces(test_repeats[1:9], n_facet_col = 3,
  #             xlim = c(100, 200),
  #             ylim = c(0,2000))


  # plot_fragments(test_repeats[1:4])
  # plot_repeat_correction_model(test_repeats)

  suppressMessages(
    suppressWarnings(
      test_metrics_grouped <- calculate_instability_metrics(
        fragments_list = test_repeats,
        grouped = TRUE,
        peak_threshold = 0.05,
        window_around_index_peak = c(-40, 40)
      )
    )
  )


  # Left join
  plot_data <- merge(test_metrics_grouped, metadata, by = "unique_id", all.x = TRUE)

  # Filter
  plot_data <- plot_data[plot_data$day > 0 & plot_data$modal_peak_height > 500, ]

  # Group by
  plot_data <- split(plot_data, plot_data$group_id)

  # Mutate
  for (i in seq_along(plot_data)) {
    plot_data[[i]]$rel_gain <- plot_data[[i]]$average_repeat_gain / median(plot_data[[i]]$average_repeat_gain[which(plot_data[[i]]$treatment == 0)])
  }

  plot_data <- do.call(rbind, plot_data)

  # Revise genotype levels

  plot_data$genotype <- factor(plot_data$genotype, levels = c("non-edited", "edited"))



  # ggplot2::ggplot(plot_data,
  #                 ggplot2::aes(as.factor(treatment), rel_gain,
  #            colour = as.factor(treatment))) +
  #   ggplot2::geom_boxplot(outlier.shape = NA) +
  #   ggplot2::geom_jitter() +
  #   ggplot2::facet_wrap(ggplot2::vars(genotype)) +
  #   ggplot2::labs(y = "Average repeat gain\n(relative to DMSO)",
  #        x = "Branaplam (nM)") +
  #   ggplot2::theme(legend.position = "none")


  medians <- aggregate(rel_gain ~ treatment + genotype, plot_data, median, na.rm = TRUE)

  testthat::expect_true(all(round(medians$rel_gain, 5) == c(1.00000, 0.85562, 0.70208, 0.56223, 1.00000, 1.18592, 1.10739, 1.00075)))
})

