testthat::test_that("find_fragments", {
  # simple option:

  file_list <- instability::cell_line_fsa_list
  suppressWarnings(
    test_ladders <- find_ladders(file_list[which(names(file_list) == "20230413_B03.fsa")],
      ladder_channel = "DATA.105",
      signal_channel = "DATA.1",
      ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )

  suppressWarnings(
    peak_list <- find_fragments(test_ladders,
      minimum_peak_signal = 20,
      min_bp_size = 100
    )
  )



  # plot_traces(peak_list[which(names(peak_list) =="20230413_B03.fsa")],
  #             xlim = c(450,650),
  #             ylim = c(0, 1000))


  extracted_fragments <- extract_fragments(peak_list[which(names(peak_list) == "20230413_B03.fsa")])


  tall_peaks <- extracted_fragments[which(extracted_fragments$size > 500 & extracted_fragments$height > 700), ]


  testthat::expect_true(all(round(tall_peaks$size, 4) == c(547.5879, 550.0142, 552.4404, 554.8666, 557.2929)))
})



# metadata transfer



testthat::test_that("metadata transfer", {
  file_list <- instability::cell_line_fsa_list


  suppressWarnings(
    test_ladders <- find_ladders(file_list[1],
      ladder_channel = "DATA.105",
      signal_channel = "DATA.1",
      ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      max_combinations = 2500000,
      ladder_selection_window = 8,
      show_progress_bar = FALSE
    )
  )
  suppressWarnings(
    metadata_added <- add_metadata(test_ladders, metadata,
      unique_id = "unique_id",
      plate_id = "plate_id",
      group_id = "cell_line",
      metrics_baseline_control = "metrics_baseline_control_TF",
      size_standard = "repeat_positive_control_TF",
      size_standard_repeat_length = "repeat_positive_control_length"
    )
  )

  peak_list <- find_fragments(metadata_added,
    minimum_peak_signal = 20,
    min_bp_size = 300
  )

  testthat::expect_true(metadata_added[[1]]$unique_id == peak_list[[1]]$unique_id)
  testthat::expect_true(metadata_added[[1]]$plate_id == peak_list[[1]]$plate_id)
  testthat::expect_true(metadata_added[[1]]$group_id == peak_list[[1]]$group_id)
  testthat::expect_true(metadata_added[[1]]$metrics_baseline_control == peak_list[[1]]$metrics_baseline_control)
  testthat::expect_false(metadata_added[[1]]$size_standard & peak_list[[1]]$size_standard)
  testthat::expect_true(is.na(metadata_added[[1]]$size_standard_repeat_length) & is.na(peak_list[[1]]$size_standard_repeat_length))
})
