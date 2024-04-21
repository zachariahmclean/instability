

testthat::test_that("find ladder peaks",{

  file_list <- instability::cell_line_fsa_list




  test_processed <- process_ladder_signal(file_list[[1]]$Data$DATA.105,
                                          scans = 0:(length(file_list[[1]]$Data$DATA.105) -1),
                                          spike_location = 1000,
                                          smoothing_window = 10)


  ladder_sizes=c(50, 75, 100, 139, 150, 160, 200, 300, 350, 400, 450, 490, 500)


  test_ladder_peaks <- find_ladder_peaks(
    test_processed,
    length(ladder_sizes)
    )

  testthat::expect_true(length(test_ladder_peaks) >= length(ladder_sizes))


  test_ladder_peaks_40 <- find_ladder_peaks(
    test_processed,
    40
  )

  testthat::expect_true(length(test_ladder_peaks_40) >= 40)

})


test_that("iterative ladder", {



  scans_162 <- c(1548 ,
                 1624,
                 1771,
                 1913,
                 2124,
                 2142,
                 2200,
                 2259,
                 2503,
                 2802,
                 3131,
                 3376,
                 3439,
                 3758,
                 4050,
                 4287,
                 4336)



  ladder_sizes=c(50, 75, 100, 139, 150, 160, 200, 300, 350, 400, 450, 490, 500)


  iteration_result <- ladder_iteration(ladder_sizes, scans_162)

  expect_true(round(mean(iteration_result$scan),3) == 2877.923)



})



test_that("fit ladder", {


  file_list <- instability::cell_line_fsa_list

  test_ladder_signal <- file_list[[1]]$Data$DATA.105
  test_scans <- 0:(length(file_list[[1]]$Data$DATA.105) -1)



  test_fit <- fit_ladder(
    ladder = test_ladder_signal,
    scans = test_scans,
    spike_location = NULL,
    ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
    hq_ladder = FALSE,
    smoothing_window = 5,
    max_combinations = 2500000,
    ladder_selection_window = 5
  )


  mod <- lm(scan ~ size, data = test_fit)

  testthat::expect_true(round(mod$coefficients[[1]],3) == 1347.69)
  testthat::expect_true(round(mod$coefficients[[2]],5) == 6.25118)

})



test_that("local southern", {


  file_list <- instability::cell_line_fsa_list

  test_ladder_signal <- file_list[[1]]$Data$DATA.105
  test_scans <- 0:(length(file_list[[1]]$Data$DATA.105) -1)



  test_fit <- fit_ladder(
    ladder = test_ladder_signal,
    scans = test_scans,
    ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
    hq_ladder = FALSE,
    spike_location = NULL,
    smoothing_window = 5,
    max_combinations = 2500000,
    ladder_selection_window = 5
  )

  #ladder all as one lm


  mod <- lm(size ~ scan, data = test_fit)
  single_rsqd <- summary(mod)$r.squared




  mod_fit <- local_southern_fit(test_fit$scan, test_fit$size)

  multi_rsqd <- sapply(mod_fit, function(x) summary(x$mod)$r.squared)

  rsq_diff <- sum(multi_rsqd - single_rsqd)


  predicted_sizes <- local_southern_predict(mod_fit, test_scans)

  #plot(predicted_sizes, file_list[[1]]$Data$DATA.1)


  testthat::expect_true(round(rsq_diff,5) == -0.00343)


})



test_that("find ladders", {

  file_list <- instability::cell_line_fsa_list

  suppressWarnings(
  test_ladders <- find_ladders(file_list[which(names(file_list) == "20230413_B03.fsa")],
               ladder_channel = "DATA.105",
               signal_channel = "DATA.1",
               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
               hq_ladder = FALSE,
               max_combinations = 2500000,
               ladder_selection_window = 8)
  )

  test_ladders_fixed <- fix_ladders_auto(test_ladders, "20230413_B03.fsa")

  testthat::expect_true(all(test_ladders$`20230413_B03.fsa`$ladder_df$scan == c(1555, 1568, 1633, 1926, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)))
  testthat::expect_true(all(test_ladders_fixed$`20230413_B03.fsa`$ladder_df$scan == c(1633, 1926, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)))



})





test_that("fix ladders", {

  file_list <- instability::cell_line_fsa_list

  suppressWarnings(
    test_ladders <- find_ladders(file_list[which(names(file_list) == "20230413_B03.fsa")],
                                 ladder_channel = "DATA.105",
                                 signal_channel = "DATA.1",
                                 ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                                 hq_ladder = FALSE,
                                 max_combinations = 2500000,
                                 ladder_selection_window = 8)
  )

  test_ladders_fixed <- fix_ladders_auto(test_ladders, "20230413_B03.fsa")



  testthat::expect_true(all(test_ladders$`20230413_B03.fsa`$ladder_df$scan == c(1555, 1568, 1633, 1926, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)))
  testthat::expect_true(all(test_ladders_fixed$`20230413_B03.fsa`$ladder_df$scan == c(1633, 1926, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)))



})

test_that("fix ladders manual", {

  example_list <- list(
    "20230413_A01.fsa" = data.frame(
      size = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      scan = c(1593, 1671, 1825, 1971, 2208, 2269, 2329, 2581, 2888, 3228, 3479, 3543, 3872, 4170, 4412, 4460)
    )
  )

  test_ladders <- find_ladders(instability::cell_line_fsa_list[1],
                               ladder_channel = "DATA.105",
                               signal_channel = "DATA.1",
                               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                               hq_ladder = FALSE,
                               max_combinations = 2500000,
                               ladder_selection_window = 8)

  test_ladders_fixed_manual <- fix_ladders_manual(
    test_ladders,
    example_list
  )

  expect_true(nrow(test_ladders_fixed_manual[[1]]$ladder_df) == 16)

})





testthat::test_that("fix_ladders_auto", {



  # simple option:
  file_list <- instability::cell_line_fsa_list

  suppressWarnings(
    test_ladders <- find_ladders(file_list[which(names(file_list) == "20230413_B03.fsa")],
                                 ladder_channel = "DATA.105",
                                 signal_channel = "DATA.1",
                                 ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                                 hq_ladder = FALSE,
                                 max_combinations = 2500000,
                                 ladder_selection_window = 8)
  )

  fixed_ladder <- fix_ladders_auto(test_ladders, "20230413_B03.fsa")

  test_df <- fixed_ladder[[1]]$trace_bp_df


  peaks_df <- find_fragment_peaks(test_df,
                      minimum_peak_signal = 20,
                      minpeakdistance = 20,
                      smoothing_window = 21)

  # ggplot2::ggplot(test_df,
  #                 aes(size, signal)) +
  #   geom_line() +
  #   geom_point(data = peaks_df,
  #              aes(size, height),
  #              shape = 1, colour = "blue") +
  #   scale_x_continuous(limits = c(475, 625)) +
  #   scale_y_continuous(limits = c(-50,1000))


  suppressWarnings(
    peak_list <- find_fragments(fixed_ladder,
                                minimum_peak_signal = 20,
                                min_bp_size = 100)
  )




  expect_true(all(fixed_ladder[[1]]$ladder_df[c(1,2,3), "scan"] == c(1633, 1926, 2159)))

  # ggplot2::ggplot(peak_list[[1]]$trace_bp_df,
  #                 aes(size, signal)) +
  #   geom_line() +
  #   geom_point(data = peak_list[[1]]$peak_table_df,
  #              aes(size, height),
  #              shape = 1, colour = "blue") +
  #   geom_hline(yintercept = 50) +
  #   scale_x_continuous(limits = c(475, 625)) +
  #   scale_y_continuous(limits = c(-50,1000))


})

testthat::test_that("find_fragments", {



  # simple option:

  file_list <- instability::cell_line_fsa_list
  suppressWarnings(
  test_ladders <- find_ladders(file_list,
                               ladder_channel = "DATA.105",
                               signal_channel = "DATA.1",
                               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                               hq_ladder = FALSE,
                               max_combinations = 2500000,
                               ladder_selection_window = 8)
)

fixed_ladder <- fix_ladders_auto(test_ladders, "20230413_B03.fsa")

suppressWarnings(
peak_list <- find_fragments(fixed_ladder,
                            minimum_peak_signal = 20,
                            min_bp_size = 100)
)



# plot_traces(peak_list[which(names(peak_list) =="20230413_B03.fsa")],
#             xlim = c(450,650),
#             ylim = c(0, 1000))


extracted_fragments <- extract_fragments(peak_list[which(names(peak_list) =="20230413_B03.fsa")])


tall_peaks <- extracted_fragments[which(extracted_fragments$size > 500 & extracted_fragments$height > 700), ]


testthat::expect_true(all(round(tall_peaks$size, 4) == c(547.5879, 550.0142, 552.4404, 554.8666, 557.2929)))

})



#metadata transfer



testthat::test_that("metadata transfer", {

  file_list <- instability::cell_line_fsa_list


  suppressWarnings(
  test_ladders <- find_ladders(file_list,
                               ladder_channel = "DATA.105",
                               signal_channel = "DATA.1",
                               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                               hq_ladder = FALSE,
                               max_combinations = 2500000,
                               ladder_selection_window = 8)
  )
  test_ladders_fixed <- fix_ladders_auto(test_ladders, "20230413_B03.fsa")
  suppressWarnings(
  metadata_added <- add_metadata(test_ladders_fixed, metadata,
                                 unique_id = "unique_id",
                                  plate_id = "plate_id",
                                  group_id = "cell_line",
                                  metrics_baseline_control = "metrics_baseline_control_TF",
                                  size_standard = "repeat_positive_control_TF",
                                  size_standard_repeat_length = "repeat_positive_control_length")
  )

  peak_list <- find_fragments(metadata_added,
                              minimum_peak_signal = 20,
                              min_bp_size = 300)

  testthat::expect_true(metadata_added[[1]]$unique_id == peak_list[[1]]$unique_id)
  testthat::expect_true(metadata_added[[1]]$plate_id == peak_list[[1]]$plate_id)
  testthat::expect_true(metadata_added[[1]]$group_id == peak_list[[1]]$group_id)
  testthat::expect_true(metadata_added[[1]]$metrics_baseline_control == peak_list[[1]]$metrics_baseline_control)
  testthat::expect_false(metadata_added[[1]]$size_standard & peak_list[[1]]$size_standard)
  testthat::expect_true(is.na(metadata_added[[1]]$size_standard_repeat_length) & is.na(peak_list[[1]]$size_standard_repeat_length))

})



testthat::test_that("full pipline", {



  # simple option:

  # files <- list.files("data/GT_Z.McLean_2023-04-14/", full.names = TRUE)



  suppressWarnings(
  test_ladders <- find_ladders(cell_line_fsa_list,
                               ladder_channel = "DATA.105",
                               signal_channel = "DATA.1",
                               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                               hq_ladder = FALSE,
                               max_combinations = 2500000,
                               ladder_selection_window = 5)
  )
  test_ladders_fixed <- fix_ladders_auto(test_ladders, "20230413_B03.fsa")

  # plot_ladders(test_ladders_fixed[1:9], n_facet_col = 3,
  #              xlim = c(1000, 4800),
  #              ylim = c(0, 15000))


  # # Start a PDF device
  # pdf(file = "C:/Users/zlm2/Downloads/ladder.pdf", width = 12, height = 6) # Set width and height as desired
  #
  # # Loop through the list of plots
  # for (i in seq_along(test_ladders_fixed)) {
  #   test_ladders_fixed[[i]]$plot_ladder(xlim = c(1400, 4500))
  # }
  #
  # # Close the PDF device
  # dev.off()


  peak_list <- find_fragments(test_ladders_fixed,
                              minimum_peak_signal = 20,
                              min_bp_size = 300)

  fragment_metadata <- add_metadata(
    fragments_list = peak_list,
    metadata_data.frame = metadata,
    unique_id = "unique_id",
    plate_id = "plate_id",
    group_id = "cell_line",
    metrics_baseline_control = "metrics_baseline_control_TF",
    size_standard = "repeat_positive_control_TF",
    size_standard_repeat_length = "repeat_positive_control_length")

  fragment_alleles <- find_alleles(
    fragments_list = fragment_metadata,
    number_of_peaks_to_return = 1)


  suppressWarnings(

  test_repeats <- call_repeats(
    fragments_list = fragment_alleles,
    repeat_length_correction = "from_metadata"
  )
  )

  # plot_traces(test_repeats[1:9], n_facet_col = 3,
  #             xlim = c(400, 550),
  #             ylim = c(0,2000))


  # plot_fragments(test_repeats[1:4])
  # plot_repeat_correction_model(test_repeats)

  suppressWarnings(
  test_metrics_grouped <- calculate_instability_metrics(
    fragments_list = test_repeats,
    grouped = TRUE,
    peak_threshold = 0.05,
    window_around_main_peak = c(-40, 40))

  )


  # Left join
  plot_data <- merge(test_metrics_grouped, metadata, by = "unique_id", all.x = TRUE)

  # Filter
  plot_data <- plot_data[plot_data$day > 0 & plot_data$modal_peak_height > 500, ]

  # Group by
  plot_data <- split(plot_data, plot_data$cell_line)

  # Mutate
  for (i in seq_along(plot_data)) {
    plot_data[[i]]$rel_gain <- plot_data[[i]]$average_repeat_gain / median(plot_data[[i]]$average_repeat_gain[which(plot_data[[i]]$treatment == 0)])
  }

  plot_data <- do.call(rbind, plot_data)

  # Revise genotype levels

  plot_data$genotype <- factor(plot_data$genotype, levels = c("non-edited", "edited"))



  # ggplot2::ggplot(plot_data,
  #        aes(as.factor(treatment), rel_gain,
  #            colour = as.factor(treatment))) +
  #   ggplot2::geom_boxplot(outlier.shape = NA) +
  #   ggplot2::geom_jitter() +
  #   facet_wrap(ggplot2::vars(genotype)) +
  #   labs(y = "Average repeat gain\n(relative to DMSO)",
  #        x = "Branaplam (nM)") +
  #   ggplot2::theme(legend.position = "none")


  medians <- aggregate(rel_gain~treatment + genotype, plot_data, median, na.rm = TRUE)

  round(medians$rel_gain, 5)
  expect_true(all(round(medians$rel_gain, 5) == c(1.00000, 0.85527, 0.70219, 0.56016, 1.00000, 1.17530, 1.10860, 1.00459)))


})
