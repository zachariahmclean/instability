

test_that("find ladder peaks",{

  files <- list.files("data/GT_Z.McLean_2023-04-14/", full.names = TRUE)

  file_list <- read_fsa(files[1])

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


  files <- list.files("data/GT_Z.McLean_2023-04-14/", full.names = TRUE)

  file_list <- read_fsa(files[1])

  test_ladder_signal <- file_list[[1]]$Data$DATA.105
  test_scans <- 0:(length(file_list[[1]]$Data$DATA.105) -1)



  test_fit <- fit_ladder(
    ladder = test_ladder_signal,
    scans = test_scans,
    ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
    hq_ladder = FALSE
  )


  mod <- lm(scan ~ size, data = test_fit)

  testthat::expect_true(round(mod$coefficients[[1]],3) == 1347.69)
  testthat::expect_true(round(mod$coefficients[[2]],5) == 6.25118)

})



test_that("local southern", {


  files <- list.files("data/GT_Z.McLean_2023-04-14/", full.names = TRUE)

  file_list <- read_fsa(files[1])

  test_ladder_signal <- file_list[[1]]$Data$DATA.105
  test_scans <- 0:(length(file_list[[1]]$Data$DATA.105) -1)



  test_fit <- fit_ladder(
    ladder = test_ladder_signal,
    scans = test_scans,
    ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
    hq_ladder = FALSE
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

  files <- list.files("data/GT_Z.McLean_2023-04-14/", full.names = TRUE)

  file_list <- read_fsa(files)

  suppressWarnings(
  test_ladders <- find_ladders(file_list[which(names(file_list) == "20230413_B03.fsa")],
               ladder_channel = "DATA.105",
               signal_channel = "DATA.1",
               sample_id_channel = 'SpNm.1',
               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
               hq_ladder = FALSE,
               max_combinations = 2500000,
               ladder_selection_window = 8)
  )

  test_ladders_fixed <- fix_ladders(test_ladders, "20230413_B03.fsa")



  test_ladders$`20230413_B03.fsa`$plot_ladder(scan_limits = c(1400, 4500))

  test_ladders_fixed$`20230413_B03.fsa`$plot_ladder(scan_limits = c(1400, 4500))

  testthat::expect_true(all(test_ladders$`20230413_B03.fsa`$ladder_df$scan == c(1555, 1568, 1633, 1926, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)))
  testthat::expect_true(all(test_ladders_fixed$`20230413_B03.fsa`$ladder_df$scan == c(1633, 1926, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)))



})


testthat::test_that("find_fragment_peaks", {



  # simple option:

  file_list <- read_fsa("data/GT_Z.McLean_2023-04-14/20230413_B03.fsa")

  test_ladders <- find_ladders(file_list,
                               ladder_channel = "DATA.105",
                               signal_channel = "DATA.1",
                               sample_id_channel = 'SpNm.1',
                               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                               hq_ladder = FALSE,
                               max_combinations = 2500000,
                               ladder_selection_window = 8)


peak_list <- vector("list", length(test_ladders))

for (i in seq_along(peak_list)) {

  peak_list[[i]]  <- test_ladders[[i]]$call_peaks(
    smoothing_window = 5,
    minumum_peak_signal = 20,
    min_bp_size = 100
  )

}

#
# peak_list[[1]]$trace_bp_df |>
#     ggplot(aes(size, signal)) +
#     geom_line() +
#     geom_point(data = peak_list[[1]]$peak_table_df,
#                aes(y = height))+
#     ggplot2::scale_x_continuous(limits = c(420,600)) +
#   ggplot2::scale_y_continuous(limits = c(NA, 5000))


peak_list <- find_fragments(test_ladders,
                            smoothing_window = 5,
                            minumum_peak_signal = 20,
                            min_bp_size = 100)


})



#metadata transfer



testthat::test_that("metadata transfer", {

  file_list <- read_fsa("data/GT_Z.McLean_2023-04-14/20230413_B03.fsa")

  test_ladders <- find_ladders(file_list,
                               ladder_channel = "DATA.105",
                               signal_channel = "DATA.1",
                               sample_id_channel = 'SpNm.1',
                               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                               hq_ladder = FALSE,
                               max_combinations = 2500000,
                               ladder_selection_window = 8)
  test_ladders_fixed <- fix_ladders(test_ladders, "20230413_B03.fsa")
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
                              smoothing_window = 5,
                              minumum_peak_signal = 20,
                              min_bp_size = 300)

  testthat::expect_true(metadata_added[[1]]$unique_id == peak_list[[1]]$unique_id)
  testthat::expect_true(metadata_added[[1]]$plate_id == peak_list[[1]]$plate_id)
  testthat::expect_true(metadata_added[[1]]$group_id == peak_list[[1]]$group_id)
  testthat::expect_true(metadata_added[[1]]$metrics_baseline_control == peak_list[[1]]$metrics_baseline_control)
  testthat::expect_true(metadata_added[[1]]$size_standard == peak_list[[1]]$size_standard)
  testthat::expect_true(metadata_added[[1]]$size_standard_repeat_length == peak_list[[1]]$size_standard_repeat_length)

})



testthat::test_that("full pipline", {



  # simple option:

  files <- list.files("data/GT_Z.McLean_2023-04-14/", full.names = TRUE)

  file_list <- read_fsa(files)


  test_ladders <- find_ladders(file_list,
                               ladder_channel = "DATA.105",
                               signal_channel = "DATA.1",
                               sample_id_channel = 'SpNm.1',
                               ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
                               hq_ladder = FALSE,
                               max_combinations = 2500000,
                               ladder_selection_window = 8)
  test_ladders_fixed <- fix_ladders(test_ladders, "20230413_B03.fsa")


  # # Start a PDF device
  # pdf(file = "C:/Users/zlm2/Downloads/ladder.pdf", width = 12, height = 6) # Set width and height as desired
  #
  # # Loop through the list of plots
  # for (i in seq_along(test_ladders_fixed)) {
  #   test_ladders_fixed[[i]]$plot_ladder(scan_limits = c(1400, 4500))
  # }
  #
  # # Close the PDF device
  # dev.off()

  peak_list <- find_fragments(test_ladders_fixed,
                              smoothing_window = 5,
                              minumum_peak_signal = 20,
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

  test_repeats <- call_repeats(
    fragments_list = fragment_alleles,
    repeat_length_correction = "from_metadata"
  )

  test_metrics_grouped <- calculate_instability_metrics(
    fragments_list = test_repeats,
    grouped = TRUE,
    peak_threshold = 0.05,
    window_around_main_peak = c(-40, 40))

  plot_data <- test_metrics_grouped |>
    dplyr::left_join(metadata) |>
    dplyr::filter(day >0,
                  modal_peak_height > 500) |>
    dplyr::group_by(cell_line) |>
    dplyr::mutate(rel_gain = average_repeat_gain / median(average_repeat_gain[which(treatment == 0)]),
                  genotype = forcats::fct_rev(genotype))

  ggplot2::ggplot(plot_data,
         aes(as.factor(treatment), rel_gain,
             colour = as.factor(treatment))) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter() +
    facet_wrap(ggplot2::vars(genotype)) +
    labs(y = "Average repeat gain\n(relative to DMSO)",
         x = "Branaplam (nM)") +
    ggplot2::theme(legend.position = "none")


})
