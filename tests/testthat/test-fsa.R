
#what is going on here?
## I think I was trying out ladder recursive as an option for fun
## is the actual code the one to use?


test_that("recursive ladder", {



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



  bench::mark(ladder_recursive(ladder_sizes, scans_162, choose = 5))
  bench::mark(ladder_recursive(ladder_sizes, scans_162, choose = 12))
  ladder_recursive(ladder_sizes, scans_162)



  # simulate 1200 liz

  ladder_1200 <- c(20, 30, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600, 614, 620, 640, 660, 680, 700, 714, 720, 740, 760, 780, 800, 820, 840, 850, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1040, 1060, 1080, 1100, 1120, 1160, 1200)
  random_noise_peaks <-  runif(n = 12, 80, 5000)
  simulated_1200_scans <- jitter((ladder_1200)) / 0.25
  simulated_1200_scans_w_random <- sort(c(simulated_1200_scans, random_noise_peaks))
  bench::mark( ladder_recursive(ladder_1200, simulated_1200_scans_w_random, choose = 10))
  found <- ladder_recursive(ladder_1200, simulated_1200_scans_w_random, choose = 10)
  actual <- data.frame(scan = simulated_1200_scans,
                       size = ladder_1200)

  library(ggplot2)
  ggplot(data = found,
         aes(scan, size)) +
    geom_point(colour = "blue",
               shape = 1) +
    geom_point(data = actual,
               colour = "red",
               shape = 2,
               size = 0.5) +
    geom_vline(xintercept = random_noise_peaks)






})


# ladder code development

test_that("ladder_dev", {


  #read in samples

  files <- list.files("data/fsa_files/2022-01-03_JR_traces/", full.names = TRUE)

  file_list <- read_fsa(files)
  ladder_1200 <- c(20, 30, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560, 580, 600, 614, 620, 640, 660, 680, 700, 714, 720, 740, 760, 780, 800, 820, 840, 850, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1040, 1060, 1080, 1100, 1120, 1160, 1200)


  test_processed <- process_ladder_signal(file_list[[1]]$Data$DATA.105,
                                          scans = 0:(length(file_list[[1]]$Data$DATA.105) -1),
                                          spike_location = 1000)

  test_ladder_peaks <- find_ladder_peaks(
    test_processed,
    ladder_peak_threshold = 0.8)



  fit_unique_chunks(
    ladder_1200,
    test_ladder_peaks$scan,
    10,
    rsq_threshold_ladder = 0.999,
    rsq_threshold_observed = 0.9999)

  # can include logic to find when ladder peaks have been droped.
  # do an lm from the chunk before,
  # predict scans based off ladder sizes,
  # if there is no observed peak within ~10 scans of where we expect one based off ladder, remove ladder peak




  #identify the edges of the

  called_scans_list <- vector('list', length(identified_peaks_over_threshold))
  for (i in seq_along(identified_peaks_over_threshold)) {
    mod <- lm(scan ~ size, identified_peaks_over_threshold[[i]])
    called_scans_list[[i]] <- as.numeric(predict.lm(mod, data.frame(size = ladder_1200)))
  }



  called_size_list <- vector('list', length(identified_peaks_over_threshold))
  for (i in seq_along(identified_peaks_over_threshold)) {
    mod <- lm(size ~ scan, identified_peaks_over_threshold[[i]])
    called_size_list[[i]] <- as.numeric(predict.lm(mod, observed_peaks))
  }




  lag_delta <- function(x){
    delta <- rep(NA_real_, length(x))
    for (i in 2:length(x)) {
      delta[i] <- x[i] - x[i - 1]

    }
    return(delta)
  }

  ladder_delta <- lag_delta(ladder_1200)
  min_ladder_delta <- min(ladder_delta, na.rm = TRUE)
  observed_delta <- lag_delta(called_size_list[[1]])

  ladder_noise_warning <- rep(NA_real_, length(observed_delta))
  warning_value = 1
  for (i in 2:length(observed_delta)) {
    if(observed_delta[i] < min_ladder_delta * 0.9){
      ladder_noise_warning[i - 1] <- warning_value
      ladder_noise_warning[i] <- warning_value
      warning_value <- warning_value + 1
    }
  }



  mod_1_end = identified_peaks[[1]]$scan[length(identified_peaks[[1]]$scan)]
  mod_2_start = identified_peaks[[2]]$scan[1]

  called_scan <- vector("numeric", length(observed_peaks$scan))
  for (i in seq_along(observed_peaks$scan)) {

    if(observed_peaks$scan[i] <= mod_1_end) {
      called_scan[i] <- called_scans_list[[1]][i]
    }
    else if(observed_peaks$scan[i] >= mod_2_start){
      called_scan[i] <- called_scans_list[[2]][i]
    }
    else{
      called_scan[i] <- mean(called_scans_list[[1]][i], called_scans_list[[2]][i])
    }
  }


  nearest_standard <- vector('numeric', nrow(observed_peaks))
  for (i in 1:nrow(observed_peaks)) {
    delta <- abs(observed_peaks$scan[i] - called_scan)

    nearest_standard[i] <- ladder_1200[which(delta == min(delta))]

  }




  # detrend
  #
  # smooth
  #
  # find_peaks
  #
  # pick peaks two std above median signal



  tictoc::tic()

  raw_ladder <- file_list[[1]]$Data$DATA.105
  test <- fit_ladder(raw_ladder,
                     scans = NULL,
                     ladder_sizes = ladder_1200,
                     hq_ladder = FALSE,
                     spike_location = 1000,
                     max_combinations = 2500000)

  tictoc::toc()

  test |>
    ggplot(aes(scan, size)) +
    geom_point()

  data.frame(scan = 0:(length(raw_ladder)-1),
             signal = raw_ladder) |>
    ggplot(aes(scan, signal)) +
    geom_point() +
    geom_text(data = test,
              aes(scan, 10000,
                  label = size),
              angle=90,
              size = 2)


  ladder_list <- identify_ladders(file_list,
                                  ladder_channel = "DATA.105",
                                  signal_channel = "DATA.1",
                                  sample_id_channel = 'SpNm.1',
                                  ladder_sizes=c(50, 75, 100, 139, 150, 160, 200, 300, 350, 400, 450, 490, 500), hq_ladder=T, spike_location = NULL, noise_left = T,
                                  ladder_peak_threshold = 0.985,
                                  ladder_filter_threshold = 0.1,
                                  max_combinations = 2500000)





  plate_combined_df <- extract_peak_table(ladder_list)


  plate_peaks_df <- call_peaks(plate_combined_df,
                               unique_id_col= "unique_id",
                               n_alleles = 2,
                               repeat_unit = 3,
                               allele_height_multiplier = 1,
                               allele_size_gap = 6,
                               bp_lower_limit=100,
                               bp_upper_limit=500,
                               height_threshold = 50,
                               local_peak_scan_window = -2:2,
                               fft_scan_window = 200)




  plot_traces(plate_peaks_df,
              ladder_list,
              names(ladder_list),
              show_called_peaks = TRUE,
              facet_scales = "free",
              zoom_in_on_peak = 2,
              zoom_range = c(-50,50))




})


