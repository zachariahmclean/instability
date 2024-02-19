
#what is going on here?
## I think I was trying out ladder recursive as an option for fun
## is the actual code the one to use?


test_that("find ladder peaks",{

  files <- list.files("data/GT_Z.McLean_2023-04-14/", full.names = TRUE)

  file_list <- read_fsa(files[1])

  test_processed <- process_ladder_signal(file_list[[1]]$Data$DATA.105,
                                          scans = 0:(length(file_list[[1]]$Data$DATA.105) -1),
                                          spike_location = 1000)


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

  testthat::expect_true(round(mod$coefficients[[1]],3) == 1348.234)
  testthat::expect_true(round(mod$coefficients[[2]],5) == 6.24925)

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


  testthat::expect_true(round(rsq_diff,5) == -0.00276)


})




test_that("fit ladder bug where the selected peaks are different length to ladder", {

  #still need to fix where it has NA in size for some reason?


  files <- list.files("data/GT_Z.McLean_2023-04-14/", full.names = TRUE)

  file_list <- read_fsa("data/GT_Z.McLean_2023-04-14/20230413_C04.fsa")

  test_position <- which(names(file_list) == "20230413_E10.fsa")


  test_ladder_signal <- file_list[[test_position]]$Data$DATA.105
  test_scans <- 0:(length(file_list[[test_position]]$Data$DATA.105) -1)

  #not fixed! too many scans are being selected


  test_fit <- fit_ladder(
    ladder = test_ladder_signal,
    scans = test_scans,
    ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
    hq_ladder = FALSE
  )


  mod <- lm(scan ~ size, data = test_fit)

  testthat::expect_true(round(mod$coefficients[[1]],3) == 1348.234)
  testthat::expect_true(round(mod$coefficients[[2]],5) == 6.24925)

})




test_that("find ladders", {

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


})



#trying out fitting a polynomial to detect peaks. Scan along whole trace and find polynomial inflections


testthat::test_that("peak detection", {

 #  library(pracma)
 #  library(ggplot2)
 #
 # # Option 1:
 #
 #
 #  detect_peaks <- function(data, polynomial_degree, peak_window_size) {
 #    require(pracma) # for polyfit() and polyder() functions
 #
 #    # Initialize variables to store detected peaks
 #    peaks <- numeric(0)
 #
 #    # Iterate through each data point
 #    for (i in 1:(length(data) - peak_window_size)) {
 #      # Define window boundaries
 #      start <- max(1, i - floor(peak_window_size / 2))
 #      end <- min(length(data), i + floor(peak_window_size / 2))
 #
 #      # Fit polynomial curve to data within the window
 #      poly_fit <- tryCatch(pracma::polyfit(start:end, data[start:end], polynomial_degree), error = function(e) NULL)
 #
 #      # Check if polynomial fit is successful
 #      if (!is.null(poly_fit)) {
 #        # Compute second derivative of the polynomial curve
 #        second_derivative <- pracma::polyder(pracma::polyder(poly_fit))
 #
 #        # Plot polynomial fit
 #        plot(x = start:end, y = data[start:end], type = "l", ylim = range(data), main = "Polynomial Fit",
 #             xlab = "Data Point Index", ylab = "Signal Value")
 #        lines(start:end, polyval(poly_fit, start:end), col = "red")
 #
 #        # Check if the second derivative changes sign at the current data point
 #        if (!is.na(second_derivative[1]) && !is.na(second_derivative[2]) && second_derivative[1] < 0 && second_derivative[2] > 0) {
 #          points(i, data[i], col = "blue", pch = 16)
 #          peaks <- c(peaks, i)
 #        }
 #      }
 #    }
 #
 #    # Return the positions of detected peaks
 #    return(peaks)
 #  }
 #
 #  #option two
 #
 #
 #  # Function to detect peaks from signal data
 #  detect_peaks <- function(signal_data) {
 #    # Calculate the first derivative
 #    first_derivative <- diff(signal_data)
 #
 #    # Find zero crossings in the first derivative
 #    zero_crossings <- which(first_derivative[-1] * first_derivative[-length(first_derivative)] < 0)
 #
 #    # Initialize list to store peak positions
 #    peak_positions <- c()
 #
 #    # Identify peak positions
 #    for (i in zero_crossings) {
 #      if (first_derivative[i] > 0) {
 #        peak_positions <- c(peak_positions, i)
 #      }
 #    }
 #
 #    # Represent peaks in a list
 #    peak_list <- list()
 #
 #    # Extract characteristics for each peak
 #    for (position in peak_positions) {
 #      peak <- list()
 #      peak$retention_time <- position # Retention time
 #      peak$amplitude <- signal_data[position] # Amplitude
 #
 #      # Calculate peak shape using cubic spline interpolation
 #      # Assuming k = 1 (3 consecutive time points on each side of the peak center)
 #      spline_points <- stats::spline(x = (position - 1):(position + 1), y = signal_data[(position - 1):(position + 1)], n = 5 * 2 + 1)
 #      peak$shape <- spline_points$y
 #      browser()
 #
 #      peak_list <- c(peak_list, list(peak))
 #    }
 #
 #    return(peak_list)
 #  }
 #
 #  # Example usage
 #  peaks <- detect_peaks(file_list[[1]]$Data$DATA.1)
 #  print(peaks)
 #
 #
 #
 #  # Example usage
 #
 #  files <- list.files("data/GT_Z.McLean_2023-04-14/", full.names = TRUE)
 #
 #  file_list <- read_fsa("data/GT_Z.McLean_2023-04-14/20230413_C04.fsa")
 #
 #
 #  test_scan <-   4200:4250
 #  test_signal <- file_list[[1]]$Data$DATA.1[test_scan]
 #
 #  # Set peak detection parameters
 #  polynomial_degree <- 3
 #  peak_window_size <- 16
 #
 #
 #  # Detect peaks
 #  detected_peaks <- detect_peaks(test_signal, polynomial_degree, peak_window_size)
 #  print(detected_peaks)
 #
 #
 #
 #  data.frame(scan = test_scan,
 #             signal = test_signal) |>
 #    ggplot(aes(scan, signal)) +
 #    geom_line() +
 #    geom_vline(xintercept = detected_peaks + 4200)
 #

})

