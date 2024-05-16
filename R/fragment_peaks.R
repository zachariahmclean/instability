# peak calling ------------------------------------------------------------


find_fragment_peaks <- function(trace_bp_df,
                                smoothing_window,
                                minimum_peak_signal,
                                ...) {
  # smoothed_signal <- moving_average(trace_bp_df$signal,
  #                                   n = smoothing_window)
  #
  smoothed_signal <- pracma::savgol(
    trace_bp_df$signal,
    smoothing_window
  )

  # deals with cases of user overriding values
  if ("peakpat" %in% ...names()) {
    peaks <- pracma::findpeaks(smoothed_signal,
      minpeakheight = minimum_peak_signal,
      ...
    )
  } else if ("minpeakheight" %in% ...names()) {
    peaks <- pracma::findpeaks(smoothed_signal,
      peakpat = "[+]{6,}[0]*[-]{6,}", # see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
      ...
    )
  } else {
    peaks <- pracma::findpeaks(smoothed_signal,
      peakpat = "[+]{6,}[0]*[-]{6,}", # see https://stackoverflow.com/questions/47914035/identify-sustained-peaks-using-pracmafindpeaks
      minpeakheight = minimum_peak_signal,
      ...
    )
  }

  n_scans <- length(trace_bp_df$signal)
  window_width <- 3

  # go through raw signal and make sure that the identified scan in the smoothed signal is still the highest
  # it will also deal with cases where the scans have the same height (which.max will chose first)
  peak_position <- numeric(nrow(peaks))
  for (i in seq_along(peak_position)) {
    if (peaks[i, 2] + window_width > 1 & peaks[i, 2] + window_width < n_scans) { # make sure that the subsetting would be in bounds when taking window into account
      max_peak <- which.max(trace_bp_df$signal[(peaks[i, 2] - window_width):(peaks[i, 2] + window_width)])

      peak_position[i] <- peaks[i, 2] - window_width - 1 + max_peak
    } else {
      peak_position[i] <- peaks[i, 2]
    }
  }

  df <- trace_bp_df[peak_position, c("scan", "size", "signal", "off_scale")]
  colnames(df) <- c("scan", "size", "height", "off_scale")

  # remove shoulder peaks
  df2 <- deshoulder(df, shoulder_window = 1.5)

  return(df2)
}
