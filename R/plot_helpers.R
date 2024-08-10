
# plot ladder helper ------------------------------------------------------

plot_ladder_helper <- function(fragments_trace,
                               xlim, ylim,
                               plot_title){
  plot(fragments_trace$trace_bp_df$scan, fragments_trace$trace_bp_df$ladder_signal,
       xlab = "Scan", ylab = "Ladder Signal",
       main = ifelse(is.null(plot_title), fragments_trace$unique_id, plot_title),
       type = "l",
       xlim = xlim,
       ylim = ylim
  )

  # Adding text
  text(fragments_trace$ladder_df$scan, rep(max(fragments_trace$trace_bp_df$ladder_signal) / 3, nrow(fragments_trace$ladder_df)),
       labels = fragments_trace$ladder_df$size,
       adj = 0.5, cex = 0.7, srt = 90
  )

  # Adding vertical lines with transparency
  for (i in 1:nrow(fragments_trace$ladder_df)) {
    abline(
      v = fragments_trace$ladder_df$scan[i],
      lty = 3,
      col = rgb(1, 0, 0, alpha = 0.3)
    )
  }
}


# plot_fragments_helper ---------------------------------------------------

plot_fragments_helper <- function(fragment_repeats,
                                  ylim,
                                  xlim,
                                  plot_title){

  if (is.null(fragment_repeats$repeat_table_df)) {
    data <- fragment_repeats$peak_table_df
    data$x <- data$size
  } else {
    data <- fragment_repeats$repeat_table_df
    data$x <- data$repeats
  }

  if (nrow(data) == 0) {
    plot.new()
    title(main = fragment_repeats$unique_id)
    return()
  }

  if (!is.null(xlim)) {
    if (length(xlim == 2) & is(xlim, "numeric")) {
      data <- data[which(data$x < xlim[2] & data$x > xlim[1]), ]
    } else {
      stop(
        call. = FALSE,
        "xlim must be a numeric vector with length of 2"
      )
    }
  }


  allele_1_mode <- ifelse(is.null(fragment_repeats$repeat_table_df), round(fragment_repeats$allele_1_size), round(fragment_repeats$allele_1_repeat))
  allele_2_mode <- ifelse(is.null(fragment_repeats$repeat_table_df), round(fragment_repeats$allele_2_size), round(fragment_repeats$allele_2_repeat))

  # Fill missing y values with zeros
  rounded_x <- round(data$x)
  all_x_values <- seq(min(rounded_x), max(rounded_x))
  y_values <- rep(0, length(all_x_values))
  for (i in seq_along(rounded_x)) {
    y_values[which(all_x_values == rounded_x[i])] <- data[which(data$x == data$x[i]), "height"]
  }



  barplot(
    names.arg = all_x_values,
    height = y_values,
    main = ifelse(is.null(plot_title), fragment_repeats$unique_id, plot_title),
    xlab = ifelse(is.null(fragment_repeats$repeat_table_df), "Size", "Repeat"),
    ylab = "Signal",
    ylim = ylim,
    beside = TRUE,
    col = sapply(all_x_values, function(x) if (!is.na(allele_1_mode) && x == allele_1_mode) "red" else if (!is.na(allele_2_mode) && x == allele_2_mode) "blue" else "gray")
  )

}


# plot_trace helper -------------------------------------------------------

plot_trace_helper <- function(fragments,
                              show_peaks ,
                              x_axis ,
                              ylim ,
                              xlim ,
                              height_color_threshold ,
                              plot_title ){
  if (is.null(fragments$trace_bp_df)) {
    stop(
      call. = FALSE,
      paste(fragments$unique_id, "This sample does not have trace data. Use fsa files as inputs to pipeline to plot trace.")
    )
  }

  #there must be a simpler way of the following if else below
  if (is.null(x_axis) && is.null(fragments$repeat_table_df)) {
    data <- fragments$trace_bp_df
    data$x <- data$size
    x_axis_label <- "Size"
  } else if (is.null(x_axis) && !is.null(fragments$repeat_table_df)) {
    data <- fragments$trace_bp_df
    data$x <- data$calculated_repeats
    x_axis_label <- "Repeats"
  } else if (x_axis == "size") {
    data <- fragments$trace_bp_df
    data$x <- data$size
    x_axis_label <- "Size"
  } else {
    data <- fragments$trace_bp_df
    data$x <- data$calculated_repeats
    x_axis_label <- "Repeats"
  }

  if (!is.null(xlim)) {
    data <- data[which(data$x < xlim[2] & data$x > xlim[1]), ]
  }

  plot(data$x,
       data$signal,
       main = ifelse(is.null(plot_title), fragments$unique_id, plot_title),
       type = "l",
       xlab = x_axis_label,
       ylab = "Signal",
       ylim = ylim
  )


  if (any(data$off_scale)) {
    abline(v = data[which(data$off_scale), "x"], col = adjustcolor("red", alpha.f = 0.3), lwd = 2.5)
  }

  # add points onto plot showing peaks
  if (!is.null(fragments$peak_table_df) && show_peaks) {
    if (is.null(x_axis) && is.null(fragments$repeat_table_df)) {
      peak_table <- fragments$peak_table_df
      peak_table$x <- peak_table$size
    } else if (is.null(x_axis) && !is.null(fragments$repeat_table_df)) {
      peak_table <- fragments$repeat_table_df
      peak_table$x <- peak_table$repeats
    } else if (x_axis == "size") {
      peak_table <- fragments$peak_table_df
      peak_table$x <- peak_table$size
    } else {
      peak_table <- fragments$repeat_table_df
      peak_table$x <- peak_table$repeats
    }

    #exit early if the peak table is empty
    if(nrow(peak_table) == 0){
      return()
    }

    if (!is.null(xlim)) {
      peak_table <- peak_table[which(peak_table$x < xlim[2] & peak_table$x > xlim[1]), ]
    }

    tallest_peak_height <- peak_table[which(peak_table$height == max(peak_table$height)), "height"]
    tallest_peak_x <- peak_table[which(peak_table$height == tallest_peak_height), "x"]
    if (!is.null(fragments$allele_1_height) && !is.na(fragments$allele_1_height)) {
      tallest_peak_height <- fragments$allele_1_height
      #find the tallest peak x axis position
      if (is.null(x_axis) && is.na(fragments$allele_1_repeat)) {
        tallest_peak_x <- fragments$allele_1_size
      } else if (is.null(x_axis) && !is.na(fragments$allele_1_repeat)) {
        tallest_peak_x <- fragments$allele_1_repeat
      } else if (x_axis == "size") {
        tallest_peak_x <- fragments$allele_1_size
      } else {
        tallest_peak_x <- fragments$allele_1_repeat
      }
    }

    peaks_above <- peak_table[which(peak_table$height > tallest_peak_height * height_color_threshold), ]
    peaks_below <- peak_table[which(peak_table$height < tallest_peak_height * height_color_threshold), ]

    # Adding peaks
    points(peaks_above$x,
           peaks_above$height,
           col = "blue"
    )
    points(peaks_below$x,
           peaks_below$height,
           col = "purple"
    )
    points(tallest_peak_x,
           tallest_peak_height,
           col = "green")

    # Draw horizontal dotted lines to connect repeats to their actual place on the plot
    if (!is.null(peak_table$repeats) && !is.null(peak_table$calculated_repeats)) {
      for (i in 1:nrow(peak_table)) {
        segments(x0 = peak_table$repeats[i],
                 y0 = peak_table$height[i],
                 x1 = peak_table$calculated_repeats[i],
                 y1 = peak_table$height[i],
                 lty = 2)
      }
    }

  }


  if(!is.null(fragments$index_repeat) && !is.na(fragments$index_repeat)){
    abline(v = fragments$index_repeat, col = "black", lwd = 2, lty = 3)
  }
}
