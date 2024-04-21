test_that("multiplication works", {

#   # fsa_list <- instability::cell_line_fsa_list[1]
#   file_list <- list.files("C:/Users/zlm2/R/Fragment analysis/2024-4-20 instability mouse test/data/mouse_fsa/", recursive = TRUE, full.names = TRUE)
#   fsa_list <- read_fsa(file_list[grep("321417m_B09_B9.fsa",file_list)])
#
#   suppressWarnings(
#     test_ladders <- find_ladders(fsa_list,
#                                  ladder_channel = "DATA.105",
#                                  signal_channel = "DATA.1",
#                                  ladder_sizes = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
#                                  hq_ladder = FALSE,
#                                  max_combinations = 2500000,
#                                  ladder_selection_window = 8)
#   )
#
#   test_fragments <- find_fragments(test_ladders,
#                  peakpat = "[+]{6,}[0]*[-]{6,}")
#
#   test_fragments |>
#     plot_traces(xlim = c(425, 430),
#                 ylim = c(0,600), n_facet_col = 1)
#
#
#
#   tmp <- extract_fragments(test_fragments)
#
#   for(i in seq_along(tmp)){
#
#   }
#
#
#
#
#   smoothed_signal <- moving_average(test_ladders[[1]]$trace_bp_df$signal,
#                                     n = 5)
#
#
#   smoothed_savgol <- pracma::savgol(test_ladders[[1]]$trace_bp_df$signal,
#                                     21)
#
#   #smoothed_whittaker <- pracma::whittaker(test_ladders[[1]]$trace_bp_df$signal)
#
#
#
# bind_rows(
#  data.frame(
#    scan = test_ladders[[1]]$trace_bp_df$scan,
#    smooth = smoothed_signal,
#    smoothing_alg = rep("zm",length(test_ladders[[1]]$trace_bp_df$scan))
#  ),
#  data.frame(
#    scan = test_ladders[[1]]$trace_bp_df$scan,
#    smooth = smoothed_savgol,
#    smoothing_alg = rep("smoothed_savgol",length(test_ladders[[1]]$trace_bp_df$scan))
#  ),
#  data.frame(
#    scan = test_ladders[[1]]$trace_bp_df$scan,
#    smooth = test_ladders[[1]]$trace_bp_df$signal,
#    smoothing_alg = rep("none",length(test_ladders[[1]]$trace_bp_df$scan))
#  )
#
#  ) |>
#   filter(between(scan, 10190, 10250)) |>
#   ggplot(aes(scan , smooth, colour = smoothing_alg)) +
#   geom_line() +
#   geom_point()
#
#
#
#
#
# #
#
#
#
#
# pracma::findpeaks(smoothed_savgol,
#                   peakpat = "[+]{5,}[0]*[-]{5,}")
#
#



})
