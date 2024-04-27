test_that("multiplication works", {

  library(shiny)
  library(ggplot2)

  # Sample Data
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


  example_list <- list(
    "20230413_B03.fsa" = data.frame(
      size = c(35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500),
      scan = c(1555, 1633, 1783, 1827, 2159, 2218, 2278, 2525, 2828, 3161, 3408, 3470, 3792, 4085, 4322, 4370)
    )
  )

  suppressWarnings(
    test_ladders_fixed_manual <- fix_ladders_manual(
      test_ladders,
      example_list
    )
  )


  # Shiny App
  ui <- fluidPage(
    titlePanel("Interactive Size Standard Ladder"),
    plotOutput("plot")
  )

  server <- function(input, output) {
    output$plot <- renderPlotly({
  #do this https://stackoverflow.com/questions/74148981/how-to-create-a-draggable-plot-in-r-shiny-using-a-reactive-dataframe
    })
  }

  shinyApp(ui = ui, server = server)



})
