# modules -----------------------------------------------------------------

plot_module_ui <- function(id) {
  tagList(
    plotlyOutput(NS(id, "plot"))
  )
}

plot_module_server  <- function(id, fragment_ladder) {
  moduleServer(id, function(input, output, session) {

    ladders <- reactiveValues()

    observe({
      ladders$scan <- fragment_ladder$ladder_df$scan
      ladders$size <- fragment_ladder$ladder_df$size
    })

    output$plot <- renderPlotly({

      shapes_with_labels <- list()
      text_annotations <- list()
      for (i in 1:length(ladders$scan)) {
        shapes_with_labels[[i]] <- list(
          type = "rect",
          x0 = ladders$scan[i] - 1,   # Adjust as needed for the positions of your shapes
          x1 = ladders$scan[i] + 1,   # Adjust as needed for the positions of your shapes
          y0 = 0,
          y1 = max(fragment_ladder$trace_bp_df$ladder_signal),
          yref = "paper",
          fillcolor = "rgba(0,0,0,0)",  # Transparent fill
          line = list(
            color = "black",
            width = 1
          ),
          editable = TRUE  # Allow shape editing
        )

        # Add text annotation
        text_annotations[[i]] <- list(
          x = ladders$scan[i] + 19,  # X-position of the text
          y = max(fragment_ladder$trace_bp_df$ladder_signal) / 2,  # Adjust Y-position as needed
          text = ladders$size[i],
          showarrow = FALSE,  # Remove arrow if not desired
          textanchor = "end",  # Horizontal text alignment
          yanchor = "middle",   # Vertical text alignment
          font = list(
            color = "black",
            size = 10
          ),
          textangle = 270
        )
      }

      plot_ly(fragment_ladder$trace_bp_df, x = ~scan, y = ~ladder_signal, type = 'scatter', mode = "lines") %>%
        layout(shapes = shapes_with_labels, annotations = text_annotations) %>%
        # allow to edit plot by dragging lines
        config(edits = list(shapePosition = TRUE))


    })

    # update ladders reactive values in response to changes in shape anchors
    observe({
      ed <- event_data("plotly_relayout")
      scan_positions <- ed[grepl("^shapes.*x.*", names(ed))]
      if (length(scan_positions) != 2) return()
      row_index <- unique(as.numeric(sub(".*\\[(.*?)\\].*", "\\1", names(scan_positions)[1])) + 1)
      new_scans <- round(as.numeric(scan_positions))
      ladders$scan[row_index] <- new_scans[1] + 1 # +1 because the shape is a box 1 above and one below scan
    })

    return(list(ladders = reactive(ladders)))

  })
}



# Shiny App ---------------------------------------------------------------

ui <-  fluidPage(
  fluidRow(column(2,h5("sample:"),
                  selectInput("unique_id_selection", "sample:", NULL)),
           column(10,
                  plot_module_ui("plot_module")

                  )),
  fluidRow(column(2,h5("Ladders:"),tableOutput('data1table')))
)



###





server_function <- function(input, output, session, fragment_trace_list) {

  unique_ids <- names(fragment_trace_list)

  observe({
    updateSelectInput(session, "unique_id_selection", choices = unique_ids)
  })

  selected_fragments_trace <- reactive({

    if(input$unique_id_selection == ""){
      fragment_trace_list[[1]]
    }
    else{
      fragment_trace_list[[which(names(fragment_trace_list) == input$unique_id_selection)]]
    }

  })

  plot_output <- plot_module_server("plot_module", selected_fragments_trace())

  output$data1table <- renderTable({
    as.data.frame(reactiveValuesToList(plot_output$ladders()))
  })
}


interactive_ladders <- function(fragment_trace_list) {
  # Launch the Shiny app with fragment_trace_list passed as a parameter
  shinyApp(
    ui = ui,
    server = function(input, output, session) {
      server_function(input, output, session, fragment_trace_list)
    }
  )
}
