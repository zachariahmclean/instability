# modules -----------------------------------------------------------------
## select sample module
sample_selection_ui <- function(id){
  # shinyWidgets::pickerInput(shiny::NS(id, "unique_id_selection"), "Sample:", NULL,
  #                           options = shinyWidgets::pickerOptions(iconBase = "fas")
  #                           )
  shiny::selectInput(shiny::NS(id, "unique_id_selection"), "Sample:", NULL)


}

sample_selection_module <- function(id, fragment_trace_list){
  shiny::moduleServer(id, function(input, output, session){


    # ladder_warning <- shiny::reactive({
    #   sapply(reactiveValuesToList(fragment_trace_list),
    #          function(x){
    #            if(is.null(tryCatch(ladder_rsq_warning_helper(x, 0.998),
    #                                warning = function(w) w))){
    #              NA
    #            }
    #            else{
    #              "fa-triangle-exclamation"
    #            }
    #          })
    # })
#
#
#     shiny::observe({
#       shinyWidgets::updatePickerInput(session, "unique_id_selection", choices = names(fragment_trace_list),
#                                       choicesOpt = list(icon = ladder_warning())
#                                       )
#     })


    shiny::observe({
      shiny::updateSelectInput(session, "unique_id_selection", choices = names(fragment_trace_list)
      )
    })


    selected_fragments_trace <- shiny::reactive({

      # if(is.null(input$unique_id_selection)){
      if(input$unique_id_selection == ""){

        fragment_trace_list[[names(fragment_trace_list)[1]]]
      }
      else{
        fragment_trace_list[[input$unique_id_selection]]
      }

    })

    return(list(sample = selected_fragments_trace,
                input_unique_id_selection = shiny::reactive(input$unique_id_selection)))

  })
}
## plot module
plot_module_ui <- function(id) {
  shiny::tagList(
    plotly::plotlyOutput(shiny::NS(id, "plot"))
  )
}

plot_module_server <- function(id, fragment_ladder, input_unique_id_selection) {
  shiny::moduleServer(id, function(input, output, session) {

    # Initialize ladders as NULL
    ladders <- shiny::reactiveValues(scan = NULL, size = NULL)
    relayout_data <- shiny::reactiveVal(NULL)  # Initialize relayout_data

    # Reset ladders and relayout_data when unique_id_selection changes
    shiny::observeEvent(input_unique_id_selection(), {
      ladders$scan <- NULL
      ladders$size <- NULL
      relayout_data(NULL)
    })

    shiny::observe({
        ladders$scan <- fragment_ladder()$ladder_df$scan
        ladders$size <- fragment_ladder()$ladder_df$size
    })

    output$plot <- plotly::renderPlotly({
      if (is.null(ladders$scan) || is.null(ladders$size)) {
        # Return a blank plot if ladders are not initialized
        return(plotly::plot_ly())
      }

      shapes_with_labels <- list()
      text_annotations <- list()
      for (i in 1:length(ladders$scan)) {
        shapes_with_labels[[i]] <- list(
          type = "rect",
          x0 = ladders$scan[i] - 1,   # Adjust as needed for the positions of your shapes
          x1 = ladders$scan[i] + 1,   # Adjust as needed for the positions of your shapes
          y0 = 0,
          y1 = max(fragment_ladder()$trace_bp_df$ladder_signal),
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
          y = max(fragment_ladder()$trace_bp_df$ladder_signal) / 2,  # Adjust Y-position as needed
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

      p <- plotly::plot_ly(fragment_ladder()$trace_bp_df, x = ~scan, y = ~ladder_signal, type = 'scatter', mode = "lines")
      p <- plotly::layout(p, shapes = shapes_with_labels, annotations = text_annotations, title = fragment_ladder()$unique_id)
        # allow to edit plot by dragging lines
      plotly::config(p, edits = list(shapePosition = TRUE))
    })

    # Reset relayout_data when plot is clicked or dragged
    shiny::observeEvent(plotly::event_data("plotly_relayout"), {
      relayout_data(plotly::event_data("plotly_relayout"))
    })

    # Capture relayout_data
    shiny::observe({
      if (!is.null(relayout_data())) {
        ed <- relayout_data()
        scan_positions <- ed[grepl("^shapes.*x.*", names(ed))]
        if (length(scan_positions) != 2) return()
        row_index <- unique(as.numeric(sub(".*\\[(.*?)\\].*", "\\1", names(scan_positions)[1])) + 1)
        new_scans <- round(as.numeric(scan_positions))
        ladders$scan[row_index] <- new_scans[1] + 1 # +1 because the shape is a box 1 above and one below scan
      }
    })

    # Reset other reactive values if needed

    return(list(ladders = shiny::reactive(ladders)))
  })
}

## export ladder fixes

ladder_export_ui <- function(id) {
  shiny::tagList(
    shiny::downloadButton(shiny::NS(id, "download"))
  )
}

ladder_export_server <- function(id, manual_ladder_list) {
  shiny::moduleServer(id, function(input, output, session) {
    output$download <- shiny::downloadHandler(
      filename = function() {
        paste0(format(Sys.time(), '%Y-%m-%d_%H%M%S'),"_ladder_df_list", ".rds")
      },
      content = function(file) {
        saveRDS(shiny::reactiveValuesToList(manual_ladder_list), file)
      }
    )
  })
}



# Shiny App ---------------------------------------------------------------

ui <-  shiny::fluidPage(
  shiny::fluidRow(shiny::column(2,
                  sample_selection_ui("sample_selection"),
                  ladder_export_ui("ladder_df_list_download")
                  ),
                  shiny::column(10,
                  plot_module_ui("plot_module")

                  )),
  shiny::fluidRow(shiny::column(2,shiny::h5("Ladders:"),shiny::tableOutput('data1table')))
)



###
server_function <- function(input, output, session, fragment_trace_list) {


  fragment_trace_list_reactive <- shiny::reactiveValues()
  for (sample_name in names(fragment_trace_list)) {
    fragment_trace_list_reactive[[sample_name]] <- fragment_trace_list[[sample_name]]
  }
  manual_ladder_list <- shiny::reactiveValues()

  selected_fragments_trace <- sample_selection_module("sample_selection", fragment_trace_list_reactive)


  plot_output <- plot_module_server("plot_module",
                                    selected_fragments_trace$sample,
                                    selected_fragments_trace$input_unique_id_selection)

  output$data1table <- shiny::renderTable({
    as.data.frame(shiny::reactiveValuesToList(plot_output$ladders()))
  })

  # have a reactive list that gets updated when you change the stuff
  shiny::observe({
    sample_unique_id <- selected_fragments_trace$sample()$unique_id

    selected_ladder_df <- selected_fragments_trace$sample()$ladder_df
    selected_sample_scans <- selected_ladder_df[which(!is.na(selected_ladder_df$size)), "scan"]

    plot_ladder_df <- as.data.frame(reactiveValuesToList(plot_output$ladders()))
    plot_scans <- plot_ladder_df[which(!is.na(plot_ladder_df$size)), "scan"]

    #skip if ladder info hasn't been updated
    if(identical(selected_sample_scans, plot_scans)){
      return()
    }
    else if(nrow(as.data.frame(reactiveValuesToList(plot_output$ladders()))) == 0 ){
      return()
    }

    manual_ladder_list[[sample_unique_id]] <- as.data.frame(reactiveValuesToList(plot_output$ladders()))
    fragment_trace_list_reactive[[sample_unique_id]] <- fix_ladders_manual(
      reactiveValuesToList(fragment_trace_list_reactive)[sample_unique_id],
      reactiveValuesToList(manual_ladder_list)[sample_unique_id]
    )[[1]]

  })

  # export data
  ladder_export_server("ladder_df_list_download", manual_ladder_list)

}




fix_ladder_interactive <- function(fragment_trace_list) {

  # Launch the Shiny app with fragment_trace_list passed as a parameter
  shiny::shinyApp(
    ui = ui,
    server = function(input, output, session) {
      server_function(input, output, session, fragment_trace_list)
    }
  )
}
