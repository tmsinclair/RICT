#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("hSSD"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            textInput("chem", "Chemical", ""),
            fileInput("tox","Known toxicity data",
                      accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv")),
            actionButton("view_tox", "View or Hide"),
            fileInput("tax_tox","Taxonomy of toxicity data",
                      accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv")),
            actionButton("view_tax_tox", "View or Hide"),
            fileInput("tax_pred","Taxonomy of assemblage to be predicted",
                      accept = c(
                          "text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv")),
            actionButton("view_tax_pred", "View or Hide"),
            sliderInput("samples", "Number of samples", 500, min = 100, max = 1000, step = 5),
            sliderInput("burn", "Number to burn", 10, min = 0, max = 100),
            actionButton("run_hSSD", "Run hSSD"),
            selectInput("dataset", "Choose a dataset:",
                        choices = c("None","Predicted toxicity", "hSSD Summary"))
            #        radioButtons("filetype", "File type:",
            #                     choices = c("csv", "tsv")),
            #        downloadButton('downloadData', 'Download')
        ),

        # Show a plot of the generated distribution
        mainPanel(
            tableOutput("tox"),
            tableOutput("tax_tox"),
            tableOutput("tax_pred"),
            plotOutput("hSSD"),
            dataTableOutput("pred_tox")
        )
    )
)

#####
# Define server logic required to draw a histogram
server <- function(input, output) {

    #####
    tox_tab <- reactive({
        if(is.null(input$tox))     return(NULL)
        tox <- read_csv(input$tox$datapath)
    })

    observeEvent(input$view_tox, {
        if(input$view_tox %% 2 == 1){
            output$tox <- renderTable(tox_tab())
        }
        else{
            output$tox <- NULL
        }
    })

    tax_tox_tab <- reactive({
        if(is.null(input$tax_tox))     return(NULL)
        tax_tox <- read_csv(input$tax_tox$datapath)
    })

    observeEvent(input$view_tax_tox, {
        if(input$view_tax_tox %% 2 == 1){
            output$tax_tox <- renderTable(tax_tox_tab())
        }
        else{
            output$tax_tox <- NULL
        }
    })

    tax_pred_tab <- reactive({
        if(is.null(input$tax_pred))     return(NULL)
        else{tax_pred <- read_csv(input$tax_pred$datapath)
        }
    })

    observeEvent(input$view_tax_pred, {
        if(input$view_tax_pred %% 2 == 1){
            output$tax_pred <- renderTable(tax_pred_tab())
        }
        else{
            output$tax_pred <- NULL
        }
    })

    #####
    observeEvent(input$run_hSSD, {

    })
        #####
        output$hSSD <- renderPlot({
            ggplot() +
                #(data = bootdat, aes(x = newxs, y = value, group = variable), col = 'steelblue', alpha = 0.05) +
                geom_point(data = pred_tox, aes(x = tox, y = frac)) +
                geom_line(data = pdat, aes(x = newxs, y = py), col = 'red') +
                geom_line(data = pdat, aes(x = newxs, y = lwr), linetype = 'dashed') +
                geom_line(data = pdat, aes(x = newxs, y = upr), linetype = 'dashed') +
                geom_text(data = pred_tox, aes(x = fit, y = frac, label = latin), hjust = 1, size = 2.5) +
                theme_bw() +
                scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000,10000,100000), limits = c(0.003, max(pred_tox$tox))) +
                labs(x = paste('Concentration of ',input$chem," (ug/L)",sep = ""),
                     y = 'Fraction of species affected')
        })
    })

    observeEvent(input$view_tax_pred, {
        if(input$view_tax_pred %% 2 == 1){
            output$table <- renderTable(tax_pred_tab())
        }
        else{
            output$table <- NULL
        }
    })

    output$downloadData <- downloadHandler(

        # This function returns a string which tells the client
        # browser what name to use when saving the file.
        filename = function() {
            paste(input$dataset, input$filetype, sep = ".")
        },

        # This function should write data to a file given to it by
        # the argument 'file'.
        content = function(file) {
            sep <- switch(input$filetype, "csv" = ",", "tsv" = "\t")

            # Write to a file specified by the 'file' argument
            write.table(datasetInput(), file, sep = sep,
                        row.names = FALSE)
        }
    )

}

# Run the application
shinyApp(ui = ui, server = server)
