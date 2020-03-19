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
library(shinydashboard)

# Define UI for application that draws a histogram
ui <-  dashboardPage(
    dashboardHeader(title = "RICT Compare Module Visualisation"),
    dashboardSidebar(
        tags$br(),
        fileInput("results","Compared site results",
                  accept = c(
                      "text/csv",
                      "text/comma-separated-values,text/plain",
                      ".csv")),
        tags$br(),
        menuItem("Plot", tabName = "plot", icon = icon("chart-line")),
        menuItem("Table", tabName = "table", icon = icon("table")),
        menuItem("Input data", tabName = "input_data", icon = icon("database")),
        tags$br(),
        #Depreciated code for showing and hiding the table
        #actionButton("view_results", "View or hide input data", width = "250px"),
        selectInput("select_sites", "Select sites compared",
                    c()),
        actionButton("make_plot", "Produce plot")
    ),
    dashboardBody(
        fluidRow(
            tabItems(
                tabItem(tabName = "plot",
                        box(width = 12,
                            tags$h2(textOutput("title_plot"), align="center"),
                            plotOutput("corrPlot")
                        )
                ),

                tabItem(tabName = "table",
                        box(width = 12,
                            tags$h2(textOutput("title_matrix"), align="center"),
                            tableOutput("quality_matrix")
                            )
                ),

                tabItem(tabName = "input_data",
                        box(width = 12,
                            tags$h2("Raw data inputted by user"),
                            tableOutput("results_tab")
                        )
                )
            )
        )
    )
)



# Define server logic required to draw a histogram
server <- function(input, output, session) {

    #create a function to turn the results into matrix format (although data frame in R)
    quality_matrix <- function(x,y){
        z <- NULL

        for(i in y){
            w <- x %>%
                select(starts_with(paste0("Probability.Result.A.in.", i)))

            colnames(w) <- paste0("Probability.Result.A.in.",y)

            z <- rbind(z,w)

        }
        z <- cbind(paste0("Result.A.in.",y),z)

        colnames(z) <- c("",paste0("Result.B.in.",y))

        return(z)
    }

    results_tab <- reactive({
        if(is.null(input$results))     return(NULL)
        results <- read_csv(input$results$datapath)
    })

    observe({
        if(!is.null(results_tab())){
        output$results_tab <- renderTable(results_tab())
        updateSelectInput(session, "select_sites",
                          label = "Select sites compared",
                          choices = paste(results_tab()$'Result B', "to", results_tab()$'Result A')
        )
        }
    })

#Depreciated code for showing and hiding the raw data table
    # observeEvent(input$view_results, {
    #     if(input$view_results %% 2 == 1){
    #         output$results_tab <- renderTable(results_tab())
    #     }
    #     else{
    #         output$results_tab <- NULL
    #     }
    # })

    observeEvent(input$make_plot, {
        if(is.null(input$results))     return(NULL)
        results_obj <- data.frame(results_tab())
        result_A <- gsub(".* to ","",input$select_sites)
        result_B <- gsub(" to .*","",input$select_sites)

        sites_compared <- results_obj[results_obj$Result.A == result_A &
                                          results_obj$Result.B == result_B,]
             # filter('Result A' == result_A &
             #            'Result B' == result_B)

        #select just probability comparisons for the difference in water quality
        comp_probs <- sites_compared %>%
            select(starts_with("Probability.Result.A.in."))

        #The water quality categories
        categories <-  c("High", "Good", "Moderate", "Poor", "Bad")

        #####

        quality_matrix <- quality_matrix(comp_probs, categories)

        output$title_matrix <- renderText({
            paste0("Probability of classifications for ", gsub(" to .*","",input$select_sites), " (Result B) compared to ", gsub(".* to ","",input$select_sites) ," (Result A)")
        })
        output$quality_matrix <- renderTable(quality_matrix)

        expanded_categories <- expand.grid(Result_B = categories, Result_A = categories)

        #attach the probabilities into the site quality combinations
        expanded_categories$probability <- as.numeric(comp_probs)

        output$title_plot <- renderText({
            paste0("Probability of classifications for ", gsub(" to .*","",input$select_sites), " (Result B) compared to ", gsub(".* to ","",input$select_sites) ," (Result A)")
            })

        output$corrPlot <- renderPlot({
            ggplot(expanded_categories, aes(x = Result_B, y = Result_A, fill= probability)) +
                geom_tile()+
                geom_text(aes(label=probability, colour = as.factor(ifelse(probability / 100 < 0.6, 0, 1))), size = 6) +
                scale_fill_gradientn(colours = c("#340044","#2C4279","#1E8644","#64CA44","#FFE51E"), limits = c(0, 100))+
                scale_colour_manual(values = c("#FFFFFF","#000000"))+
                guides(size = "legend", colour = "none")+
                theme(aspect.ratio=1)
        })

        })
}

# Run the application
shinyApp(ui = ui, server = server)

