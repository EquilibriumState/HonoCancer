#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(ggpubr)
library(DT)
library(dplyr)
load("shortdata.RData")
col_names = colnames(dataClin)
vars <- colnames(cpm_fixed)
# Define UI for application that draws a histogram
ui <- tagList(
    shinythemes::themeSelector(),
    navbarPage(
        theme = "cerulean",  # <--- To use a theme, uncomment this
        "Breast Cancer",
        tabPanel("Gene expression correlation",
    # Application title
    fluidRow(
        sidebarPanel(
            selectInput('xcol', 'X gene', vars),
            selectInput('ycol', 'Y gene', vars, selected = vars[[2]]),
            actionButton("exclude_toggle", "Select"),
            actionButton("exclude_reset", "Reset"),
            helpText("Choose part on image which you want to focus on, and click select",
                     "to present it on the other graph",
                     "If you want to choose specific patients group, go to tab \"Information about patientsp\"",
                     "\n\nCreated by Marcin Radziszewski")
            
        ),
        mainPanel(
            fluidRow(
                column(width = 6,
                       plotOutput("plot1", height = 300,
                                  brush = brushOpts(
                                      id = "plot1_brush",
                                      resetOnNew = TRUE
                                  )
                       )
                ),
                column(width = 6,
                       plotOutput("plot2", height = 300)
                )
            )
            )),
            
                     
        fluidRow(
            column(4,h4("Brushed points"),verbatimTextOutput("brush_info")),
        column(4,h4('All data info'),verbatimTextOutput("summary")),
        column(4,h4('Selected data info'),verbatimTextOutput("summary2")))),
    tabPanel("Information about patients", 
        column(4,
               selectInput("submitter_id",
                           "submitter_id:",
                           c("All",
                             unique(as.character(dataClin$submitter_id))))
        ),
        column(4,
               selectInput("tumor_stage",
                           "Tumor stage:",
                           c("All",
                             unique(as.character(dataClin$tumor_stage))))
        ),
        column(4,
               selectInput("prior_malignancy",
                           "Prior malignancy:",
                           c("All",
                             unique(as.character(dataClin$prior_malignancy))))
        ),
        column(4,
               selectInput("vital_status",
                           "Vital Status:",
                           c("All",
                             unique(as.character(dataClin$vital_status))))
        ),
        column(4,
               selectInput("treatments_radiation_treatment_or_therapy",
                           "Radiation, treatment or therapy:",
                           c("All",
                             unique(as.character(dataClin$treatments_radiation_treatment_or_therapy))))
        ),
        DT::dataTableOutput("dt")
    )
    
)
)
# Define server logic required to draw a histogram
server <- function(input, output) {
    ranges2 <- reactiveValues(x = NULL, y = NULL)
    # For storing which rows have been excluded
    vals <- reactiveValues(
        keeprows = rep(TRUE, nrow(cpm_fixed) ))
    output$plot2 <- renderPlot({
        # Plot the kept and excluded points as two separate data sets
        data <- dataClin
        if (input$submitter_id != "All") {
            data <- data[data$submitter_id == input$submitter_id,]
        }
        if (input$tumor_stage != "All") {
            data <- data[data$tumor_stage == input$tumor_stage,]
        }
        if (input$prior_malignancy != "All") {
            data <- data[data$prior_malignancy == input$prior_malignancy,]
        }
        if (input$vital_status != "All") {
            data <- data[data$vital_status == input$vital_status,]
        }
        if (input$treatments_radiation_treatment_or_therapy != "All") {
            data <- data[data$treatments_radiation_treatment_or_therapy == input$treatments_radiation_treatment_or_therapy,]
        }
        
        for_brush = cpm_fixed[data$submitter_id,c(input$xcol,input$ycol)]
        keep    <- for_brush[ !vals$keeprows, , drop = FALSE]
        
        ggplot(keep, aes_string(x=input$xcol, y=input$ycol)) + 
            geom_point()+
            geom_smooth(method=lm, color="black")+
            labs(title="",
                 x="cpm for gene x", y = "cpm for gene y")+
            theme_classic()
    })
    output$plot1 <- renderPlot({
        data <- dataClin
        if (input$submitter_id != "All") {
            data <- data[data$submitter_id == input$submitter_id,]
        }
        if (input$tumor_stage != "All") {
            data <- data[data$tumor_stage == input$tumor_stage,]
        }
        if (input$prior_malignancy != "All") {
            data <- data[data$prior_malignancy == input$prior_malignancy,]
        }
        if (input$vital_status != "All") {
            data <- data[data$vital_status == input$vital_status,]
        }
        if (input$treatments_radiation_treatment_or_therapy != "All") {
            data <- data[data$treatments_radiation_treatment_or_therapy == input$treatments_radiation_treatment_or_therapy,]
        }
        
        for_brush = cpm_fixed[data$submitter_id,c(input$xcol,input$ycol)]
        ggplot(for_brush, aes_string(x=input$xcol, y=input$ycol)) + 
            geom_point()+
            geom_smooth(method=lm, color="black")+
            labs(title="",
                 x="cpm for gene x", y = "cpm for gene y")+
            theme_classic() 
    })
    observe({
        brush <- input$plot1_brush
        if (!is.null(brush)) {
            ranges2$x <- c(brush$xmin, brush$xmax)
            ranges2$y <- c(brush$ymin, brush$ymax)
            
        } else {
            ranges2$x <- NULL
            ranges2$y <- NULL
        }
    })
    output$summary <- renderPrint({
        cor.test(cpm_fixed[,input$xcol], cpm_fixed[,input$ycol], method="pearson")
    })
    output$summary2 <- renderPrint({
        keep    <- cpm_fixed[ !vals$keeprows, , drop = FALSE]
        cor.test(keep[,input$xcol], keep[,input$ycol], method="pearson")
    })
    # Toggle points that are brushed, when button is clicked
    observeEvent(input$exclude_toggle, {
        vals$keeprows <- rep(TRUE, nrow(cpm_fixed))
        res <- brushedPoints(cpm_fixed, input$plot1_brush, allRows = TRUE)
        
        vals$keeprows <- xor(vals$keeprows, res$selected_)
    })
    # Reset all points
    observeEvent(input$exclude_reset, {
        vals$keeprows <- rep(TRUE, nrow(cpm_fixed))
    })
    
    output$dt <- DT::renderDataTable(DT::datatable(filter = 'top',{
        data <- dataClin
        if (input$submitter_id != "All") {
            data <- data[data$submitter_id == input$submitter_id,]
        }
        if (input$tumor_stage != "All") {
            data <- data[data$tumor_stage == input$tumor_stage,]
        }
        if (input$prior_malignancy != "All") {
            data <- data[data$prior_malignancy == input$prior_malignancy,]
        }
        if (input$vital_status != "All") {
            data <- data[data$vital_status == input$vital_status,]
        }
        if (input$treatments_radiation_treatment_or_therapy != "All") {
            data <- data[data$treatments_radiation_treatment_or_therapy == input$treatments_radiation_treatment_or_therapy,]
        }
        data
    }))
    output$brush_info <- renderPrint({
        
        data <- dataClin
        if (input$submitter_id != "All") {
            data <- data[data$submitter_id == input$submitter_id,]
        }
        if (input$tumor_stage != "All") {
            data <- data[data$tumor_stage == input$tumor_stage,]
        }
        if (input$prior_malignancy != "All") {
            data <- data[data$prior_malignancy == input$prior_malignancy,]
        }
        if (input$vital_status != "All") {
            data <- data[data$vital_status == input$vital_status,]
        }
        if (input$treatments_radiation_treatment_or_therapy != "All") {
            data <- data[data$treatments_radiation_treatment_or_therapy == input$treatments_radiation_treatment_or_therapy,]
        }
        
        for_brush = cpm_fixed[data$submitter_id,c(input$xcol,input$ycol)]
        brushedPoints(for_brush, input$plot1_brush)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

