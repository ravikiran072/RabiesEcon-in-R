####################################
# Ravikiran Keshavamurthy                   #
# Shiny app for RabiesEcon (CDC) #
# 7/11/2023
####################################


# Import libraries
library(shiny)
library(shinythemes)
library(data.table)
library(RCurl)


# # Read data
# weather <- read.csv(text = getURL("https://raw.githubusercontent.com/dataprofessor/data/master/weather-weka.csv"))
# weather$play <- as.factor(weather$play)
# weather$outlook <- as.factor(weather$outlook)
setwd("C:/Users/mfl3/OneDrive - CDC/Rabies/Rabies_R0/RabiesEcon_in_R/REcon_shiny")

# Build model
source("Input.R")
source("initial_run.R")

# Save model to RDS file
# saveRDS(model, "model.rds")

# Read in the RF model
#model <- readRDS("model.rds")

####################################
# User interface                   #
####################################

ui <- fluidPage(theme = shinytheme("united"),
                
                # Page header
                headerPanel('RabiesEcon'),
                
                # Input values
                sidebarPanel(
                  tags$label(h3('Input parameters')),
                  numericInput("Km2_of_program_area", 
                               label = "Square kilometers (km2) of program area", 
                               value = 14000),
                  actionButton("submitbutton", "Submit", class = "btn btn-primary")
                ),
                
                mainPanel(
                  tags$label(h3('Results-Fig')), # Status/Output Text Box
                  verbatimTextOutput('contents'),
                  tableOutput('tabledata') # Prediction results table
                  
                )
)

####################################
# Server                           #
####################################

server <- function(input, output, session) {
  
  # Input Data
  datasetInput <- reactive({  
    Km2_of_program_area <- input$Km2_of_program_area 
    # outlook,temperature,humidity,windy,play
    source("Input.R")
    source("initial_run.R")
    
    
    #Output <- data.frame(R0_dog_to_dog)
    
    
  })
  
  # Status/Output Text Box
  output$contents <- renderPrint({
    if (input$submitbutton>0) { 
      isolate("Calculation complete.") 
    } else {
      return("Server is ready for calculation.")
    }
  })
  
  # Prediction results table
  output$tabledata <- renderTable({
    print(head(initial_run))
  })
  
}

####################################
# Create the shiny app             #
####################################
shinyApp(ui = ui, server = server)
