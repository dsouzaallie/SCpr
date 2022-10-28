library(shiny)
library("shinythemes")

shinyUI(fluidPage(theme=shinytheme("darkly"),
                  tabsetPanel(type="tabs",
                              #Output Plots onto actual Tabs.
                              tabPanel("File Format.",
                                       fileInput('input$filename$datapath', 'Single Cell Filename'),
                                       actionButton('runAnnotation', 'Format Single Cell File.')),
                              
                              tabPanel("PCAness.",
                                       actionButton('runPCA', 'PCAness.')),
                              
                              tabPanel("Compute Normalization.",
                                       actionButton('normMatrix', 'Normalization.')),
                              tabPanel("Heatmap.",
                                       actionButton('createHeat', 'Heatmap.')),
                              tabPanel("Frequency.",
                                       actionButton('createFrequency', 'Frequency.')),
                              tabPanel("Venn.",
                                       actionButton('createVenn', 'Venn.')),
                              tabPanel("Correlation",
                                       actionButton('createCorrelation', 'Correlation.')),
                         
                  )
)

)