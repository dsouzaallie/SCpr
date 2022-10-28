library(shiny)
library("shinythemes")

shinyUI(fluidPage(theme=shinytheme("darkly"),
                  tabsetPanel(type="tabs",
                              #Output Plots onto actual Tabs.
                              tabPanel("Jurkat File Format.",
                                       fileInput('input$Jurkatfilename$datapath', 'Jurkat Filename'),
                                       actionButton('runJurkatAnnotation', 'Format Jurkat File.')),
                              
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