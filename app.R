library(shinyalert) # pop-up error messages
library(shinycssloaders) # spinning plot loading icon
library(icesTAF) # contains mkdir
library(tidyverse)
library(shiny)
source("dye_calling_support_script.R")

library(shiny)

facet_grid_validation <- function(df_melt, title) {
    p <- ggplot(df_melt, aes(x = Temperature, # temperature on X
                             y = value, # RFU on y
                             color = conc, # colored by the state
                             linetype = protein,
                             group = dye_conc_protein_channel # group means series, as in, this defines the unique data sets
    )) +
        geom_line(size = 0.3, alpha = 0.8) + # change the line type depending on the dye concentration # linetype = df_melt$conc #
        facet_grid(dye~channel_f, scales = "free") +
        labs(title = title, color = "[Dye] (µM)") +
        scale_color_viridis_c(end = 0.2, begin = 1) +
        theme_bw() +
        
        scale_linetype_manual(values = c("dashed", "solid")) +
        facet_no_y_theme
    
    p # return the plot
}

facet_wrap_validation <- function(df_melt, title) {
    p <- ggplot(df_melt, aes(x = Temperature, # temperature on X
                             y = value, # RFU on y
                             color = conc, # colored by the state
                             linetype = protein,
                             group = dye_conc_protein_channel # group means series, as in, this defines the unique data sets
    )) +
        geom_line(size = 0.3, alpha = 0.8) + # change the line type depending on the dye concentration # linetype = df_melt$conc #
        facet_wrap(~dye_channel_f, scales = "free",  ncol = 6) +
        labs(title = title, color = "[Dye] (µM)") +
        scale_color_viridis_c(end = 0.2, begin = 1) +
        theme_bw() +
        
        scale_linetype_manual(values = c("dashed", "solid")) +
        facet_no_y_theme
    
    p # return the plot
}



# Define server logic required to draw a histogram
server <- function(input, output) {
    
    values <- reactiveValues()
    
    ###### when the analyze button is clicked, this is when almost everything happens
    observeEvent(input$analyze_button, {
    #### read in the layout
    tryCatch({
            protein_layout_path <- input$protein_layout$datapath  ##### NEEDS UI ELEMENT
            print(protein_layout_path)
            values$protein_layout <-  make_layout(input$protein_layout$datapath) %>%
                                        rename(dye = compound) %>%
                rename(conc = concentration)
            
    }, # end tryCatch function
    error = function(e){
        shinyalert("Protein layout file failed to upload!", "Please check file and try again.")
        values$protein_layout <- NULL
    }
    )
    
    ##### read in the raw data files
   tryCatch({
        protein_path <- input$protein_file$datapath ##### NEEDS UI ELEMENT
        print(protein_path)
        values$protein_file <-  read_qTower(protein_path) %>% 
            mutate(screen = rep(input$exp_num, times = nrow(.)), ##### NEEDS UI ELEMENT
                   Temperature = .$Temperature + (input$start_T-1))

        values$protein_file %>% head()

    }, # end tryCatch function
    error = function(e){
        shinyalert("Protein data file failed to upload!", "Please check file and try again.")
        values$protein_file <- NULL
    }
    )

    # join layouts with data files, and join buffer and protein files
    tryCatch({
        
        values$df_all <- left_join(values$protein_file, values$protein_layout, by = "well") %>%
                            filter(dye != "Empty") %>%
                            unite(dye_conc_protein_channel, c(dye, conc, protein, channel), sep = "-", remove = FALSE) %>%
                            mutate(conc = as.numeric(conc)) %>%
                            unite("dye_channel", c(dye, channel_f), remove = FALSE) %>% 
                            arrange( dye, channel_f) %>% 
                            mutate(dye_channel_f = factor(dye_channel, levels = unique(dye_channel)))
        # %>% filter(dye == "SYPRO") # use this when testing the app
        
        values$new_vals <- tibble(dye_channel_f = values$df_all$dye_channel_f, 
                            assignment = "none") %>%
                            distinct(dye_channel_f, .keep_all = TRUE)
        
    }, # end tryCatch function
    error = function(e){
        shinyalert("Protein, buffer, and layout files don't match!", "Please check files and try again.")
        values$protein_file_layout <- NULL
        values$buffer_file_layout <- NULL
        values$df_all <- NULL
    })
    
        
    #### create the interactive hit-calling plot 
    output$plot <- renderPlot( values$plot ) 
    
    observeEvent({input$plot_dblclick$panelvar1}, {
        print("plot double click")
        values$new_vals <- values$new_vals %>%
            mutate( assignment = replace(assignment, dye_channel_f == input$plot_dblclick$panelvar1, input$double_click_meaning))
        print(values$new_vals)
    })
    
    observeEvent({input$plot_click$panelvar1}, {
        print("plot  click")
        values$new_vals <- values$new_vals %>%
            mutate(assignment = replace(assignment, dye_channel_f == input$plot_click$panelvar1, input$single_click_meaning)) %>%
            mutate( assignment = replace(assignment, dye_channel_f == input$plot_brush$panelvar1, "none"))
        print(values$new_vals)
    })
    
    observeEvent(input$plot_brush$panelvar1, {
        print("plot brush")
        values$new_vals <- values$new_vals %>%
            mutate( assignment = replace(assignment, dye_channel_f == input$plot_brush$panelvar1, "none"))
        print(values$new_vals)
    })
    
    output$data_table <-  DT::renderDataTable({  
        values$new_vals  %>%
            filter(assignment != "none") %>%
            arrange(assignment)
        })
    
    ###### download the plot
    
    output$download_files <- downloadHandler(
        filename = function() {
            
            paste0(input$exp_num, "_", as.character(base::Sys.Date()), "_df_all.rds")
        },
        content = function(file) {
            # write the RDS
            write_rds(values$df_all, file)
        }
    )
    
    output$download_plot <- downloadHandler(
        filename = function() {
            
            paste0(input$exp_num, "_", as.character(base::Sys.Date()), "all_screen_plot.pdf")
        },
        content = function(file) {
            # save the primary plot
            values$plot_ratio <- values$df_all %>%
                select(dye) %>%
                distinct() %>%
                nrow()
            
            ggsave(file, values$plot, width = 10,  height =  1.5*(1/1.618 * values$plot_ratio + 1))  # ,
            
        }
    )
    
    ### download the hit list
    output$save_dt_button <- downloadHandler(
        filename = function() {
            paste0(input$exp_num, "_", as.character(base::Sys.Date()), "_manually_called_hits.csv")
        },
        content = function(file) {
            write.csv(values$new_vals %>% arrange(assignment), file)
        }
    )
    
    values$plot <- values$df_all %>%
        facet_wrap_validation( . , title = paste0(input$exp_num, ": validation"))
    
    }) ## end the observeEvent actions for the analyze button
    
    
    output$download_summary_files <- downloadHandler(
        filename = function() {
            
            paste0(input$exp_num, "_", as.character(base::Sys.Date()), "_notes_and_next_steps.txt")
        },
        content = function(file) {
            writeLines(c(paste("Experiment", paste0(input$exp_num, input$protein), sep = ": "),
                         paste("Now what?", input$now_what, sep = ": "),
                         paste("Dye daughter plates", input$daughter_num, sep = ": "),
                         paste("Buffer used for screen", input$buffer_type, sep = ":"),
                         paste("Hit calling script version", input$script_version, sep = ": "),
                         paste("Instrument used", input$instrument, sep = ": "),
                         paste("Plate type", input$plate_type, sep = ": "),
                         paste("Plate lot", input$plate_lot, sep = ": "),
                         paste("Analysis completed on", as.character(base::Sys.Date()), sep = ": "),
                         paste("Additional notes", input$additional_notes, sep = ": "),
                         paste("Session info", capture.output(sessionInfo()))

            ), 
            file)
        }
    )
    
    
    output$download_example_layout <- downloadHandler(
        filename = "example_validation_layout.csv",
        content = function(fname) {
            write.csv( read_csv("example_validation_layout.csv"), fname, row.names = FALSE, quote = FALSE)
        }
    )
    
    
    
}


# Define UI for application that draws a histogram
ui <- navbarPage(useShinyalert(),
                 tabPanel("Instructions",
                          fluidRow(
                              # left panel: interactive sliders
                              column(3),
                              column(6,
                                     tags$h1("Instructions for validation plotting"), 
                                     tags$p("Sorry I haven't written full instructions yet. It works a lot like the hit-calling app, though . . ."),
                                     tags$br(),
                                     tags$br(),
                                     tags$p("In the 'Enter Experiment Information' tab:") %>% strong(),
                                     tags$p("1. Upload your raw qTower data file."),
                                     tags$p("2. Upload your validation layout (you can use the same one that you uploaded to the Echo instructions app.)"),
                                     tags$p("The layout has to have some specific names!") %>% strong(),
                                     tags$p("Layout must be named in the same way as it is for Echo uploads. That is:"),
                                     tags$p("Dyes are listed under the variable name 'compound'"),
                                     tags$p("Final dye concentration is listed under the variable name 'concentration'"),
                                     tags$p("Whether something is buffer or protein is listed under the variable name 'protein'"),
                                     tags$br(),
                                     tags$br(),
                                     tags$p("In the 'Examine results, write next steps' tab:") %>% strong(),
                                     tags$p("You can save notes on specific dyes using double-clicks and clicks on that panel of the plot. If you make a mistake and want to remove an assignment, click and drag anywhere in that panel. You can set your own designations for clicking and double-clicking; they default to the same as the hit-calling app: 'hit' and 'sensitive'."),
                                     downloadButton("download_example_layout", "Download example layout file", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff"),
                                     tags$br(),
                                     tags$br(),
                                     ### download button ####
                                     
                                     tags$p("Lmk if you need anything: 360 927 0262. Good luck!")
                              ),
                              
                              column(3)
                              
                              
                          )),
                 #"UCSF Dye Screen Processing",
                 tabPanel("Enter experiment information",
                          column(3,
                                 p("Basic screen information", style = "font-size: 16px;",align = "center"),
                                 textInput("protein", "Protein", "SP0XX"),
                                 textInput("exp_num", "Experiment number", "Ds00XX"),
                                 textInput("buffer", "Buffer used", "Buffer componets here"),
                                 textAreaInput("additional_notes", "Additional notes","Notes on the screen.", height = 100)
                          ),
                          
                          column(3,
                                 p("Technical notes", style = "font-size: 16px;",align = "center"),
                                 numericInput("start_T", "Starting temperature (C)", 25),
                                 numericInput("end_T", "Ending temperature (C)", 94),
                                 textInput("daughter_num", "Validation daughter number", "ENTER VALIDATION DAUGHTER NUMBER!!!"),
                                 textInput("script_version", "Analysis script version", "dye_calling_support_script.R"),
                                 textInput("instrument", "qPCR instrument used", "qTower384g quentin Towerntino"),
                                 textInput("plate_type", "Type of PCR plate used", "Axygen PCR-284-LC480WNFBC"),
                                 textInput("plate_lot", "PCR plate lot number", "lot 23517000")
                                 # ,
                                 # checkboxInput("save_outputs", "Save outputs", value = TRUE, width = NULL)
                          ),
                          
                          column(3,
                                 p("Screening files", style = "font-size: 16px;",align = "center"),
                                 #textAreaInput("path_in", "Path to data", "~/Box Sync/data/dye_screens/dye_screening_code/mock_screen/"),
                                 
                                 tags$hr(),
                                 p("Select raw data"),
                                 fileInput("protein_file", "Upload raw qTower data (csv)",
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 
                                 fileInput("protein_layout", "Upload validation layout (csv)",
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 
                                 
                                 
                                 actionButton("analyze_button", "Analyze")
                          )),
                 
                 tabPanel("Examine results, write next steps",
                          fluidRow(
                              # left panel: interactive sliders
                              column(9,
                                     wellPanel(id = "facet_plot",
                                               #plotOutput("plot1", brush = "plot_brush", height = "5000px"),
                                               p("Plot is ready when this grey box appears white. You'll probably have to scroll down in this window to see your results."),
                                               plotOutput("plot",
                                                          click = "plot_click",
                                                          dblclick = "plot_dblclick",
                                                          brush = "plot_brush",
                                                          height = "5000px"
                                               ) %>% withSpinner(color="#525252"),
                                               style = "overflow-y:scroll; max-height: 600px")
                              ),
                              column(3,
                                     textAreaInput("now_what", "Notes on next steps","What's the next step? Gets written to 'experimental notes'.", height = 150),
                                     textInput("double_click_meaning", "Double-click designation", "hit"),
                                     textInput("single_click_meaning", "Single-click designation", "sensitive"),
                                     
                                    
                                     DT::dataTableOutput("data_table"),
                                     downloadButton("download_files", "Download processed data (rds)"), 
                                     downloadButton("download_summary_files", "Download experimental notes"), 
                                     downloadButton("download_plot", "Download plot (takes a minute or so!)"),
                                     downloadButton("save_dt_button", "Save hit list to final folder")
                                     
                                     ))
                 )
                 
)

shinyApp(ui, server)
