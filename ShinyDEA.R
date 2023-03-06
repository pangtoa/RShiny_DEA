
library(DT)
library(plotly)
library(shiny)
library(readr)
library(shinydashboard)
library(shinymanager)
library(nonparaeff)
library(productivity)

#tfp_df <- read.csv("https://www.dropbox.com/s/3basp4lokyvpetg/global_ota_final_v2.csv?dl=1")
#mpi_df <- read.csv("https://www.dropbox.com/s/ymfhw6udwpudg3p/europe_tourism_data.csv?dl=1")

# credentials <- data.frame(
#   user = c("htm"),
#   password = c("htm123"), stringsAsFactors = FALSE)

header <- 
  dashboardHeader(
    title = "ShinyDEA",
    tags$li(a(href = 'https://www.shinyapps.io/',
              title = "DEA Analyzer"),
            class = "dropdown")
  )

sidebar <- 
  dashboardSidebar(
    sidebarMenu(
      menuItem("About ShinyDEA",
               tabName = "tab_about",
               icon = icon("info")),
      hr(),
      menuItem(
        "Static DEA",
        icon = icon("code"),
        menuItem(
          "BCC/CCR/SBM",
          tabName = "StaticDEA", 
          selected = FALSE
        ),
        fileInput("static_data", "Choose CSV File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        expandedName = "childfullMenuItem",
        startExpanded = FALSE
      ),
      
      hr(),
      
      menuItem(
        "Dynamic DEA",
        icon = icon("code"),
        menuItem(
          "Malmquist Index",
          tabName = "DynamicDEA", 
          selected = FALSE
        ),
        
        fileInput("dynamic_data", "Choose CSV File",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        expandedName = "childfullMenuItem",
        startExpanded = FALSE
      ),
      hr(),
      div(style="margin-left:10px; text-align:center","Created by Dohyung Bang", 
          br(), 
          "Email: dohyung.bang@gmail.com")

    )
  )

body <- 
  dashboardBody(
    tabItems(
      tabItem(tabName = "tab_about",
              fluidRow(
                box(width=12, 
                    h3("About ShinyDEA", align = "Left", style = "font-family: 'times'"),
                    div(style = "margin-bottom:-17em", imageOutput("image")),
                    h4("ShinyDEA is a", span("FREE DEA software", style = "color:blue"), "build by R Shiny.", style = "font-family: 'times'; font-si18pt"),
                    br(),
                    h4("This provides basic static models, such as", em("CCR, "),  em("BCC, "), "and ", em("slack-based model"), "(SBM),", style = "font-family: 'times'; font-si18pt"), 
                    h4("and also provides a dynamic model.", em("Malmquist Productivity Index (MPI)."), style = "font-family: 'times'; font-si18pt"),
                    br(),
                    h4("You can find example files below:", style = "font-family: 'times'; font-si18pt"),
                    h4("- ", tags$a(href="https://www.dropbox.com/s/ymfhw6udwpudg3p/europe_tourism_data.csv?dl=1", "Europe tourism productivity (for Dynamic DEA)")),
                    h4("- ", tags$a(href="https://www.dropbox.com/s/3basp4lokyvpetg/global_ota.csv?dl=1", "OTA total factor productivity (for Static DEA)")),
                    h4("- ", tags$a(href="https://www.dropbox.com/s/g72l13490paj25y/hotel_data.csv?dl=1", "Hotel productivity (for Static DEA)")),
                    br(),
                    h4("Extended models", em("Two-stage Network Model, "), em("General Network Model,"), "and", em("alternative Malmquist indexes"), ", will be added in the future.", 
                       style = "font-family: 'times'; font-si18pt"), 
                    br(),
                    h5("Contact: dohyung.bang@gmail.com", style = "font-family: 'times'")
                    #imageOutput("image"),
                    )

                
              )
      ),
      tabItem(tabName = "StaticDEA",
              navbarPage("Static DEA",
                         inverse = TRUE,
                         tabPanel(strong("Data Description"),
                                  fluidRow(
                                    box(width = 12,
                                      column(3,
                                             fluidRow(uiOutput("StaticDMU")),
                                             fluidRow(uiOutput("StaticOrientation"))
                                             ),
                                      column(3,uiOutput("StaticInputs")
                                             ),
                                      column(3,uiOutput("StaticOutputs")
                                             ),
                                      column(3, actionButton("StaticRun", "Run!")
                                      )
                                    )
                                    ),
                                    
                                    fluidRow(
                                      box(width = 12,
                                          DT::dataTableOutput(outputId = "StaticDataDescription"))
                                    )
                                           

                         ),
                         tabPanel(strong("CCR/BCC"),
                                  fluidRow(
                                    box(width = 12,
                                        DT::dataTableOutput(outputId = "CCRBCC_Result"))
                                  )),
                         tabPanel(strong("SBM"),
                                  fluidRow(
                                    box(width = 12,
                                        DT::dataTableOutput(outputId = "SBM_Result"))
                                  ))
                         
                        )
              
              ),
       tabItem(tabName = "DynamicDEA",
               navbarPage("Dynamic DEA",
                          inverse = TRUE,
                          tabPanel(strong("Data Description"),
                                   fluidRow(
                                     box(width = 12,
                                         column(3,
                                                fluidRow(uiOutput("DynamicDMU")),
                                                fluidRow(uiOutput("DynamicTimeID"))
                                         ),
                                         column(3,uiOutput("DynamicInputs")
                                         ),
                                         column(3,uiOutput("DynamicOutputs")
                                         ),
                                         column(3,
                                                fluidRow(uiOutput("DynamicOrientation")),
                                                fluidRow(uiOutput("DynamicRTS")),
                                                fluidRow(actionButton("DynamicRun", "Run!"))
                                         )
                                     )
                                   ),

                                   fluidRow(
                                     box(width = 12,
                                         DT::dataTableOutput(outputId = "DynamicDataDescription"))
                                   )


                          ),
                          tabPanel(strong("MPI"),
                                   fluidRow(
                                     box(width = 12,
                                         DT::dataTableOutput(outputId = "MPI_Result"))
                                   ))
                          
               )
               
       )
          
    
    )
  )

ui <- 
  dashboardPage(
    title = "ShinyDEA",
    header,
    sidebar,
    body,
    skin = "yellow")

# ui <- secure_app(ui)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # res_auth <- secure_server(check_credentials = check_credentials(credentials))
  # 
  # output$auth_output <- renderPrint({
  #   reactiveValuesToList(res_auth)
  # })
  
  output$image <- renderImage({
    filename <- normalizePath(file.path("./landing_page.png"))
    list(
      src = filename, 
      height = 150,
      margin = 0
    )
  }, deleteFile = FALSE)
  
  StaticData <- eventReactive(input$static_data, {
    
    input_file <- input$static_data
    df <- read.csv(input_file$datapath)

    return(df)
  })
  
  output$StaticDataDescription <- renderDataTable({
    
    df <- StaticData()
    df
  }, options = list(pageLength = 100))
  
  output$StaticDMU <- renderUI({
    df <- StaticData()
    selectInput("selected_static_dmu",
                "Select DMU ID Variable",
                choices = names(df))
  })
  
  output$StaticOrientation <- renderUI({
    df <- StaticData()
    selectInput("selected_static_orientation",
                "Select Orientation",
                choices = c("Input-oriented", "Output-oriented"))
  })
  
  output$StaticInputs <- renderUI({
    df <- StaticData()
    checkboxGroupInput("selected_static_inputs",
                       "Select Input Variables",
                       choices = names(df))
  })
  
  output$StaticOutputs <- renderUI({
    df <- StaticData()
    checkboxGroupInput("selected_static_outputs",
                       "Select Output Variables",
                       choices = names(df))
  })
  
  
  CCRBCC_RESULT <- eventReactive(input$StaticRun, {
    
    df <- 
      StaticData() %>% 
      rename(dmu = input$selected_static_dmu)
    
    n_output <- length(input$selected_static_outputs)
    dea_set <- 
      df %>% 
      select(input$selected_static_outputs, input$selected_static_inputs) %>% 
      as.data.frame()
    
    orientation <- ifelse(input$selected_static_orientation == "Input-oriented", 1, 2)
    
    te <- dea(dea_set, noutput = n_output, orientation = orientation, rts = 1, onlytheta = T) 
    pte <- dea(dea_set, noutput = n_output, orientation = orientation, rts = 2, onlytheta = T) 

    if (orientation == 2){
      te <- 1/te 
      pte <- 1/pte 
    }
    
    
    eff_crs <- dea(dea_set, noutput = 1, orientation = orientation, rts = 1)
    
    lambdas <- eff_crs[,2:(nrow(eff_crs)+1)]
    eff_crs$scale_eff <- apply(lambdas, 1, sum)
    
    eff_crs$scale_eff_cate <- ifelse(eff_crs$scale_eff > 1, "DRS",
                                     ifelse(1 - eff_crs$scale_eff < 0.000001, "CRS", "IRS"))
    eff_crs$scale_eff_cate <- ifelse(eff_crs$scale_eff > 1.00001, "DRS",
                                     ifelse(1 - eff_crs$scale_eff < 0.000001, "CRS", "IRS"))
    dmu <- df$dmu
    scale_eff <- eff_crs$scale_eff_cate
    
    #merge
    idx <- 1:length(dmu)
    dea_result <- cbind(idx, dmu, te, pte, scale_eff)
    names(dea_result) <- c("idx", "dmu", "eff_ccr", "eff_bcc", "scale_eff")
    dea_result <-
      dea_result %>% 
      mutate(eff_ccr = eff_ccr %>% round(digit =3),
             eff_bcc = eff_bcc %>% round(digit =3))
    
    benchmark <- cbind(dmu, eff_crs[,c(2:(nrow(eff_crs)+1))])
    
    # lambda 
    for (i in 1:nrow(eff_crs)){
      a <- benchmark[, i+1]
      a[a <= 0] <- NA
      a[a > 0] <- benchmark[i,1]
      benchmark[, i+1] <- a
    }
    
    # concatenate benchmark
    benchmark_lists <- benchmark[,2]
    
    for (i in 2:length(benchmark[,1])){
      benchmark_lists <- paste0(benchmark_lists, ", ", benchmark[,i+1])
      benchmark_lists <- gsub('NA, ', '', benchmark_lists)
    }
    
    benchmark_lists <- gsub(', NA', '', benchmark_lists)
    dea_result$benchmark <- benchmark_lists
    
    return(dea_result)
  })
  
  output$CCRBCC_Result <- renderDataTable({
    
    ccrbcc_result <- CCRBCC_RESULT()
    ccrbcc_result
    
    datatable(ccrbcc_result, extensions = 'Buttons', 
              options = list(paging = TRUE,
                             scrollX=TRUE, 
                             searching = TRUE,
                             ordering = TRUE,
                             dom = 'Bfrtip',
                             buttons = list(
                               list(extend = 'csv', 
                                    filename = paste0("ccr_bcc_result_", Sys.Date())
                                    )
                             ),
                             pageLength = 100000))
    
  })
  
  
  SBM_RESULT <- eventReactive(input$StaticRun, {
    
    df <- 
      StaticData() %>% 
      rename(dmu = input$selected_static_dmu)
    
    dea_inputs <- df %>% select(input$selected_static_inputs)
    dea_outputs <- df %>% select(input$selected_static_outputs)
    
    dmu <- df$dmu
    n_output <- length(input$selected_static_outputs)
    
    dea_set <- 
      df %>% 
      select(input$selected_static_outputs, input$selected_static_inputs) %>% 
      as.data.frame()
    
    sbm_result <- 
      sbm.vrs(base = dea_set, noutput = n_output)
    
    sbm_eff <- sbm_result$eff %>% round(digit = 3)
    sbm_input_slacks <- 
      sbm_result %>% 
      select(contains("x")) / dea_inputs
    
    names(sbm_input_slacks) <- paste0(names(dea_inputs), "_inefficiency(%)")
    
    sbm_output_slacks <- 
      sbm_result %>% 
      select(contains("y")) / dea_outputs
    
    names(sbm_output_slacks) <- paste0(names(dea_outputs), "_inefficiency(%)")
    
    idx <- 1:length(dmu)
    sbm_result_final <- cbind(idx, dmu, sbm_eff, 
                              sbm_input_slacks %>% round(digit = 3), 
                              sbm_output_slacks %>% round(digit = 3))
    
    benchmark <- cbind(dmu, sbm_result[,2:(nrow(dea_set)+1)])
    
    # lambda 
    for (i in 1:nrow(dea_set)){
      a <- benchmark[, i+1]
      a[a <= 0] <- NA
      a[a > 0] <- benchmark[i,1]
      benchmark[, i+1] <- a
    }
    
    # concatenate benchmark
    benchmark_list <- benchmark[,2]
    
    for (i in 2:length(benchmark[,1])){
      benchmark_list <- paste0(benchmark_list, ", ", benchmark[,i+1])
      benchmark_list <- gsub('NA, ', '', benchmark_list)
    }
    
    benchmark_list <- gsub(', NA', '', benchmark_list)
    sbm_result_final$benchmark <- benchmark_list
    
    return(sbm_result_final)
  })
  
  output$SBM_Result <- renderDataTable({
    
    sbm_result <- SBM_RESULT()
    sbm_result
    
    datatable(sbm_result, extensions = 'Buttons', 
              options = list(paging = TRUE,
                             scrollX=TRUE, 
                             searching = TRUE,
                             ordering = TRUE,
                             dom = 'Bfrtip',
                             buttons = list(
                               list(extend = 'csv', 
                                    filename = paste0("ccr_bcc_result_", Sys.Date())
                               )
                             ),
                             pageLength = 100000))
  })
  
  ########################################### Dynamic ####################################
  
  DynamicData <- eventReactive(input$dynamic_data, {
    
    input_file <- input$dynamic_data
    df <- read.csv(input_file$datapath)
    
    return(df)
  })
  
  output$DynamicDataDescription <- renderDataTable({
    
    df <- DynamicData()
    df
  }, options = list(pageLength = 100000))
  
  output$DynamicDMU <- renderUI({
    
    df <- DynamicData()
    selectInput("selected_dynamic_dmu",
                "Select DMU ID Variable",
                choices = names(df))
  })
  
  output$DynamicTimeID <- renderUI({
    
    df <- DynamicData()
    selectInput("selected_dynamic_time_id",
                "Select Time ID Variable",
                choices = names(df))
  })
  
  output$DynamicOrientation <- renderUI({
    
    df <- DynamicData()
    selectInput("selected_dynamic_orientation",
                "Select Orientation",
                choices = c("Input-oriented", "Output-oriented"))
  })
  
  output$DynamicRTS <- renderUI({
    
    df <- DynamicData()
    selectInput("selected_dynamic_rts",
                "Select Return-to-Scale",
                choices = c("Variable RTS (VRS)", "Constant RTS (CRS)"))
  })
  
  output$DynamicInputs <- renderUI({
    
    df <- DynamicData()
    checkboxGroupInput("selected_dynamic_inputs",
                       "Select Input Variables",
                       choices = names(df))
  })
  
  output$DynamicOutputs <- renderUI({
    
    df <- DynamicData()
    checkboxGroupInput("selected_dynamic_outputs",
                       "Select Output Variables",
                       choices = names(df))
  })
  
  MPI_RESULT <- eventReactive(input$DynamicRun, {
    
    df <- DynamicData()
    
    mpi_inputs <- df %>% select(input$selected_dynamic_inputs)
    mpi_outputs <- df %>% select(input$selected_dynamic_outputs)
    
    
    orientation <- ifelse(input$selected_dynamic_orientation == "Input-oriented", "in", "out")
    rts <- ifelse(input$selected_dynamic_rts == "Variable RTS (VRS)", "vrs", "crs")
    
    mpi_result <- 
      malm(data = df, 
           id.var = input$selected_dynamic_dmu, 
           time.var = input$selected_dynamic_time_id, 
           x.vars = names(mpi_inputs),
           y.vars = names(mpi_outputs),
           rts = rts, 
           orientation = orientation)
    
    # mpi_result <- 
    #   malm(data = mpi_df, 
    #        id.var = "Full_Name", 
    #        time.var = "year", 
    #        x.vars = c("TTCAPITAL" ,"TTLABOR"),
    #        y.vars = c("TTGDP_DIR", "ITA"),
    #        rts = "crs", 
    #        orientation = "out")
    
    mpi_result <- 
      mpi_result$Changes %>% 
      select(-obtech, -ibtech, -matech)
    
    if (rts == "vrs"){
      
      names(mpi_result) <- c("dmu", "year.to", "year.from", "mpi", "eff.ch", "tech.ch", "pure.eff.ch", "scale.eff.ch")
      mpi_result <- 
        mpi_result %>% 
        mutate(mpi = mpi %>% round(digit = 3),
               eff.ch = eff.ch %>% round(digit = 3),
               tech.ch = tech.ch %>% round(digit = 3),
               pure.eff.ch = pure.eff.ch %>% round(digit = 3),
               scale.eff.ch = scale.eff.ch %>% round(digit = 3))
      
    } else {
        
      names(mpi_result) <- c("dmu", "year.to", "year.from", "mpi", "eff.ch", "tech.ch")
      mpi_result <- 
        mpi_result %>% 
        mutate(mpi = mpi %>% round(digit = 3),
               eff.ch = eff.ch %>% round(digit = 3),
               tech.ch = tech.ch %>% round(digit = 3))
      
      }
    
    
    return(mpi_result)
    
    })
  
  output$MPI_Result <- renderDataTable({
    
    mpi_result <- MPI_RESULT()
    mpi_result
    
    datatable(mpi_result, extensions = 'Buttons', 
              options = list(paging = TRUE,
                             scrollX=TRUE, 
                             searching = TRUE,
                             ordering = TRUE,
                             dom = 'Bfrtip',
                             buttons = list(
                               list(extend = 'csv', 
                                    filename = paste0("mpi_result_", Sys.Date())
                               )
                             ),
                             pageLength = 1000000))
  })

}

# Run the application 
shinyApp(ui = ui, server = server)
