library(shiny)
library(shinyFiles)
library(stringi)
library(lubridate)
library(dplyr)
library(tidyr)
library(patientProfilesVis)
library(RColorBrewer)
library(ggplot2)
library(DT)
library(shinyjs)


ui <- fluidPage(
  useShinyjs(),
  titlePanel("Hand Grip Assessment Calibration Analysis"),
  
  # Add custom CSS for better visuals
  tags$style(HTML("
    .container-fluid {
      background-color: #f8f9fa;
      padding: 20px;
    }
    .panel-heading {
      background-color: #007bff;
      color: white;
      padding: 10px;
    }
    .btn-secondary {
      background-color: #007bff;
      border-color: #007bff;
    }
    .btn-secondary:hover {
      background-color: #0056b3;
      border-color: #0056b3;
    }
    .btn-primary {
      background-color: #007bff;
      border-color: #007bff;
    }
    .btn-primary:hover {
      background-color: #0056b3;
      border-color: #0056b3;
    }
    .reference-table td, .reference-table th {
      width: 50%;
      text-align: center;
    }
    .highlight-message-container {
      display: flex;
      justify-content: center;
      margin-top: 5px;
    }
    .highlight-message {
      background-color: #d9e6fe;
      padding: 10px 15px;
      font-size: 16px;
      text-align: left;
      border: 2px solid #cccccc;
      line-height: 1.6;
      word-wrap: break-word;
      font-family: Arial, sans-serif;
      max-width: 600px;
      width: 100%;
      box-shadow: 2px 2px 5px rgba(0, 0, 0, 0.1);
      margin: 0 auto; /* Center the box */
    }
    .highlight-message h3 {
      text-align: center;
      font-weight: bold;
      margin: 0 0 10px 0; /* Minimal spacing around header */
      font-size: 18px;
      color: #333333;
    }
    .highlight-message ul {
      padding-left: 20px;
      margin: 0; /* Remove extra margins */
      list-style-type: disc;
    }
    .highlight-message li {
      margin-bottom: 5px;
      color: #555555;
    }
    .sidebar-panel {
      background-color: #ffffff;
      padding: 15px;
      border-radius: 8px;
      box-shadow: 0px 0px 10px rgba(0, 0, 0, 0.1);
    }
    .wellPanel {
      padding: 5px;
      background-color: #f9f9f9;
      border-radius: 8px;
    }
    .shiny-input-container {
      margin-bottom: 5px;
    }
    .btn-primary {
      width: 100%;
      padding: 10px;
      font-size: 16px;
    }
    #logFolder {
      width: 100%;
      background-color: blue;
      color: white;
      border: none;
      padding: 10px 15px;
      border-radius: 4px;
      cursor: pointer;
    }
    #logFolder:hover {
      background-color: darkblue;
    }
    .select-input-container {
      margin-bottom: 10px;
    }
    .file-input-container {
      margin-bottom: 10px;
      font-size: 14px;
    }
     .file-input-container {
      margin-bottom: 10px;
      font-size: 14px;
    }
    .highlight-message-container {
      margin-bottom: 20px; /* Add space after the message box */
    }
    .reference-table-container {
      margin-top: 20px; /* Add space before the reference table */
      margin-bottom: 10px;
    }
    .button-container {
      display: flex;
      justify-content: space-between; /* Add space between buttons */
      gap: 10px; /* Optional: Add spacing between buttons */
      width: 100%; /* Ensure container spans full width */
    }
    .btn-primary, .btn-secondary {
      flex: 1; /* Each button takes up equal space */
      max-width: 100%; /* Ensure buttons can expand fully */
      padding: 10px;
      font-size: 16px;
      background-color: #007bff; /* Button background color */
      border-color: #007bff; /* Button border color */
      color: white; /* Make button text white */
    }
    
    .btn-primary:hover, .btn-secondary:hover {
      background-color: #0056b3; /* Darker blue on hover */
      border-color: #0056b3; /* Border changes on hover */
      color: white; /* Ensure text remains white on hover */
  }
  ")),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      # Highlight message container inside the wellPanel
      div(class = "highlight-message-container", 
          div(class = "highlight-message", uiOutput("referenceMessage"))
      ),
      
      # File Input: Assessment File
      div(class = "file-input-container", 
          fileInput("assessmentFile", "Upload Assessment CSV File", accept = ".csv")
      ),
      
      # File Input: Log Files
      div(class = "file-input-container", 
          fileInput("logFiles", "Upload Log Files" , multiple = FALSE)
      ),
      
      # Text Input: Log Folder Number
      div(class = "text-input-container",
          selectInput("log_folder",
                      label= "Select Log Folder Number",
                      choices = c("","001", "003", "005", "007", "008", "009", "010", "012", "013"),
                      selected = "")
      ),
      
      # Site Number Selection
      div(class = "select-input-container", 
          selectInput(
            inputId = "siteNumber",
            label = "Site Number",
            choices = c("","2301", "2504", "2103", "2503", "2304", "2202", "2505", "2701", "3001"),
            selected = "" # Ensures no option is selected by default
          )
      ),
      
      # Assesment Number Selection
      div(class = "select-input-container", 
          selectInput(
            inputId = "assessmentID",
            label = "Assessment Type",
            choices = c("Hand Grip Relaxation Time","vHOT Trial 2","vHOT Trial 1","9-HPT","Stairs Ascend/Descend","5xSTS","10-MWRT","QMT","MMT"),
            selected = "Hand Grip Relaxation Time" # Ensures no option is selected by default
          )
      ),
      
      # Process Button
      div(class = "shiny-input-container", 
          actionButton("process", "Process Data", class = "btn-primary")
      ),
      div(class = "reference-table-container",
          DT::dataTableOutput("referenceTable")),
    ),
    
    mainPanel(
      width = 9, # Adjust this width to fill the remaining space
      tabsetPanel(
        tabPanel("Calibration Plot",
                 div(class = "plot-container", plotOutput("calibrationPlot", height = "600px")),
                 uiOutput("downloadButtonUI", class = "button-container")
        )
      )
    )
  )
)



server <- function(input, output, session) {
  
  reference_data <- data.frame(
    Log_Folder = c("001", "003", "005", "007", "008", "009", "010", "012", "013"),
    Site = c("3001", "2701", "2505", "2202", "2304", "2503", "2103", "2504", "2301")
  )
  
  # Reactive value to control reference table
  #reference_table_data <- reactiveVal(reference_data)
  
  # Render the Reference List table
  output$referenceTable <- DT::renderDataTable({
    #req(reference_table_data())
    datatable(reference_data, options = list(pageLength = 10))
  })
  
  output$referenceMessage <- renderUI({
    HTML("
            <h3>Instructions</h3>
            <ul>
              <li>Please select the combination of Site and files from the Log Folder as per the reference list given below.</li>
              <li>The <strong>Calibration Plot Tab</strong> will show the plots on the right once the process is completed.</li>
            </ul>
          ")
  })
  
  site_specific_data <- reactiveVal(NULL)
  
  # Use eventReactive to trigger processing when the actionButton is clicked
  process_data <- eventReactive(input$process, {
    
    # Reference table data
    print("Process button clicked")  # Add a print statement here
    start_time <- Sys.time()
    print(paste("Start Time:", start_time))  # Log the start time
    
    withProgress(message = 'Processing file...', value = 0, {
      req(input$assessmentFile, input$logFiles, input$siteNumber)
      
      # Define the folder-to-site mapping
      file_site_mapping <- data.frame(
        Log_Folder = c("001", "003", "005", "007", "008", "009", "010", "012", "013"),
        Site = c("3001", "2701", "2505", "2202", "2304", "2503", "2103", "2504", "2301")
      )
      
      log_file_zip <- input$logFiles$datapath  # Path of the uploaded zip file
      zip_file_dir <- dirname(log_file_zip)  # Get the directory of the uploaded file
      
      # Unzip the ZIP file to the extraction directory
      unzip_start_time <- Sys.time()
      print(paste("Unzip Start Time:", unzip_start_time))
      unzip(log_file_zip, exdir = zip_file_dir)
      unzip_end_time <- Sys.time()
      print(paste("Unzip End Time:", unzip_end_time))  # Log end time for unzip
      print(paste("Unzip Duration:", unzip_end_time - unzip_start_time))  # Calculate duration
      
      # List the extracted files and directories
      extracted_files <- list.files(zip_file_dir, full.names = TRUE)
      
      # Extract the first directory inside the zip (assumed to be the subfolder with .txt files)
      subfolder_name <- basename(extracted_files[2])
      #print(subfolder_name)
      
      # Check if the first item is a folder (not a file)
      if (dir.exists(file.path(zip_file_dir, subfolder_name))) {
        # If it’s a folder, proceed to list the .txt files inside
        subfolder_path <- file.path(zip_file_dir, subfolder_name)
        log_files <- list.files(subfolder_path, pattern = "\\.txt$", full.names = TRUE)
        #print("Log files found in subfolder:")
        #print(log_files)
      } else {
        # If no subfolder, look for .txt files in the main directory
        log_files <- list.files(zip_file_dir, pattern = "\\.txt$", full.names = TRUE)
        #print("Log files found in main directory:")
        #print(log_files)
      }
      
      # If no .txt files are found, stop with an error
      if (length(log_files) == 0) {
        stop("No .txt files found in the ZIP archive.")
      }
      
      file_names <- input$logFiles$name  # Get the names of the uploaded files
      
      site_number <- input$siteNumber
      log_folder <- input$log_folder  # User input for log folder number
      
      # Define the folder-to-site mapping again (redundant)
      file_site_mapping <- data.frame(
        Log_Folder = c("001", "003", "005", "007", "008", "009", "010", "012", "013"),
        Site = c("3001", "2701", "2505", "2202", "2304", "2503", "2103", "2504", "2301")
      )
      
      # Check if the folder number corresponds to the selected site number
      expected_site_number <- file_site_mapping$Site[file_site_mapping$Log_Folder == log_folder]
      
      if (length(expected_site_number) == 0 || expected_site_number != site_number) {
        showModal(modalDialog(
          title = "Error",
          paste("Please ensure that the selected Log Folder number ('", log_folder, "') corresponds to the selected Site Number (", site_number, ").", sep = ""),
          easyClose = TRUE,
          footer = NULL
        ))
        return(NULL)
      }
      
      assessment_data <- read.csv(input$assessmentFile$datapath, sep = ",", header = TRUE)
      assessment_data$Subject <- substr(trimws(assessment_data$Subject.ID), nchar(trimws(assessment_data$Subject.ID)) - 3, nchar(trimws(assessment_data$Subject.ID)))
      col <- ncol(assessment_data)
      
      #print("PRINTING NMBER COLLL")
      #print(col)
      #print(head(assessment_data))
      #print(colnames(assessment_data))
      
      incProgress(0.1, detail = "Reading the file")
      
      extract_data <- function(file) {
        content <- readLines(file)
        content <- paste(content, collapse = "\n")
        split_content <- stri_split_regex(content, "\\nHIDNAME:", simplify = TRUE)
        
        data_between_dates <- vector("list", length = ncol(split_content) - 1)
        for (i in 1:(ncol(split_content) - 1)) {
          data_between_dates[[i]] <- split_content[1, i + 1]
          data_between_dates[[i]] <- paste(data_between_dates[[i]], basename(file), sep = "\n")
        }
        return(data_between_dates)
      }
      
      # Apply the function to each selected log file
      extract_data_start_time <- Sys.time()
      all_data <- lapply(log_files, extract_data)
      extract_data_end_time <- Sys.time()
      print(paste("Data extraction Duration:", extract_data_end_time - extract_data_start_time))
      #print(all_data)
      
      parse_dates_from_data <- function(data) {
        date_string <- sub(".*DATE:\\s*([A-Za-z]+, [A-Za-z]+ \\d{1,2}, \\d{4}).*", "\\1", data)
        parsed_date <- mdy(date_string)
        formatted_date <- format(parsed_date, "%m/%d/%y")
        return(formatted_date)
      }
      
      find_closest_data <- function(assessment_date, all_data) {
        closest <- NULL
        min_diff <- -1000
        
        for (file_data in all_data) {
          for (data_chunk in file_data) {
            dates <- parse_dates_from_data(data_chunk)
            dates_un <- unique(dates)
            
            for (date in dates) {
              diff <- difftime(as.Date(date,format = "%m/%d/%y"), as.Date(assessment_dates,format = "%m/%d/%y"), units = "days")
              diff_find <- Sys.time()
              
              days <- as.numeric(gsub("Time difference of ([0-9]+) days", "\\1", diff))
              days_time <- Sys.time()
              print(paste("Day calculation Duration:", diff_find - days_time))
              if (days <= 0 && !is.na(days)) {
                if (days > min_diff && grepl("Hand Dyno 100", data_chunk)) {
                  min_diff <- days
                  closest <- data_chunk
                  closest <- paste(closest, days, sep = "\n")
                }
              }
            }
          }
        }
        
        return(closest)
      }
      
      sortDays <- function(days) {
        convert_to_sortable <- function(day) {
          if (grepl("Screening", day)) {
            return(as.numeric(gsub("Screening ", "", day)) - 1) 
          } else if (day == "Unscheduled") {
            return(Inf)
          } else {
            return(as.numeric(gsub("Day ", "", day)))
          }
        }
        sortable_days <- sapply(days, convert_to_sortable)
        sorted_days <- days[order(sortable_days)]
        return(sorted_days)
      }
      
      incProgress(0.5, detail = "Sorting the file")
      
      # Populate assessment_data with closest calibration data
      for (i in 1:nrow(assessment_data)) {
        assessment_dates <- as.Date(assessment_data[i,]$Date.Assessment.was.conducted)
        if (grepl(input$assessmentID, assessment_data[i,]$Assessment) && grepl(site_number, assessment_data[i,]$Site.ID)) { 
          closest_data_for_assessment <- find_closest_data(assessment_dates, all_data)
          
          if (!is.null(closest_data_for_assessment)) {
            str_with_tabs <- gsub("\t", "", closest_data_for_assessment)
            
            # Split the string into lines
            split_str <- strsplit(str_with_tabs, "\n")[[1]]
            
            assessment_data[i, col+1] <- paste(split_str, collapse = "\\t")
            assessment_data[i, col+2] <- split_str[length(split_str)]
            
            colnames(assessment_data)[col + 1] <- "split_structure_col"
            colnames(assessment_data)[col + 2] <- "Calculated_Days"
          }
        }
      }
      
      upload_dir <- dirname(input$assessmentFile$datapath) 
      output_file_name <- paste0(input$siteNumber, "_", input$log_folder, "_output.csv")
      output_file <- file.path(zip_file_dir, output_file_name)
      
      site_specific <- assessment_data[grepl(site_number, assessment_data$Site.ID),]
      site_specific <- site_specific[grepl("Hand Grip Relaxation Time", site_specific$Assessment),]
      site_specific$Calculated_Days <- as.numeric(site_specific$Calculated_Days)
      site_specific$calibration <- round(abs(site_specific$Calculated_Days))
      
      # Check if the file already exists, and remove it if necessary
      if (file.exists(output_file)) {
        file.remove(output_file)
      }
      
      write.csv(site_specific, output_file, row.names = FALSE)
      
      end_time <- Sys.time()
      print(paste("End Time:", end_time))  # Log the end time
      print(paste("Total Duration:", end_time - start_time))  # Total duration of the process
      
      return(site_specific)
    })
  })
  
  # Render Summary Table
  output$summaryTable <- renderTable({
    req(process_data())
    process_data()
  })
  
  lbPlotsColorShape <- reactive({
    req(process_data())
    site_specific <- process_data()
    
    site_specific <- site_specific %>%
      mutate(
        calibration = round(calibration), 
        color = case_when(
          calibration == 0 ~ "green",
          calibration == 1 ~ "blue",
          calibration == 2 ~ "orange",
          calibration == 3 ~ "yellow",
          calibration == 4 ~ "purple",
          calibration == 5 ~ "brown",
          calibration == 6 ~ "pink",
          calibration == 7 ~ "red",
          TRUE ~ "default_color"
        )
      )
    
    list_of_color <- unique(unlist(site_specific$color))
    order_color <- c("green", "blue", "orange", "yellow", "purple", "brown", "pink", "red")
    filtered_order <- order_color[order_color %in% list_of_color]
    
    subjectProfileEventPlot(
      data = site_specific,
      paramVar = "Subject",
      timeVar = "Visit",
      colorLab = "Calibration vs. Assessment (Days)",
      colorVar = "calibration",
      colorPalette = filtered_order,
      subjectVar = "Assessment",
      shapePalette = rep("square", 7),
      title = "Assessment Calibration Indicator",
      yLab = "Subject"
    )
  })
  
  # Render the plot in UI
  output$calibrationPlot <- renderPlot({
    req(lbPlotsColorShape())
    lbPlotsColorShape()
  })
  
  # Render download button UI conditionally after plot is shown
  output$downloadButtonUI <- renderUI({
    req(lbPlotsColorShape()) # Ensure process button is clicked before showing the buttons
    
    div(
    class = "button-container",
    tagList(
      downloadButton("downloadPlot", "Download Plot", class = "btn-primary"),
      downloadButton("downloadData", "Download Data", class = "btn-secondary")
    )
    )
  })
  
  # Download the plot as PNG
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0("calibration_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 10, height = 6)  # Open a PDF device
      print(lbPlotsColorShape())        # Use the reactive plot object
      dev.off()                         # Close the PDF device
    }
  )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      # Use user inputs to dynamically generate the file name
      site <- input$siteNumber
      assessment <- input$assessmentID
      paste0("calibration_data_", site, "_", assessment, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      # Write the processed data to the file
      write.csv(site_specific_data(), file, row.names = FALSE)
    }
  )
  
  
}



shinyApp(ui, server)

