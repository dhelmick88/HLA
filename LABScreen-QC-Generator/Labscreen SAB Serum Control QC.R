library(dplyr)
library(openxlsx)
library(ggplot2)
library(patchwork)
library(gridExtra)


#################################################
# SETTINGS, update Line 12
#################################################

input_folder <- "Select Input Folder Path"

output_file <- file.path(
  input_folder,
  "Luminex_Bead_Statistics_By_Control.xlsx"
)

graph_dir <- file.path(
  input_folder,
  "Luminex QC_Graphs"
)
pdf_dir <- file.path(
  input_folder,
  "Luminex QC_PDFs"
)

if(!dir.exists(pdf_dir)){
  dir.create(pdf_dir, recursive = TRUE)
}


if(!dir.exists(graph_dir)){
  dir.create(graph_dir, recursive = TRUE)
}

#################################################
# QC CONTROL DEFINITIONS
#################################################

qc_controls <- data.frame(
  Control = c(
    "INTCTL",
    "LS_NEG",
    "LS_POLY",
    "TRU_NEG"
  ),
  Pattern = c(
    "SERA CTR_LUM_INTCTL",
    "CTR_LUM_LS_NEG",
    "CTR_LUM_LS_POLY",
    "CTR_LUM_TRU_NEG"
  ),
  stringsAsFactors = FALSE
)

#################################################
# GET FILES
#################################################

files <- list.files(
  input_folder,
  pattern = "\\.csv$",
  full.names = TRUE
)

all_results <- data.frame()

#################################################
# PROCESS FILES
#################################################

for(file in files){
  
  cat("Reading:", basename(file), "\n")
  
  lines <- readLines(file, warn = FALSE)
  
  #################################################
  # EXTRACT PROTOCOL DESCRIPTION
  #################################################
  
  protocol_line <- grep(
    "ProtocolDescription",
    lines,
    value = TRUE,
    ignore.case = TRUE
  )
  
  if(length(protocol_line) > 0){
    
    protocol_description <- tryCatch({
      
      parts <- strsplit(protocol_line[1], ",")[[1]]
      
      if(length(parts) >= 2){
        trimws(gsub('"', "", parts[2]))
      } else {
        NA_character_
      }
      
    }, error = function(e) NA_character_)
    
  } else {
    
    protocol_description <- NA_character_
    
  }
  
  #################################################
  # DETERMINE HLA CLASS
  #################################################
  
  hla_class <- case_when(
    grepl("Class\\s*(1|I)\\b", protocol_description, ignore.case = TRUE) ~ "Class I",
    grepl("Class\\s*(2|II)\\b", protocol_description, ignore.case = TRUE) ~ "Class II",
    TRUE ~ "Unknown"
  )
  
  #################################################
  # RUN DATE
  #################################################
  
  run_date <- tryCatch({
    
    date_line <- strsplit(lines[3], ",")[[1]]
    date_text <- gsub('"', "", date_line[2])
    as.Date(date_text, format = "%m/%d/%Y")
    
  }, error = function(e) NA)
  
  #################################################
  # READ DATA TABLE
  #################################################
  
  dat <- tryCatch(
    
    read.csv(
      file,
      skip = 54,
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    
    error = function(e) NULL
    
  )
  
  if(is.null(dat)){
    cat("Skipping file because it could not be read:", basename(file), "\n")
    next
  }
  
  if(!"Sample" %in% names(dat)){
    cat("Skipping file because Sample column was not found:", basename(file), "\n")
    next
  }
  
  dat$Sample <- trimws(dat$Sample)
  
  bead_cols <- grep(
    "^\\d{3}$",
    names(dat),
    value = TRUE
  )
  
  if(length(bead_cols) == 0){
    cat("Skipping file because no bead columns were found:", basename(file), "\n")
    next
  }
  
  #################################################
  # FIND QC CONTROLS
  #################################################
  
  for(i in 1:nrow(qc_controls)){
    
    control_name <- qc_controls$Control[i]
    control_pattern <- qc_controls$Pattern[i]
    
    qc_rows <- grepl(
      control_pattern,
      dat$Sample,
      ignore.case = TRUE
    )
    
    qc_dat <- dat[qc_rows, ]
    
    if(nrow(qc_dat) == 0) next
    
    for(b in bead_cols){
      
      mfi_values <- suppressWarnings(
        as.numeric(
          gsub(",", "", qc_dat[[b]])
        )
      )
      
      temp <- data.frame(
        Date = run_date,
        RunID = tools::file_path_sans_ext(basename(file)),
        HLA_Class = hla_class,
        Control = control_name,
        Sample = qc_dat$Sample,
        Bead = b,
        MFI = mfi_values,
        File = basename(file),
        stringsAsFactors = FALSE
      )
      
      all_results <- rbind(
        all_results,
        temp
      )
    }
  }
}

#################################################
# STOP IF NO DATA FOUND
#################################################

if(nrow(all_results) == 0){
  stop("No QC data found. Check input folder, CSV files, QC sample names, and skip row.")
}

#################################################
# CREATE RUN SEQUENCE NUMBER
#################################################

run_lookup <- all_results %>%
  distinct(Date, RunID, File, HLA_Class) %>%
  arrange(Date, RunID, File) %>%
  mutate(
    Run_Number = row_number()
  )

all_results <- all_results %>%
  left_join(
    run_lookup %>% select(RunID, File, Run_Number),
    by = c("RunID", "File")
  )

#################################################
# SUMMARY STATS
#################################################

stats <- all_results %>%
  group_by(
    Control,
    Sample,
    Bead,
    HLA_Class
  ) %>%
  summarise(
    N = n(),
    Mean = mean(MFI, na.rm = TRUE),
    SD = sd(MFI, na.rm = TRUE),
    CV = round((SD / Mean) * 100, 2),
    Min = min(MFI, na.rm = TRUE),
    Max = max(MFI, na.rm = TRUE),
    Range = Max - Min,
    .groups = "drop"
  )

#################################################
# MEDIAN BY RUN
#################################################

median_data <- all_results %>%
  group_by(
    Control,
    HLA_Class,
    Bead,
    RunID,
    File,
    Run_Number
  ) %>%
  summarise(
    Median_MFI = median(MFI, na.rm = TRUE),
    .groups = "drop"
  )

#################################################
# EXPORT EXCEL
#################################################

wb <- createWorkbook()

for(control_name in unique(all_results$Control)){
  
  for(class_name in unique(all_results$HLA_Class)){
    
    raw_data <- all_results %>%
      filter(
        Control == control_name,
        HLA_Class == class_name
      ) %>%
      arrange(Run_Number, Sample, Bead)
    
    stat_data <- stats %>%
      filter(
        Control == control_name,
        HLA_Class == class_name
      ) %>%
      arrange(Sample, Bead)
    
    if(nrow(raw_data) == 0) next
    
    class_short <- ifelse(
      class_name == "Class I",
      "C1",
      ifelse(
        class_name == "Class II",
        "C2",
        "UNK"
      )
    )
    
    raw_sheet <- paste0(
      "Raw_",
      control_name,
      "_",
      class_short
    )
    
    stat_sheet <- paste0(
      "Stats_",
      control_name,
      "_",
      class_short
    )
    
    addWorksheet(wb, raw_sheet)
    writeData(wb, raw_sheet, raw_data)
    
    addWorksheet(wb, stat_sheet)
    writeData(wb, stat_sheet, stat_data)
    
  }
}

addWorksheet(wb, "Run_Index")
writeData(wb, "Run_Index", run_lookup)

saveWorkbook(
  wb,
  output_file,
  overwrite = TRUE
)

cat("Excel file saved:\n")
cat(output_file, "\n\n")

#################################################
# DIAGNOSTIC OUTPUT BEFORE GRAPHS
#################################################

cat("Rows in all_results:", nrow(all_results), "\n")
cat("Rows in median_data:", nrow(median_data), "\n")
cat("Graph directory:", graph_dir, "\n\n")

#################################################
# GENERATE PDF REPORTS
#################################################

for(control_name in unique(median_data$Control)){
  
  control_subset <- median_data %>%
    filter(Control == control_name)
  
  for(class_name in unique(control_subset$HLA_Class)){
    
    class_subset <- control_subset %>%
      filter(HLA_Class == class_name)
    
    plot_list <- list()
    
    for(bead_id in unique(class_subset$Bead)){
      
      bead_data <- class_subset %>%
        filter(Bead == bead_id) %>%
        arrange(Run_Number)
      
      if(nrow(bead_data) == 0) next
      
      mean_mfi <- mean(
        bead_data$Median_MFI,
        na.rm = TRUE
      )
      
      sd_mfi <- sd(
        bead_data$Median_MFI,
        na.rm = TRUE
      )
      
      if(is.na(sd_mfi)){
        sd_mfi <- 0
      }
      
      plot_color <- ifelse(
        class_name == "Class I",
        "blue",
        "red"
      )
      
      p <- ggplot(
        bead_data,
        aes(
          x = Run_Number,
          y = Median_MFI
        )
      ) +
        geom_line(
          color = plot_color,
          linewidth = 0.8
        ) +
        geom_point(
          color = plot_color,
          size = 1.5
        ) +
        
        geom_hline(
          yintercept = mean_mfi,
          color = "darkgreen"
        ) +
        
        geom_hline(
          yintercept = mean_mfi + (2 * sd_mfi),
          color = "red",
          linetype = "dashed"
        ) +
        
        geom_hline(
          yintercept = mean_mfi - (2 * sd_mfi),
          color = "red",
          linetype = "dashed"
        ) +
        
        labs(
          title = paste(
            class_name,
            "Bead",
            bead_id
          ),
          x = "Run Number",
          y = "Median MFI"
        ) +
        
        theme_bw() +
        
        theme(
          plot.title = element_text(
            size = 9
          ),
          axis.text = element_text(
            size = 7
          ),
          axis.title = element_text(
            size = 8
          )
        )
      
      plot_list[[length(plot_list) + 1]] <- p
    }
    
    pdf_file <- file.path(
      pdf_dir,
      paste0(
        control_name,
        "_",
        gsub(" ", "_", class_name),
        ".pdf"
      )
    )
    
    pdf(
      pdf_file,
      width = 11,
      height = 8.5
    )
    
    plots_per_page <- 6
    
    for(i in seq(
      1,
      length(plot_list),
      by = plots_per_page
    )){
      
      grid.arrange(
        grobs = plot_list[
          i:min(
            i + plots_per_page - 1,
            length(plot_list)
          )
        ],
        ncol = 2,
        nrow = 3
      )
      
    }
    
    dev.off()
    
    cat(
      "Created PDF:",
      basename(pdf_file),
      "\n"
    )
    
  }
}

cat(
  "\nFinished generating PDF QC reports.\n"
)

