#```R
library(shiny)
library(glmnet)
library(data.table)
library(sva) # for ComBat
library(preprocessCore) # for quantile normalization

# Increase maximum upload size to 30MB
options(shiny.maxRequestSize = 40 * 1024^2)

# UI Definition
ui <- fluidPage(
  titlePanel("ImmunoResponse Predictor"),
  sidebarLayout(
    sidebarPanel(
      br(), br(), fileInput("testFile", "Upload Test.csv (For mUC or mRCC)", accept = c(".csv")),
      br(),
      selectInput("modelFile", "Select trained model", 
                  choices = c("select model" = "",
                              "mUC model" = "logistic-Model-train-muc-test-muc.rds", 
                              "mRCC model" = "logistic-Model-train-rcc-test-rcc.rds")),
      br(), br(),
      actionButton("predictButton", "Generate predictions"),
      br(), br(), br(), br(), br(),
      uiOutput("downloadUI")
    ),
    mainPanel(
      h4("Instructions:"),
      p("1. Upload your Test.csv file."),
      p("2. Select a pre-trained model from the dropdown."),
      p("3. Click 'Generate Predictions' to process the data."),
      p("4. Download the predictions as a CSV file."),
      br(),
      textOutput("status"),
      textOutput("predictionsCount")
    )
  )
)

# Server Definition
server <- function(input, output, session) {
  predictions <- reactiveVal(NULL)
  
  testData <- reactive({
    req(input$testFile)
    expr.test0 <- fread(input$testFile$datapath, header = TRUE, sep = ",")
    sampleIDs <- expr.test0[[1]]  # First column as sample IDs
    expr.test0 <- expr.test0[, -1, with = FALSE]  # Remove sample ID column
    expr.test0_matrix <- as.matrix(expr.test0)
    rownames(expr.test0_matrix) <- sampleIDs
    return(list(data = expr.test0_matrix, sampleIDs = sampleIDs))
  })
  
  # Function to compute average cosine distances for a single test sample
  compute_average_cosine_distances <- function(y1, y2, x) {
    cosine_distance <- function(a, b) {
      norm_a <- sqrt(sum(a^2))
      norm_b <- sqrt(sum(b^2))
      if (norm_a == 0 || norm_b == 0) return(NA)
      1 - sum(a * b) / (norm_a * norm_b)
    }
    
    avg_distance_R <- if (nrow(y1) > 0) mean(apply(y1, 1, function(sample) cosine_distance(x, sample)), na.rm = TRUE) else NA
    avg_distance_NR <- if (nrow(y2) > 0) mean(apply(y2, 1, function(sample) cosine_distance(x, sample)), na.rm = TRUE) else NA
    
    return(c(avg_distance_R, avg_distance_NR))
  }
  
  # Function to apply ORR-based prior and compute matching statistics
  apply_orr_prior <- function(results_df, orr = 0.2) {
    results_df <- results_df[order(results_df$CosDist_2_Rs), ]
    total_samples <- nrow(results_df)
    x <- round(total_samples * orr)
    results_df$CosineDist_prior <- c(rep("R", x), rep("NR", total_samples - x))
    
    rs_matching <- sum(results_df$LogitDA_pred_label[1:x] == results_df$CosineDist_prior[1:x])
    nrs_matching <- sum(results_df$LogitDA_pred_label[(x + 1):total_samples] == 
                          results_df$CosineDist_prior[(x + 1):total_samples])
    total_matching <- rs_matching + nrs_matching
    prior_final_percentage <- round((total_matching / total_samples) * 100)

    list(
      results_df = results_df,
      prior_final_percentage = prior_final_percentage
    )
  }
  
  # Function to save results and display percentage
  save_and_report_results <- function(results_df, prior_final_percentage) {
    results_df$`% supporting by CosineDist_prior` <- ""
    results_df$`% supporting by CosineDist_prior`[1] <- prior_final_percentage
    return(results_df)
  }
  
  observeEvent(input$predictButton, {
    req(input$testFile, input$modelFile)
    tryCatch({
      withProgress(message = "Generating Predictions", value = 0.1, {
        model_path <- file.path("models", input$modelFile)
        if (!file.exists(model_path)) stop("Model file not found: ", model_path)
        bestModel <- readRDS(model_path)
        
        test_data <- testData()
        expr.test0_matrix <- test_data$data
        sampleIDs <- test_data$sampleIDs
        
        gene_ids_clean <- sub("^X", "", bestModel$beta@Dimnames[[1]])
        test_colnames_clean <- sub("^X", "", colnames(expr.test0_matrix))
        test_cols <- match(gene_ids_clean, test_colnames_clean)
        
        if (any(is.na(test_cols))) {
          missing_genes <- gene_ids_clean[is.na(test_cols)]
          stop("Missing genes in test data: ", paste(missing_genes, collapse = ", "))
        }
        
        if (input$modelFile == "logistic-Model-train-muc-test-muc.rds") {
          # --- mUC Model Block ---
          if (length(test_cols) != 49) {
            stop("mUC model expects 49 features, but test data has ", length(test_cols), " features.")
          }
          
          train_zip_path <- file.path("data/mUC_log2TPMp1_train.zip")
          if (!file.exists(train_zip_path)) stop("Training zip file not found: ", train_zip_path)
          train_data_raw_mUC <- loadTable(file = unz(train_zip_path, "log2TPMp1_train.csv"), 
                                          transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE)
          if (!is.matrix(train_data_raw_mUC)) {
            stop("mUC training data is not a matrix")
          }
          
          sampleIDs_train <- rownames(train_data_raw_mUC)
          train_data_raw_mUC <- scale(train_data_raw_mUC)
          
          train_data <- train_data_raw_mUC[, test_cols, drop = FALSE]
          x.test <- scale(expr.test0_matrix[, test_cols, drop = FALSE])

          
          annot_zip_path <- file.path("data/mUC_response_train.zip")
          if (!file.exists(annot_zip_path)) stop("Annotation zip file not found: ", annot_zip_path)
          sampleAnnot.train <- read.csv(unz(annot_zip_path, "response_train.csv"))
          
          # Use RNASEQ_SAMPLE_ID for mUC (assuming same structure)
          if (!"RNASEQ_SAMPLE_ID" %in% colnames(sampleAnnot.train)) {
            stop("No RNASEQ_SAMPLE_ID column in mUC annotations")
          }
          
          if (!any(sampleIDs_train %in% sampleAnnot.train$RNASEQ_SAMPLE_ID)) {

            stop("No common sample IDs between mUC training data and annotations")
          }
          
          common_samples <- intersect(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID)
          if (length(common_samples) == 0) {
            stop("No common sample IDs between mUC training data and annotations")
          }
          
          # Subset to common samples
          sample_idx <- which(sampleIDs_train %in% common_samples)
          sampleIDs_train <- sampleIDs_train[sample_idx]
          train_data <- train_data[sampleIDs_train, , drop = FALSE]
          sampleAnnot.train <- sampleAnnot.train[match(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID), ]
  
          train_data_ordered <- train_data[match(sampleAnnot.train$RNASEQ_SAMPLE_ID, rownames(train_data)), , drop = FALSE]
          y1 <- train_data_ordered[sampleAnnot.train$RESPONSE == 1 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
          y2 <- train_data_ordered[sampleAnnot.train$RESPONSE == 0 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]

          
        } else if (input$modelFile == "logistic-Model-train-rcc-test-rcc.rds") {
          # --- RCC Model Block ---
          if (length(test_cols) != 27) {
            stop("RCC model expects 27 features, but test data has ", length(test_cols), " features.")
          }
          
          train_zip_path <- file.path("data/RCC_log2TPMp1_train.zip")
          if (!file.exists(train_zip_path)) stop("Training zip file not found: ", train_zip_path)
          train_data_raw_mRCC <- loadTable(file = unz(train_zip_path, "log2TPMp1_train.csv"), 
                                           transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE)
          if (!is.matrix(train_data_raw_mRCC)) {
            stop("RCC training data is not a matrix")
          }
          
          sampleIDs_train <- rownames(train_data_raw_mRCC)
          x.test_full <- expr.test0_matrix

          
          # ComBat + Quantile Normalization on full gene set
          expr.train0.t <- t(train_data_raw_mRCC)
          expr.test0.t <- t(x.test_full)
          
          train_len <- ncol(expr.train0.t)
          test_len <- ncol(expr.test0.t)
          bat <- c(rep("train", train_len), rep("test", test_len))
          
          data_all <- cbind(expr.train0.t, expr.test0.t)
          
          combat_data <- tryCatch({
            ComBat(dat = t(t(data_all)), batch = bat, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
          }, error = function(e) {
            stop("ComBat failed: ", e$message)
          })
          
          expr.all <- normalize.quantiles(combat_data, copy = TRUE)
          colnames(expr.all) <- colnames(combat_data)
          rownames(expr.all) <- gsub("^X[-|_]", "", rownames(combat_data))
          
          expr.all <- t(expr.all)
          
          # Validate gene IDs for subsetting
          expr_all_colnames_clean <- gsub("^X[-|_]", "", colnames(expr.all))
          test_cols_updated <- match(gene_ids_clean, expr_all_colnames_clean)
          if (any(is.na(test_cols_updated))) {
            stop("Gene IDs not found in expr.all: ", paste(gene_ids_clean[is.na(test_cols_updated)], collapse = ", "))
          }
          
          # Subset to model genes after normalization
          train_data <- expr.all[rownames(expr.all) %in% sampleIDs_train, test_cols_updated, drop = FALSE]
          x.test <- expr.all[rownames(expr.all) %in% sampleIDs, test_cols_updated, drop = FALSE]
          
          if (nrow(train_data) == 0 || nrow(x.test) == 0) {
            stop("Empty train_data or x.test after subsetting. Check sample IDs.")
          }

          stopifnot(all(colnames(x.test) == colnames(train_data)))

          
          # Check annotation file path
          annot_zip_path <- file.path("data/RCC_response_train.zip")
          if (!file.exists(annot_zip_path)) stop("Annotation zip file not found: ", annot_zip_path)
          sampleAnnot.train <- read.csv(unz(annot_zip_path, "response_train.csv"))

          
          # Use RNASEQ_SAMPLE_ID
          if (!"RNASEQ_SAMPLE_ID" %in% colnames(sampleAnnot.train)) {
            stop("No RNASEQ_SAMPLE_ID column in RCC annotations")
          }
          
          # Log mismatches
          if (!any(sampleIDs_train %in% sampleAnnot.train$RNASEQ_SAMPLE_ID)) {
            stop("No common sample IDs between RCC training data and annotations")
          }
          
          common_samples <- intersect(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID)
          if (length(common_samples) == 0) {
            stop("No common sample IDs between RCC training data and annotations")
          }
          
          # Subset to common samples
          sample_idx <- which(sampleIDs_train %in% common_samples)
          sampleIDs_train <- sampleIDs_train[sample_idx]
          train_data <- train_data[sampleIDs_train, , drop = FALSE]
          sampleAnnot.train <- sampleAnnot.train[match(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID), ]

          
          train_data_ordered <- train_data[match(sampleAnnot.train$RNASEQ_SAMPLE_ID, rownames(train_data)), , drop = FALSE]
          y1 <- train_data_ordered[sampleAnnot.train$RESPONSE == 1 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
          y2 <- train_data_ordered[sampleAnnot.train$RESPONSE == 0 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
        }
        
        # Manual prediction
        x.test <- as.matrix(x.test)
        if (any(is.na(x.test))) stop("NA values detected in x.test before prediction")

        beta <- as.matrix(bestModel$beta)
        a0 <- bestModel$a0
        if (ncol(x.test) != nrow(beta)) {
          stop("Dimension mismatch: x.test has ", ncol(x.test), " columns, beta has ", nrow(beta), " rows")
        }
        
        pred_prob <- 1 / (1 + exp(-(x.test %*% beta + a0)))
        pred_prob <- as.vector(pred_prob)
        pred_class <- ifelse(pred_prob > 0.5, 1, 0)
        pred_labels <- ifelse(pred_class == 1, "R", "NR")
        
        results_list <- list()
        for (i in 1:nrow(x.test)) {
          x_sample <- x.test[i, , drop = FALSE]
          sample_id <- rownames(x.test)[i]
          
          avg_distances <- compute_average_cosine_distances(y1, y2, x_sample)
          
          results_list[[i]] <- data.frame(
            CosDist_2_Rs = avg_distances[1],
            CosDist_2_NRs = avg_distances[2],
            LogitDA_pred_label = pred_labels[i]
          )
        }
        
        results_df <- do.call(rbind, results_list)
        results_df$sampleID <- sampleIDs
        results_df <- results_df[match(sampleIDs, results_df$sampleID), ]
        results_df <- results_df[, c("sampleID", setdiff(names(results_df), "sampleID"))]
        
        stopifnot(all(results_df$sampleID == sampleIDs))
        
        orr_results <- apply_orr_prior(results_df, orr = 0.2)
        results_df <- save_and_report_results(
          results_df = orr_results$results_df,
          prior_final_percentage = orr_results$prior_final_percentage
        )

        
        predictions(results_df)
        output$downloadUI <- renderUI({ downloadButton("downloadPredictions", "Download predictions") })
        output$status <- renderText(sprintf("Predictions generated successfully! %% supporting by CosineDist_prior: %.2f%%", orr_results$prior_final_percentage))
      })
    }, error = function(e) {
      output$status <- renderText(paste("Error generating predictions:", e$message))
    })
  })
  
  output$predictionsCount <- renderText({
    preds <- predictions()
    if (is.null(preds)) {
      return("Rows in predictions: 0")
    }
    paste("Rows in predictions:", nrow(preds))
  })
  
  output$downloadPredictions <- downloadHandler(
    filename = function() {
      paste("predictions-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      preds <- predictions()
      if (is.null(preds)) {
        stop("No predictions available to download. Please generate predictions first.")
      }
      fwrite(preds, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  output$status <- renderText("")
}

# Run the app
shinyApp(ui = ui, server = server)