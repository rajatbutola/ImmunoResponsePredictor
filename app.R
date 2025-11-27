#```R

library(shiny)
library(glmnet)
library(data.table)

## --- Bioconductor annotation libs for ID mapping ---
library(org.Hs.eg.db)
library(AnnotationDbi)

options(shiny.maxRequestSize = 40 * 1024^2)

###############################################################
# Helpers
###############################################################

loadTable <- function(file, transpose = FALSE, convertToMatrix = TRUE,
                      sep = ",", header = TRUE) {
  data <- read.csv(file, sep = sep, header = header,
                   row.names = 1, check.names = FALSE)
  if (transpose) data <- t(data)
  if (convertToMatrix) data <- as.matrix(data)
  return(data)
}

# last-resort NA scrubber
scrub_na <- function(M, fill = 0) {
  M[is.na(M)] <- fill
  M
}

# Generic cleaning of gene IDs:
# - trim, strip transcript version, X1234 -> 1234, uppercase
clean_gene_ids <- function(ids) {
  ids2 <- trimws(ids)
  ids2 <- sub("\\.\\d+$", "", ids2)          # strip .version at end
  ids2 <- sub("^X([0-9]+)$", "\\1", ids2)    # X1234 -> 1234
  toupper(ids2)
}

###############################################################
# Build synonym-to-model mapping tables (for renaming cols)
###############################################################

# Map everything → mUC Entrez panel (X<entrez>)
build_mUC_mapping <- function(mUC_entrez, mUC_model_ids) {
  keys <- as.character(mUC_entrez)  # Entrez IDs as character
  
  sym   <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = keys,
    column = "SYMBOL", keytype = "ENTREZID",
    multiVals = "first"
  )
  ensg  <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = keys,
    column = "ENSEMBL", keytype = "ENTREZID",
    multiVals = "first"
  )
  alias_list <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = keys,
    column = "ALIAS", keytype = "ENTREZID",
    multiVals = "list"
  )
  
  syn_list <- vector("list", length(keys))
  names(syn_list) <- mUC_model_ids
  
  for (i in seq_along(keys)) {
    ent <- keys[i]
    mid <- mUC_model_ids[i]   # e.g. "X9700"
    
    s <- c(
      mid,           # exact model ID (X9700)
      ent            # plain Entrez ("9700")
    )
    if (!is.na(sym[i]))  s <- c(s, sym[i])
    if (!is.na(ensg[i])) s <- c(s, ensg[i])
    
    if (!is.null(alias_list[[i]]) && length(alias_list[[i]]) > 0) {
      s <- c(s, alias_list[[i]])
    }
    
    syn_list[[i]] <- unique(clean_gene_ids(s))
  }
  
  all_syn <- unlist(syn_list, use.names = FALSE)
  all_mid <- rep(names(syn_list), times = lengths(syn_list))
  keep    <- !duplicated(all_syn)
  
  synonym_to_model <- all_mid[keep]
  names(synonym_to_model) <- all_syn[keep]
  
  list(
    synonym_to_model = synonym_to_model,   # cleaned synonym -> "X9700"
    model_ids        = mUC_model_ids
  )
}

# Map everything → mRCC ENST panel (ENSTxxx.y) using org.Hs.eg.db
build_mRCC_mapping <- function(mRCC_model_ids) {
  base <- sub("\\.\\d+$", "", mRCC_model_ids)  
  
  # Use ENSEMBLTRANS to go transcript → gene
  sym <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = base,
    column = "SYMBOL", keytype = "ENSEMBLTRANS",
    multiVals = "first"
  )
  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = base,
    column = "ENTREZID", keytype = "ENSEMBLTRANS",
    multiVals = "first"
  )
  ensg <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = base,
    column = "ENSEMBL", keytype = "ENSEMBLTRANS",
    multiVals = "first"
  )
  
  # From ENTREZ, pull ALIAS symbols
  alias_list <- rep(list(NULL), length(base))
  names(alias_list) <- base
  valid_entrez <- !is.na(entrez)
  if (any(valid_entrez)) {
    alias_tmp <- AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = entrez[valid_entrez],
      column = "ALIAS", keytype = "ENTREZID",
      multiVals = "list"
    )
    alias_list[valid_entrez] <- alias_tmp
  }
  
  syn_list <- vector("list", length(mRCC_model_ids))
  names(syn_list) <- mRCC_model_ids
  
  for (i in seq_along(mRCC_model_ids)) {
    mid <- mRCC_model_ids[i]   
    b   <- base[i]           
    
    s <- c(mid, b)
    
    if (!is.na(sym[i]))   s <- c(s, sym[i])
    if (!is.na(ensg[i]))  s <- c(s, ensg[i])
    if (!is.na(entrez[i])) s <- c(s, as.character(entrez[i]))
    
    if (!is.null(alias_list[[i]]) && length(alias_list[[i]]) > 0) {
      s <- c(s, alias_list[[i]])
    }
    
    syn_list[[i]] <- unique(clean_gene_ids(s))
  }
  
  all_syn <- unlist(syn_list, use.names = FALSE)
  all_mid <- rep(names(syn_list), times = lengths(syn_list))
  keep    <- !duplicated(all_syn)
  
  synonym_to_model <- all_mid[keep]
  names(synonym_to_model) <- all_syn[keep]
  
  list(
    synonym_to_model = synonym_to_model,
    model_ids        = mRCC_model_ids
  )
}

# Generic colname normalizer that only RENAMEs columns, does not touch values
normalize_colnames_generic <- function(cols, mapping) {
  cleaned <- clean_gene_ids(cols)
  out     <- cols
  hits    <- cleaned %in% names(mapping$synonym_to_model)
  out[hits] <- mapping$synonym_to_model[ cleaned[hits] ]
  out
}

###############################################################
# ===================== UI =========================
###############################################################
ui <- fluidPage(
  titlePanel("ImmunoResponse Predictor"),
  sidebarLayout(
    sidebarPanel(
      br(), br(),
      fileInput("testFile", "Upload Test.csv (For mUC or mRCC)", accept = c(".csv")),
      br(),
      selectInput(
        "modelFile", "Select trained model", 
        choices = c(
          "select model" = "",
          "mUC model"  = "logistic-Model-train-muc-test-muc.rds", 
          "mRCC model" = "logistic-Model-train-rcc-test-rcc.rds"
        )
      ),
      br(),
      actionButton("predictButton", "Make predictions"),
      br(), br(), br(), br(), br(),
      uiOutput("downloadUI")
    ),
    mainPanel(
      h4("Instructions:"),
      p("1. ", strong("Upload Test Data:"), "Upload a CSV file containing gene-expression matrix."),
      p("2. ", strong("File Structure:"), "The rows of the file represent individual samples, while the columns correspond to gene expression data."),
      p("3. ", strong("Gene IDs:"), "User can use one of the following gene identifiers in the columns:"),
      tags$ul(
        tags$li(strong("Gene Symbols"), "(e.g., TP53)"),
        tags$li(strong("Entrez Gene IDs"), "(e.g., 7157)"),
        tags$li(strong("Ensembl Gene IDs"), "(e.g., ENSG00000141510)"),
        tags$li(strong("Ensembl Transcript IDs"), "(e.g., ENST00000269305 or ENST00000269305.3)")
      ),
      p("4. ", strong("Select Model:"), "Choose a pre-trained model (either mUC or mRCC) from the dropdown."),
      p("5. ", strong("Generate Predictions:"), "Click the 'Generate predictions' button to process the data."),
      p("6. ", strong("Applicability Metric:"), "After predictions, the system calculates a score based on how closely your data matches the model’s trained biological structure."),
      p("7. ", strong("Download Predictions:"), "After predictions are made, you can download a CSV file containing:"),
      tags$ul(
        tags$li(strong("Sample ID")),
        tags$li(strong("Cosine Distances from Rs and NRs groups")),
        tags$li(strong("Predicted Response"), "(R or NR)"),
        tags$li(strong("Applicability Score"), "(Cosine Distance)")
      ),
      
      br(),
      textOutput("status"),
      textOutput("predictionsCount")
    )
  )
)

###############################################################
# ===================== SERVER =========================
###############################################################
server <- function(input, output, session) {
  predictions <- reactiveVal(NULL)
  
  # Read uploaded test file
  testData <- reactive({
    req(input$testFile)
    expr.test0 <- fread(
      input$testFile$datapath,
      header = TRUE,
      sep = ",",
      na.strings = c("", "NA")
    )
    sampleIDs <- expr.test0[[1]]                # First column as sample IDs
    expr.test0 <- expr.test0[, -1, with = FALSE]  # Remove sample ID column
    expr.test0_matrix <- as.matrix(expr.test0)
    rownames(expr.test0_matrix) <- sampleIDs
    expr.test0_matrix[is.na(expr.test0_matrix)] <- 0
    list(data = expr.test0_matrix, sampleIDs = sampleIDs)
  })
  
  # Cosine-distance helpers
  compute_average_cosine_distances <- function(y1, y2, x) {
    cosine_distance <- function(a, b) {
      norm_a <- sqrt(sum(a^2, na.rm = TRUE))
      norm_b <- sqrt(sum(b^2, na.rm = TRUE))
      if (is.na(norm_a) || is.na(norm_b) || norm_a == 0 || norm_b == 0) return(NA_real_)
      1 - sum(a * b, na.rm = TRUE) / (norm_a * norm_b)
    }
    avg_distance_R  <- if (nrow(y1) > 0) mean(apply(y1, 1, function(sample) cosine_distance(x, sample)), na.rm = TRUE) else NA
    avg_distance_NR <- if (nrow(y2) > 0) mean(apply(y2, 1, function(sample) cosine_distance(x, sample)), na.rm = TRUE) else NA
    c(avg_distance_R, avg_distance_NR)
  }
  
  # ORR prior
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
    list(results_df = results_df, prior_final_percentage = prior_final_percentage)
  }
  
  save_and_report_results <- function(results_df, prior_final_percentage) {
    results_df$`% supporting by CosineDist_prior` <- ""
    results_df$`% supporting by CosineDist_prior`[1] <- prior_final_percentage
    results_df
  }
  
  observeEvent(input$predictButton, {
    req(input$testFile, input$modelFile)
    tryCatch({
      withProgress(message = "Generating predictions", value = 0.1, {
        # ---------------- Load model ----------------
        model_path <- file.path("models", input$modelFile)
        if (!file.exists(model_path)) stop("Model file not found: ", model_path)
        
        bestModel <- readRDS(model_path)
        model_gene_ids <- bestModel$beta@Dimnames[[1]]
        
        # Build mapping depending on model
        if (input$modelFile == "logistic-Model-train-muc-test-muc.rds") {
          mUC_entrez <- as.numeric(sub("^X", "", model_gene_ids))
          mapping <- build_mUC_mapping(mUC_entrez, model_gene_ids)
        } else if (input$modelFile == "logistic-Model-train-rcc-test-rcc.rds") {
          mapping <- build_mRCC_mapping(model_gene_ids)
        } else {
          stop("Please select a valid model.")
        }
        
        # ---------------- Test data ----------------
        test_data <- testData()
        expr.test0_matrix <- test_data$data
        sampleIDs <- test_data$sampleIDs
        
        # Normalize uploaded gene IDs to model IDs
        colnames(expr.test0_matrix) <- normalize_colnames_generic(colnames(expr.test0_matrix), mapping)
        
        # Model gene panel (strip any X prefixes)
        gene_ids_clean       <- sub("^X", "", model_gene_ids)
        test_colnames_clean  <- sub("^X", "", colnames(expr.test0_matrix))
        match_idx            <- match(gene_ids_clean, test_colnames_clean)
        missing_mask         <- is.na(match_idx)
        
        # ------------- Branch by model: load TRAIN + align ----------------
        if (input$modelFile == "logistic-Model-train-muc-test-muc.rds") {
          # ----- mUC block -----
          train_zip_path <- file.path("models/mUC_log2TPMp1_train.zip")
          if (!file.exists(train_zip_path)) stop("Training zip file not found: ", train_zip_path)
          train_data_raw_mUC <- loadTable(
            file = unz(train_zip_path, "log2TPMp1_train.csv"),
            transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
          )
          if (!is.matrix(train_data_raw_mUC)) stop("mUC training data is not a matrix")
          sampleIDs_train <- rownames(train_data_raw_mUC)
          
          # If some genes missing in test → add zero columns
          if (any(missing_mask)) {
            add_mat <- matrix(0,
                              nrow = nrow(expr.test0_matrix),
                              ncol = sum(missing_mask))
            colnames(add_mat) <- gene_ids_clean[missing_mask]
            expr.test0_matrix <- cbind(expr.test0_matrix, add_mat)
            test_colnames_clean <- c(test_colnames_clean, gene_ids_clean[missing_mask])
            match_idx <- match(gene_ids_clean, test_colnames_clean)
          }
          
          # Build test matrix in model-gene order
          x.test <- expr.test0_matrix[, match_idx, drop = FALSE]
          
          # Align train columns with test order and scale
          train_cols_order <- match(
            sub("^X", "", colnames(x.test)),
            sub("^X", "", colnames(train_data_raw_mUC))
          )
          if (any(is.na(train_cols_order))) {
            missing_in_train <- colnames(x.test)[is.na(train_cols_order)]
            stop("Model genes missing in mUC training data: ", paste(missing_in_train, collapse = ", "))
          }
          trainM <- train_data_raw_mUC[, train_cols_order, drop = FALSE]
          testM  <- x.test
          
          train_data <- scale(trainM)
          x.test     <- scale(testM)
          train_data[is.na(train_data)] <- 0
          x.test[is.na(x.test)]         <- 0
          
          # TRAIN annotations
          annot_zip_path <- file.path("models/mUC_response_train.zip")
          if (!file.exists(annot_zip_path)) stop("Annotation zip file not found: ", annot_zip_path)
          sampleAnnot.train <- read.csv(unz(annot_zip_path, "response_train.csv"))
          
          if (!"RNASEQ_SAMPLE_ID" %in% colnames(sampleAnnot.train)) {
            stop("No RNASEQ_SAMPLE_ID column in mUC annotations")
          }
          if (!any(sampleIDs_train %in% sampleAnnot.train$RNASEQ_SAMPLE_ID)) {
            stop("No common sample IDs between mUC training data and annotations")
          }
          common_samples <- intersect(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID)
          sample_idx <- which(sampleIDs_train %in% common_samples)
          sampleIDs_train <- sampleIDs_train[sample_idx]
          train_data <- train_data[sampleIDs_train, , drop = FALSE]
          sampleAnnot.train <- sampleAnnot.train[match(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID), ]
          
          train_data_ordered <- train_data[
            match(sampleAnnot.train$RNASEQ_SAMPLE_ID, rownames(train_data)),
            , drop = FALSE
          ]
          y1 <- train_data_ordered[sampleAnnot.train$RESPONSE == 1 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
          y2 <- train_data_ordered[sampleAnnot.train$RESPONSE == 0 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
          
        } else if (input$modelFile == "logistic-Model-train-rcc-test-rcc.rds") {
          # ----- mRCC block -----
          train_zip_path <- file.path("models/standardized_QN_TPM_train.csv.gz")
          if (!file.exists(train_zip_path)) stop("Training csv.gz file not found: ", train_zip_path)
          train_data_raw_mRCC <- loadTable(
            file = gzfile(train_zip_path, "rt"),
            transpose = FALSE, convertToMatrix = TRUE, sep = ",", header = TRUE
          )
          if (!is.matrix(train_data_raw_mRCC)) stop("RCC training data is not a matrix")
          sampleIDs_train <- rownames(train_data_raw_mRCC)
          
          # Per-gene training means (for imputation)
          train_means_vec <- colMeans(
            train_data_raw_mRCC[, intersect(colnames(train_data_raw_mRCC), gene_ids_clean), drop = FALSE],
            na.rm = TRUE
          )
          train_means_full <- setNames(rep(NA_real_, length(gene_ids_clean)), gene_ids_clean)
          common_genes <- intersect(names(train_means_vec), gene_ids_clean)
          train_means_full[common_genes] <- train_means_vec[common_genes]
          
          # Impute missing genes in TEST using training means (fallback 0)
          if (any(missing_mask)) {
            add_vals <- train_means_full[missing_mask]
            add_vals[is.na(add_vals)] <- 0
            add_mat <- matrix(
              rep(add_vals, each = nrow(expr.test0_matrix)),
              nrow = nrow(expr.test0_matrix),
              byrow = FALSE
            )
            colnames(add_mat) <- gene_ids_clean[missing_mask]
            expr.test0_matrix <- cbind(expr.test0_matrix, add_mat)
            test_colnames_clean <- c(test_colnames_clean, gene_ids_clean[missing_mask])
            match_idx <- match(gene_ids_clean, test_colnames_clean)
          }
          
          # Build matrices in model-gene order (no extra scaling: QN+standardized already)
          x.test <- expr.test0_matrix[, match_idx, drop = FALSE]
          train_cols_order <- match(
            sub("^X", "", colnames(x.test)),
            sub("^X", "", colnames(train_data_raw_mRCC))
          )
          if (any(is.na(train_cols_order))) {
            missing_in_train <- colnames(x.test)[is.na(train_cols_order)]
            stop("Model genes missing in RCC training data: ", paste(missing_in_train, collapse = ", "))
          }
          train_data <- train_data_raw_mRCC[, train_cols_order, drop = FALSE]
          
          # TRAIN annotations
          annot_zip_path <- file.path("models/RCC_response_train.zip")
          if (!file.exists(annot_zip_path)) stop("Annotation zip file not found: ", annot_zip_path)
          sampleAnnot.train <- read.csv(unz(annot_zip_path, "response_train.csv"))
          
          if (!"RNASEQ_SAMPLE_ID" %in% colnames(sampleAnnot.train)) {
            stop("No RNASEQ_SAMPLE_ID column in RCC annotations")
          }
          if (!any(sampleIDs_train %in% sampleAnnot.train$RNASEQ_SAMPLE_ID)) {
            stop("No common sample IDs between RCC training data and annotations")
          }
          common_samples <- intersect(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID)
          sample_idx <- which(sampleIDs_train %in% common_samples)
          sampleIDs_train <- sampleIDs_train[sample_idx]
          train_data <- train_data[sampleIDs_train, , drop = FALSE]
          sampleAnnot.train <- sampleAnnot.train[match(sampleIDs_train, sampleAnnot.train$RNASEQ_SAMPLE_ID), ]
          
          train_data_ordered <- train_data[
            match(sampleAnnot.train$RNASEQ_SAMPLE_ID, rownames(train_data)),
            , drop = FALSE
          ]
          y1 <- train_data_ordered[sampleAnnot.train$RESPONSE == 1 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
          y2 <- train_data_ordered[sampleAnnot.train$RESPONSE == 0 & !is.na(sampleAnnot.train$RESPONSE), , drop = FALSE]
        } else {
          stop("Please select a valid model.")
        }
        
        # Final guard: ensure no NA in x.test before prediction
        if (anyNA(x.test)) {
          x.test <- scrub_na(x.test, fill = 0)
          output$status <- renderText("Some NA values detected after preprocessing; replaced with 0 to proceed.")
        }
        
        # ---------------- Manual logistic prediction ----------------
        x.test <- as.matrix(x.test)
        if (any(is.na(x.test))) stop("NA values detected in x.test before prediction")
        beta <- as.matrix(bestModel$beta)
        a0   <- bestModel$a0
        if (ncol(x.test) != nrow(beta)) {
          stop("Dimension mismatch: x.test has ", ncol(x.test), " columns, beta has ", nrow(beta), " rows")
        }
        pred_prob  <- 1 / (1 + exp(-(x.test %*% beta + a0)))
        pred_prob  <- as.vector(pred_prob)
        pred_class <- ifelse(pred_prob > 0.5, 1, 0)
        pred_labels <- ifelse(pred_class == 1, "R", "NR")
        
        # ---------------- Cosine distance block ----------------
        results_list <- vector("list", nrow(x.test))
        for (i in seq_len(nrow(x.test))) {
          x_sample <- x.test[i, , drop = FALSE]
          avg_distances <- compute_average_cosine_distances(y1, y2, x_sample)
          results_list[[i]] <- data.frame(
            CosDist_2_Rs        = avg_distances[1],
            CosDist_2_NRs       = avg_distances[2],
            LogitDA_Score       = pred_prob[i],
            LogitDA_pred_label  = pred_labels[i]
          )
        }
        results_df <- do.call(rbind, results_list)
        results_df$sampleID <- sampleIDs
        results_df <- results_df[match(sampleIDs, results_df$sampleID), ]
        results_df <- results_df[, c("sampleID", setdiff(names(results_df), "sampleID"))]
        stopifnot(all(results_df$sampleID == sampleIDs))
        
        # ---------------- ORR prior + final table ----------------
        orr_results <- apply_orr_prior(results_df, orr = 0.2)
        results_df  <- save_and_report_results(
          results_df = orr_results$results_df,
          prior_final_percentage = orr_results$prior_final_percentage
        )
        predictions(results_df)
        
        # ---------------- UI feedback ----------------
        msg_bits <- c("Predictions generated successfully!")
        if (any(missing_mask)) {
          listed <- paste(head(gene_ids_clean[missing_mask], 8), collapse = ", ")
          suffix <- if (sum(missing_mask) > 8) " ..." else ""
          msg_bits <- c(msg_bits)
        } else {
          msg_bits <- c(msg_bits, "No missing genes; exact match to model panel.")
        }
        msg_bits <- c(msg_bits, sprintf("%% supporting by CosineDist_prior: %.2f%%",
                                        orr_results$prior_final_percentage))
        output$status <- renderText(paste(msg_bits, collapse = " | "))
        
        output$downloadUI <- renderUI({ downloadButton("downloadPredictions", "Download predictions") })
      })
    }, error = function(e) {
      output$status <- renderText(paste("Error generating predictions:", e$message))
    })
  })
  
  output$predictionsCount <- renderText({
    preds <- predictions()
    if (is.null(preds)) return("Rows in predictions: 0")
    paste("Rows in predictions:", nrow(preds))
  })
  
  output$downloadPredictions <- downloadHandler(
    filename = function() paste("predictions-", Sys.Date(), ".csv", sep = ""),
    content = function(file) {
      preds <- predictions()
      if (is.null(preds)) stop("No predictions available to download. Please generate predictions first.")
      fwrite(preds, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  output$status <- renderText("")
}

shinyApp(ui = ui, server = server)
