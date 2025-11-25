# ImmunoResonsePredictor - A Shiny based GUI

**ImmunoResonsePredictor** is an interactive R-based Shiny web application that allows users to generate predictions using pre-trained predictors on uploaded test data. The app supports logistic regression models for two different types of cancer: mUC (metastasis urothelial carcinoma) and mRCC (metastasis renal cell carcinoma). Users can upload a test dataset, select a model, and generate predictions, which can then be downloaded as a CSV file for further analysis.

![image](https://github.com/user-attachments/assets/3c6a115c-b703-48e4-a0d8-50cefd5ff24a)


## Features

- **Upload Test Data**: Upload a CSV file containing gene-expression matrix.
- **File Structure**: The rows of the file represent individual samples, while the columns correspond to gene expression data.
- **Gene IDs**: You can use one of the following gene identifiers in the columns:

   - Gene Symbols (e.g., TP53)
   - Entrez Gene IDs (e.g., 7157)
   - Ensembl Gene IDs (e.g., ENSG00000141510)
   - Ensembl Transcript IDs (e.g., ENST00000269305)

- **Gene Panel**: For **mUC model**, the file must include **49 specific signature genes**; for **mRCC model**, the file must have **27 signature genes**.
- **Select a Pre-Trained Model**: Choose from two pre-trained LogitDA models:
  - mUC Model
  - mRCC Model
- **Generate Predictions**: After uploading your file and selecting the model, click on the ‘Generate predictions’ button to process the data.
- The system will:
  - Match your uploaded gene expression data with the required signature genes for the selected model.
  - Handle missing genes by either imputing data from training means or replacing with zeros.
  - Normalize gene IDs to ensure consistency across the dataset.
  - Run predictions using the LogitDA algorithm to classify each sample as Responsive (R) or Non-Responsive (NR).
- **Applicability Metric**:
  - Alongside predictions, the system computes an **applicability score** for each sample, based on cosine-distance similarity between your dataset and the model’s trained R/NR profiles.
  - This metric gives you an interpretability score that indicates how closely your samples align with the biological structure learned by the model.
  - The higher the applicability score, the better the model's predictions align with the expected biological pattern.
- **Download Predictions**: Once predictions are made, you can download the results as a CSV file containing:
  - Sample ID
  - Cosine Distances from Rs and NRs groups
  - Predicted Response (R or NR)
  - Applicability Score (cosine distance)
- **Model Validation**: Ensures that the test dataset has the required features (genes) for the selected model and notifies the user if there are any discrepancies.


## How It Works

1. **Upload Your Test Data**: The test data should be in CSV format with features (genes) as columns and sample IDs as the first column.
2. **Select a Pre-trained Model**: Choose either the mUC or mRCC model.
3. **Generate Predictions**: Once the model and test data are uploaded, click the **"Generate Predictions"** button. The app will process the data and display results.
4. **Download the Results**: After generating the predictions, a **download button** will appear, allowing you to save the predictions as a CSV file.

## Requirements

- **R Packages**:  
  This project uses [renv](https://rstudio.github.io/renv/) to manage all required packages and their versions.  
  The recommended way to install all dependencies is to run in your R console (from the project directory):

  ```r
  install.packages("renv")
  renv::restore()

or the packages can be installed manually as given below:
- **R**: Version 4.0.0 or higher.
- **Libraries**: The following libraries are required:
  ```r
  install.packages(c(
  "shiny",
  "glmnet",
  "data.table"
  ))

 **Install Bioconductor packages**
  ```r
      if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
      BiocManager::install(c("sva", "preprocessCore"))
```
- A user just needs to run in their R console:
  ```r
  source("install.R")

- **Trained Models**: The app requires pre-trained logistic regression models saved as `.rds` files. The models for mUC and mRCC are included in this project inside models folder.


---
# The GUI can either be downloaded directly from this repository or accessed online via following link:

    https://logitda.shinyapps.io/immunoresponsepredictor/

---
# 1.  Running the ImmunoResponsePredictor locally:

## Detailed steps to run this code locally
- Running the ImmunoResonsePredictor on local system. 
- To download and run the Shiny app from the repository at https://github.com/rajatbutola/ImmunoResponsePredictor on your local system, follow these step-by-step instructions.

## Prerequisites

Before starting, ensure the following tools are installed:

- [**R**](https://cran.r-project.org/): Version **4.3.3** or higher  
- [**RStudio**](https://posit.co/download/rstudio-desktop/): Recommended IDE for running Shiny apps

## Step 1: Clone the Repository 
- First download the source code of ImmunoResponcePredictor from Github to your local system. Open your terminal or command prompt and navigate to your project directory:

  ```r
  cd /path/to/your/projects


- Run the following command to clone the repository from GitHub:

  ```r
  git clone https://github.com/rajatbutola/ImmunoResponsePredictor.git

- After this you will successfully download the code into a folder named ImmunoResponsePredictor containing the repository files.

## Requirements

## Step 2. Install Required R Packages:

### R Packages:
- This project uses renv to manage all required packages and their versions.
- The recommended way to install all dependencies is to run in your R console (from the project directory):

  ```r
  install.packages("renv")
  renv::restore()

- or the packages can be installed manually as given below:

  ```r
  Libraries: The following libraries are required:
  install.packages(c(
  "shiny",
  "glmnet",
  "data.table"
  ))

- Install Bioconductor packages
  ```r
    if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
    BiocManager::install(c("sva", "preprocessCore"))

- A user just needs to run in their R console:

  ```r
  source("install.R")


## Step 3: Confirm Directory Structure 
- After installing the packages please check the Shiny App Structure:

- The folder should be arranged as below:


ImmunoResponsePredictor/

 ├── app.R
 
 ├── renv.lock
 
 ├── README.md
 
 ├── Dockerfile
 
 ├── models/
 
 │   ├── logistic-Model-train-muc-test-muc.rds
 
 │   └── logistic-Model-train-rcc-test-rcc.rds


    
- Apart from this you will also need train and test datasets (count matrix) which are not provided here with this Github repository due to copyright and sensitivity issues. The data can be provided upon request. 

## Step 4. Run the Shiny App Locally

A.	Open the main app file:
In RStudio, open app.R.

B.	Run the app:
Click the “Run App” button in RStudio (top-right of the script editor).
Or, in the R console, run:
  ```r
  shiny::runApp()
```

(Ensure the working directory is set to the ImmunoResponsePredictor folder.)

Step 5. Access the app:
-	The Shiny app should open in a browser window.
-	To interact with the app to verify it works use a “Browse” button for uploading .csv files containing gene expression data, 
-	a status indicator will confirm that the “Upload complete”, 
-	then on a dropdown menu select the appropriate pre-trained LogitDA model (mUC or mRCC), and click “Make predictions” button to initiate response prediction on the uploaded data, and 
-	at last, click on “Download predictions” button to export the results in .csv format. 
-	In addition, the interface provides real-time feedback, including the number of rows in the output predictions and confirmation messages upon successful prediction generation.


# 2.  Running the ImmunoResponsePredictor Online:
You can use the GUI directly in your browser without downloading or installing anything locally. Simply visit the following link:

    https://logitda.shinyapps.io/immunoresponsepredictor/

## Important Note for mRCC Model Users
The mRCC model requires intensive preprocessing using "ComBat + Quantile" normalization due to batch effect correction and distribution alignment. These steps are computationally heavy and memory-intensive, making it difficult to run them directly on the hosting server due to RAM limitations.

To address this, we recommend users preprocess their test data locally using the provided preprocess.R script included in the repository.

### Please preprocess your test gene expression matrix using ComBat + Quantile normalization before uploading it to the web app.

## Instructions for Using the Online App
- Upload your preprocessed ```Test.csv``` file.

- Select the desired pre-trained model (```mUC``` or ```mRCC```) from the dropdown.

- Click the “Make predictions” button to initiate processing.

- Download the resulting predictions as a ```.csv``` file.

The app will provide real-time feedback during each step, confirming successful uploads and prediction generation.              


Clone the repository:

```bash
git clone https://github.com/rajatbutola/ImmunoResponsePredictor.git

 
