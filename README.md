# ImmunoResonsePredictor- A Shiny based GUI

**ImmunoResonsePredictor** is an interactive R-based Shiny web application that allows users to generate predictions using pre-trained predictors on uploaded test data. The app supports logistic regression models for two different types of cancer: mUC (metastasis urothelial carcinoma) and mRCC (metastasis renal cell carcinoma). Users can upload a test dataset, select a model, and generate predictions, which can then be downloaded as a CSV file for further analysis.

## Features

- **Upload Test Data**: Upload a CSV file containing the test dataset.
- **Select Trained Model**: Choose from two pre-trained LogitDA models:
  - mUC Model
  - mRCC Model
- **Generate Predictions**: Click a button to process the data and generate predictions.
- **Download Predictions**: After generating the predictions, users can download the results as a CSV file.
- **Model Validation**: Ensures that the test dataset has the required features (genes) for the selected model and notifies the user if there are any discrepancies.

## How It Works

1. **Upload Your Test Data**: The test data should be in CSV format with features (genes) as columns and sample IDs as the first column.
2. **Select a Pre-trained Model**: Choose either the mUC or mRCC model.
3. **Generate Predictions**: Once the model and test data are uploaded, click the **"Generate Predictions"** button. The app will process the data and display results.
4. **Download the Results**: After generating the predictions, a **download button** will appear, allowing you to save the predictions as a CSV file.

## Requirements

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

## 

**ImmunoResponsePredictor** is an interactive Shiny-based application designed to predict immunotherapy response using pre-trained logistic models on gene expression data.

---

## 📌 Features

- Upload gene expression `.csv` files
- Choose between pre-trained mUC and mRCC models
- Generate and download immunotherapy response predictions
- Real-time feedback on upload and prediction status

---

## 🧰 Prerequisites

Before starting, ensure the following tools are installed:

- [**R**](https://cran.r-project.org/): Version **4.3.3** or higher  
- [**RStudio**](https://posit.co/download/rstudio-desktop/): Recommended IDE for running Shiny apps

---

## 📥 Step 1: Clone the Repository

Open your terminal or command prompt and navigate to your project directory:

```bash
cd /path/to/your/projects


**Detailed Steps to run this code**
- Running the ImmunoResonsePredictor on local system. 

- To download and run the Shiny app from the repository at https://github.com/rajatbutola/ImmunoResponsePredictor on your local system, follow these step-by-step instructions.

Prerequisites
Before starting, ensure you have:

R: Version 4.3.3 or higher (based on your package list, which includes base 4.3.3).
RStudio: Recommended for running Shiny apps (download from rstudio.com).

Step 1. First download the source code of ImmunoResponcePredictor from Github to your local system. For that open the command line interface and navigate to the directory you wanted GUI to save.

Example: cd /path/to/your/projects

Run the following command to clone the repository from GitHub:

git clone https://github.com/rajatbutola/ImmunoResponsePredictor.git

After this you will successfully download the code into a folder named ImmunoResponsePredictor containing the repository files.

Requirements
Step 2. Install Required R Packages:

R Packages:
This project uses renv to manage all required packages and their versions.
The recommended way to install all dependencies is to run in your R console (from the project directory):

install.packages("renv")
renv::restore()

or the packages can be installed manually as given below:

Libraries: The following libraries are required:
install.packages(c(
"shiny",
"glmnet",
"data.table"
))

Install Bioconductor packages

    if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
    BiocManager::install(c("sva", "preprocessCore"))

A user just needs to run in their R console:

source("install.R")


Step 3. After installing the packages please check the Shiny App Structure:

The folder should be arranged as below:

app.R
renv.lock
README.md
Dockerfile
models\
    logistic-Model-train-muc-test-muc.rds
    logistic-Model-train-rcc-test-rcc.rds
    
Apart from this you will also need train and test datasets (count matrix) which are not provided here with this Github repository due to copyright and sensitivity issues. The data can be provided upon request. 

Step 4. Run the Shiny App Locally

A.	Open the main app file:
In RStudio, open app.R.

B.	Run the app:
Click the “Run App” button in RStudio (top-right of the script editor).
Or, in the R console, run:
shiny::runApp()

(Ensure the working directory is set to the ImmunoResponsePredictor folder.)

Step 5. Access the app:
•	The Shiny app should open in a browser window.
•	To interact with the app to verify it works use a “Browse” button for uploading .csv files containing gene expression data, 
•	a status indicator will confirm that the “Upload complete”, 
•	then on a dropdown menu select the appropriate pre-trained LogitDA model (mUC or mRCC), and click “Make predictions” button to initiate response prediction on the uploaded data, and 
•	at last, click on “Download predictions” button to export the results in .csv format. 
•	In addition, the interface provides real-time feedback, including the number of rows in the output predictions and confirmation messages upon successful prediction generation.





Clone the repository:

```bash
git clone https://github.com/rajatbutola/ImmunoResponsePredictor.git

 
