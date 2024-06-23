# High-dimensional Expected Shortfall Regression

## Overview

This project implements the estimation and inference of the two-step approach for high-dimensional expected shortfall regression proposed in Zhang et. al (2024+). 

## Directory Structure

- **src/**: Contains the core implementation of the method.
  - `highD_2step.R`: Main R script with the `highdim_2step` and `highdim_inf` functions.
- **Simulation/**: Includes scripts and results for simulations that demonstrate the usage of the implemented functions.
- **DataApplication/**: Contains scripts and results for a data application, health disparity research, illustrating the practical utility of the method.

## Key Functions

### highdim_2step
This function performs the estimation procedure for the high-dimensional expected shortfall regression. It involves computing an $\ell_1$-penalized quantile regression estimator and a lasso-type estimator with some adjusted response variable. Specifically, we employ the `R` package `conquer` with the default tuning parameters to obtain an $\ell_1$-penalized smoothed quantile regression estimator. Given a quantile regression estimator, we compute the expected shortfall regression estimator via the `glmnet` package. The sparsity tuning parameters for both steps are selected using ten-fold cross-validation.

### highdim_inf
This function carries out the inference procedure based on the estimations obtained from `highdim_2step`. Based on Zhang et. al (2024+), we compute the debiased estimator in (3.6) using `glmnet` package, where the tuning parameter is chosen by ten-fold cross-validation with the "one standard error rule" (see Section~7.10 in Hastie et al. 2009), that is, we choose the largest lambda whose error is no more than one standard error above the minimum mean cross-validated error. We then estimate the asymptotic variance of the debiased estimator using refitted cross-validation as described in Section 2.4, and calculate confidence intervals. 


## Workflow

This section outlines the steps to use the `highdim_2step` and `highdim_inf` functions for high-dimensional data analysis.

### Step 1: Setup Environment

1. **Clone the repository**:
    ```{bash}
    git clone https://github.com/shushuzh/ES_highD.git
    cd ES_highD
    ```

2. **Install necessary packages**: Ensure that you have the required R packages installed. You can install them using the following R command:
    ```{r,eval = FALSE}
    install.packages(c("glmnet", "conquer", "dplyr","Matrix","MultiRNG","janitor"))
    ```

### Step 2: Load Functions

**Source the R script**: Load the functions from the `src/highD_2step.R` file.
    ```{r,eval = FALSE}
    source("src/highD_2step.R")
    ```

### Step 3: Prepare Data

**Load your data**: Import your dataset into R. Ensure your data is preprocessed with a design matrix (`x`) and response (`y`).
    ```{r,eval = FALSE}
    data <- read.csv("path/to/your/data.csv")  # Replace with your actual data file
    ```

### Step 4: Estimation Procedure

**Run the estimation procedure**: Use the `highdim_2step` function to estimate the `alpha`-th conditional expected shortfall. You can choose not to standardize the covariates by setting "standardize=F". 
    ```{r,eval = FALSE}
    result_estimation <- highdim_2step(x,y,alpha,standardize=F)
    ```

### Step 5: Inference Procedure


**Run the inference procedure**: Use the `highdim_inf` function to perform inference on the covariates of interest. Specify the column number of the desired covariate using the `col` parameter. This function provides confidence intervals at the `conf.level` confidence level.
    ```{r,eval = FALSE}
    result_inference <- highdim_inf(x,y,alpha,col = col ,res_est=result_estimation, conf.level = 0.95, standardize=F)
    ```

### Step 6: Analyze Results

**Review results**: Examine the output from the estimation and inference procedures to draw conclusions.
    ```{r,eval = FALSE}
    print(result_estimation$theta_hat) # ES estimator (slope)
    print(result_estimation$theta0_hat) # ES estimator (intercept)
    print(result_inference$theta_debias) # Debiased ES estimator
    print(result_inference$Conf.int) # Confidence interval for the parameter of interest
    ```
---

This workflow should help you utilize the `highdim_2step` and `highdim_inf` functions effectively for your high-dimensional data analysis. For detailed examples, refer to the scripts in the `Simulation` and `DataApplication` folders.


## Demonstrations
### Simulations
To understand how to use the functions for simulations, refer to "Simulation/simulation_final.R". This folder also contains the results obtained from 500 replications of the simulations, which can serve as a reference for expected outcomes.

### Data Applications
For practical applications, we analyze data from the National Health and Nutrition Examination Survey(NHANES) for the years 2017 to 2020 (https:
//wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?Cycle=2017-2020). Our study focuses on examining health disparities in high cotinine values among different ethnic groups.

The data preprocessing is handled in "DataApplication/get_data.R" and "DataApplication/disparity_ES.R". The analysis itself is conducted in "DataApplication/disparity_ES.R". Results at various quantile levels are recorded in "DataApplication/results.csv".

## References

- Zhang S., He X., Tan, K. M. & Zhou, W.-X. (2024+), High-Dimensional expected shortfall regression. https://arxiv.org/abs/2307.02695. 

- Man, R., Pan, X., Tan, K. M. & Zhou, W.-X. (2024), ‘A unified algorithm for penalized convolution smoothed quantile regression’, Journal of Computational and Graphical Statistics 33(2), 625–637.



