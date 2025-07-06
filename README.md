
# Fast and Robust PRESS for PCA

## Overview

This project explores robust methods for determining the optimal number of components to retain in **Principal Component Analysis (PCA)** using **Prediction Error Sum of Squares (PRESS)** and its robust variant **R-PRESS**. Standard methods like the elbow method or percentage of variance explained are not ideal when prediction accuracy or robustness is of primary interest. This work implements and evaluates a fast, robust approach to computing PRESS values, especially for high-dimensional data.

This project was conducted as part of the **Robust Statistics (G0B16)** course during Spring 2022.

## Contents

- 📄 `group7 - report_comments.pdf` — Full project report with theoretical background, simulations, and application.
- 🧪 R scripts — Code for simulations, method implementation, and data analysis.
- 📊 Plots — Visual comparisons of PRESS vs. R-PRESS and application results.
- 📁 `data/` — (Optional) Folder containing datasets used in the project.

## Objectives

1. **Understand the PRESS criterion** for component selection in PCA and its advantages over traditional methods.
2. **Explore robustness** of the PRESS method and introduce a robust version (R-PRESS).
3. **Implement the fast CV algorithm** for robust PCA based on Hubert & Engelen (2007).
4. **Reproduce key experiments** from the paper to compare PRESS and R-PRESS.
5. **Apply the method** to a high-dimensional real-world dataset and interpret the results.

## Methodology

- Literature review of the fast and robust CV method from:
  > Hubert, M., & Engelen, S. (2007). *Fast cross-validation of high-breakdown resampling methods for PCA*. Computational Statistics & Data Analysis, 51(10), 5013–5024. [DOI](https://doi.org/10.1016/j.csda.2006.08.031)

- Implementation of the robust R-PRESS procedure in R.
- Simulation studies to compare PRESS and R-PRESS under various conditions.
- Application of R-PRESS to a high-dimensional dataset (`p > n`) with exploratory analysis and interpretation.

## Requirements

- R (version ≥ 4.0)
- R packages:
  - `rrcov`
  - `robustbase`
  - `ggplot2`
  - `MASS`

## Usage

To run the analysis:

1. Open R or RStudio.
2. Source the main script (e.g., `main.R`) to run simulations and data analysis:
   ```R
   source("main.R")
   ```
3. Results including PRESS plots and diagnostics will be saved to the `plots/` folder.

## Results

- **R-PRESS** provides more stable and robust estimates for selecting the number of components in PCA, especially under contamination or when `p > n`.
- The fast CV algorithm significantly reduces computational cost compared to brute-force LOO-CV.

## Contributors

- Group 7, Robust Statistics (Spring 2022)
- Course Instructors: Peter Rousseeuw, Ioannis Kalogridis, Jordy Menvouta

## License

This project is for academic purposes only.
