V Gene Chi-Square Analysis

Overview

This project contains a Python script (v_gene_disease_analysis.py) that analyzes the distribution of V gene usage across different disease stages. The analysis aims to understand whether there is a significant association between V gene usage and disease progression. The insights gained from this analysis may provide valuable clues for targeted therapeutic interventions and enhance our understanding of immune response characteristics in different disease states.

Features

The script performs the following tasks:

- Counts V Genes by Disease Stage: It counts the occurrence of each V gene in different samples, associating these counts with the disease stages based on metadata provided.
- Grouped Bar Plot Visualization: Generates a grouped bar plot to visualize the distribution of V gene usage across different disease stages.
- Chi-Square Test: Performs a Chi-Square test to determine if there is a statistically significant association between V gene usage and disease stages.
- Standardized Residuals Calculation: Calculates standardized residuals to identify which V genes and disease stages contribute most to any observed statistical association.
- Heatmap Visualization of Residuals: Creates a heatmap of the residuals to identify patterns of over- or under-representation of specific V genes in each disease stage.

Input Files:
- A CSV containing metadata for the samples, with the following columns (example: metadata_Disease_ETC.csv):
  - 'sample_id': A unique identifier for each sample.
  - 'disease_stage': The disease stage for each sample.
  - SEA.txt: This file should be a tab-separated file containing V gene usage data for each sample. It should include:
  - 'sample_id': A unique identifier for each sample (must match with the metadata file).
  - 'V': The V gene used by each sequence in the sample.

Output Files:
- A CSV file containing the counts of V genes for each disease stage, which is generated as an intermediate result.
- A CSV file containing standardized residuals, useful for identifying significant deviations from expected values.

How to Run:
Ensure you have Python 3 and the required libraries installed. The script relies on the following Python packages:
pandas
scipy
seaborn
matplotlib
numpy

Results and Insights
The grouped bar plot provides a visualization of how V gene usage is distributed across different disease stages, helping to identify any visible trends.
The Chi-Square test results indicate whether there is a statistically significant relationship between V gene usage and disease stages.
The heatmap of standardized residuals highlights which specific V genes are over- or under-represented in certain disease stages, pointing to potentially interesting immunological patterns.

Notes
Significance Threshold: The threshold of |2| for standardized residuals is used to identify significant deviations from expected values, as values beyond ±2 are typically considered noteworthy in statistical analysis.
Error Handling: The script checks if sample IDs from the V gene file are found in the metadata file, and warns if they are missing to ensure proper matching.

Contributing

Contributions are welcome! If you would like to improve the code or add more features, feel free to submit a pull request. Please ensure that your contributions are well-documented.
