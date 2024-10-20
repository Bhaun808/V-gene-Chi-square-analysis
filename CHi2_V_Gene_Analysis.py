# V Gene Chi-Square Analysis

"""
This script is designed to analyze the distribution of V gene usage across different disease stages. It performs the following tasks:
1. Counts V genes in various samples and associates them with disease stages based on metadata.
2. Visualizes the distribution of V genes by disease stage using a grouped bar plot.
3. Performs a Chi-Square test to determine if there is a statistically significant association between V gene usage and disease stage.
4. Calculates standardized residuals to identify specific V genes and disease stages that contribute most to any observed association.
5. Visualizes these residuals in a heatmap to help identify patterns of over- or under-representation of V genes in specific disease stages.

The insights gained from this analysis are important for understanding how immune response characteristics, such as V gene usage, vary with disease progression, and may provide clues for targeted therapeutic interventions.

Input File Descriptions:
1. metadata_Disease_ETC.csv: This file should be a CSV containing metadata for the samples. It must include the following columns:
   - 'sample_id': A unique identifier for each sample.
   - 'disease_stage': The disease stage for each sample.

2. SEA.txt: This file should be a tab-separated file containing V gene usage data for each sample. It must include the following columns:
   - 'sample_id': A unique identifier for each sample (matching the metadata file).
   - 'V': The V gene used by each sequence in the sample.

These input files are essential for linking V gene usage to disease progression and enabling the statistical analysis performed in this script.
"""

import pandas as pd
import os
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Step 1: Count V Genes Based on Disease Stages

# Load metadata file
# NOTE: This metadata file contains sample information, including disease stages for each sample.
metadata_df = pd.read_csv('metadata_Disease_ETC.csv')

# Path to the V gene file to analyze
# NOTE: This file contains sample IDs and V gene usage details, which will be analyzed to understand the distribution of V genes across different disease stages.
v_gene_file = 'SEA.txt'

# Initialize a dictionary to store the counts per disease stage
# NOTE: Storing counts per disease stage is crucial for understanding how V gene usage varies with disease progression, which may help in identifying key immunological patterns.
v_gene_counts = {}

# Function to count V genes in a file and associate them with the disease stage
def count_v_genes(file_path):
    """
    Counts V genes from the provided file and associates them with corresponding disease stages.

    Parameters:
    file_path (str): The path to the file containing V gene data for various samples.

    Expected Output:
    Updates the global v_gene_counts dictionary with the counts of each V gene per disease stage, allowing for subsequent statistical analysis.
    """
    try:
        df = pd.read_csv(file_path, sep='\t')
        for sample_id in df['sample_id'].unique():
            if sample_id in metadata_df['sample_id'].values:
                disease_stage = metadata_df.loc[metadata_df['sample_id'] == sample_id, 'disease_stage'].values[0]
                v_gene_count = df[df['sample_id'] == sample_id]['V'].value_counts().to_dict()
                if disease_stage not in v_gene_counts:
                    v_gene_counts[disease_stage] = {}
                for v_gene, count in v_gene_count.items():
                    if v_gene not in v_gene_counts[disease_stage]:
                        v_gene_counts[disease_stage][v_gene] = 0
                    v_gene_counts[disease_stage][v_gene] += count
            else:
                print(f"Warning: sample_id {sample_id} not found in metadata.")
    except Exception as e:
        print(f"Error reading {file_path}: {e}")

# Process the specified V gene file
print(f"Processing file: {v_gene_file}")
if os.path.exists(v_gene_file):
    count_v_genes(v_gene_file)
else:
    print(f"File not found: {v_gene_file}")

# Prepare data for CSV
data = []
for disease_stage, v_genes in v_gene_counts.items():
    for v_gene, count in v_genes.items():
        data.append({'disease_stage': disease_stage, 'V': v_gene, 'count': count})

results_df_v_gene = pd.DataFrame(data)

# Export the results
# NOTE: Saving this intermediate result is important for ease of future analysis or record-keeping.
results_df_v_gene.to_csv('SEA_v_gene_counts_by_disease_stage.csv', index=False)
print("V gene counts by disease stage have been exported to SEA_v_gene_counts_by_disease_stage.csv")

# Step 2: Grouped Bar Plot

# Load the results from the CSV file
results_df = pd.read_csv('SEA_v_gene_counts_by_disease_stage.csv')

# Create a grouped bar plot
# NOTE: The grouped bar plot helps visualize the distribution of V gene usage across different disease stages, providing insights into how V gene expression varies with disease progression.
plt.figure(figsize=(14, 8))
sns.barplot(x='V', y='count', hue='disease_stage', data=results_df)
plt.xticks(rotation=90)
plt.title('V Gene Counts by Disease Stage')
plt.show()

# Step 3: Chi-Square Test for V Gene Counts

# Pivot the DataFrame to create a contingency table
contingency_table = results_df.pivot_table(values='count', index='V', columns='disease_stage', aggfunc='sum', fill_value=0)

# Perform Chi-Square Test
# NOTE: The Chi-Square test is used here to determine if there is a statistically significant association between V gene usage and disease stage.
chi2, p, dof, expected = stats.chi2_contingency(contingency_table)

# Print the results
print(f"Chi-Square Statistic: {chi2}")
print(f"P-value: {p}")
print(f"Degrees of Freedom: {dof}")
print("Expected Frequencies:")
print(expected)

# Determine if the p-value is significant (commonly, p < 0.05)
if p < 0.05:
    print("There is a significant association between V genes and disease stages.")
else:
    print("There is no significant association between V genes and disease stages.")

# Step 4: Plot Significant Residuals and Export to CSV

# Calculate standardized residuals
# NOTE: Standardized residuals are calculated to understand how observed values deviate from expected values. They help identify which specific V genes and disease stages contribute most to the Chi-Square statistic. Significant residuals (absolute value > 2) indicate cells where the observed frequency is much higher or lower than expected, suggesting potential areas of interest.
residuals = (contingency_table - expected) / np.sqrt(expected)

# Print standardized residuals
print("Standardized Residuals:")
print(residuals)

# Highlight significant residuals (e.g., residuals with absolute value > 2)
# NOTE: A threshold of 2 is chosen because, in standard normal distribution, residuals beyond ±2 are considered significantly different from the expected values. This helps in identifying which specific observations deviate significantly from the expected frequency, providing insights into patterns that may be biologically relevant.
significant_residuals = residuals[(residuals > 2) | (residuals < -2)]
print("Significant Standardized Residuals (abs > 2):")
print(significant_residuals)

# Export standardized residuals to CSV for further analysis if needed
residuals.to_csv('SEA_standardized_residuals.csv')

# Visualize significant residuals using a heatmap
# NOTE: The heatmap helps visualize the magnitude and direction of the residuals across different V genes and disease stages. It is particularly useful for identifying patterns where certain V genes are over- or under-represented in specific disease stages, which can provide valuable biological insights.
plt.figure(figsize=(12, 8))
sns.heatmap(residuals, annot=True, fmt=".2f", cmap="coolwarm", center=0)
plt.title('Standardized Residuals of V Gene Counts by Disease Stage')
plt.show()
