# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 09:46:23 2024

@author: Jason
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

# Load the CSV file
file_path = r'C:/Users/Jason/Desktop/stochICE/examples/Chateauguay/MonteCarlo/SimulationTifs/Stats/largest_polygon_areas.csv'  # Replace with your actual file path
data = pd.read_csv(file_path)

# Define a function to extract the needed parts from the filenames
def extract_parts(filename):
    parts = filename.split('_')
    Q = float(parts[1])
    thickness = float(parts[2])
    phi = float(parts[3])
    toe_location = float(parts[4].replace('.tif', ''))  # Assuming the fourth part is 'toe location' and removing '.tif'
    return Q, thickness, phi, toe_location

# Apply the function to extract 'Q', 'thickness', 'phi', and 'toe location'
data[['Q', 'thickness', 'phi', 'toe_location']] = data['file'].apply(lambda x: pd.Series(extract_parts(x)))

# Save the new dataframe to a CSV file
output_file_path = r'C:\Users\Jason\Desktop\stochICE\examples\Chateauguay\MonteCarlo\SimulationTifs\Stats\largest_polygon_areas_formatted.csv'  # Replace with the desired output path
data.to_csv(output_file_path, index=False)

# Function to calculate p-values for the correlation matrix
def calculate_pvalues(df):
    df = df.dropna()
    pvalues = pd.DataFrame(np.ones((df.shape[1], df.shape[1])), columns=df.columns, index=df.columns)
    for row in df.columns:
        for col in df.columns:
            if row != col:
                _, pvalue = stats.pearsonr(df[row], df[col])
                pvalues.loc[row, col] = pvalue
    return pvalues

# Function to create and plot the correlation matrix with significance
def plot_correlation_matrix(data):
    # Select only the numerical columns (area, Q, thickness, phi, toe_location)
    numerical_data = data[['area', 'Q', 'thickness', 'phi', 'toe_location']]

    # Create the correlation matrix
    corr_matrix = numerical_data.corr()

    # Calculate p-values
    pvalues = calculate_pvalues(numerical_data)
    
    # Mask for the upper triangle (to hide repeated values)
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))

    # Plot the correlation matrix using a heatmap with French text
    plt.figure(figsize=(8, 6))
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', fmt='.2f', mask=mask, cbar=True, 
                annot_kws={"size": 12}, vmin=-1, vmax=1)

    # Add titles and labels in French
    plt.title("Matrice de corrélation", fontsize=16)
    plt.xticks(ticks=[0.5, 1.5, 2.5, 3.5, 4.5], labels=["aire d'inondation", "Q", "Épaisseur", "Phi", "position pied"], fontsize=12)
    plt.yticks(ticks=[0.5, 1.5, 2.5, 3.5, 4.5], labels=["aire d'inondation", "Q", "Épaisseur", "Phi", "position pied"], fontsize=12)
    plt.xlabel("Variables", fontsize=14)
    plt.ylabel("Variables", fontsize=14)
    
    # Show the plot
    plt.show()

    # Print the p-values for the correlations
    print("P-values for the correlations:")
    print(pvalues)

# Call the function to plot the correlation matrix
plot_correlation_matrix(data)
