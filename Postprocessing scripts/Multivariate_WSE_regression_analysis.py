
folder_path = r"C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder\StochICE_data_Chateauguay_final_Q_93a650_base_fixed_no\Results\Individual_WSE_profiles" # Update this to your directory


# -*- coding: utf-8 -*-
"""
User Guide
----------
This script performs a multivariate regression analysis of water surface elevation (WSE)
at specified chainages, using input variables (Q, thickness, phi, porosity, pied, tete)
extracted from HEC-RAS CSV filenames. For each chainage, it generates:
- A CSV file with the variables and WSE data.
- Scatter plots of WSE vs. each variable, saved as PDFs.
- A text file with the regression summary.

Parameters to modify:
- `folder_path`: Directory containing the WSE profile CSV files.
- `river_name`: River name in the CSV files (e.g., "River 1").
- `chainages`: List of chainages (in meters) to analyze.

Outputs:
- Folder `regression_analysis` within `folder_path`, with subfolders:
  - `chainage_<chainage>` for each chainage, containing:
    - CSV: `data_chainage_<chainage>.csv`.
    - PDFs: `WSE_vs_<variable>_chainage_<chainage>.pdf`.
    - Text: `regression_summary_chainage_<chainage>.txt`.

Requirements:
- Python libraries: pandas, numpy, statsmodels, matplotlib
"""

import os
import glob
import pandas as pd

import statsmodels.api as sm
import matplotlib.pyplot as plt

# User-Specified Parameters
# -------------------------
folder_path = r"C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder\StochICE_data_Chateauguay_final_Q_93a650_base_fixed_no\Results\Individual_WSE_profiles" # Update this to your directory
river_name = "River 1"  # Update this to your river name
chainages = [1266, 2446, 4479]  # Update this with your chainages (m)

# Create Main Output Directory
# ----------------------------
main_output_folder = os.path.join(folder_path, "regression_analysis")
os.makedirs(main_output_folder, exist_ok=True)

# Function to Parse Filename
# -------------------------
def parse_filename(filename, river_name="River 1"):
    """
    Parse a filename like 'River 1_WSE_117.46_0.55_56.0_0.31_1000.0_3000.0_5.csv'
    to extract Q, thickness, phi, porosity, pied, and tete.

    Args:
        filename (str): Path or name of the file.
        river_name (str): Expected river name.

    Returns:
        dict: Variables mapped to their float values.

    Raises:
        ValueError: If the filename format is invalid.
    """
    base_name = os.path.basename(filename)
    parts = base_name.split('_')
    if len(parts) != 9:
        raise ValueError(f"Filename '{base_name}' must have 9 parts, got {len(parts)}.")
    if parts[0] != river_name or parts[1] != "WSE":
        raise ValueError(f"Filename '{base_name}' must start with '{river_name}_WSE_'.")
    variables = ['Q', 'thickness', 'phi', 'porosity', 'pied', 'tete']
    try:
        values = [float(parts[i]) for i in range(2, 8)]
    except ValueError as e:
        raise ValueError(f"Error parsing numbers in '{base_name}': {e}")
    return dict(zip(variables, values))

# Process Each Chainage
# ---------------------
for chainage in chainages:
    print(f"\nProcessing chainage {chainage} m...")
    
    # Create Subfolder for This Chainage
    chainage_folder = os.path.join(main_output_folder, f"chainage_{chainage}")
    os.makedirs(chainage_folder, exist_ok=True)
    
    # Find CSV Files
    file_pattern = os.path.join(folder_path, f"{river_name}_WSE_*.csv")
    file_list = glob.glob(file_pattern)
    
    if not file_list:
        print(f"No files found for chainage {chainage}.")
        continue
    
    # Collect Data
    data = []
    for file in file_list:
        try:
            # Parse variables from filename
            params = parse_filename(file, river_name)
            # Read WSE from CSV
            df = pd.read_csv(file)
            wse_values = df.loc[df["Chainage (m)"] == chainage, "wse (m)"].values
            if len(wse_values) == 1:
                data.append({**params, 'WSE': wse_values[0]})
            else:
                print(f"Warning: Chainage {chainage} not found in {os.path.basename(file)}")
        except Exception as e:
            print(f"Error processing {os.path.basename(file)}: {e}")
    
    # Check Data Availability
    if not data:
        print(f"No valid data for chainage {chainage}.")
        continue
    
    # Create DataFrame
    df_data = pd.DataFrame(data)
    print(f"Collected {len(df_data)} simulations for chainage {chainage}.")
    
    # Save Data
    data_csv_path = os.path.join(chainage_folder, f"data_chainage_{chainage}.csv")
    df_data.to_csv(data_csv_path, index=False)
    print(f"Data saved to {data_csv_path}")
    
    # Plot WSE vs. Each Variable
    variables = ['Q', 'thickness', 'phi', 'porosity', 'pied', 'tete']
    for variable in variables:
        plt.figure()
        plt.scatter(df_data[variable], df_data['WSE'])
        plt.xlabel(variable)
        plt.ylabel('WSE (m)')
        plt.title(f'WSE vs. {variable} at Chainage {chainage} m')
        plt.grid(True)
        plot_filename = f'WSE_vs_{variable}_chainage_{chainage}.pdf'
        plot_path = os.path.join(chainage_folder, plot_filename)
        plt.savefig(plot_path, format='pdf')
        plt.close()
    
    # Perform Regression
    X = df_data[variables]
    X = sm.add_constant(X)  # Add intercept
    y = df_data['WSE']
    
    try:
        model = sm.OLS(y, X).fit()
        summary_path = os.path.join(chainage_folder, f"regression_summary_chainage_{chainage}.txt")
        with open(summary_path, 'w') as f:
            f.write(model.summary().as_text())
        print(f"Regression summary saved to {summary_path}")
    except Exception as e:
        print(f"Regression failed for chainage {chainage}: {e}")