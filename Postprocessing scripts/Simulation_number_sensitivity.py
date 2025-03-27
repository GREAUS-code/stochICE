
# -*- coding: utf-8 -*-
"""
Guide de l'utilisateur (Français)
--------------------------------
Ce script analyse les niveaux maximums cumulatifs d'eau (WSE) à partir de fichiers
CSV générés par HEC-RAS, en sélectionnant les simulations une par une de manière
aléatoire sans remplacement. Pour chaque chaînage spécifié, il produit :
- Un fichier CSV avec le nombre de simulations, le WSE maximum cumulatif et le
  nom du fichier correspondant.
- Un graphique PDF en escalier montrant l'évolution du WSE maximum cumulatif.

Paramètres à modifier :
- `folder_path` : Chemin vers les fichiers CSV de profils WSE.
- `river_name` : Nom de la rivière dans les fichiers CSV (ex. "River 1").
- `chainages` : Liste des chaînages (en mètres) à analyser.
- `y_limits` : Dictionnaire des limites y (min, max) pour chaque chaînage.

Sorties :
- Dossier `sensibilite_nombre_simulations` dans `folder_path`.
- Fichiers CSV : `max_wse_cumule_chainage_<chainage>.csv`.
- Graphiques PDF : `max_wse_cumule_vs_simulations_chainage_<chainage>.pdf`.

Utilisation :
1. Ajustez les paramètres ci-dessus selon vos données.
2. Exécutez le script. Vérifiez les messages pour les avertissements ou erreurs.

User Guide (English)
-------------------
This script analyzes cumulative maximum water surface elevations (WSE) from
HEC-RAS-generated CSV files, selecting simulations one by one randomly without
replacement. For each specified chainage, it generates:
- A CSV file with the number of simulations, cumulative maximum WSE, and the
  corresponding filename.
- A step-plot PDF showing the evolution of the cumulative maximum WSE.

Parameters to modify:
- `folder_path`: Path to the WSE profile CSV files.
- `river_name`: River name in the CSV files (e.g., "River 1").
- `chainages`: List of chainages (in meters) to analyze.
- `y_limits`: Dictionary of y-axis limits (min, max) for each chainage.

Outputs:
- Folder `sensibilite_nombre_simulations` within `folder_path`.
- CSV files: `max_wse_cumule_chainage_<chainage>.csv`.
- PDF plots: `max_wse_cumule_vs_simulations_chainage_<chainage>.pdf`.

Usage:
1. Adjust the parameters above based on your data.
2. Run the script. Check messages for warnings or errors.
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Paramètres spécifiés par l'utilisateur
# --------------------------------------
folder_path = (
    r"C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder"
    r"\StochICE_data_Chateauguay_final_Q_93a650_base_fixed_no"
    r"\Results\Individual_WSE_profiles"
)  # Chemin vers les fichiers CSV
river_name = "River 1"  # Nom de la rivière
chainages = [1266, 2446, 4479]  # Liste des chaînages à analyser (m)
y_limits = {  # Limites y pour chaque chaînage (m) : {chainage: (min, max)}
    1266: (20.0, 25.0),  # Limites pour chaînage 1266
    2446: (25.0, 30.0),  # Limites pour chaînage 2446
    4479: (25.0, 30.0),  # Limites pour chaînage 4479
}


# Création du dossier de sortie
# -----------------------------
sim_sensitivity_folder = os.path.join(folder_path, "sensibilite_nombre_simulations")
os.makedirs(sim_sensitivity_folder, exist_ok=True)


# Traitement des chaînages
# ------------------------
for chainage in chainages:
    # Vérification des limites y pour le chaînage
    if chainage not in y_limits:
        raise ValueError(
            f"Aucune limite y spécifiée pour le chaînage {chainage}. "
            "Ajoutez-la dans y_limits."
        )
    y_min, y_max = y_limits[chainage]

    # Recherche des fichiers CSV
    file_pattern = os.path.join(folder_path, f"{river_name}_WSE_*.csv")
    file_list = glob.glob(file_pattern)

    # Extraction des WSE pour le chaînage
    wse_values = []
    file_names = []
    for file in file_list:
        df = pd.read_csv(file)
        wse = df.loc[df["Chainage (m)"] == chainage, "wse (m)"].values
        if len(wse) == 1:
            wse_values.append(wse[0])
            file_names.append(os.path.basename(file))
        else:
            print(f"Avertissement : Chaînage {chainage} non trouvé dans {file}")

    # Vérification des données disponibles
    n_files = len(wse_values)
    if n_files == 0:
        print(f"Aucune donnée trouvée pour le chaînage {chainage}. Passage à la suivante.")
        continue

    print(f"{n_files} fichiers de simulation trouvés pour le chaînage {chainage} m.")

    # Mélange aléatoire des indices
    indices = list(range(n_files))
    np.random.shuffle(indices)

    # Initialisation des variables
    current_max = -np.inf
    max_file = None
    data = []

    # Calcul du maximum cumulatif
    for i in range(1, n_files + 1):
        selected_index = indices[i - 1]
        wse = wse_values[selected_index]
        file_name = file_names[selected_index]

        if wse > current_max:
            current_max = wse
            max_file = file_name

        data.append({
            "Nombre de simulations": i,
            "Maximum cumulatif WSE": current_max,
            "Nom du fichier": max_file,
        })

    # Création du DataFrame
    df = pd.DataFrame(data)

    # Enregistrement du CSV
    csv_filename = f"max_wse_cumule_chainage_{chainage}.csv"
    csv_path = os.path.join(sim_sensitivity_folder, csv_filename)
    df.to_csv(csv_path, index=False)
    print(f"Fichier CSV enregistré à {csv_path}")

    # Création du graphique
    plt.figure()
    plt.step(
        df["Nombre de simulations"],
        df["Maximum cumulatif WSE"],
        where="post",
        linewidth=0.5,
        color="black",
    )
    plt.xlabel("Nombre de simulations incluses")
    plt.ylabel("niveau maximum cumulatif (m)")
    plt.ylim(y_min, y_max)

    # Enregistrement du PDF
    plot_filename = f"max_wse_cumule_vs_simulations_chainage_{chainage}.pdf"
    plot_path = os.path.join(sim_sensitivity_folder, plot_filename)
    plt.savefig(plot_path, format="pdf")
    plt.close()
    print(f"Graphique PDF enregistré à {plot_path}")


# # -*- coding: utf-8 -*-
# """
# Script to analyze cumulative maximum WSE values as simulations are included one by one,
# include the corresponding filename, and generate both CSV files and PDF plots.
# """

# import os
# import glob
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

# # User-specified parameters (modify these as needed)
# folder_path = r'C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder\StochICE_data_Chateauguay_final_Q_93a650_base_fixed_no\Results\Individual_WSE_profiles'  # Path to CSV files
# river_name = 'River 1'                  # Name of the river
# chainages = [1266, 2446, 4479]         # List of chainages to analyze

# # Create a subfolder for outputs
# sim_sensitivity_folder = os.path.join(folder_path, 'Sim_nmb_sensitivity')
# os.makedirs(sim_sensitivity_folder, exist_ok=True)

# # Process each chainage
# for chainage in chainages:
#     # Find all CSV files for the specified river
#     file_pattern = os.path.join(folder_path, f"{river_name}_WSE_*.csv")
#     file_list = glob.glob(file_pattern)
    
#     # Extract WSE at the specified chainage from each file
#     wse_values = []
#     file_names = []  # Store filenames corresponding to WSE values
#     for file in file_list:
#         df = pd.read_csv(file)
#         wse = df.loc[df['Chainage (m)'] == chainage, 'wse (m)'].values
#         if len(wse) == 1:
#             wse_values.append(wse[0])
#             file_names.append(os.path.basename(file))  # Store just the filename
#         else:
#             print(f"Warning: Chainage {chainage} not found in {file}")
    
#     # Total number of files with valid data
#     N = len(wse_values)
#     if N == 0:
#         print(f"No data found for chainage {chainage}. Skipping.")
#         continue
    
#     print(f"Found {N} simulation files with data at chainage {chainage} m.")
    
#     # Shuffle the indices to select files randomly without replacement
#     indices = list(range(N))
#     np.random.shuffle(indices)
    
#     # Initialize variables
#     current_max = -np.inf
#     max_file = None
#     data = []
    
#     # Accumulate files one by one and track the cumulative maximum WSE
#     for i in range(1, N + 1):
#         selected_index = indices[i - 1]
#         wse = wse_values[selected_index]
#         file_name = file_names[selected_index]
        
#         if wse > current_max:
#             current_max = wse
#             max_file = file_name
        
#         # Record the current cumulative maximum and the file it came from
#         data.append({
#             'Number of Simulations': i,
#             'Cumulative Maximum WSE': current_max,
#             'File Name': max_file
#         })
    
#     # Create a DataFrame
#     df = pd.DataFrame(data)
    
#     # Save to CSV with the specified filename
#     csv_filename = f'cumulative_max_wse_chainage_{chainage}.csv'
#     csv_path = os.path.join(sim_sensitivity_folder, csv_filename)
#     df.to_csv(csv_path, index=False)
#     print(f"CSV saved to {csv_path}")
    
#     # Create and save the PDF plot
#     plt.figure()
#     plt.step(df['Number of Simulations'], df['Cumulative Maximum WSE'], where='post', linewidth=1,color='black')
#     plt.xlabel('Number of Simulations Included')
#     plt.ylabel('Cumulative Maximum WSE (m)')
#     plt.title(f'Cumulative Maximum WSE at Chainage {chainage} m for {river_name}')
    
    
#     # Save the plot as a PDF file
#     plot_filename = f'cumulative_max_wse_vs_simulations_chainage_{chainage}.pdf'
#     plot_path = os.path.join(sim_sensitivity_folder, plot_filename)
#     plt.savefig(plot_path, format='pdf')
#     plt.close()  # Close the figure to free memory
    
#     print(f"PDF plot saved to {plot_path}")



# # -*- coding: utf-8 -*-
# """
# Script to analyze maximum WSE values across different sample sizes, include the corresponding filename,
# and generate both CSV files and PDF plots.
# """

# import os
# import glob
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt

# # User-specified parameters (modify these as needed)
# folder_path = r'C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder\StochICE_data_Chateauguay_final_Q_93a650_base\Results\Individual_WSE_profiles'  # Path to CSV files
# river_name = 'River 1'                  # Name of the river
# chainages = [1266, 2446, 4479]                # List of chainages to analyze
# X = 5                                   # Sample size increment (e.g., 5, 10, 15, ...)

# # Create a subfolder for outputs
# sim_sensitivity_folder = os.path.join(folder_path, 'Sim_nmb_sensitivity')
# os.makedirs(sim_sensitivity_folder, exist_ok=True)

# # Process each chainage
# for chainage in chainages:
#     # Find all CSV files for the specified river
#     file_pattern = os.path.join(folder_path, f"{river_name}_WSE_*.csv")
#     file_list = glob.glob(file_pattern)
    
#     # Extract WSE at the specified chainage from each file
#     wse_values = []
#     file_names = []  # Store filenames corresponding to WSE values
#     for file in file_list:
#         df = pd.read_csv(file)
#         wse = df.loc[df['Chainage (m)'] == chainage, 'wse (m)'].values
#         if len(wse) == 1:
#             wse_values.append(wse[0])
#             file_names.append(os.path.basename(file))  # Store just the filename
#         else:
#             print(f"Warning: Chainage {chainage} not found in {file}")
    
#     # Total number of files with valid data
#     N = len(wse_values)
#     if N == 0:
#         print(f"No data found for chainage {chainage}. Skipping.")
#         continue
    
#     print(f"Found {N} simulation files with data at chainage {chainage} m.")
    
#     # Generate sample sizes: X, 2X, 3X, ..., up to N
#     sample_sizes = list(range(X, N + 1, X))
#     if sample_sizes[-1] < N:
#         sample_sizes.append(N)
    
#     # Store results for the CSV
#     data = []
    
#     # Calculate maximum WSE and corresponding filename for each sample size
#     for S in sample_sizes:
#         # Randomly select S files without replacement
#         selected_indices = np.random.choice(N, S, replace=False)
#         selected_wse = [wse_values[i] for i in selected_indices]
#         selected_files = [file_names[i] for i in selected_indices]
        
#         # Find the maximum WSE and its corresponding file
#         max_wse = max(selected_wse)
#         max_index = selected_wse.index(max_wse)
#         max_file = selected_files[max_index]
        
#         # Add to data list
#         data.append({
#             'Sample Size': S,
#             'Maximum WSE': max_wse,
#             'File Name': max_file
#         })
    
#     # Create a DataFrame
#     df = pd.DataFrame(data)
    
#     # Save to CSV with the specified filename
#     csv_filename = f'max_wse_data_chainage_{chainage}.csv'
#     csv_path = os.path.join(sim_sensitivity_folder, csv_filename)
#     df.to_csv(csv_path, index=False)
#     print(f"CSV saved to {csv_path}")
    
#     # Create and save the PDF plot
#     plt.figure()
#     plt.plot(df['Sample Size'], df['Maximum WSE'], marker='o', linewidth=2)
#     plt.xlabel('Sim. Nmb. sample size')
#     plt.ylabel('Maximum WSE (m)')
#     plt.title(f'Maximum WSE at Chainage {chainage} m for {river_name}')
#     plt.grid(True)
    
#     # Save the plot as a PDF file
#     plot_filename = f'max_wse_vs_sample_size_chainage_{chainage}.pdf'
#     plot_path = os.path.join(sim_sensitivity_folder, plot_filename)
#     plt.savefig(plot_path, format='pdf')
#     plt.close()  # Close the figure to free memory
    
#     print(f"PDF plot saved to {plot_path}")