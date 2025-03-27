
# -*- coding: utf-8 -*-
"""
Guide de l'utilisateur (Français)
--------------------------------
Ce script génère une enveloppe des niveaux d'eau (WSE) à partir de fichiers CSV
générés par HEC-RAS, en calculant les percentiles et en superposant le profil de
lit et des rectangles représentant des ponts ou seuils. Les rectangles sont
placés derrière les enveloppes pour une meilleure visibilité.

Paramètres à modifier :
- `directory` : Chemin vers les fichiers CSV de profils WSE.
- `bed_profile_path` : Chemin vers le fichier CSV du profil de lit.
- `pdf_export_path` : Chemin pour exporter le graphique PDF.
- `y_min`, `y_max` : Limites de l'axe y (en mètres).
- `shift` : Décalage du chaînage pour aligner avec le profil de lit.
- `rectangles_data` : Liste des données pour les rectangles (chaînage, haut, bas).

Sorties :
- Un graphique PDF des enveloppes WSE exporté à `pdf_export_path`.

Utilisation :
1. Ajustez les paramètres ci-dessus selon vos données.
2. Exécutez le script. Vérifiez le message de confirmation ou les erreurs.

User Guide (English)
-------------------
This script generates a water surface elevation (WSE) envelope from HEC-RAS-
generated CSV files, computing percentiles and overlaying the bed profile and
rectangles representing bridges or weirs. Rectangles are plotted behind the
envelopes for better visibility.

Parameters to modify:
- `directory`: Path to the WSE profile CSV files.
- `bed_profile_path`: Path to the bed profile CSV file.
- `pdf_export_path`: Path to export the PDF plot.
- `y_min`, `y_max`: Y-axis limits (in meters).
- `shift`: Chainage shift to align with the bed profile.
- `rectangles_data`: List of rectangle data (chainage, top, bottom).

Outputs:
- A PDF plot of WSE envelopes saved at `pdf_export_path`.

Usage:
1. Adjust the parameters above based on your data.
2. Run the script. Check the confirmation message or errors.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


# Paramètres spécifiés par l'utilisateur
# --------------------------------------
directory = (
    r"C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder"
    r"\StochICE_data_Chateauguay_final_Q_93a650_base_fixed_no"
    r"\Results\Individual_WSE_profiles"
)  # Chemin vers les fichiers CSV

bed_profile_path = (
    r"C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder"
    r"\Bed profile.csv"
)  # Chemin vers le profil de lit

pdf_export_path = (
    r"C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder"
    r"\StochICE_data_Chateauguay_final_Q_93a650_base_fixed_no"
    r"\Results\WSE_envelopes\wse_envelope.pdf"
)  # Chemin pour exporter le PDF

y_min = 15.0  # Limite inférieure de l'axe y (m)
y_max = 40.0  # Limite supérieure de l'axe y (m)
shift = 2377.30  # Décalage du chaînage (m)

rectangles_data = [  # Données des rectangles : chaînage, haut, bas (m)
    {"chainage": 7784, "top": 38.75, "bottom": 24.71},
    {"chainage": 7460, "top": 35.21, "bottom": 23.84},
    {"chainage": 6490, "top": 24.80, "bottom": 24.05},
    {"chainage": 4770, "top": 28.71, "bottom": 19.02},
    {"chainage": 2680, "top": 25.3, "bottom": 16.7},
]


# Lecture du profil de lit
# ------------------------
bed_df = pd.read_csv(bed_profile_path)
bed_df.iloc[:, 0] += shift  # Ajustement du chaînage
chainage_black_line = bed_df.iloc[:, 0]
elevation_black_line = bed_df.iloc[:, 1]


# Collecte des données WSE
# ------------------------
filter_str = "River 1"
files = [
    f for f in os.listdir(directory)
    if filter_str in f and f.endswith(".csv")
]
if not files:
    print(f"Aucun fichier correspondant trouvé dans {directory}.")
    exit()

dfs = []
for file in files:
    try:
        df = pd.read_csv(os.path.join(directory, file))
        dfs.append(df)
    except Exception as e:
        print(f"Erreur lors de la lecture de {file} : {e}")

if not dfs:
    print("Aucun fichier CSV valide à tracer.")
    exit()

chainage = dfs[0]["Chainage (m)"]
chainage_shifted = chainage + shift
wse_series = [df["wse (m)"] for df in dfs]
wse_df = pd.concat(wse_series, axis=1)


# Calcul des percentiles
# ----------------------
q10 = wse_df.quantile(0.10, axis=1)
q25 = wse_df.quantile(0.25, axis=1)
median_wse = wse_df.quantile(0.50, axis=1)
q75 = wse_df.quantile(0.75, axis=1)
q90 = wse_df.quantile(0.90, axis=1)
min_wse = wse_df.min(axis=1)
max_wse = wse_df.max(axis=1)


# Création du graphique
# ---------------------
plt.rcParams.update({"font.size": 8})
fig, ax = plt.subplots(figsize=(7, 3))

# Ajout des rectangles (en arrière-plan avec zorder=0)
for rect in rectangles_data:
    x_start = rect["chainage"]
    y_start = rect["bottom"]
    height = rect["top"] - rect["bottom"]
    rect_patch = Rectangle(
        (x_start, y_start),
        40,
        height,
        facecolor="gray",
        edgecolor="none",
        alpha=1,
        zorder=0,  # Derrière tout
    )
    ax.add_patch(rect_patch)

# Ajout des enveloppes WSE (zorder=1) et profil de lit (zorder=2)
plt.fill_between(
    chainage_shifted,
    min_wse,
    max_wse,
    color="blue",
    alpha=0.1,
    label="0 à 100 percentile",
    zorder=1,
)
plt.fill_between(
    chainage_shifted,
    q10,
    q90,
    color="blue",
    alpha=0.2,
    label="10e à 90e percentile",
    zorder=1,
)
plt.fill_between(
    chainage_shifted,
    q25,
    q75,
    color="blue",
    alpha=0.4,
    label="25e à 75e percentile",
    zorder=1,
)
plt.plot(
    chainage_shifted,
    median_wse,
    color="blue",
    linewidth=2,
    label="WSE médiane",
    zorder=1,
)
plt.plot(
    chainage_black_line,
    elevation_black_line,
    color="black",
    linewidth=1.5,
    zorder=2,
)

# Configuration du graphique
plt.xlabel("Chaînage (m)", fontsize=8)
plt.ylabel("Niveau d'eau (m)", fontsize=8)
plt.ylim(y_min, y_max)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.tight_layout(pad=0.1)

# Exportation du PDF
plt.savefig(pdf_export_path, format="pdf", bbox_inches="tight")
plt.show()
print(f"PDF exporté à : {pdf_export_path}")

