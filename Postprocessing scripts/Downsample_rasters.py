
import os
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import Affine
import numpy as np

# === USER CONFIGURATION ===
input_folder = r"C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder\StochICE_data_Chateauguay_final_Q_93a650_base_fixed_no\Results\Individual_depth_tifs"
output_folder = r"C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder\StochICE_data_Chateauguay_final_Q_93a650_base_fixed_no\Results\Individual_depth_tifs_downsampled"
os.makedirs(output_folder, exist_ok=True)

# === Function to downsample by factor of 2
def downsample_raster(input_path, output_path):
    with rasterio.open(input_path) as src:
        scale_factor = 2
        new_height = src.height // scale_factor
        new_width = src.width // scale_factor

        # Compute new transform
        transform = src.transform
        new_transform = Affine(transform.a * scale_factor, transform.b, transform.c,
                                transform.d, transform.e * scale_factor, transform.f)

        profile = src.profile.copy()
        profile.update({
            'height': new_height,
            'width': new_width,
            'transform': new_transform,
            'compress': 'deflate',
            'tiled': True
        })

        with rasterio.open(output_path, 'w', **profile) as dst:
            for i in range(1, src.count + 1):
                data = src.read(i, out_shape=(new_height, new_width), resampling=Resampling.average)
                dst.write(data, i)

# === Process all GeoTIFFs in the folder
for filename in os.listdir(input_folder):
    if filename.lower().endswith(('.tif', '.tiff')):
        input_path = os.path.join(input_folder, filename)
        output_path = os.path.join(output_folder, filename)
        print(f"Downsampling {filename}...")
        downsample_raster(input_path, output_path)

print("All rasters downsampled and saved.")







import os
import time
import numpy as np
import rasterio
from rasterio.windows import Window
import csv

# Parameters to configure
input_folder = r"C:\Users\dugj2403\Desktop\Chateauguay_final_5mars2025\HECRAS_project_folder\StochICE_data_Chateauguay_final_Q_93a650_base_fixed_no\Results\Individual_depth_tifs_downsampled"
tile_size = 512
coordinates = [(285449.8742,5024607.9990)]  # list of (x, y) coordinates in projected CRS

# Gather input GeoTIFF files
file_list = [os.path.join(input_folder, f) for f in os.listdir(input_folder) 
              if f.lower().endswith(('.tif', '.tiff'))]
file_list.sort()
num_files = len(file_list)
if num_files == 0:
    raise FileNotFoundError("No GeoTIFF files found in the input folder.")

# Open all datasets
datasets = [rasterio.open(fp) for fp in file_list]
ref_profile = datasets[0].profile
height, width = ref_profile['height'], ref_profile['width']
crs = datasets[0].crs
transform = datasets[0].transform

# Validate that all rasters have identical shape and georeferencing
for ds in datasets[1:]:
    if ds.width != width or ds.height != height:
        raise ValueError("All input rasters must have the same dimensions.")
    if ds.crs != crs or ds.transform != transform:
        raise ValueError("All input rasters must have the same CRS and alignment.")

# Prepare output raster files for each percentile
percentiles = [25, 50, 75, 90]
out_profile = ref_profile.copy()
out_profile.update(count=1, dtype=np.float32)
output_files = {}
for p in percentiles:
    out_path = os.path.join(input_folder, f"percentile_{p}_{num_files}maps_new.tif")
    output_files[p] = rasterio.open(out_path, 'w', **out_profile)

# Process rasters in tiles
start_time = time.time()
tiles_y = (height + tile_size - 1) // tile_size
tiles_x = (width  + tile_size - 1) // tile_size
total_tiles = tiles_y * tiles_x
tile_count = 0

print(f"Computing percentiles for {num_files} maps in tiles of {tile_size}x{tile_size}...")
for row_off in range(0, height, tile_size):
    for col_off in range(0, width, tile_size):
        win_h = min(tile_size, height - row_off)
        win_w = min(tile_size, width - col_off)
        window = Window(col_off=col_off, row_off=row_off, width=win_w, height=win_h)
        # Read tile from all rasters
        tile_stack = []
        for ds in datasets:
            tile_data = ds.read(1, window=window)
            tile_stack.append(tile_data)
        tile_stack = np.array(tile_stack)  # shape: (num_files, win_h, win_w)
        # Compute percentiles along axis=0 (time dimension)
        p25 = np.percentile(tile_stack, 25, axis=0)
        p50 = np.percentile(tile_stack, 50, axis=0)
        p75 = np.percentile(tile_stack, 75, axis=0)
        p90 = np.percentile(tile_stack, 90, axis=0)
        # Write results to output files
        output_files[25].write(p25.astype(np.float32), 1, window=window)
        output_files[50].write(p50.astype(np.float32), 1, window=window)
        output_files[75].write(p75.astype(np.float32), 1, window=window)
        output_files[90].write(p90.astype(np.float32), 1, window=window)
        # Progress update
        tile_count += 1
        print(f"Tile {tile_count}/{total_tiles} done (row_off={row_off}, col_off={col_off})")

# Close output files
for dst in output_files.values():
    dst.close()
elapsed = time.time() - start_time
print(f"Percentile raster computation complete in {elapsed:.1f} seconds.")

# Extract values at specified coordinates
coord_headers = [f"{x}_{y}" for (x, y) in coordinates]
all_coords_data = []  # list of [val_at_coord1, val_at_coord2, ...] for each map
for i, ds in enumerate(datasets):
    row_vals = []
    for (x, y) in coordinates:
        try:
            row, col = ds.index(x, y)
        except Exception:
            row, col = None, None
        if row is None or col is None or row < 0 or row >= height or col < 0 or col >= width:
            # Coordinate is out of bounds
            val = np.nan
        else:
            # Read the pixel value at (row, col)
            pixel_window = Window(col_off=col, row_off=row, width=1, height=1)
            arr = ds.read(1, window=pixel_window)
            val = float(arr[0, 0]) if arr.size > 0 else np.nan
        row_vals.append(val)
    all_coords_data.append(row_vals)

# Compute percentile values for each coordinate across all maps
all_data_array = np.array(all_coords_data, dtype=float)  # shape: (num_files, num_coords)
coord_values = {25: [], 50: [], 75: [], 90: []}
for j in range(len(coordinates)):
    col_data = all_data_array[:, j]
    # Use nanpercentile to ignore missing data (nan) if any
    coord_values[25].append(float(np.nanpercentile(col_data, 25)))
    coord_values[50].append(float(np.nanpercentile(col_data, 50)))
    coord_values[75].append(float(np.nanpercentile(col_data, 75)))
    coord_values[90].append(float(np.nanpercentile(col_data, 90)))

# Write CSV files for each percentile
for p in percentiles:
    csv_path = os.path.join(input_folder, f"percentile_{p}_values_at_coordinates.csv")
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Header: optional "Map" label + coordinate headers
        writer.writerow(["Map"] + coord_headers)
        # One row per map
        for idx, row_vals in enumerate(all_coords_data):
            map_label = os.path.basename(file_list[idx])
            writer.writerow([map_label] + [row_vals[j] for j in range(len(coordinates))])
        # Final row: percentile value
        writer.writerow([f"P{p}"] + [coord_values[p][j] for j in range(len(coordinates))])
    print(f"Wrote {csv_path}")

# Close input datasets
for ds in datasets:
    ds.close()
    
print("Processing complete. Percentile GeoTIFFs and CSV files have been saved.")




























