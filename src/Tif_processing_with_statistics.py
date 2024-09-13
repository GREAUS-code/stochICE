import os
import numpy as np
import pandas as pd
from osgeo import gdal, ogr, gdalconst
from qgis.core import (
    QgsRasterLayer, QgsVectorLayer, QgsField, QgsProject, QgsProcessingFeedback, QgsGeometry
)
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry
from qgis.PyQt.QtCore import QVariant
from qgis import processing

from qgis.analysis import QgsNativeAlgorithms
from processing.core.Processing import Processing




# Initialize processing algorithms
Processing.initialize()

# Define input and output directories
input_directory = 'C:/Users/Jason/Desktop/stochICE/examples/Chateauguay/MonteCarlo/SimulationTifs'
output_raster_directory = os.path.join(input_directory, 'ProcessedTifs')
output_shapefile_directory = os.path.join(input_directory, 'ProcessedShapefiles')
output_stats_directory = os.path.join(input_directory, 'Stats')  # New directory for stats (max pixel values)
output_areas_file = os.path.join(output_stats_directory, 'largest_polygon_areas.csv')  # File for polygon areas

# Create output directories if they do not exist
os.makedirs(output_raster_directory, exist_ok=True)
os.makedirs(output_shapefile_directory, exist_ok=True)
os.makedirs(output_stats_directory, exist_ok=True)

# Step 1: Find all .tif files in the specified directory
tif_files = [f for f in os.listdir(input_directory) if f.endswith('.tif')]

# Array to hold maximum pixel values across all rasters
max_array = None
count = 0

# List to store the areas of the largest polygons
polygon_areas = []

# Step 2: Loop through each .tif file
for tif_file in tif_files:
    print(f"Processing file: {tif_file}")
    
    # Load the raster layer
    input_raster_path = os.path.join(input_directory, tif_file)
    raster_layer = QgsRasterLayer(input_raster_path, tif_file)
    
    # Define output paths for raster and shapefile
    output_raster_path = os.path.join(output_raster_directory, f'processed_{tif_file}')
    output_shapefile_path = os.path.join(output_shapefile_directory, f'{os.path.splitext(tif_file)[0]}_largest_polygon.shp')
    intermediate_shapefile_path = os.path.join(output_shapefile_directory, f'{os.path.splitext(tif_file)[0]}_polygon_intermediate.shp')

    # Step 3: Prepare the raster for reclassification using QgsRasterCalculatorEntry
    entries = []
    ras = QgsRasterCalculatorEntry()
    ras.ref = 'ras@1'
    ras.raster = raster_layer
    ras.bandNumber = 1
    entries.append(ras)
    
    # Perform raster reclassification
    calc = QgsRasterCalculator('("ras@1" > -9999)*1', output_raster_path, 'GTiff', raster_layer.extent(), raster_layer.width(), raster_layer.height(), entries)
    calc.processCalculation()

    # Step 4: Polygonize the reclassified raster
    processing.run("gdal:polygonize", {
        'INPUT': output_raster_path,  
        'BAND': 1,
        'FIELD': 'DN',
        'EIGHT_CONNECTEDNESS': False,
        'EXTRA': '',
        'OUTPUT': intermediate_shapefile_path  
    })

    # Step 5: Find and keep only the largest polygon
    polygon_layer = QgsVectorLayer(intermediate_shapefile_path, 'Polygons', 'ogr')

    # Add a new field to store the area of each polygon
    polygon_layer.startEditing()
    polygon_layer.dataProvider().addAttributes([QgsField('area', QVariant.Double)])
    polygon_layer.updateFields()

    largest_area = 0
    largest_feature_id = None

    # Calculate area for each polygon and track the largest one
    for feature in polygon_layer.getFeatures():
        geometry = feature.geometry()
        area = geometry.area()
        polygon_layer.changeAttributeValue(feature.id(), feature.fieldNameIndex('area'), area)
        
        if area > largest_area:
            largest_area = area
            largest_feature_id = feature.id()

    # Record the area of the largest polygon
    polygon_areas.append({'file': tif_file, 'area': largest_area})

    # Remove all polygons except the largest one
    for feature in polygon_layer.getFeatures():
        if feature.id() != largest_feature_id:
            polygon_layer.deleteFeature(feature.id())

    polygon_layer.commitChanges()

    # Step 6: Save the largest polygon to a shapefile
    options = QgsVectorFileWriter.SaveVectorOptions()
    options.driverName = "ESRI Shapefile"
    QgsVectorFileWriter.writeAsVectorFormat(polygon_layer, output_shapefile_path, options)

    # Step 7: Mask the original raster using the largest polygon (Use gdal:cliprasterbymasklayer)
    masked_raster_path = os.path.join(output_raster_directory, f'masked_{tif_file}')
    
    processing.run("gdal:cliprasterbymasklayer", {
        'INPUT': input_raster_path,  # The original raster file
        'MASK': output_shapefile_path,  # The largest polygon shapefile
        'SOURCE_CRS': None,
        'TARGET_CRS': None,
        'TARGET_EXTENT': None,
        'NODATA': None,  # Can specify NODATA value here, e.g., -9999 if required
        'ALPHA_BAND': False,
        'CROP_TO_CUTLINE': False,  # Crop the raster to the mask
        'KEEP_RESOLUTION': True,  # Keep the original resolution
        'SET_RESOLUTION': False,
        'X_RESOLUTION': None,
        'Y_RESOLUTION': None,
        'MULTITHREADING': False,
        'OPTIONS': 'COMPRESS=DEFLATE',
        'DATA_TYPE': 0,  # Default data type
        'EXTRA': '',
        'OUTPUT': masked_raster_path  # Path to the masked output raster
    })

    # Step 8: Read the masked raster and update the maximum pixel values
    masked_ds = gdal.Open(masked_raster_path)
    masked_data = masked_ds.GetRasterBand(1).ReadAsArray()

    # Initialize max_array for the first time
    if max_array is None:
        max_array = np.full_like(masked_data, -np.inf, dtype=np.float64)  # Start with -inf to track maximums

    # Update max_array with maximum values, ignoring NoData (-9999)
    max_array = np.where(masked_data != -9999, np.maximum(max_array, masked_data), max_array)

    count += 1

    # Step 9: Delete the uncompressed 'processed_' TIFF file
    if os.path.exists(output_raster_path):
        os.remove(output_raster_path)
        print(f"Deleted temporary file: {output_raster_path}")

# Step 10: Save the maximum value raster
if count > 0:
    max_raster_path = os.path.join(output_stats_directory, 'max_raster.tif')

    # Use the first raster's metadata as a reference for output
    reference_raster = gdal.Open(os.path.join(input_directory, tif_files[0]))
    geo_transform = reference_raster.GetGeoTransform()
    projection = reference_raster.GetProjection()
    x_size = reference_raster.RasterXSize
    y_size = reference_raster.RasterYSize

    # Save maximum raster
    driver = gdal.GetDriverByName('GTiff')
    max_ds = driver.Create(max_raster_path, x_size, y_size, 1, gdal.GDT_Float32, options=['COMPRESS=DEFLATE'])
    max_ds.SetGeoTransform(geo_transform)
    max_ds.SetProjection(projection)
    max_band = max_ds.GetRasterBand(1)
    max_band.WriteArray(max_array)
    max_band.FlushCache()
    max_ds = None  # Close the file

    print(f"Maximum value raster saved at {max_raster_path}")

# Step 11: Export polygon areas to CSV
polygon_df = pd.DataFrame(polygon_areas)
polygon_df.to_csv(output_areas_file, index=False)
print(f"Polygon areas saved to {output_areas_file}")


