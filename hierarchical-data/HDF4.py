# Import packages
import os
import warnings

import matplotlib.pyplot as plt
import numpy.ma as ma
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
from shapely.geometry import mapping, box
import earthpy as et
import earthpy.spatial as es
import earthpy.plot as ep

warnings.simplefilter("ignore")

# Get MODIS data
et.data.get_data("cold-springs-modis-h4")

# Get boundary data
et.data.get_data("cold-springs-fire")

# Set working directory
os.chdir(os.path.join(et.io.HOME, "earth-analytics", "data"))

# Create path to MODIS data
modis_pre_path = os.path.join("cold-springs-modis-h4", "07_july_2016", "MOD09GA.A2016189.h09v05.006.2016191073856.hdf")

# Open data 
modis_pre = rxr.open_rasterio(modis_pre_path, masked=True)

# Explore the data 
modis_pre_qa = modis_pre[0]
modis_pre_qa["num_observations_1km"]
modis_pre_bands = modis_pre[1]

# View all the subdatasets
import rasterio as rio
with rio.open(modis_pre_path) as groups:
    for name in groups.subdatasets:
        print(name)

# Process MODIS bands stored in a HDF4 file

# View the first band
ep.plot_bands(modis_pre_bands.sur_refl_b01_1)
plt.show()
#or
# ep.plot_bands(modis_pre_bands['surf_refl_b01_1'])
# plt.show()

# Graficar todas las bandas
ep.plot_bands(modis_pre_bands.to_array().values, figsize=(10, 6))
plt.show()

# Graficar RGB
ep.plot_rgb(modis_pre_bands.to_array().values, 
            rgb=[0,2,1], 
            title = 'MODIS RGB')
plt.show()

# Open for groups and clip

# List of desired bands to open
desired_bands = ['sur_refl_b01_1', 
               'sur_refl_b02_1', 
               'sur_refl_b03_1',
               'sur_refl_b04_1',
               'sur_refl_b07_1']

# Read in the fire boundary shapefile
fire_boundary_path = os.path.join("cold-springs-fire",
                                  "vector_layers",
                                  "fire-boundary-geomac",
                                  "co_cold_springs_20160711_2200_dd83.shp")
fire_boundary = gpd.read_file(fire_boundary_path)

# Reproject the fire boundary to match the MODIS data
if not fire_boundary.crs == modis_pre_bands.rio.crs:
    fire_bound_sin = fire_boundary.to_crs(modis_pre_bands.rio.crs)

# Create the crop box
crop_bound_box = [box(*fire_bound_sin.total_bounds)]

# Open and clip with the box
modis_bands_clip = rxr.open_rasterio(modis_pre_path, 
                                     masked = True, 
                                     variable = desired_bands).rio.clip(crop_bound_box,
                                                                        from_disk=True,
                                                                        all_touched=True).squeeze()

# Open and clip with boundary geometry
modis_pre_clip_geom = rxr.open_rasterio(modis_pre_path, 
                                        masked = True, 
                                        ariable = desired_bands),rio.clip(fire_bound_sin.geometry.apply(mapping),
                                                                          all_touched=True,
                                                                          from_disk=True).squeeze()

# View clipped data

from rasterio.plot import plotting_extent

modis_extent = plotting_extent(modis_bands_clip.to_array().values[0])

# With crop box
f, ax = plt.subplots()
ep.plot_bands(modis_bands_clip.to_array().values[0],
              ax=ax,
              extent=modis_extent,
              title='MODIS clipped with crop box')
fire_bound_sin.plot(ax=ax,
                    color='green')
plt.show()

# With fire boundary
f, ax = plt.subplots()
ep.plot_bands(modis_pre_clip_geom.to_array().values[0],
              ax=ax,
              extent=modis_extent,
              title='MODIS clipped with crop grometry')
fire_bound_sin.plot(ax=ax,
                    color='green')
plt.show()

# Join everything with a function

def open_clean_bands(band_path, crop_layer, variable=None):
    
    """Open, subset and crop a MODIS h4 file.

    Parameters
    -----------
    band_path : string 
        A path to the array to be opened.
    crop_layer : geopandas GeoDataFrame
        A geopandas dataframe to be used to crop the raster data using rioxarray clip().
    variable : List
        A list of variables to be opened from the raster.

    Returns
    -----------
    band : xarray DataArray
        Cropped xarray DataArray
    """
    
    crop_bound_box = [box(*crop_layer.total_bounds)]
    band = rxr.open_rasterio(band_path,
                             masked=True,
                             variable=variable).rio.clip(crop_bound_box,
                                                         from_disk=True,
                                                         all_touched=True).squeeze()
    return band

# Use the function to open and crop
clean_bands = open_clean_bands(modis_pre_path, fire_boundary, desired_bands)

# View the data
ep.plot_bands(clean_bands.to_array().values,
              title=desired_bands,
              figsize=(10, 6))
plt.show()

# Export the cropped data to GeoTIFF

stacked_file_path = os.path.join(os.path.dirname(modis_pre_path),
                                 'final_output',
                                 'modis_band_crop.tif')
if not os.path.exists(os.path.dirname(stacked_file_path)):
    os.makedirs(os.path.dirname(stacked_file_path))

clean_bands.rio.to_raster(stacked_file_path)