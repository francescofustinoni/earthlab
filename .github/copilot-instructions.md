## Purpose
This project contains exploratory geospatial analysis (mostly Jupyter notebooks) focused on Earth-observation workflows (NDVI, trend analysis) using Google Earth Engine (GEE), geemap, rasterio and geopandas. These instructions give agents targeted, actionable guidance for making safe, correct edits.

## Quick start (developer environment — Windows PowerShell)
- Create and activate a venv and install deps:
  - python -m venv .venv
  - .\\.venv\\Scripts\\Activate.ps1
  - pip install -r requirements.txt
- Notebooks are the primary artifacts; open them in VS Code or Jupyter/Lab.

## Key integration points & external deps
- Google Earth Engine Python API (`ee`) and `geemap` — many notebooks call `ee.Authenticate()` and `ee.Initialize()`; expect interactive browser auth when running locally.
- Raster IO: `rasterio` for reading/writing GeoTIFFs.
- Vector IO: `geopandas` for shapefiles (see `Ecuador/data/ec_shp/ec.shp`).
- Statistical: `pymannkendall` is used for trend tests.

## Project structure & important files (examples)
- `Ecuador/BolivarProvince/VegetationChangeDetection.ipynb` — canonical example of patterns: GEE image collections, masking functions, NDVI workflows, exports to local paths under `Ecuador/data/Bolivar`.
- `data/` — contains inputs and outputs (shapefiles, LAZ lidar, and exported GeoTIFFs). Many code paths use absolute Windows paths; prefer converting these to workspace-relative paths during edits.

## Conventions and patterns agents should follow
- Notebooks use functions for masking and NDVI (e.g., `mask_landsat_clouds`, `mask_sentinel2_clouds`, `add_ndvi`). When refactoring, keep function names and signatures stable unless all call sites are updated.
- When creating or changing exports, prefer smaller regions or export-to-Drive/Assets because GEE will block exports >50MB (noted in notebook comments).
- File naming patterns: `medianNDVIBolivar{YEAR}.tif` and `yearly_mean_ndvi_bolivar_{START}_{END}.csv`. Follow existing naming for continuity.
- Use `geemap.shp_to_ee(out_path)` to convert local shapefiles to EE objects; keep `clip(region)` usage to limit processing.

## Code examples agents may use or modify (follow these idioms)
- Initialize GEE (interactive):
  - try:
      ee.Initialize()
    except Exception:
      ee.Authenticate(); ee.Initialize()
- Export small raster via geemap:
  - geemap.ee_export_image(img, filename=out_tif, scale=30, region=ee_region.geometry(), file_per_band=False)

## Export & performance notes
- Large AOIs may cause GEE export limits. Alternatives: increase `scale`, clip the AOI to subregions, export to Google Drive or GEE Asset instead of local download.
- Local runs that call `getInfo()` or `.get(...).getInfo()` are blocking and can be slow. Use `.reduceRegion(..., bestEffort=True)` carefully and prefer asynchronous exports for large reductions.

## Safety & operational constraints
- Do not attempt to embed or exfiltrate API keys or credentials into the repo. GEE auth is interactive — prompt the user before automating.
- Avoid long-running exports without user confirmation. Notify the user and explain options (drive/asset, tiling, scale changes).

## Tests, formatting, and commits
- There is no automated test suite in the repo; changes should be validated by running the relevant notebook and a small local export/read cycle (e.g., open exported GeoTIFF with rasterio).
- If you add new Python dependencies, update `requirements.txt`.

## When in doubt
- Prefer minimal, reversible edits: make small changes, run the notebook cells locally, and ask the user before performing big exports or reorganizing data folders.

If anything above is unclear or you want more examples from a specific notebook, tell me which area to expand.
