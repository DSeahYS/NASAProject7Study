NIR GRB/GW Counterpart Pipeline
Project Overview
This pipeline is designed for processing near-infrared (NIR) observations to detect and analyze gamma-ray burst (GRB) afterglows and gravitational wave (GW) electromagnetic counterparts. It specifically aligns with the requirements of NASA's Project #7, which focuses on "Gamma-ray burst and Gravitational Wave counterpart searches in the near-infrared using ground-based methods."

The pipeline is optimized for data from instruments similar to RIMAS (Rapid infrared IMAger Spectrometer) at the Lowell Discovery Telescope and the PRIME (PRime-focus Infrared Microlensing Experiment) telescope in South Africa.

Installation
Requirements
Python 3.8 or higher

PyVO for Virtual Observatory access

Astropy for astronomical data handling

Photutils for source detection and photometry

Matplotlib for visualization

Jupyter for example notebooks

Setup
bash
# Clone the repository
git clone [repository-url]
cd NasaProject7

# Create and activate virtual environment (optional but recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install pyvo astropy photutils matplotlib jupyter
Pipeline Workflow
The pipeline follows this general workflow:

Data Acquisition: Query astronomical archives for NIR observations or monitor GRB/GW alerts

NIR Processing: Apply NIR-specific calibration and extract photometry

Analysis: Generate light curves and search for GRB/kilonova signatures

Visualization: Create plots and visualizations of results

Module Descriptions
Data Acquisition
query.py: Interface with Virtual Observatory services to find NIR observations

alerts.py: Monitor GCN notices and LIGO/Virgo alerts for new events

NIR Processing
calibration.py: NIR-specific image calibration, including thermal background correction

photometry.py: Source detection and photometric measurements optimized for NIR

Analysis
lightcurves.py: Extract and analyze light curves from time-series observations

counterparts.py: Identify and characterize potential GRB afterglows and kilonovae

Visualization
plots.py: Create publication-quality plots of light curves, color evolution, and sky maps

Example Usage
Basic Query for NIR Data
python
from data_acquisition import query
from astropy.coordinates import SkyCoord
import astropy.units as u

# Define GRB position
grb_position = SkyCoord(ra=10.625, dec=-20.88, unit="deg")

# Query for NIR observations within 5 arcminutes
results = query.query_vo_archives(position=grb_position, 
                                 radius=5*u.arcmin, 
                                 bands=['J', 'H', 'K'])

# Download the data
fits_files = query.download_fits_data(results, destination_dir="./data")
Process NIR Image
python
from nir_processing import calibration, photometry
from astropy.nddata import CCDData

# Load and calibrate image
raw_image = CCDData.read("raw_image.fits")
calibrated = calibration.calibrate_nir_image(raw_image)

# Remove sky background
sky_subtracted = calibration.subtract_sky_background(calibrated)

# Detect sources
sources = photometry.detect_sources(sky_subtracted, threshold=5.0)

# Perform photometry
photometry_results = photometry.perform_aperture_photometry(sky_subtracted, sources)
Analyze Light Curve
python
from analysis import lightcurves
import numpy as np

# Time series data
times = np.array([1.0, 2.0, 3.0, 4.0, 5.0])  # Days since trigger
fluxes = np.array([100.0, 60.0, 40.0, 30.0, 24.0])  # μJy
errors = np.array([5.0, 4.0, 3.0, 2.5, 2.0])

# Fit power-law decay
result = lightcurves.fit_powerlaw(times, fluxes, errors)
print(f"Decay index: {result['index']} ± {result['index_err']}")
Create Visualization
python
from visualization import plots
import matplotlib.pyplot as plt

# Create light curve plot
fig, ax = plt.subplots(figsize=(10, 6))
plots.plot_lightcurve(times, fluxes, errors, ax=ax)
plt.savefig("grb_lightcurve.png")
Examples
The examples/ directory contains Jupyter notebooks demonstrating the pipeline:

data_acquisition.ipynb: Querying for NIR observations and monitoring alerts

nir_analysis.ipynb: Processing NIR images and extracting photometry

full_pipeline.ipynb: End-to-end example from data acquisition to analysis

Development Roadmap
Phase 1: Framework and Data Acquisition
Implement PyVO queries to major NIR archives

Create GCN notice parser for real-time alerts

Develop FITS standardization tools

Phase 2: NIR Processing
Implement NIR-specific calibration routines

Create background modeling and subtraction

Develop photometry optimized for NIR observations

Phase 3: Analysis
Implement light curve extraction and fitting

Create GRB afterglow and kilonova models

Develop coincidence checking with GW events

Phase 4: Visualization and Examples
Create publication-quality plotting functions

Develop example notebooks

Complete full pipeline demonstration

Contributions
This project is designed to showcase skills relevant to NASA's Project #7 position. Contributions and improvements are welcome to enhance its functionality and demonstrate proficiency in NIR astronomy and multi-messenger analysis.