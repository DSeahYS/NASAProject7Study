# NASA Project 7 Study Pipeline
## A focused implementation of a Near-Infrared (NIR) Gamma-Ray Burst (GRB) and Gravitational Wave (GW) counterpart detection pipeline, developed to demonstrate skills for NASA's Astrophysics Postbaccalaureate Position #7.

# Core Implementation Files
## 1. Configuration and Setup
requirements.txt - Dependencies for the project

setup.py - Package installation configuration

config/pipeline_config.yaml - Global settings for the pipeline

## 2. Data Acquisition
data_acquisition/grb_loader.py - Load GRB catalog positions

data_acquisition/wise_interface.py - Query WISE NIR observations

data_acquisition/gw_skymap.py - Handle gravitational wave localization maps

data_acquisition/crossmatch.py - Match GRBs with NIR observations

## 3. NIR Processing
nir_processing/calibration.py - NIR-specific image calibration

nir_processing/background.py - Thermal background modeling and subtraction

nir_processing/registration.py - Image alignment tools

nir_processing/photometry.py - Source detection and flux measurement

## 4. Analysis
analysis/lightcurve.py - Extract and analyze light curves

analysis/models/afterglow.py - GRB afterglow models

analysis/models/kilonova.py - Kilonova emission models

analysis/fitting.py - Parameter estimation for transient models

analysis/classification.py - Distinguish different transient types

## 5. Visualization
visualization/lightcurve_plot.py - Create light curve visualizations

visualization/skymap.py - Plot GW probability contours

visualization/multi_wavelength.py - Compare NIR with other bands

# Project Goals
This implementation demonstrates:

Proficiency with astronomical data processing

Understanding of multi-messenger source identification

Skills in NIR-specific analysis techniques

Familiarity with time-domain astronomy workflows

Note: This is an independent study project for skill demonstration and is not affiliated with NASA.
