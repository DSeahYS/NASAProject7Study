# NIR GRB/GW Counterpart Pipeline
A comprehensive Near-Infrared (NIR) pipeline for detecting and analyzing counterparts to Gamma-Ray Bursts (GRBs) and Gravitational Wave (GW) events
This multi-messenger astronomy tool combines data from various observatories and wavelength regimes to facilitate the study of high-energy transient phenomena.

# Project Overview
The NIR GRB/GW Counterpart Pipeline provides an integrated environment for:
Real-time monitoring of GRB and GW alerts
Astronomical archive queries with hybrid data acquisition (real/simulated)
NIR data processing and analysis
Visualization of light curves, finding charts, and color-magnitude diagrams
Spatial and temporal coincidence detection for multi-messenger events
Observatory visibility calculations and follow-up planning

# Core Implementation Files
## 1. Data Acquisition
data_acquisition/data_main.py - Main application with GUI interface
data_acquisition/query.py - Hybrid data acquisition system
data_acquisition/alerts.py - Multi-messenger alert processing

## 2. Configuration
config/alerts.json - Alert system configuration
config/settings.json - Global application settings

## 3. Outputs
outputs/alerts/ - Alert outputs and coincidence reports
data/fits/ - Downloaded FITS data
data/simulated/ - Simulated observation data

# Key Features
## Data Acquisition Module
Query archival NIR data from multiple sources (MAST, ESO, CADC)
Intelligent fallback to scientifically valid simulated data
GRB-specific data retrieval with coordinates lookup
FITS file standardization and metadata extraction

## Alert Processing System
Monitoring of GRB/GW event streams
O4-era gravitational wave alert handling
Spatial and temporal coincidence detection
Observatory visibility calculations for follow-up planning
Detailed reporting for coincident events

## Analysis & Visualization
Light curve generation and visualization
Finding chart creation
Color-magnitude diagram plotting
Interactive GUI for data exploration

# Application Interface
## 1. Query Tab
Enter object ID or coordinates (e.g., GRB121128A)
Select instruments and NIR bands
Execute queries to retrieve results
Download data for further analysis
![image](https://github.com/user-attachments/assets/bffcc0f9-6391-4961-982e-d7318152d61e)

## 2. Alerts Tab
Start monitoring for real-time GRB and GW alerts
Generate test alerts for development
View detailed information about received alerts
Assess follow-up observability
![image](https://github.com/user-attachments/assets/6f44cc85-39fb-4919-9d69-bd1437e102ee)

## 3. Analysis Tab
Generate light curves for GRB afterglows
Create finding charts for transient identification
Plot color-magnitude diagrams
Analyze multi-band observations
![image](https://github.com/user-attachments/assets/7d17ff06-9dbe-49c4-8725-6c5459b04b3b)
![image](https://github.com/user-attachments/assets/057c9099-16f1-4120-973d-b511176e5717)
![image](https://github.com/user-attachments/assets/e0479f5c-59b2-40c9-ba58-1b7ea434f665)

## 4. Settings Tab
Configure data directories
Set log levels
Choose alert sources
Save and load configurations
![image](https://github.com/user-attachments/assets/1902c7e1-ad12-4bf9-b7a0-b060730468bc)

# Technical Implementation
## Hybrid Data Acquisition Approach
Primary attempt to query astronomical archives (MAST, ESO, CADC)
Automatic fallback to scientifically valid simulated data
Standardized FITS headers across instruments
Support for various NIR instruments (UKIRT, VISTA, WFCAM, HAWK-I)

## Alert Processing Framework
Multi-messenger correlation with spatial-temporal coincidence detection
Real-time observatory visibility computation
Follow-up planning with prioritization algorithms
Skymap processing for gravitational wave events

## Visualization Engine
Matplotlib-based plotting framework
Interactive visualizations with Tkinter integration
Publication-quality output figures
Support for multiple visualization types

Required packages: astropy, matplotlib, tkinter, numpy, requests, pyvo

# Project Goals
## This implementation demonstrates:
Proficiency with astronomical data processing
Implementation of multi-messenger source identification algorithms
Skills in NIR-specific analysis techniques
Expertise in time-domain astronomy workflows
GUI development for scientific applications
Hybrid real/simulated data handling for robust pipeline testing

## Future Enhancements
Enhanced image differencing capabilities
Machine learning for transient classification
Direct telescope control for follow-up observations
Web-based interface for remote access
Integration with additional alert streams (Rubin, LSST)

# This project was developed as a demonstration of skills in astronomical software development, multi-messenger astronomy, and time-domain data analysis.
