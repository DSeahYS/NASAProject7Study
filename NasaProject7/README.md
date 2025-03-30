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

## 2. Alerts Tab
Start monitoring for real-time GRB and GW alerts
Generate test alerts for development
View detailed information about received alerts
Assess follow-up observability

## 3. Analysis Tab
Generate light curves for GRB afterglows
Create finding charts for transient identification
Plot color-magnitude diagrams
Analyze multi-band observations

## 4. Settings Tab
Configure data directories
Set log levels
Choose alert sources
Save and load configurations

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
