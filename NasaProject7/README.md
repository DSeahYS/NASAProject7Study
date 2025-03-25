# NASA Project 7 Study Pipeline
## A focused implementation of a Near-Infrared (NIR) Gamma-Ray Burst (GRB) and Gravitational Wave (GW) counterpart detection pipeline, developed to demonstrate skills for NASA's Astrophysics Postbaccalaureate Position #7.

## Navigation Guide
Start Here
data_acquisition/grb_loader.py

Purpose: Core functionality for loading GRB position data

Why First: Establishes the foundation for all subsequent analysis

examples/quickstart.py

Purpose: Simple end-to-end demonstration of the pipeline

Why Second: Shows how the components work together

## Core Modules (In Order)
### nir_processing/calibration.py

Purpose: NIR-specific image calibration routines

Why Important: Handles thermal background effects unique to NIR

### analysis/lightcurve.py

Purpose: Extracts and analyzes temporal evolution of sources

Why Important: Core scientific functionality for afterglow analysis

### visualization/skymap.py

Purpose: Visualizes GW probability maps with candidate counterparts

Why Important: Critical for rapid source identification

## Data Files
The data/ directory contains:

grb_sample.fits: Curated sample of GRBs with NIR counterparts

GW170817_skymap.fits: Reference probability map from landmark GW event

wise_*.fits: NIR observations from WISE

# Run quickstart example
python examples/quickstart.py
Learning Progression
Understanding Data Sources: Start with GRB positions and GW skymaps

NIR Processing: Learn infrared-specific calibration techniques

Scientific Analysis: Implement temporal and spectral analysis

Multi-messenger Integration: Combine GRB and GW data streams

Project Goals
This implementation demonstrates:

Proficiency with astronomical data processing

Understanding of multi-messenger source identification

Skills in NIR-specific analysis techniques

Familiarity with time-domain astronomy workflows

Note: This is an independent study project for skill demonstration and is not affiliated with NASA.
