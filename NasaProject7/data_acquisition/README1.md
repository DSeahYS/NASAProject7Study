# Overview
## The Data Acquisition module is the foundation of the NIR GRB/GW Counterpart Pipeline. It is responsible for gathering astronomical data from various sources, monitoring real-time alerts, and preparing data for further analysis. This module provides a complete graphical user interface for astronomers to query, monitor, and process data related to gamma-ray bursts (GRBs) and gravitational wave (GW) events.

## Current Status (March 2025)
Hybrid Data Acquisition Approach
The system implements a hybrid strategy that:

Real-Time Queries: Attempts to query major astronomical archives like MAST, ESO, and CADC.
Simulated Data Fallback: Automatically generates scientifically valid simulated data when real observations are unavailable.
Alert Monitoring: Tracks GRB and GW alerts using infrastructure inspired by LVK alert systems.
Simulated Results: Currently relies on simulated results 95% of the time due to API failures or limited archive access.

## Features
### 1. Virtual Observatory Queries
Search astronomical archives for NIR observations using VO protocols.
Automatic fallback to simulated data when real queries fail.
Instrument-specific filters (J/H/K/Y bands).

### 2. GRB/GW Alert Monitoring
Real-time alert tracking for GRBs and GWs.
Simulated skymap generation for GW events.
Visibility calculations for observatories like LDT/PRIME.

### 3. Data Export
Save results as CSV or Excel files with automatic validation.
Standardized FITS headers for compatibility across instruments.

### 4. Visualization
Light curve generation and visualization.
Finding chart creation for transient detection.
Color-magnitude diagram plotting.

### 1. Query Tab

Enter object ID or coordinates (e.g., GRB121128A).
Select instruments and NIR bands.
Execute queries to retrieve results.

### 2. Alerts Tab

Start monitoring to track alerts in real-time.
Generate test alerts for development purposes.
View alert details including RA, Dec, and confidence levels.

### 3. Analysis Tab

Generate light curves, finding charts, and color-magnitude diagrams.
Analyze transient behavior in NIR bands.

### 4. Settings Tab

Set data directories and log levels.
Choose alert sources (GCN, SWIFT, Fermi, Simulation).
Save or reset settings as needed.

## Data Sources
Attempted Real Services:
MAST
ESO TAP
CADC
Swift Archive

## Simulation Frameworks:
GW event timelines from LVK alert infrastructure.
Kilonova models matching GRB observations (GRB 230307A kilonova).
Instrument characteristics matching THESEUS mission requirements.

## Validation & Performance
Latency Metrics:
Process	Current Performance	Goal (O4 Standards)
Data Acquisition	3–5 min	<1 min
Alert Processing	2–3 min	<30 sec
Simulation Fallback	Instant	N/A
Science Validation:
Validated against recent multi-messenger observations:

## GRB 211106A/211227A observations
Compatible with THESEUS mission requirements (ESA documentation).

## Expected Limitations:
Real archive queries may fail due to service changes or network issues.
Simulated data contains placeholder values marked as "SIMULATED."
GW skymaps use simplified HEALPix implementations.

## Recommended Workflow:
Use simulated data for pipeline development.
Attempt real queries for specific targets.
Validate against public GWTC catalogs.

## Development Roadmap
Current Phase: Data Acquisition
✅ Query interface with multiple archive support
✅ Alert processing system with coincidence detection
✅ Basic visualization capabilities

## Next Phase: Data Processing
⬜ Implement NIR data calibration routines
⬜ Create background subtraction algorithms
⬜ Develop astrometric solution validation

## Notes
Time-domain simulations use March 2025 as the reference epoch.
Coordinate systems match 7DT standards.
Export formats are compatible with A&A publication requirements.

# This project is part of ongoing research in multi-messenger astronomy. This is just a showcase for NASA. This project is not conducted under NASA Supervision

