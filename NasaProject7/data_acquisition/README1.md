Data Acquisition Module for NIR GRB/GW Counterpart Pipeline
Overview
The Data Acquisition module forms the foundation of the NIR GRB/GW Counterpart Pipeline, responsible for gathering astronomical data from various sources, monitoring real-time alerts, and preparing data for further analysis. This module consists of three main components:

data_main.py: The main entry point for data acquisition functionality

query.py: Handles querying astronomical archives and databases

alerts.py: Monitors and processes GCN notices and gravitational wave alerts

Features
Query Virtual Observatory (VO) archives for NIR observations

Monitor GCN (GRB Coordinates Network) notices for gamma-ray burst alerts

Receive and process LIGO/Virgo gravitational wave alerts

Check target visibility from observatories (e.g., Lowell Discovery Telescope)

Download and standardize FITS files

Extract metadata from astronomical observations

Save data in multiple formats (CSV, Excel)

Installation
Prerequisites
Python 3.8 or higher

Required Python packages:

pyvo

astropy

requests

pandas

matplotlib

tkinter (for GUI components)

Setup
bash
# Install required packages
pip install pyvo astropy requests pandas matplotlib

# If you're running on Linux, you may need to install tkinter separately
# For Ubuntu/Debian:
# sudo apt-get install python3-tk
Usage
Querying Virtual Observatory Archives
python
from data_acquisition.query import query_vo_archives
from astropy.coordinates import SkyCoord
import astropy.units as u

# Define target position
target_position = SkyCoord(ra=123.456, dec=78.90, unit="deg")

# Query for NIR observations within 10 arcminutes
results = query_vo_archives(
    position=target_position,
    radius=10*u.arcmin,
    bands=['J', 'H', 'K']
)

# Display results
print(f"Found {len(results)} observations")
Downloading FITS Data
python
from data_acquisition.query import download_fits_data

# Download the FITS files from query results
fits_files = download_fits_data(
    results,
    destination_dir="./data/fits_data"
)

# Process the downloaded files
for fits_file in fits_files:
    print(f"Downloaded: {fits_file}")
Monitoring GCN Notices
python
from data_acquisition.alerts import parse_gcn_notice, check_visibility

# Parse a GCN notice
notice = parse_gcn_notice("path/to/notice.xml")

# Check if the GRB is observable from your telescope
visibility = check_visibility(
    notice['ra'],
    notice['dec'],
    observatory="LDT"
)

if visibility['observable']:
    print(f"Target is observable for {visibility['hours_visible']} hours tonight")
    print(f"Transit time: {visibility['transit_time']}")
else:
    print(f"Target not observable: {visibility['reason']}")
Simulating LIGO/Virgo Alerts
python
from data_acquisition.alerts import retrieve_gw_alerts

# Generate simulated gravitational wave alerts
gw_alerts = retrieve_gw_alerts(simulated=True)

for alert in gw_alerts:
    print(f"GW Alert: {alert['id']}")
    print(f"Localization area: {alert['area']} square degrees")
    print(f"Distance: {alert['distance']} ± {alert['distance_error']} Mpc")
Saving Data
python
from data_acquisition.query import save_to_csv, save_to_excel
import pandas as pd

# Create a DataFrame with results
df = pd.DataFrame({
    'target': ['GRB 230101A', 'GRB 230215B'],
    'ra': [123.456, 234.567],
    'dec': [45.678, 56.789],
    'observation_date': ['2023-01-01', '2023-02-15'],
    'filter': ['J', 'K'],
    'exposure_time': [300, 600]
})

# Save to CSV
csv_path = save_to_csv(df, output_dir="results", filename="observations.csv")

# Save to Excel
excel_path = save_to_excel(df, output_dir="results", filename="observations.xlsx")
Module Components
data_main.py
The main execution script that coordinates the data acquisition process. It provides:

Command-line interface for initiating data queries

Graphical user interface for interactive use

Workflow coordination between query and alert modules

query.py
Handles interaction with astronomical databases and data processing:

query_vo_archives(): Searches Virtual Observatory services for observations

download_fits_data(): Downloads FITS files from query results

standardize_fits_headers(): Normalizes header information across different data sources

extract_metadata(): Extracts key information from FITS headers

save_data_to_csv() / save_to_csv(): Exports data to CSV format

save_to_excel(): Exports data to Excel format

alerts.py
Monitors and processes astronomical alert systems:

parse_gcn_notice(): Extracts information from GCN notices

check_visibility(): Determines if a target is observable from a specific location

retrieve_gw_alerts(): Gets gravitational wave alerts from LIGO/Virgo or simulates them

calculate_transit_time(): Computes when a target will transit at a given observatory

Troubleshooting
Connection Issues: If you encounter 404 errors when connecting to VO services, check your internet connection and verify the service endpoints are correct.

Missing Dependencies: Ensure all required Python packages are installed.

Import Errors: Verify your Python path includes the project directory.

FITS Processing Errors: The system will create placeholder files if download or processing fails.

Data Sources
The module connects to multiple astronomical data services:

MAST (Mikulski Archive for Space Telescopes)

ESO (European Southern Observatory)

CADC (Canadian Astronomy Data Centre)

GCN (GRB Coordinates Network)

LIGO/Virgo Gravitational Wave Network

Notes
This module is designed to work with both real and simulated data, allowing for testing and development when actual astronomical events are not available.