NIR GRB/GW Pipeline: Data Acquisition Module
Overview
The Data Acquisition module forms the foundation of the NIR GRB/GW Counterpart Pipeline, responsible for gathering astronomical data from various sources, monitoring real-time alerts, and preparing data for further analysis. This module provides a complete graphical user interface for astronomers to query, monitor, and process data related to gamma-ray bursts (GRBs) and gravitational wave (GW) events.

Features
Virtual Observatory Queries: Search astronomical archives for NIR observations

GRB/GW Alert Monitoring: Track and process GCN notices and LIGO/Virgo alerts

Data Management: Download, standardize, and export observation data

User-friendly Interface: Tabbed GUI with dedicated sections for different functions

Multi-format Export: Save results as CSV or Excel files

Customizable Settings: Configure directory paths, observatories, and logging

Installation
Prerequisites
Python 3.8 or higher

Required Python packages:

text
pip install pyvo astropy requests pandas matplotlib tkinter numpy
Application Interface
The application is divided into four main sections:

1. VO Queries Tab
Query Virtual Observatory archives for astronomical data:
![image](https://github.com/user-attachments/assets/9e81526f-44ef-4336-93de-cd3eb375a7fb)

Query Parameters:

Target Name field (e.g., "M31")

Coordinates (RA, Dec) input fields

Search Radius in arcminutes

NIR Band selection (J, H, K, Y checkboxes)

Control Buttons:

Resolve Name: Convert object name to coordinates

Execute Query: Send query to VO archives

Download Data: Retrieve FITS files for selected observations

Results Table:

Displays query results with ID, RA, DEC, INSTRUMENT, FILTER, DATE_OBS columns

Shows data from various instruments (RIMAS, PRIME, SOFI, HAWK-I, ISAAC)

Export buttons for CSV and Excel formats

2. GRB/GW Alerts Tab
Monitor and process astronomical alerts:
![image](https://github.com/user-attachments/assets/977f768b-2bae-4e48-a573-58f7116f8ecd)

Alert Controls:

GCN Monitoring: Start Monitor and Query Recent buttons

LIGO/Virgo Alerts: Get Alerts and Download Skymap buttons

Visibility Check: Determine if targets are observable

Recent Alerts:

Table showing alert time, type, source, coordinates, and significance

Displays both GW events (H1, L1, V1 sources) and GRB events (Swift-BAT)

Alert Details:

Detailed information panel for the selected alert

Shows event specifics like time, instruments, significance, and ID

3. Data Export Tab
Export and manage query results:
![image](https://github.com/user-attachments/assets/2b477156-d530-4579-8a44-0b55cba87f68)

Export Controls:

Format selection (CSV/Excel)

Dataset Size configuration

Output Directory selection with Browse button

Generate Large Dataset button for simulations

Export History:

Records of previous exports with timestamp, type, row count, and file path

4. Settings Tab
Configure application parameters:
![image](https://github.com/user-attachments/assets/31d659b3-0096-403b-88ab-6be246c49270)

General Settings:

Data Directory path configuration

Observatory Settings:

Default Observatory selection (e.g., LDT - Lowell Discovery Telescope)

Logging Settings:

Log Level selection (INFO, DEBUG, WARNING, ERROR)

Console Output:

Live log of application activities and operations

Usage Examples
Searching for NIR Observations
Navigate to the VO Queries tab

Enter a target name (e.g., "M31") or coordinates

Set search radius (e.g., 10 arcmin)

Select desired NIR bands (J, H, K, Y)

Click "Resolve Name" to convert target to coordinates

Click "Execute Query" to search archives

Review results in the table

Click "Download Data" to retrieve FITS files

Use "Export to CSV" or "Export to Excel" to save results

Monitoring GRB/GW Alerts
Navigate to the GRB/GW Alerts tab

Click "Start Monitor" to begin listening for GCN notices

Click "Get Alerts" to retrieve LIGO/Virgo alerts

Select an alert from the Recent Alerts table

View detailed information in the Alert Details panel

Click "Check Visibility" to determine if the target is observable

Click "Download Skymap" for gravitational wave events

Exporting and Managing Data
Navigate to the Data Export tab

Select export format (CSV or Excel)

Set dataset size or click "Generate Large Dataset" for simulations

Specify output directory

Track export history in the table below

Configuring Settings
Navigate to the Settings tab

Set data directory path

Select default observatory

Choose appropriate log level

Monitor operations in the Console Output

Click "Apply Settings" to save changes

Data Sources
The module connects to multiple astronomical data services:

MAST (Mikulski Archive for Space Telescopes)

ESO (European Southern Observatory)

CADC (Canadian Astronomy Data Centre)

GCN (GRB Coordinates Network)

LIGO/Virgo Gravitational Wave Network

Supported Instruments
The application is optimized for NIR instruments including:

RIMAS (Rapid infrared IMAger Spectrometer)

PRIME (PRime-focus Infrared Microlensing Experiment)

SOFI (Son OF ISAAC)

HAWK-I (High Acuity Wide-field K-band Imager)

ISAAC (Infrared Spectrometer And Array Camera)

Troubleshooting
No Results Found: The application will generate simulated data when real observations are not available

Connection Issues: Check internet connectivity or VPN settings if archive queries fail

File Download Errors: The system creates placeholder files when downloads fail

Alert Monitoring: Requires stable internet connection for real-time notifications

Notes
The current timestamp displayed in the application shows the real-time date and time

Sample data shown is simulated with future dates (2025-03-01) for testing purposes

The application supports both real and simulated data for development and testing
