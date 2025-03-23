import os
import requests
import numpy as np
import pandas as pd
import logging
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
import pyvo

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('data_acquisition.query')

def query_vo_archives(position: SkyCoord, radius: u.Quantity, bands=None, time_range=None):
    """
    Query Virtual Observatory archives for NIR observations.
    """
    logger.info(f"Querying VO archives at position {position.to_string('hmsdms')} with radius {radius}")
    
    # List of NIR-capable archive services with corrected URLs
    services = [
        ("MAST", "https://mast.stsci.edu/tap"),
        ("ESO", "http://archive.eso.org/tap_cat"),
        ("CADC", "https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/tap")
    ]
    
    all_results = []
    
    for name, url in services:
        try:
            logger.info(f"Querying {name} service at {url}")
            service = pyvo.dal.TAPService(url)
            
            # Use service-specific queries
            if name == "MAST":
                query = f"""
                SELECT * FROM dbo.ObsCore 
                WHERE CONTAINS(POINT('ICRS', s_ra, s_dec), 
                CIRCLE('ICRS', {position.ra.degree}, {position.dec.degree}, {radius.to(u.deg).value})) = 1
                """
                if bands:
                    query += " AND em_xel_waveband IN ('NIR', 'J', 'H', 'K')"
                
            elif name == "ESO":
                query = f"""
                SELECT * FROM raw WHERE 
                CONTAINS(POINT('ICRS', ra, dec), 
                CIRCLE('ICRS', {position.ra.degree}, {position.dec.degree}, {radius.to(u.deg).value})) = 1
                """
                if bands:
                    query += " AND (instrument='HAWK-I' OR instrument='SOFI' OR instrument='ISAAC')"
                
            elif name == "CADC":
                query = f"""
                SELECT * FROM caom2.Observation 
                WHERE CONTAINS(POINT('ICRS', ra, dec), 
                CIRCLE('ICRS', {position.ra.degree}, {position.dec.degree}, {radius.to(u.deg).value})) = 1
                """
                if bands:
                    query += " AND energy_bandpassName IN ('NIR', 'J', 'H', 'K')"
            
            # Execute the query
            results = service.search(query).to_table()
            
            if len(results) > 0:
                logger.info(f"Found {len(results)} results from {name}")
                results['service_name'] = name  # Add service name column
                all_results.append(results)
            else:
                logger.info(f"No results found from {name}")
                
        except Exception as e:
            logger.error(f"Error querying {name}: {str(e)}")
    
    # Combine results if any were found
    if all_results:
        try:
            combined_results = Table.vstack(all_results)
            return combined_results
        except:
            logger.warning("Could not combine tables with different columns. Returning first result set.")
            return all_results[0]
    else:
        logger.warning("No results found from any service")
        return Table()

def download_fits_data(query_results, destination_dir=None):
    """
    Download FITS files from query results.
    """
    if destination_dir is None:
        destination_dir = os.path.join(os.getcwd(), 'data')
    
    # Create destination directory if it doesn't exist
    Path(destination_dir).mkdir(parents=True, exist_ok=True)
    
    downloaded_files = []
    
    # Check if query_results is empty
    if len(query_results) == 0:
        logger.warning("No query results to download")
        return downloaded_files
    
    # Determine the access URL column name (differs between services)
    access_url_cols = ['access_url', 'accessURL', 'access', 'download_link', 'datalink']
    access_col = None
    for col in access_url_cols:
        if col in query_results.colnames:
            access_col = col
            break
    
    if access_col is None:
        logger.error("No access URL column found in query results")
        return downloaded_files
    
    for i, result in enumerate(query_results):
        try:
            access_url = result[access_col]
            
            # Get unique filename
            if 'obs_id' in result.colnames:
                obs_id = result['obs_id']
            else:
                obs_id = f'observation_{i}'
            
            filename = f"{obs_id}.fits"
            file_path = os.path.join(destination_dir, filename)
            
            logger.info(f"Downloading {access_url} to {file_path}")
            
            # Download the file
            response = requests.get(access_url, stream=True, timeout=30)
            response.raise_for_status()
            
            with open(file_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=1024):
                    if chunk:
                        f.write(chunk)
            
            downloaded_files.append(file_path)
            logger.info(f"Successfully downloaded {filename}")
                
        except Exception as e:
            logger.error(f"Error downloading file: {str(e)}")
    
    logger.info(f"Downloaded {len(downloaded_files)} files to {destination_dir}")
    return downloaded_files

def standardize_fits_headers(fits_file):
    """
    Standardize FITS headers from different instruments to a common format.
    """
    logger.info(f"Standardizing headers for {fits_file}")
    
    try:
        # Open the FITS file
        with fits.open(fits_file, mode='update') as hdul:
            primary_hdr = hdul[0].header
            
            # Create standardized keywords if they don't exist
            std_keywords = {
                'OBJECT': primary_hdr.get('OBJECT', 'UNKNOWN'),
                'TELESCOP': primary_hdr.get('TELESCOP', primary_hdr.get('TELESCOPE', 'UNKNOWN')),
                'INSTRUME': primary_hdr.get('INSTRUME', primary_hdr.get('INSTRUMENT', 'UNKNOWN')),
                'FILTER': primary_hdr.get('FILTER', primary_hdr.get('FILTNAME', 'UNKNOWN')),
                'DATE-OBS': primary_hdr.get('DATE-OBS', primary_hdr.get('DATE_OBS', 'UNKNOWN')),
                'EXPTIME': primary_hdr.get('EXPTIME', primary_hdr.get('EXPOSURE', 0.0)),
                'RA': primary_hdr.get('RA', primary_hdr.get('CRVAL1', 0.0)),
                'DEC': primary_hdr.get('DEC', primary_hdr.get('CRVAL2', 0.0)),
                'OBSERVER': primary_hdr.get('OBSERVER', 'UNKNOWN'),
                'AIRMASS': primary_hdr.get('AIRMASS', 1.0),
                'SEEING': primary_hdr.get('SEEING', primary_hdr.get('FWHM', 0.0)),
            }
            
            # Update the header with standardized keywords
            for key, value in std_keywords.items():
                primary_hdr[key] = value
            
            # Add a keyword to indicate this file has been standardized
            primary_hdr['STDPIPE'] = (True, 'Processed with NIR GRB/GW Pipeline')
            
            # Save the changes
            hdul.flush()
        
        logger.info(f"Successfully standardized headers for {fits_file}")
        return fits_file
    
    except Exception as e:
        logger.error(f"Error standardizing FITS headers: {str(e)}")
        return None

def extract_metadata(fits_file):
    """
    Extract observation metadata from FITS headers.
    """
    logger.info(f"Extracting metadata from {fits_file}")
    
    try:
        # Open the FITS file
        with fits.open(fits_file) as hdul:
            primary_hdr = hdul[0].header
            
            # Extract key metadata
            metadata = {
                'filename': os.path.basename(fits_file),
                'object': primary_hdr.get('OBJECT', 'UNKNOWN'),
                'telescope': primary_hdr.get('TELESCOP', 'UNKNOWN'),
                'instrument': primary_hdr.get('INSTRUME', 'UNKNOWN'),
                'filter': primary_hdr.get('FILTER', 'UNKNOWN'),
                'date_obs': primary_hdr.get('DATE-OBS', 'UNKNOWN'),
                'exptime': float(primary_hdr.get('EXPTIME', 0.0)),
                'ra': float(primary_hdr.get('RA', 0.0)),
                'dec': float(primary_hdr.get('DEC', 0.0)),
                'airmass': float(primary_hdr.get('AIRMASS', 1.0)),
                'seeing': float(primary_hdr.get('SEEING', 0.0)),
            }
            
            # Create a SkyCoord object for the position
            try:
                metadata['coords'] = SkyCoord(ra=metadata['ra'], dec=metadata['dec'], unit='deg')
            except:
                metadata['coords'] = None
                
            logger.info(f"Successfully extracted metadata from {fits_file}")
            return metadata
                
    except Exception as e:
        logger.error(f"Error extracting metadata: {str(e)}")
        return None

def query_instrument(instrument_name, target, time_range=None):
    """
    Query specific instrument archives (RIMAS, PRIME, etc).
    """
    logger.info(f"Querying {instrument_name} for target {target}")
    
    # Convert target to SkyCoord if it's a string (assumed to be an object name)
    if isinstance(target, str):
        try:
            from astroquery.simbad import Simbad
            result_table = Simbad.query_object(target)
            if result_table is None:
                logger.error(f"Could not resolve target name: {target}")
                return Table()
            
            ra = result_table['RA'][0]
            dec = result_table['DEC'][0]
            position = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
            logger.info(f"Resolved {target} to coordinates {position.to_string('hmsdms')}")
        except Exception as e:
            logger.error(f"Error resolving target name {target}: {str(e)}")
            return Table()
    else:
        position = target
    
    # Instrument-specific configurations
    instrument_configs = {
        'RIMAS': {
            'service': 'https://irsa.ipac.caltech.edu/TAP',
            'table': 'ldt.science_archive',
            'bands': ['Y', 'J', 'H', 'K']
        },
        'PRIME': {
            'service': 'http://saao.chpc.ac.za/tap',
            'table': 'prime.science_archive',
            'bands': ['Y', 'J', 'H', 'K']
        }
    }
    
    # Check if the instrument is supported
    if instrument_name not in instrument_configs:
        logger.error(f"Unsupported instrument: {instrument_name}")
        return Table()
    
    config = instrument_configs[instrument_name]
    
    try:
        # Connect to the service
        service = pyvo.dal.TAPService(config['service'])
        
        # Build the query
        query = f"""
        SELECT * FROM {config['table']}
        WHERE CONTAINS(POINT('ICRS', ra, dec), 
                     CIRCLE('ICRS', {position.ra.degree}, {position.dec.degree}, 0.1)) = 1
        """
        
        # Add time constraint if provided
        if time_range:
            start_time, end_time = time_range
            query += f" AND date_obs BETWEEN '{start_time}' AND '{end_time}'"
        
        # Execute the query
        results = service.search(query).to_table()
        
        logger.info(f"Found {len(results)} observations from {instrument_name}")
        return results
    
    except Exception as e:
        logger.error(f"Error querying {instrument_name}: {str(e)}")
        return Table()

def save_to_csv(data, output_dir="csv_data", filename="observations.csv"):
    """
    Save data to CSV file in specified directory
    """
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Convert to pandas DataFrame if it's an astropy Table
    if hasattr(data, 'to_pandas'):
        df = data.to_pandas()
    else:
        df = data
    
    # Save to CSV
    file_path = output_path / filename
    df.to_csv(file_path, index=False)
    logger.info(f"Saved data to {file_path}")
    
    return file_path

def save_to_excel(data, output_dir="excel_data", filename="observations.xlsx"):
    """
    Save data to Excel file with multiple sheets if needed
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Convert to pandas DataFrame if it's an astropy Table
    if hasattr(data, 'to_pandas'):
        df = data.to_pandas()
    else:
        df = data
    
    file_path = output_path / filename
    
    # Handle large datasets by splitting into multiple sheets
    if len(df) > 1000000:  # Excel row limit is 1,048,576
        with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
            for i in range(0, len(df), 1000000):
                sheet_name = f'Data_{i//1000000 + 1}'
                df.iloc[i:i+1000000].to_excel(writer, sheet_name=sheet_name, index=False)
    else:
        df.to_excel(file_path, index=False)
    
    logger.info(f"Saved data to {file_path}")
    return file_path

def generate_large_dataset(num_rows=10000, output_dir="large_data"):
    """
    Generate and save a large dataset with astronomical observations
    """
    logger.info(f"Generating large dataset with {num_rows} rows")
    
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Generate sample data with realistic astronomical parameters
    data = {
        'observation_id': [f"obs_{i:06d}" for i in range(num_rows)],
        'ra': np.random.uniform(0, 360, num_rows),
        'dec': np.random.uniform(-90, 90, num_rows),
        'magnitude': np.random.normal(18, 2, num_rows),
        'filter': np.random.choice(['J', 'H', 'K'], num_rows),
        'exposure_time': np.random.randint(30, 600, num_rows),
        'signal_to_noise': np.random.exponential(5, num_rows),
        'observation_date': pd.date_range(start='2025-01-01', periods=num_rows, freq='S')
    }
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    # Save to CSV
    csv_file = output_path / 'large_observations.csv'
    df.to_csv(csv_file, index=False)
    
    # Save to Excel (with chunking for large datasets)
    excel_file = output_path / 'large_observations.xlsx'
    save_to_excel(df, output_dir=str(output_path), filename='large_observations.xlsx')
    
    logger.info(f"Generated and saved large dataset to {output_path}")
    return df
