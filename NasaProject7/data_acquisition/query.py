# query.py with hybrid real/simulated approach

import os
import time
import json
import requests
from pathlib import Path
import logging
import pyvo
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
from astropy import units as u
import datetime
from tenacity import retry, stop_after_attempt, wait_exponential

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('data_acquisition.query')

def standardize_fits_headers(fits_file):
    """Standardize FITS headers from different instruments to a common format."""
    logger.info(f"Standardizing headers for {fits_file}")
    try:
        # Check if this is a real FITS file or a placeholder
        if not str(fits_file).endswith('.fits') or os.path.getsize(fits_file) < 2880:
            logger.warning(f"{fits_file} appears to be a placeholder, not a real FITS file")
            return fits_file

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
    """Extract observation metadata from FITS headers."""
    logger.info(f"Extracting metadata from {fits_file}")
    try:
        # Check if this is a placeholder text file
        if not str(fits_file).endswith('.fits') or os.path.getsize(fits_file) < 2880:
            # Parse the placeholder text file
            metadata = {
                'filename': os.path.basename(fits_file),
                'object': 'SIMULATED',
                'telescope': 'SIMULATED',
                'instrument': 'SIMULATED',
                'filter': 'SIMULATED',
                'date_obs': '2025-03-24T00:00:00',
                'exptime': 0.0,
                'ra': 0.0,
                'dec': 0.0,
                'airmass': 1.0,
                'seeing': 1.0,
            }
            
            # Try to extract values from the simulated file
            try:
                with open(fits_file, 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if 'RA:' in line:
                            metadata['ra'] = float(line.split('RA:')[1].strip())
                        if 'Dec:' in line:
                            metadata['dec'] = float(line.split('Dec:')[1].strip())
                        if 'Filter:' in line:
                            metadata['filter'] = line.split('Filter:')[1].strip()
                        if 'Instrument:' in line:
                            metadata['instrument'] = line.split('Instrument:')[1].strip()
                        if 'Date:' in line:
                            metadata['date_obs'] = line.split('Date:')[1].strip()
            except:
                pass
                
            logger.info(f"Extracted metadata from placeholder file {fits_file}")
            return metadata
            
        # Real FITS file processing
        with fits.open(fits_file) as hdul:
            primary_hdr = hdul[0].header
            
            # Extract basic metadata
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
            
            # Create a SkyCoord object for convenient coordinate transformations
            if metadata['ra'] != 0.0 or metadata['dec'] != 0.0:
                metadata['coords'] = SkyCoord(metadata['ra'], metadata['dec'], unit='deg')
                
            logger.info(f"Successfully extracted metadata from {fits_file}")
            return metadata
            
    except Exception as e:
        logger.error(f"Error extracting metadata: {str(e)}")
        return None

def save_data_to_csv(data, output_dir="csv_outputs", filename="observations.csv"):
    """
    Save data to CSV file in specified directory
    
    Parameters
    ----------
    data : pandas DataFrame or astropy Table
        Data to save
    output_dir : str, optional
        Directory to save CSV file
    filename : str, optional
        Name of CSV file
        
    Returns
    -------
    output_path : Path
        Path to saved CSV file
    """
    import pandas as pd
    from pathlib import Path
    
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

def validate_mast_request(request):
    """Validate MAST API request structure."""
    required_keys = ['service', 'params', 'format']
    if not all(k in request for k in required_keys):
        logger.warning("Invalid MAST request structure: missing required keys")
        return False
    
    required_params = ['coordinates', 'radius']
    if not all(k in request['params'] for k in required_params):
        logger.warning("Invalid MAST params: missing coordinates or radius")
        return False
        
    valid_services = ['Mast.Caom.Cone', 'Mast.Caom.Filtered']
    if request['service'] not in valid_services:
        logger.warning(f"Unsupported MAST service: {request['service']}")
        return False
        
    return True

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=4, max=10))
def mast_query_with_retry(request):
    """Execute MAST API query with retry logic."""
    logger.debug(f"Executing MAST query with retry: {json.dumps(request, indent=2)}")
    response = requests.post(
        'https://mast.stsci.edu/api/v0/invoke', 
        json=request,
        timeout=30
    )
    response.raise_for_status()
    return response

def query_vo_archives(position: SkyCoord, radius: u.Quantity, bands=None, time_range=None):
    """
    Query Virtual Observatory archives for NIR observations.
    Falls back to simulated data if real queries fail.
    """
    logger.info(f"Querying VO archives at position {position.to_string('hmsdms')} with radius {radius}")
    
    all_results = []
    success = False
    
    # Try MAST API first (using the corrected API approach)
    try:
        logger.info(f"Querying MAST service via API")
        
        # Convert SkyCoord to decimal degrees
        ra = position.ra.degree
        dec = position.dec.degree
        radius_deg = radius.to(u.deg).value
        
        # Build MAST API request with the correct format
        mast_request = {
            "request": {
                "service": "Mast.Caom.Cone",
                "params": {
                    "coordinates": f"{ra},{dec}",
                    "radius": f"{radius_deg} deg",
                    "pagesize": 2000,
                    "page": 1,
                    "removenullcolumns": True,
                    "removecache": True
                },
                "format": "json",
                "cachebreaker": str(int(datetime.datetime.now().timestamp()))
            }
        }
        
        # Add NIR filter constraint
        if bands:
            if "J" in bands:
                mast_request["request"]["params"]["filters"] = [
                    {"paramName": "wavelength_range", "values": [{"min": 1.0, "max": 1.4}]}
                ]
            elif "H" in bands:
                mast_request["request"]["params"]["filters"] = [
                    {"paramName": "wavelength_range", "values": [{"min": 1.4, "max": 1.8}]}
                ]
            elif "K" in bands:
                mast_request["request"]["params"]["filters"] = [
                    {"paramName": "wavelength_range", "values": [{"min": 1.8, "max": 2.5}]}
                ]
        
        # Log the request for debugging
        logger.debug(f"MAST API Request: {json.dumps(mast_request, indent=2)}")
        
        # Query MAST API
        try:
            response = requests.post(
                'https://mast.stsci.edu/api/v0/invoke', 
                json=mast_request,
                timeout=30
            )
            
            if response.status_code == 200:
                mast_data = response.json()
                if 'data' in mast_data and len(mast_data['data']) > 0:
                    # Convert to astropy table
                    mast_results = Table(rows=mast_data['data'], names=mast_data['fields'])
                    mast_results['service_name'] = 'MAST'
                    all_results.append(mast_results)
                    logger.info(f"Found {len(mast_results)} results from MAST")
                    success = True
                else:
                    logger.info("No results found from MAST")
            else:
                logger.warning(f"Error from MAST API: {response.status_code}, {response.text[:500]}")
        except Exception as e:
            logger.warning(f"MAST API request failed: {str(e)}")
    
    except Exception as e:
        logger.warning(f"Error in MAST query processing: {str(e)}")
    
    # Try ESO TAP service with improved query
    try:
        logger.info(f"Querying ESO service at http://archive.eso.org/tap_obs")
        service = pyvo.dal.TAPService("http://archive.eso.org/tap_obs")
        
        # Build ADQL query for ESO with correct table and column names
        # Focus specifically on NIR instruments and appropriate wavelength ranges
        eso_query = f"""
        SELECT * FROM ivoa.obscore
        WHERE 
        1=CONTAINS(POINT('ICRS', s_ra, s_dec),
                  CIRCLE('ICRS', {position.ra.degree}, {position.dec.degree}, {radius.to(u.deg).value}))
        AND dataproduct_type = 'image'
        AND (em_min BETWEEN 1.0 AND 2.5)  -- JHK bands (NIR range)
        AND instrument_name IN ('HAWK-I', 'ISAAC', 'SOFI', 'VIRCAM', 'KMOS')  -- NIR instruments
        """
        
        # Log the query for debugging
        logger.debug(f"ESO Query: {eso_query}")
        
        # Execute query with a reasonable timeout
        results = service.search(eso_query, maxrec=100, timeout=30)
        
        if results and len(results) > 0:
            results_table = results.to_table()
            logger.info(f"Found {len(results_table)} results from ESO")
            results_table['service_name'] = 'ESO'  # Add service name column
            all_results.append(results_table)
            success = True
        else:
            logger.info(f"No results found from ESO")
            
    except Exception as e:
        logger.warning(f"Error querying ESO: {str(e)}")
    
    # Try CADC service with proper async handling
    try:
        logger.info(f"Querying CADC service")
        
        # Use the correct CADC TAP service URL
        cadc_service = pyvo.dal.TAPService("https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/tap")
        
        # Build ADQL query for CADC with focus on NIR observations
        cadc_query = f"""
        SELECT * FROM caom2.Observation
        WHERE 1=CONTAINS(POINT('ICRS', position_bounds_center_ra, position_bounds_center_dec),
                         CIRCLE('ICRS', {position.ra.degree}, {position.dec.degree}, {radius.to(u.deg).value}))
        """
        
        # Add filter for NIR data if bands specified
        if bands:
            band_conditions = []
            for band in bands:
                band_conditions.append(f"energy_bandpassName = '{band}'")
            
            if band_conditions:
                cadc_query += " AND (" + " OR ".join(band_conditions) + ")"
        
        # Log the query for debugging
        logger.debug(f"CADC Query: {cadc_query}")
        logger.debug(f"CADC Service URL: {cadc_service.baseurl}")
        
        # Use sync query first (simpler, less prone to errors)
        try:
            results = cadc_service.search(cadc_query, maxrec=100)
            if results and len(results) > 0:
                results_table = results.to_table()
                logger.info(f"Found {len(results_table)} results from CADC")
                results_table['service_name'] = 'CADC'
                all_results.append(results_table)
                success = True
            else:
                logger.info(f"No results found from CADC sync query")
                
                # Fall back to async query if the sync query returns no results
                logger.info("Trying CADC async query")
                job = cadc_service.submit_job(cadc_query, maxrec=100)
                job.run()
                
                # Wait for the job to complete with a timeout
                timeout = 60  # seconds
                start_time = time.time()
                while job.phase not in ['COMPLETED', 'ERROR', 'ABORTED']:
                    if time.time() - start_time > timeout:
                        logger.warning("CADC async query timed out, aborting")
                        job.abort()
                        break
                    time.sleep(2)
                    
                if job.phase == 'COMPLETED':
                    results_table = job.fetch_result().to_table()
                    if len(results_table) > 0:
                        logger.info(f"Found {len(results_table)} results from CADC async")
                        results_table['service_name'] = 'CADC'
                        all_results.append(results_table)
                        success = True
                    else:
                        logger.info(f"No results found from CADC async query")
                else:
                    logger.warning(f"CADC async query failed: {job.phase}")
        except Exception as e:
            logger.warning(f"CADC sync query failed, error: {str(e)}")
            
    except Exception as e:
        logger.warning(f"Error setting up CADC query: {str(e)}")
    
    # If real queries failed or returned no results, fall back to simulated data
    if not success:
        logger.warning("No results from real archives. Falling back to simulated data.")
        simulated_results = _create_simulated_results(position, radius, bands)
        return simulated_results
    
    # Combine results if any were found from real queries
    try:
        combined_results = vstack(all_results)
        logger.info(f"Combined {len(combined_results)} total results from all services")
        return combined_results
    except Exception as e:
        logger.warning(f"Could not combine tables: {str(e)}. Returning first result set.")
        return all_results[0]

def _create_simulated_results(position, radius, bands=None):
    """Create simulated observation results."""
    logger.info("Generating simulated results")
    
    # Number of simulated observations
    n_results = 20
    
    # Create table with realistic columns
    results = Table()
    results['obs_id'] = [f'sim_obs_{i:04d}' for i in range(n_results)]
    results['ra'] = [position.ra.degree + (i*0.01 - 0.1) for i in range(n_results)]
    results['dec'] = [position.dec.degree + (i*0.01 - 0.1) for i in range(n_results)]
    
    # Random selection of realistic instruments
    instruments = ['RIMAS', 'PRIME', 'HAWK-I', 'SOFI', 'ISAAC']
    results['instrument'] = [instruments[i % len(instruments)] for i in range(n_results)]
    
    # Use requested bands or default to J, H, K
    if not bands:
        bands = ['J', 'H', 'K']
    results['filter'] = [bands[i % len(bands)] for i in range(n_results)]
    
    # Realistic exposure times (300-600s)
    import numpy as np
    results['exptime'] = np.random.uniform(300, 600, n_results)
    
    # Observation dates spread over the past month
    from datetime import datetime, timedelta
    now = datetime.now()
    results['date_obs'] = [(now - timedelta(days=i % 30)).strftime('%Y-%m-%dT%H:%M:%S') 
                           for i in range(n_results)]
    
    # Simulated access URLs
    results['access_url'] = [f'https://example.com/simulated/{i}.fits' for i in range(n_results)]
    results['service_name'] = ['SIMULATED'] * n_results
    
    logger.info(f"Created {n_results} simulated results")
    return results

def download_fits_data(query_results, destination_dir=None):
    """
    Download FITS files from query results with fallback to creating placeholder files.
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
    access_url_cols = ['access_url', 'accessURL', 'access', 'download_link', 'datalink', 
                       's_access_url', 'dataproduct_access_url', 'access_estsize']
    access_col = None
    for col in access_url_cols:
        if col in query_results.colnames:
            access_col = col
            break
    
    for i, result in enumerate(query_results):
        try:
            # Get unique filename
            if 'obs_id' in result.colnames:
                obs_id = result['obs_id']
            else:
                obs_id = f'observation_{i}'
            
            filename = f"{obs_id}.fits"
            file_path = os.path.join(destination_dir, filename)
            
            # Try to download if access URL is available
            if access_col and result[access_col] and not result[access_col].startswith('https://example.com'):
                try:
                    logger.info(f"Downloading {result[access_col]} to {file_path}")
                    response = requests.get(result[access_col], stream=True, timeout=30)
                    response.raise_for_status()
                    
                    with open(file_path, 'wb') as f:
                        for chunk in response.iter_content(chunk_size=1024):
                            if chunk:
                                f.write(chunk)
                    
                    downloaded_files.append(file_path)
                    logger.info(f"Successfully downloaded {filename}")
                    continue
                except requests.exceptions.RequestException as e:
                    logger.warning(f"Failed to download {result[access_col]}: {str(e)}")
                    # Will fall through to create placeholder
            
            # Create a placeholder FITS file
            logger.info(f"Creating placeholder FITS file for {obs_id}")
            
            # Get coordinates and metadata
            ra = result['ra'] if 'ra' in result.colnames else 0
            dec = result['dec'] if 'dec' in result.colnames else 0
            
            placeholder_path = create_placeholder_fits(
                obs_id=obs_id,
                ra=ra,
                dec=dec,
                filter_name=result.get('filter', 'J'),
                instrument=result.get('instrument', 'SIMULATED'),
                destination_dir=destination_dir
            )
            
            downloaded_files.append(placeholder_path)
        
        except Exception as e:
            logger.error(f"Error processing result {i}: {str(e)}")
    
    logger.info(f"Downloaded/created {len(downloaded_files)} files to {destination_dir}")
    return downloaded_files

def create_placeholder_fits(obs_id, ra, dec, filter_name='J', instrument='SIMULATED', destination_dir=None):
    """Create a placeholder FITS file with basic header info."""
    if destination_dir is None:
        destination_dir = os.path.join(os.getcwd(), 'data')
    
    # Create a simple 10x10 array of random data
    import numpy as np
    data = np.random.normal(100, 10, (10, 10))
    
    # Create a primary HDU
    hdu = fits.PrimaryHDU(data)
    
    # Add basic header info
    hdu.header['OBJECT'] = f'SIM_{obs_id}'
    hdu.header['RA'] = ra
    hdu.header['DEC'] = dec
    hdu.header['DATE-OBS'] = '2025-03-24T00:00:00'
    hdu.header['EXPTIME'] = 300.0
    hdu.header['INSTRUME'] = instrument
    hdu.header['TELESCOP'] = 'SIMULATED'
    hdu.header['FILTER'] = filter_name
    
    # Create HDU list
    hdul = fits.HDUList([hdu])
    
    # Save file
    filename = f"{obs_id}.fits"
    file_path = os.path.join(destination_dir, filename)
    hdul.writeto(file_path, overwrite=True)
    
    return file_path

def save_to_excel(data, output_dir="excel_outputs", filename="observations.xlsx"):
    """
    Save data to Excel file in specified directory
    
    Parameters
    ----------
    data : pandas DataFrame or astropy Table
        Data to save
    output_dir : str, optional
        Directory to save Excel file
    filename : str, optional
        Name of Excel file
        
    Returns
    -------
    output_path : Path
        Path to saved Excel file
    """
    import pandas as pd
    from pathlib import Path
    
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Convert to pandas DataFrame if it's an astropy Table
    if hasattr(data, 'to_pandas'):
        df = data.to_pandas()
    else:
        df = data
    
    # Save to Excel
    file_path = output_path / filename
    df.to_excel(file_path, index=False)
    logger.info(f"Saved data to {file_path}")
    return file_path

def save_to_csv(data, output_dir="csv_outputs", filename="observations.csv"):
    """Alias for save_data_to_csv for consistency with save_to_excel."""
    return save_data_to_csv(data, output_dir, filename)
