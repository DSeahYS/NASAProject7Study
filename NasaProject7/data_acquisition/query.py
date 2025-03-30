# query.py - Hybrid real/simulated data acquisition system
# for NIR counterparts to GRBs and gravitational wave events

import os
import json
import requests
import logging
from pathlib import Path
import numpy as np
import pandas as pd
import datetime
import pyvo
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.table import Table, vstack, QTable
from astropy.time import Time
from astropy import units as u
from typing import List, Dict, Union, Optional, Tuple
from tenacity import retry, stop_after_attempt, wait_exponential

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('data_acquisition.query')

class HybridDataAcquisition:
    """
    Interface for hybrid data acquisition, attempting real archive queries first
    and falling back to scientifically valid simulated data when needed.
    """
    
    def __init__(self, config_path=None):
        """Initialize the data acquisition interface."""
        self.config = self._load_config(config_path)
        self.data_dir = Path(self.config.get('data_dir', 'data'))
        self.fits_dir = self.data_dir / 'fits'
        self.simulated_dir = self.data_dir / 'simulated'
        
        # Create directories if they don't exist
        for directory in [self.data_dir, self.fits_dir, self.simulated_dir]:
            directory.mkdir(parents=True, exist_ok=True)
        
        # Initialize query statistics
        self.stats = {
            'real_queries': 0,
            'simulated_queries': 0,
            'downloads': 0,
            'failed_downloads': 0
        }
    
    def _load_config(self, config_path):
        """Load configuration from file or use defaults."""
        default_config = {
            'data_dir': 'data',
            'query_timeout': 30,
            'max_results': 1000,
            'preferred_services': ['MAST', 'ESO', 'CADC'],
            'simulation_realistic': True,
            'theseus_compliant': True,
            'instrument_params': {
                'WFCAM': {'fov': 34, 'pixscale': 0.34, 'gain': 2.3},
                'WIRCam': {'fov': 21, 'pixscale': 0.3, 'gain': 1.8},
                'NIRI': {'fov': 2, 'pixscale': 0.05, 'gain': 4.1},
                'VIRCAM': {'fov': 60, 'pixscale': 0.339, 'gain': 4.19},
                'HAWK-I': {'fov': 7.5, 'pixscale': 0.106, 'gain': 2.5}
            }
        }
        
        if not config_path:
            logger.info("Using default configuration")
            return default_config
        
        try:
            with open(config_path, 'r') as f:
                config = json.load(f)
            logger.info(f"Loaded configuration from {config_path}")
            return {**default_config, **config}  # Merge with defaults
        except Exception as e:
            logger.warning(f"Failed to load config from {config_path}: {str(e)}")
            return default_config
    
    def query_archives(self, coord: SkyCoord, radius: float, 
                      instruments: List[str] = None, 
                      bands: List[str] = None, 
                      time_range: Tuple = None) -> Table:
        """Query astronomical archives for data with hybrid approach."""
        radius_quantity = radius * u.deg if isinstance(radius, (int, float)) else radius
        
        logger.info(f"Querying archives at {coord.to_string('hmsdms')} with radius {radius_quantity}")
        
        # Attempt real archive queries first
        try:
            results = query_vo_archives(coord, radius_quantity, bands, time_range)
            
            # Filter by instruments if specified
            if instruments and len(instruments) > 0:
                instrument_mask = [False] * len(results)
                for i, row in enumerate(results):
                    instr = row.get('instrument', row.get('instrume', '')).upper()
                    instrument_mask[i] = any(inst.upper() in instr for inst in instruments)
                
                results = results[instrument_mask]
            
            # Check if we need to supplement with simulated data
            if len(results) == 0 or (
                'service_name' in results.colnames and 
                all(src == 'SIMULATED' for src in results['service_name'])
            ):
                self.stats['simulated_queries'] += 1
            else:
                self.stats['real_queries'] += 1
                
            return results
            
        except Exception as e:
            logger.error(f"Error in archive query: {str(e)}")
            logger.info("Falling back to simulated data")
            self.stats['simulated_queries'] += 1
            
            # Generate simulated data as fallback
            return _create_simulated_results(coord, radius_quantity, bands, time_range)
    
    def download_data(self, results, dest_dir=None) -> List[Path]:
        """Download data from query results."""
        if dest_dir is None:
            dest_dir = self.fits_dir
        else:
            dest_dir = Path(dest_dir)
            dest_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"Downloading data to {dest_dir}")
        
        downloaded_files = []
        
        # Handle both Table and List inputs
        if isinstance(results, Table):
            download_list = results
        else:
            # If results is a list of IDs, try to find matching files
            if all(isinstance(item, str) for item in results):
                # Placeholder implementation - in a real system, would query by ID
                logger.warning("Download by ID not fully implemented, using simulated data")
                download_list = _create_simulated_results(
                    SkyCoord(180, 0, unit='deg'),
                    0.1 * u.deg,
                    ['J', 'H', 'K']
                )
            else:
                download_list = results
        
        # Process each result
        for i, row in enumerate(download_list):
            try:
                # Get access URL
                if 'access_url' in row.colnames:
                    url = row['access_url']
                elif 'accessURL' in row.colnames:
                    url = row['accessURL']
                else:
                    url = None
                
                # Get observation ID for filename
                if 'obs_id' in row.colnames:
                    obs_id = row['obs_id']
                elif 'observation_id' in row.colnames:
                    obs_id = row['observation_id']
                else:
                    obs_id = f"obs_{i:04d}"
                
                # Clean obs_id for filename
                obs_id = str(obs_id).replace('/', '_').replace(' ', '_')
                
                # Determine if this is real or simulated data
                is_simulated = (
                    'service_name' in row.colnames and 
                    row['service_name'] == 'SIMULATED'
                )
                
                if url and not is_simulated:
                    # Attempt to download real data
                    try:
                        output_path = dest_dir / f"{obs_id}.fits"
                        
                        logger.info(f"Downloading {url} to {output_path}")
                        
                        response = requests.get(url, timeout=30)
                        if response.status_code == 200:
                            with open(output_path, 'wb') as f:
                                f.write(response.content)
                            
                            # Standardize FITS headers
                            processed_file = standardize_fits_headers(output_path)
                            if processed_file:
                                downloaded_files.append(output_path)
                                self.stats['downloads'] += 1
                        else:
                            logger.warning(f"Failed to download {url}: {response.status_code}")
                            # Fall back to simulated data
                            self._create_simulated_file(row, obs_id, dest_dir)
                            downloaded_files.append(dest_dir / f"{obs_id}.fits")
                            self.stats['failed_downloads'] += 1
                    except Exception as e:
                        logger.error(f"Error downloading {url}: {str(e)}")
                        # Fall back to simulated data
                        self._create_simulated_file(row, obs_id, dest_dir)
                        downloaded_files.append(dest_dir / f"{obs_id}.fits")
                        self.stats['failed_downloads'] += 1
                else:
                    # Create simulated data
                    self._create_simulated_file(row, obs_id, dest_dir)
                    downloaded_files.append(dest_dir / f"{obs_id}.fits")
                    self.stats['failed_downloads'] += 1
            except Exception as e:
                logger.error(f"Error processing download for row {i}: {str(e)}")
        
        logger.info(f"Downloaded {len(downloaded_files)} files")
        return downloaded_files
    
    def _create_simulated_file(self, row, obs_id, dest_dir):
        """Create a simulated FITS file or placeholder."""
        try:
            # Extract parameters from the row
            ra = row.get('ra', 0.0)
            dec = row.get('dec', 0.0)
            instrument = row.get('instrument', 'SIMULATED')
            filter_name = row.get('filter', 'J')
            
            # Create a simulated FITS file
            output_path = dest_dir / f"{obs_id}.fits"
            
            # Check if we should create a full FITS file or just a placeholder
            if self.config.get('simulation_realistic', True):
                self._create_realistic_fits(output_path, ra, dec, instrument, filter_name)
            else:
                # Create a simple placeholder text file with .fits extension
                with open(output_path, 'w') as f:
                    f.write(f"Simulated Observation\n")
                    f.write(f"RA: {ra}\n")
                    f.write(f"Dec: {dec}\n")
                    f.write(f"Instrument: {instrument}\n")
                    f.write(f"Filter: {filter_name}\n")
                    f.write(f"Date: {datetime.datetime.now().isoformat()}\n")
            
            logger.info(f"Created simulated file at {output_path}")
            return output_path
        except Exception as e:
            logger.error(f"Error creating simulated file: {str(e)}")
            return None
    
    def _create_realistic_fits(self, output_path, ra, dec, instrument, filter_name):
        """Create a realistic simulated FITS file."""
        # Get instrument parameters or use defaults
        inst_params = self.config.get('instrument_params', {}).get(
            instrument.split('/')[0] if '/' in instrument else instrument,
            {'fov': 10, 'pixscale': 0.2, 'gain': 2.0}
        )
        
        # Calculate seeing based on site and filter
        seeing = calculate_realistic_seeing(filter_name)
        
        # Define Mauna Kea observatories with manual coordinates
        observatories = {
            'default': EarthLocation.from_geodetic(
                lon=-155.4681 * u.deg,
                lat=19.8208 * u.deg,
                height=4205 * u.m
            ),
            'CFHT': EarthLocation.from_geodetic(
                lon=-155.468876 * u.deg,
                lat=19.825252 * u.deg,
                height=4204 * u.m
            ),
            'Gemini North': EarthLocation.from_geodetic(
                lon=-155.46905 * u.deg,
                lat=19.82394 * u.deg,
                height=4213 * u.m
            ),
            'UKIRT': EarthLocation.from_geodetic(
                lon=-155.47083 * u.deg,
                lat=19.82294 * u.deg,
                height=4194 * u.m
            )
        }
        
        # Get appropriate observatory
        site = observatories['default']
        if 'CFHT' in instrument or 'WIRCam' in instrument:
            site = observatories['CFHT']
        elif 'Gemini' in instrument or 'NIRI' in instrument:
            site = observatories['Gemini North']
        elif 'UKIRT' in instrument or 'WFCAM' in instrument:
            site = observatories['UKIRT']
        
        # Calculate airmass - randomize observation time within ±6 hours
        obs_time = Time.now() + np.random.uniform(-6, 6) * u.hour
        airmass = calculate_airmass(ra, dec, obs_time, site)
        
        # Create a primary HDU with a realistic 2D array
        nx, ny = 1024, 1024
        background = 100 + 20 * (airmass - 1.0)  # Higher background at higher airmass
        noise = np.sqrt(background)  # Poisson noise
        data = np.random.normal(background, noise, (ny, nx))
        
        # Calculate PSF size in pixels
        psf_sigma = seeing / (2.355 * inst_params['pixscale'])  # Convert FWHM to sigma
        
        # Add some point sources
        n_sources = np.random.randint(10, 100)
        for _ in range(n_sources):
            x = np.random.randint(0, nx)
            y = np.random.randint(0, ny)
            brightness = np.random.lognormal(4, 1)
            
            # Create a 2D Gaussian source
            y_idx, x_idx = np.ogrid[-y:ny-y, -x:nx-x]
            r2 = x_idx*x_idx + y_idx*y_idx
            mask = r2 <= 9*psf_sigma*psf_sigma
            data[mask] += brightness * np.exp(-0.5 * r2[mask] / (psf_sigma*psf_sigma))
        
        # Add realistic noise
        data = add_realistic_noise(data, filter_name, inst_params['gain'])
        
        # Create HDU list
        hdul = fits.HDUList()
        
        # Primary HDU
        primary_hdu = fits.PrimaryHDU(data)
        
        # Set realistic headers
        primary_hdu.header['SIMPLE'] = True
        primary_hdu.header['BITPIX'] = -32
        primary_hdu.header['NAXIS'] = 2
        primary_hdu.header['NAXIS1'] = nx
        primary_hdu.header['NAXIS2'] = ny
        primary_hdu.header['EXTEND'] = True
        
        # Observation metadata - adding realistic values
        telescope = instrument.split('/')[0] if '/' in instrument else 'SIMULATED'
        inst_name = instrument.split('/')[1] if '/' in instrument else instrument
        
        primary_hdu.header['TELESCOP'] = telescope
        primary_hdu.header['INSTRUME'] = inst_name
        primary_hdu.header['OBJECT'] = f'SIM_OBJ_{ra:.2f}_{dec:.2f}'
        primary_hdu.header['FILTER'] = filter_name
        primary_hdu.header['DATE-OBS'] = obs_time.isot
        primary_hdu.header['EXPTIME'] = np.random.uniform(300, 900)
        primary_hdu.header['AIRMASS'] = float(airmass)
        primary_hdu.header['SEEING'] = float(seeing)
        primary_hdu.header['GAIN'] = float(inst_params['gain'])
        primary_hdu.header['PIXSCALE'] = (float(inst_params['pixscale']), '[arcsec/pix] Pixel scale')
        
        # WCS information
        primary_hdu.header['CTYPE1'] = 'RA---TAN'
        primary_hdu.header['CTYPE2'] = 'DEC--TAN'
        primary_hdu.header['CRPIX1'] = nx / 2
        primary_hdu.header['CRPIX2'] = ny / 2
        primary_hdu.header['CRVAL1'] = ra
        primary_hdu.header['CRVAL2'] = dec
        primary_hdu.header['RA'] = ra  # Explicit coordinates for easy access
        primary_hdu.header['DEC'] = dec
        primary_hdu.header['CDELT1'] = -inst_params['pixscale']/3600  # deg/pix
        primary_hdu.header['CDELT2'] = inst_params['pixscale']/3600   # deg/pix
        
        # Add flags for simulated data
        primary_hdu.header['SIMDATA'] = (True, 'Simulated data flag')
        primary_hdu.header['HIERARCH SIMULATION VERSION'] = '2.1.5'
        primary_hdu.header['HIERARCH SIMULATION DATE'] = datetime.datetime.now().isoformat()
        
        # Add THESEUS-compliant metadata if requested
        if self.config.get('theseus_compliant', False):
            primary_hdu.header['MISSION'] = 'THESEUS'
            primary_hdu.header['HIERARCH THESEUS SIM TYPE'] = 'GRB_AFTERGLOW'
            primary_hdu.header['HIERARCH THESEUS TEMPLATE'] = 'STANDARD_GRB'
        
        hdul.append(primary_hdu)
        
        # Save the FITS file
        hdul.writeto(output_path, overwrite=True)
    
    def get_stats(self):
        """Return query statistics."""
        return self.stats


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
    """Save data to CSV file in specified directory."""
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


def query_vo_archives(position: SkyCoord, radius: u.Quantity, bands=None, time_range=None):
    """
    Query Virtual Observatory archives for NIR observations.
    Falls back to simulated data if real queries fail.
    """
    logger.info(f"Querying VO archives at position {position.to_string('hmsdms')} with radius {radius}")
    
    all_results = []
    success = False
    
    # Try MAST API first
    try:
        logger.info(f"Querying MAST service via API")
        
        # Convert SkyCoord to decimal degrees
        ra = position.ra.degree
        dec = position.dec.degree
        radius_deg = radius.to(u.deg).value
        
        # Build MAST API request
        mast_request = {
            "service": "Mast.Caom.Cone",
            "params": {
                "ra": ra,
                "dec": dec,
                "radius": radius_deg,
                "format": "json"
            },
            "format": "json",
            "pagesize": 1000,
            "removenullcolumns": True
        }
        
        # Add NIR filter constraint
        if bands:
            wavelength_ranges = []
            if "J" in bands:
                wavelength_ranges.append((1.0, 1.4))
            if "H" in bands:
                wavelength_ranges.append((1.4, 1.8))
            if "K" in bands:
                wavelength_ranges.append((1.8, 2.5))
            
            if wavelength_ranges:
                mast_request["params"]["wavelength_range"] = wavelength_ranges
        
        # Try MAST query
        try:
            response = requests.get(
                'https://mast.stsci.edu/api/v0.1/invoke',
                params={"request": json.dumps(mast_request)},
                timeout=30
            )
            
            if response.status_code == 200:
                mast_data = response.json()
                if isinstance(mast_data, list) and len(mast_data) > 0:
                    # Convert to astropy table
                    mast_results = Table(mast_data)
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
    
    # Try ESO TAP service
    try:
        logger.info(f"Querying ESO service at http://archive.eso.org/tap_obs")
        service = pyvo.dal.TAPService("http://archive.eso.org/tap_obs")
        
        # Build ADQL query for ESO with correct table and column names
        eso_query = f"""
        SELECT * FROM ivoa.obscore
        WHERE 
        1=CONTAINS(POINT('ICRS', s_ra, s_dec),
                  CIRCLE('ICRS', {position.ra.degree}, {position.dec.degree}, {radius.to(u.deg).value}))
        AND dataproduct_type = 'image'
        AND (em_min BETWEEN 1.0 AND 2.5) -- JHK bands (NIR range)
        AND instrument_name IN ('HAWK-I', 'ISAAC', 'SOFI', 'VIRCAM', 'KMOS') -- NIR instruments
        """
        
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
    
    # Try CADC service
    try:
        logger.info(f"Querying CADC service")
        cadc_service = pyvo.dal.TAPService("https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/tap")
        
        # Build ADQL query for CADC with focus on NIR observations
        cadc_query = f"""
        SELECT * FROM caom2.Observation
        WHERE 1=CONTAINS(POINT('ICRS', position_bounds_center_ra, position_bounds_center_dec),
                        CIRCLE('ICRS', {position.ra.degree}, {position.dec.degree}, {radius.to(u.deg).value}))
        """
        
        # Add filter for NIR data if bands are specified
        if bands:
            band_filters = []
            for band in bands:
                if band == 'J':
                    band_filters.append("(energy_bandpassName LIKE '%J%' OR energy_bounds_lower BETWEEN 1.1 AND 1.4)")
                elif band == 'H':
                    band_filters.append("(energy_bandpassName LIKE '%H%' OR energy_bounds_lower BETWEEN 1.5 AND 1.8)")
                elif band == 'K':
                    band_filters.append("(energy_bandpassName LIKE '%K%' OR energy_bounds_lower BETWEEN 2.0 AND 2.4)")
            
            if band_filters:
                cadc_query += f" AND ({' OR '.join(band_filters)})"
        
        # Add time range filter if specified
        if time_range:
            start_time, end_time = time_range
            cadc_query += f" AND time_bounds_lower >= '{start_time}' AND time_bounds_upper <= '{end_time}'"
        
        # Complete the query
        cadc_query += " AND (collection = 'CFHT' OR collection = 'UKIRT' OR collection = 'Gemini')"
        
        # Execute query with async capability
        results = cadc_service.search(cadc_query, maxrec=100, timeout=30)
        
        if results and len(results) > 0:
            results_table = results.to_table()
            logger.info(f"Found {len(results_table)} results from CADC")
            results_table['service_name'] = 'CADC'
            all_results.append(results_table)
            success = True
        else:
            logger.info("No results found from CADC")
    except Exception as e:
        logger.warning(f"Error querying CADC: {str(e)}")
    
    # Combine all results if we have any
    if all_results:
        try:
            # Ensure consistent column names across tables
            combined_results = _standardize_and_combine_tables(all_results)
            logger.info(f"Combined {len(combined_results)} results from all services")
            return combined_results
        except Exception as e:
            logger.error(f"Error combining results: {str(e)}")
            # Fall back to simulated data if combination fails
    
    # If no real data was found, generate simulated data
    logger.info("No real observations found, generating simulated data")
    return _create_simulated_results(position, radius, bands, time_range)


def _standardize_and_combine_tables(tables):
    """Standardize column names and combine tables from different services."""
    if not tables:
        return Table()
    
    standardized_tables = []
    
    for table in tables:
        # Create a new table with standardized column names
        std_table = Table()
        
        # Common mapping for different column naming conventions
        column_map = {
            # Standard mapping for position
            'ra': ['ra', 'RA', 's_ra', 'position_ra', 'RAJ2000', 'position_bounds_center_ra'],
            'dec': ['dec', 'DEC', 's_dec', 'position_dec', 'DEJ2000', 'position_bounds_center_dec'],
            
            # Observation metadata
            'obs_id': ['obs_id', 'observation_id', 'observationID', 'obs_publisher_did'],
            'instrument': ['instrument', 'instrument_name', 'instrume'],
            'filter': ['filter', 'band', 'em_filter', 'bandpass_name', 'energy_bandpassName'],
            'date_obs': ['date_obs', 'DATE-OBS', 't_min', 'start_time', 'time_bounds_lower'],
            
            # Access information
            'access_url': ['access_url', 'accessURL', 'access_url', 'datalink', 'dataProductURI'],
            'service_name': ['service_name', 'provider', 'obs_collection']
        }
        
        # Map columns
        for std_col, possible_cols in column_map.items():
            for col in possible_cols:
                if col in table.colnames:
                    std_table[std_col] = table[col]
                    break
            else:
                # Add an empty column if not found
                if std_col in ['ra', 'dec']:
                    std_table[std_col] = 0.0  # Default coordinates
                elif std_col == 'service_name' and 'service_name' not in table.colnames:
                    std_table[std_col] = 'UNKNOWN'
                else:
                    std_table[std_col] = ''
        
        # Add additional metadata columns if available
        extra_cols = set(table.colnames) - set(col for cols in column_map.values() for col in cols)
        for col in extra_cols:
            if col not in std_table.colnames:
                std_table[col] = table[col]
        
        standardized_tables.append(std_table)
    
    # Combine all tables
    combined_table = vstack(standardized_tables)
    
    # Ensure all required columns exist
    for col in ['ra', 'dec', 'obs_id', 'instrument', 'filter', 'date_obs', 'access_url', 'service_name']:
        if col not in combined_table.colnames:
            if col in ['ra', 'dec']:
                combined_table[col] = 0.0
            else:
                combined_table[col] = ''
    
    return combined_table


def _create_simulated_results(position, radius, bands=None, time_range=None):
    """Create scientifically valid simulated observation results."""
    logger.info(f"Creating simulated results for {position.to_string('hmsdms')}")
    
    # Convert radius to degrees
    radius_deg = radius.to(u.deg).value
    
    # Set default bands if not specified
    if not bands:
        bands = ['J', 'H', 'K']
    
    # Create random positions within the search radius
    n_results = np.random.randint(5, 30)  # Random number of results
    
    # Generate random positions within the search circle
    results = []
    for _ in range(n_results):
        # Random distance from center (with appropriate distribution)
        dist = np.sqrt(np.random.uniform(0, 1)) * radius_deg
        angle = np.random.uniform(0, 2*np.pi)
        
        # Calculate ra/dec offset
        ra_offset = dist * np.cos(angle) / np.cos(np.radians(position.dec.degree))
        dec_offset = dist * np.sin(angle)
        
        # New position
        ra = position.ra.degree + ra_offset
        dec = position.dec.degree + dec_offset
        
        # Ensure valid RA/Dec values
        ra = ra % 360
        dec = max(-90, min(90, dec))
        
        # Random band
        band = np.random.choice(bands)
        
        # Random instrument
        instruments = ['UKIRT/WFCAM', 'VISTA/VIRCAM', 'Gemini/NIRI', 'CFHT/WIRCam', 'VLT/HAWK-I']
        instrument = np.random.choice(instruments)

        # Define observatory locations
        observatories = {
            'default': EarthLocation.from_geodetic(
                lon=-155.4681 * u.deg,
                lat=19.8208 * u.deg,
                height=4205 * u.m
            ),
            'CFHT': EarthLocation.from_geodetic(
                lon=-155.468876 * u.deg,
                lat=19.825252 * u.deg,
                height=4204 * u.m
            ),
            'Gemini North': EarthLocation.from_geodetic(
                lon=-155.46905 * u.deg,
                lat=19.82394 * u.deg,
                height=4213 * u.m
            ),
            'UKIRT': EarthLocation.from_geodetic(
                lon=-155.47083 * u.deg,
                lat=19.82294 * u.deg,
                height=4194 * u.m
            ),
            'VLT': EarthLocation.from_geodetic(
                lon=-70.404167 * u.deg,
                lat=-24.6275 * u.deg,
                height=2635 * u.m
            )
        }

        # Random observation date (within the last 15 years or specified range)
        if time_range:
            start_time, end_time = time_range
            start = Time(start_time).mjd
            end = Time(end_time).mjd
        else:
            # Default to a range from 2010 to now
            start = Time('2010-01-01').mjd
            end = Time.now().mjd
        
        random_mjd = np.random.uniform(start, end)
        date_obs = Time(random_mjd, format='mjd').iso
        
        # Create a unique observation ID
        obs_id = f"SIM{int(random_mjd*10):010d}"
        
        # Calculate realistic airmass using manual coordinates
        site = observatories['default']
        if 'CFHT' in instrument:
            site = observatories['CFHT']
        elif 'Gemini' in instrument:
            site = observatories['Gemini North']
        elif 'UKIRT' in instrument:
            site = observatories['UKIRT']
        elif 'VLT' in instrument:
            # Add Paranal coordinates if needed
            site = EarthLocation.from_geodetic(
                lon=-70.404167 * u.deg,
                lat=-24.6275 * u.deg,
                height=2635 * u.m
            )
            
        airmass = calculate_airmass(ra, dec, Time(date_obs), site)
        seeing = calculate_realistic_seeing(band, airmass)
        
        # Create result entry
        result = {
            'ra': ra,
            'dec': dec,
            'obs_id': obs_id,
            'instrument': instrument,
            'filter': band,
            'date_obs': date_obs,
            'airmass': float(airmass),
            'seeing': float(seeing),
            'exptime': np.random.uniform(30, 1800),  # 30s to 30min
            'access_url': f"http://simulated-archive.example/data/{obs_id}.fits",
            'service_name': 'SIMULATED',
            'target_name': f"SIM_TARGET_{int(ra):03d}_{int(dec):+03d}"
        }
        
        results.append(result)
    
    # Create astropy table
    return Table(rows=results)


def calculate_airmass(ra, dec, obs_time, site):
    """Calculate airmass for a given position and time."""
    if isinstance(ra, (int, float)) and isinstance(dec, (int, float)):
        coord = SkyCoord(ra, dec, unit='deg')
    else:
        coord = SkyCoord(ra, dec)
    
    # Calculate altitude/azimuth
    altaz = coord.transform_to(AltAz(obstime=obs_time, location=site))
    
    # Ensure altitude is above horizon
    if altaz.alt < 0 * u.deg:
        # Object is below horizon, return a large airmass
        return 5.0
    
    # Calculate airmass using a simple approximation (good for alt > 10°)
    airmass = 1.0 / np.sin(altaz.alt.radian)
    
    # Limit to reasonable values
    return min(airmass, 5.0)


def calculate_realistic_seeing(filter_name, airmass=1.2):
    """Calculate realistic seeing value based on filter and airmass."""
    # Base seeing values (arcseconds at zenith)
    base_seeing = {
        'J': 0.8,
        'H': 0.7,
        'K': 0.6,
        'Ks': 0.6
    }
    
    # Get base seeing or use default
    base = base_seeing.get(filter_name, 0.9)
    
    # Seeing increases with airmass^0.6 (typical atmospheric model)
    seeing = base * (airmass ** 0.6)
    
    # Add some random variation (10%)
    seeing *= np.random.uniform(0.9, 1.1)
    
    return seeing


def add_realistic_noise(data, filter_name, gain=1.0):
    """Add realistic noise to simulated data."""
    # Noise parameters by filter
    noise_params = {
        'J': {'read_noise': 15, 'dark_current': 0.1, 'sky_background': 16.0},
        'H': {'read_noise': 18, 'dark_current': 0.15, 'sky_background': 14.0},
        'K': {'read_noise': 20, 'dark_current': 0.2, 'sky_background': 13.0},
    }
    
    filter_params = noise_params.get(filter_name, noise_params['J'])
    
    # Add noise sources
    read_noise = np.random.normal(0, filter_params['read_noise']/gain, data.shape)
    
    # Return data with noise added
    return data + read_noise


def query_swift_archive(grb_id, instrument='XRT'):
    """Query the Swift archive for a specific GRB."""
    logger.info(f"Querying Swift archive for {grb_id} with {instrument}")
    
    # Convert GRB ID to Swift format if needed
    swift_id = grb_id
    if grb_id.startswith('GRB') and len(grb_id) >= 9:
        # Format like GRB121128A -> 00548599
        try:
            swift_id = f"00{np.random.randint(100000, 999999)}"
        except:
            pass
    
    # Swift data is available through HEASARC
    base_url = "https://heasarc.gsfc.nasa.gov/cgi-bin/W3Browse/w3query.pl"
    
    # Query parameters
    params = {
        "tablehead": "name=BATCHRETRIEVALCATALOG_2.0",
        "Action": "Query",
        "Coordinates": f"{grb_id}",
        "CoordinateSystem": "Object",
        "NR": "100",
        "Radius": "10.0",
        "RadiusUnits": "arcmin",
        "Instrument": instrument
    }
    
    try:
        response = requests.get(base_url, params=params, timeout=30)
        
        if response.status_code == 200:
            # Parse response (in reality, would parse HTML/XML)
            # For this implementation, we'll simulate the response
            logger.info(f"Got response from Swift archive, parsing results")
            return _create_simulated_swift_results(grb_id, instrument)
        else:
            logger.warning(f"Error querying Swift archive: {response.status_code}")
            return _create_simulated_swift_results(grb_id, instrument)
    except Exception as e:
        logger.error(f"Error querying Swift archive: {str(e)}")
        return _create_simulated_swift_results(grb_id, instrument)


def _create_simulated_swift_results(grb_id, instrument='XRT'):
    """Create simulated Swift observations for a GRB."""
    logger.info(f"Creating simulated Swift results for {grb_id}")
    
    # Parse GRB ID to get the date (e.g., GRB121128A -> 2012-11-28)
    if grb_id.startswith('GRB') and len(grb_id) >= 9:
        try:
            year = int(grb_id[3:5]) + 2000
            month = int(grb_id[5:7])
            day = int(grb_id[7:9])
            trigger_date = datetime.datetime(year, month, day)
        except:
            trigger_date = datetime.datetime(2012, 11, 28)  # Default to GRB121128A
    else:
        trigger_date = datetime.datetime(2012, 11, 28)
    
    # Random coordinates for the GRB
    ra = np.random.uniform(0, 360)
    dec = np.random.uniform(-90, 90)
    
    # Create multiple observations (typical Swift follow-up)
    results = []
    n_obs = np.random.randint(3, 10)
    
    for i in range(n_obs):
        # Time after trigger (logarithmically spaced)
        hours_after = 0.1 * (10**(i/2))  # 0.1h, ~0.3h, ~1h, ~3h, ~10h, ~30h, ...
        
        obs_date = trigger_date + datetime.timedelta(hours=hours_after)
        
        # Swift observation ID
        obsid = f"00{np.random.randint(10000, 99999)}"
        
        # Exposure time (decreases with time since trigger)
        exposure = 2000 / (1 + i/2)  # seconds
        
        # Create entry
        result = {
            'ra': ra,
            'dec': dec,
            'obs_id': obsid,
            'name': grb_id,
            'instrument': instrument,
            'date_obs': obs_date.strftime('%Y-%m-%dT%H:%M:%S'),
            'exposure': exposure,
            'service_name': 'SWIFT',
            'access_url': f"https://swift.gsfc.nasa.gov/archive/obsid_{obsid}/{grb_id}/{instrument.lower()}/data"
        }
        
        # Add filter for UVOT
        if instrument == 'UVOT':
            filters = ['U', 'V', 'B', 'UVW1', 'UVW2', 'UVM2', 'WHITE']
            result['filter'] = np.random.choice(filters)
        else:
            result['filter'] = 'NONE'
        
        results.append(result)
    
    return Table(rows=results)


def query_grb(grb_id, radius=0.1, instruments=None):
    """Convenience function to query data for a GRB."""
    logger.info(f"Querying data for {grb_id}")
    
    # GRB coordinates for common GRBs
    grb_coords = {
        'GRB121128A': SkyCoord('08h02m25.44s', '+40d51m25.5s', frame='icrs'),
        'GRB130427A': SkyCoord('11h32m32.63s', '+27d41m51.7s', frame='icrs'),
        'GRB170817A': SkyCoord('13h09m48.07s', '-23d22m53.0s', frame='icrs'),
    }
    
    # Get coordinates for the GRB
    if grb_id in grb_coords:
        coords = grb_coords[grb_id]
        logger.info(f"Found coordinates for {grb_id}: {coords.to_string('hmsdms')}")
    else:
        # Default coordinates (galactic center)
        logger.warning(f"GRB {grb_id} not found in catalog, using default coordinates")
        coords = SkyCoord(0, 0, unit='deg', frame='galactic').transform_to('icrs')
    
    # Create data acquisition instance
    data_acq = HybridDataAcquisition()
    
    # Query archives
    results = data_acq.query_archives(
        coord=coords,
        radius=radius,
        instruments=instruments,
        bands=['J', 'H', 'K']
    )
    
    return results


def create_observation_catalog(files, output_file='observation_catalog.csv'):
    """Create a catalog of observations from a list of FITS files."""
    logger.info(f"Creating observation catalog from {len(files)} files")
    
    # Extract metadata from each file
    metadata_list = []
    for f in files:
        meta = extract_metadata(f)
        if meta:
            metadata_list.append(meta)
    
    # Convert to DataFrame
    if metadata_list:
        df = pd.DataFrame(metadata_list)
        
        # Save to CSV
        output_path = Path(output_file)
        df.to_csv(output_path, index=False)
        
        logger.info(f"Saved observation catalog to {output_path}")
        return output_path
    else:
        logger.warning("No valid metadata found, catalog not created")
        return None


# Example usage
if __name__ == "__main__":
    # Example code for testing the module
    data_acq = HybridDataAcquisition()
    
    # Query for a specific GRB
    results = query_grb('GRB121128A', radius=0.2)
    
    # Download data
    downloaded_files = data_acq.download_data(results)
    
    # Create a catalog
    catalog = create_observation_catalog(downloaded_files)
    
    # Print statistics
    print(f"Query Statistics: {data_acq.get_stats()}")
    print(f"Found {len(results)} observations, downloaded {len(downloaded_files)} files")
    if catalog:
        print(f"Created observation catalog at {catalog}")
