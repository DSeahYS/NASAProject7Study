import os
import time
import requests
import json
import logging
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.table import Table, QTable
from astropy import units as u
import astropy_healpix as ah  # Using astropy_healpix instead of healpy

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('data_acquisition.alerts')

class AlertMonitor:
    def __init__(self, output_dir='alerts_data'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.gcn_endpoint = "https://gcn.nasa.gov/api/notices/latest"
        self.gracedb_endpoint = "https://gracedb.ligo.org/api/superevents/public/"
        self.headers = {
            "User-Agent": "NIR_GRB_Pipeline/1.0",
            "Accept": "application/json"
        }
    
    def monitor_gcn_notices(self, callback=None, poll_interval=60):
        """
        Monitor GCN notices for new GRB alerts.
        
        Parameters
        ----------
        callback : function, optional
            Function to call when a new notice is received
        poll_interval : int, optional
            Interval in seconds between polling for new notices
        """
        logger.info(f"Starting GCN notice monitoring with interval of {poll_interval}s")
        
        # Keep track of the latest notice ID we've seen
        latest_notice_id = None
        
        try:
            while True:
                try:
                    # Request the latest notices
                    response = requests.get(self.gcn_endpoint, headers=self.headers)
                    response.raise_for_status()
                    notices = response.json()
                    
                    # Check if there are new notices
                    if notices and notices[0]['id'] != latest_notice_id:
                        for notice in notices:
                            if latest_notice_id is None or notice['id'] > latest_notice_id:
                                # Parse the notice
                                parsed_notice = self.parse_gcn_notice(notice)
                                
                                # Save the notice
                                self.save_notice(parsed_notice)
                                
                                # If it's a GRB notice and we have a callback, call it
                                if parsed_notice.get('type') == 'GRB' and callback is not None:
                                    callback(parsed_notice)
                                
                                logger.info(f"Processed new GCN notice: {parsed_notice.get('type')} from {parsed_notice.get('source')}")
                        
                        # Update the latest notice ID
                        latest_notice_id = notices[0]['id']
                    
                    # Wait for the next poll
                    time.sleep(poll_interval)
                    
                except Exception as e:
                    logger.error(f"Error polling GCN notices: {str(e)}")
                    time.sleep(poll_interval)  # Wait before trying again
                    
        except KeyboardInterrupt:
            logger.info("GCN monitoring stopped by user")
    
    def parse_gcn_notice(self, notice_content):
        """
        Parse GCN notice to extract GRB information.
        """
        logger.info("Parsing GCN notice")
        
        # Convert string to dictionary if needed
        if isinstance(notice_content, str):
            try:
                notice_content = json.loads(notice_content)
            except json.JSONDecodeError:
                logger.error("Failed to parse notice content as JSON")
                return {}
        
        # Initialize the parsed notice with basic information
        parsed_notice = {
            'type': None,
            'source': None,
            'time': None,
            'position': None,
            'error': None,
            'energy': None,
            'duration': None,
            'flux': None,
            'raw_notice': notice_content
        }
        
        try:
            # Extract the notice type
            if 'type' in notice_content:
                parsed_notice['type'] = notice_content['type']
            elif 'subject' in notice_content:
                # Try to determine the type from the subject
                subject = notice_content['subject'].upper()
                if 'GRB' in subject:
                    parsed_notice['type'] = 'GRB'
                elif 'LVC' in subject or 'LIGO' in subject:
                    parsed_notice['type'] = 'GW'
            
            # Extract the notice source
            if 'source' in notice_content:
                parsed_notice['source'] = notice_content['source']
            elif 'instrument' in notice_content:
                parsed_notice['source'] = notice_content['instrument']
            
            # Extract the trigger time
            if 'trigger_time' in notice_content:
                parsed_notice['time'] = Time(notice_content['trigger_time'])
            elif 'time' in notice_content:
                parsed_notice['time'] = Time(notice_content['time'])
            
            # Extract the position (RA, Dec)
            if all(k in notice_content for k in ('ra', 'dec')):
                ra = float(notice_content['ra'])
                dec = float(notice_content['dec'])
                parsed_notice['position'] = SkyCoord(ra=ra, dec=dec, unit='deg')
            
            # Extract the position error
            if 'error_radius' in notice_content:
                parsed_notice['error'] = float(notice_content['error_radius']) * u.deg
            elif 'error' in notice_content:
                parsed_notice['error'] = float(notice_content['error']) * u.deg
            
            # Extract GRB-specific information
            if parsed_notice['type'] == 'GRB':
                # Energy range
                if 'energy_range' in notice_content:
                    parsed_notice['energy'] = notice_content['energy_range']
                
                # Duration (T90)
                if 't90' in notice_content:
                    parsed_notice['duration'] = float(notice_content['t90']) * u.s
                
                # Flux
                if 'flux' in notice_content:
                    parsed_notice['flux'] = float(notice_content['flux'])
            
            logger.info(f"Successfully parsed notice: {parsed_notice['type']} from {parsed_notice['source']}")
            return parsed_notice
            
        except Exception as e:
            logger.error(f"Error parsing GCN notice: {str(e)}")
            return parsed_notice
    
    def save_notice(self, notice):
        """
        Save a parsed notice to disk.
        """
        try:
            # Create a unique filename based on the notice time and type
            if notice.get('time'):
                time_str = notice['time'].isot.replace(':', '-').replace('.', '-')
            else:
                time_str = datetime.now().isoformat().replace(':', '-').replace('.', '-')
            
            notice_type = notice.get('type', 'UNKNOWN')
            filename = f"{time_str}_{notice_type}.json"
            
            # Prepare notice for serialization
            serialized_notice = {}
            for k, v in notice.items():
                if k == 'raw_notice':
                    continue  # Skip raw notice
                elif isinstance(v, Time):
                    serialized_notice[k] = v.isot
                elif isinstance(v, SkyCoord):
                    serialized_notice[k] = {'ra': v.ra.degree, 'dec': v.dec.degree}
                elif isinstance(v, u.Quantity):
                    serialized_notice[k] = v.value
                else:
                    serialized_notice[k] = v
            
            # Save to file
            file_path = self.output_dir / filename
            with open(file_path, 'w') as f:
                json.dump(serialized_notice, f, indent=2)
            
            logger.info(f"Saved notice to {file_path}")
            return file_path
            
        except Exception as e:
            logger.error(f"Error saving notice: {str(e)}")
            return None
    
    def get_ligo_alerts(self, time_range=None, max_results=50):
        """
        Retrieve LIGO/Virgo alerts from GraceDB.
        
        Parameters
        ----------
        time_range : tuple, optional
            (start_time, end_time) as strings in ISO format
        max_results : int, optional
            Maximum number of results to return
            
        Returns
        -------
        alerts : list
            List of GW alert dictionaries
        """
        logger.info("Retrieving LIGO/Virgo alerts from GraceDB")
        
        # Build URL with parameters
        url = f"{self.gracedb_endpoint}?max_results={max_results}"
        
        try:
            # Request the alerts
            response = requests.get(url, headers=self.headers)
            response.raise_for_status()
            alerts_data = response.json()
            
            if 'superevents' not in alerts_data:
                logger.warning("No alerts found in response")
                return []
            
            alerts = alerts_data['superevents']
            
            # Filter by time range if provided
            if time_range:
                start_time, end_time = Time(time_range[0]), Time(time_range[1])
                filtered_alerts = []
                
                for alert in alerts:
                    if 'created' in alert:
                        alert_time = Time(alert['created'])
                        if start_time <= alert_time <= end_time:
                            filtered_alerts.append(alert)
                
                alerts = filtered_alerts
                logger.info(f"Found {len(alerts)} alerts within specified time range")
            else:
                logger.info(f"Retrieved {len(alerts)} alerts")
            
            # Save alerts to disk
            self.save_ligo_alerts(alerts)
            
            return alerts
            
        except Exception as e:
            logger.error(f"Error retrieving LIGO/Virgo alerts: {str(e)}")
            return []
    
    def save_ligo_alerts(self, alerts):
        """
        Save LIGO/Virgo alerts to disk.
        """
        try:
            # Create alerts directory if it doesn't exist
            ligo_dir = self.output_dir / 'ligo'
            ligo_dir.mkdir(exist_ok=True)
            
            # Save each alert to a separate file
            for alert in alerts:
                event_id = alert.get('superevent_id', 'unknown')
                file_path = ligo_dir / f"{event_id}.json"
                
                with open(file_path, 'w') as f:
                    json.dump(alert, f, indent=2)
                
                logger.info(f"Saved LIGO alert {event_id} to {file_path}")
            
            # Also save a summary file with all alerts
            summary_path = ligo_dir / 'all_alerts.json'
            with open(summary_path, 'w') as f:
                json.dump(alerts, f, indent=2)
            
            logger.info(f"Saved summary of {len(alerts)} LIGO alerts to {summary_path}")
            
        except Exception as e:
            logger.error(f"Error saving LIGO alerts: {str(e)}")
    
    def download_ligo_skymap(self, event_id, filename=None):
        """
        Download LIGO skymap for a given event ID.
        
        Parameters
        ----------
        event_id : str
            GraceDB superevent ID
        filename : str, optional
            Filename to save the skymap to
            
        Returns
        -------
        skymap_file : str
            Path to downloaded skymap file
        """
        logger.info(f"Downloading skymap for LIGO event {event_id}")
        
        # Create skymaps directory if it doesn't exist
        skymap_dir = self.output_dir / 'skymaps'
        skymap_dir.mkdir(exist_ok=True)
        
        if filename is None:
            filename = f"{event_id}_skymap.fits.gz"
        
        skymap_file = skymap_dir / filename
        
        # URL for the skymap
        skymap_url = f"https://gracedb.ligo.org/api/superevents/{event_id}/files/bayestar.fits.gz"
        
        try:
            # Download the skymap
            response = requests.get(skymap_url, headers=self.headers, stream=True)
            response.raise_for_status()
            
            with open(skymap_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
            
            logger.info(f"Downloaded skymap to {skymap_file}")
            return skymap_file
            
        except Exception as e:
            logger.error(f"Error downloading skymap for {event_id}: {str(e)}")
            return None
    
    def parse_ligo_skymap(self, skymap_file):
        """
        Parse LIGO probability skymap using astropy_healpix.
        """
        logger.info(f"Parsing LIGO skymap: {skymap_file}")
        
        try:
            # Read the skymap file
            skymap = QTable.read(skymap_file)
            
            # Check if this is a multi-order skymap
            is_multiorder = 'UNIQ' in skymap.colnames
            
            if is_multiorder:
                # Handle multi-order skymap
                probdensity = skymap['PROBDENSITY']
                level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
                nside = 2**level
                
                # Create a dictionary with the skymap data
                skymap_data = {
                    'skymap': skymap,
                    'is_multiorder': True,
                    'probdensity': probdensity,
                    'level': level,
                    'ipix': ipix,
                    'nside': nside,
                    'path': skymap_file
                }
            else:
                # Handle standard HEALPix map
                prob = skymap['PROB'] if 'PROB' in skymap.colnames else skymap[0].data
                
                # Try to determine NSIDE from the length of the probability array
                npix = len(prob)
                nside = int(np.sqrt(npix/12))
                
                skymap_data = {
                    'skymap': skymap,
                    'is_multiorder': False,
                    'prob': prob,
                    'nside': nside,
                    'path': skymap_file
                }
                
            # Find the highest probability pixel
            if is_multiorder:
                max_prob_idx = np.argmax(probdensity)
                ipix_max = ipix[max_prob_idx]
                level_max = level[max_prob_idx]
                # Convert the HEALPix index to sky coordinates
                lon, lat = ah.healpix_to_lonlat(ipix_max, nside=2**level_max, order='nested')
            else:
                max_prob_idx = np.argmax(prob)
                # Convert the HEALPix index to sky coordinates
                lon, lat = ah.healpix_to_lonlat(max_prob_idx, nside=nside, order='nested')
                
            max_prob_ra, max_prob_dec = lon.to(u.deg).value, lat.to(u.deg).value
            
            skymap_data['max_prob_position'] = SkyCoord(ra=max_prob_ra, dec=max_prob_dec, unit='deg')
            
            logger.info(f"Successfully parsed skymap")
            return skymap_data
            
        except Exception as e:
            logger.error(f"Error parsing LIGO skymap: {str(e)}")
            return None
            
    def check_visibility(self, coordinates: SkyCoord, observatory, time=None):
        """
        Check if target is visible from specific observatory.
        
        Parameters
        ----------
        coordinates : SkyCoord
            Sky coordinates of the target
        observatory : str or EarthLocation
            Observatory name or EarthLocation object
        time : Time or str, optional
            Time to check visibility (default: now)
            
        Returns
        -------
        visibility : dict
            Dictionary with visibility information
        """
        logger.info(f"Checking visibility of {coordinates.to_string('hmsdms')} from {observatory}")
        
        # Define observatory location if string is provided
        if isinstance(observatory, str):
            # Dictionary of observatory locations
            observatories = {
                'LDT': EarthLocation.from_geodetic(-111.4231, 34.7443, 2360),  # Lowell Discovery Telescope
                'PRIME': EarthLocation.from_geodetic(20.8112, -32.3795, 1760),  # PRIME at SAAO
                'SOAR': EarthLocation.from_geodetic(-70.7040, -30.2380, 2738),
                'Keck': EarthLocation.from_geodetic(-155.4747, 19.8260, 4160),
                'Gemini-N': EarthLocation.from_geodetic(-155.4690, 19.8238, 4213),
                'Gemini-S': EarthLocation.from_geodetic(-70.7234, -30.2407, 2722),
                'CTIO': EarthLocation.from_geodetic(-70.8150, -30.1650, 2215),
                'CFHT': EarthLocation.from_geodetic(-155.4681, 19.8250, 4204)
            }
            
            if observatory not in observatories:
                logger.error(f"Unknown observatory: {observatory}")
                return None
            
            location = observatories[observatory]
        else:
            location = observatory
        
        # Set time to now if not provided
        if time is None:
            time = Time.now()
        elif isinstance(time, str):
            time = Time(time)
        
        # Calculate target's altitude and azimuth
        altaz = coordinates.transform_to(AltAz(obstime=time, location=location))
        
        # Define visibility criteria
        sun_altaz = SkyCoord(0, 0, unit='deg').transform_to(AltAz(obstime=time, location=location))
        
        # Check if sun is below horizon (astronomical twilight)
        sun_below_horizon = sun_altaz.alt < -18 * u.deg
        
        # Check if target is above minimum altitude (30 degrees is a common threshold)
        target_above_horizon = altaz.alt > 30 * u.deg
        
        # Calculate airmass (use 999 for targets below horizon)
        if altaz.alt > 0 * u.deg:
            airmass = 1.0 / np.cos(np.pi/2 - altaz.alt.radian)
        else:
            airmass = 999.0
        
        # Calculate when the target transits (reaches highest altitude)
        time_array = Time(time) + np.linspace(0, 24, 49) * u.hour
        altaz_array = coordinates.transform_to(AltAz(obstime=time_array, location=location))
        transit_idx = np.argmax(altaz_array.alt)
        transit_time = time_array[transit_idx]
        
        # Calculate how long the target is observable tonight
        sun_altaz_array = SkyCoord(0, 0, unit='deg').transform_to(AltAz(obstime=time_array, location=location))
        is_night = sun_altaz_array.alt < -18 * u.deg
        
        if not any(is_night):
            observable_hours = 0.0
        else:
            # Find contiguous segments of nighttime
            night_segments = np.where(np.diff(is_night.astype(int)) != 0)[0] + 1
            if is_night[0]:
                night_segments = np.insert(night_segments, 0, 0)
            if is_night[-1]:
                night_segments = np.append(night_segments, len(is_night))
            
            # Reshape into pairs of start/end indices
            night_segments = night_segments.reshape(-1, 2)
            
            # Find when target is above horizon during each night segment
            observable_hours = 0.0
            for start, end in night_segments:
                target_up = altaz_array.alt[start:end] > 30 * u.deg
                if any(target_up):
                    # Calculate duration of observability in hours
                    n_steps = sum(target_up)
                    observable_hours += n_steps * 24 / 48  # 48 steps over 24 hours
        
        visibility = {
            'coordinates': coordinates,
            'observatory': observatory,
            'time': time,
            'altitude': altaz.alt,
            'azimuth': altaz.az,
            'airmass': airmass,
            'sun_altitude': sun_altaz.alt,
            'is_night': sun_below_horizon,
            'is_observable': sun_below_horizon and target_above_horizon,
            'transit_time': transit_time,
            'observable_hours': observable_hours
        }
        
        if visibility['is_observable']:
            logger.info(f"Target is currently observable from {observatory} (altitude: {altaz.alt:.1f}, airmass: {airmass:.2f})")
        else:
            if not sun_below_horizon:
                logger.info(f"Target not observable from {observatory}: daytime (sun altitude: {sun_altaz.alt:.1f})")
            else:
                logger.info(f"Target not observable from {observatory}: below altitude limit (altitude: {altaz.alt:.1f})")
        
        logger.info(f"Target will transit at {transit_time.iso} and is observable for {observable_hours:.1f} hours tonight")
        return visibility
