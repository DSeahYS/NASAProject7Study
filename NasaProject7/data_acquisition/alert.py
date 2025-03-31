# alerts.py - Enhanced Multi-Messenger Alert Processing

import os
import json
import logging
import gzip
import io
import threading
import time
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun, get_moon
from astropy import units as u
from astropy.table import QTable, Table
import astropy_healpix as ah
from astropy_healpix import HEALPix

# Optional imports for real GCN processing
try:
    from gcn_kafka import Consumer
    from voeventparse import load
    HAS_GCN = True
except ImportError:
    # Mock imports for development
    HAS_GCN = False
    
    class Consumer:
        def __init__(self, *args, **kwargs): pass
        def consume(self, *args, **kwargs): return []
    
    def load(*args, **kwargs): return None

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('data_acquisition.alerts')

class MultiMessengerAlertHandler:
    """O4-era alert processor with low-latency GW/GRB correlation"""
    
    def __init__(self, config_path: str = 'config/alerts.json', alert_callback=None):
        """Initialize the alert handler"""
        self.config = self._load_config(config_path)
        self.locations = self._init_observatories()
        self.skymap_cache = {}
        self.grb_cache = []
        self.gw_cache = []
        self.coincidences = []
        self.monitoring_active = False
        self.alert_thread = None
        self.alert_lock = threading.Lock()  # Thread-safe access to alerts
        self.alert_callback = alert_callback  # Callback to notify UI
        
        # Alert statistics
        self.metrics = {
            'processed': 0,
            'gw_events': 0,
            'grb_events': 0,
            'coincident': 0,
            'false_alarms': 0
        }
        
        # Create output directories if they don't exist
        self.output_dir = Path(self.config.get('output_dir', 'outputs/alerts'))
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def _load_config(self, path: str) -> Dict:
        """Load observatory and analysis parameters"""
        try:
            with open(path) as f:
                config = json.load(f)
                
            # Validate configuration
            required = ['grace_dbs', 'gcn_streams', 'observatories']
            for key in required:
                if key not in config:
                    raise ValueError(f"Missing required config key: {key}")
            return config
        except FileNotFoundError:
            logger.warning(f"Config file {path} not found, using default configuration")
            # Provide default configuration
            return {
                'grace_dbs': {
                    'url': 'https://gracedb.ligo.org/api/',
                    'timeout': 30
                },
                'gcn_streams': {
                    'domain': 'gcn.nasa.gov',
                    'secret': 'demo_secret'
                },
                'observatories': {
                    'LDT': {
                        'lon': -111.421,
                        'lat': 34.744,
                        'elevation': 2360
                    },
                    'PRIME': {
                        'lon': 149.066,
                        'lat': -31.276,
                        'elevation': 1165
                    }
                },
                'output_dir': 'outputs/alerts',
                'skymap_cache_size': 10,
                'coincidence_window': 600,  # seconds
                'probability_threshold': 0.1,
                'notification_email': 'alerts@example.com'
            }
    
    def _init_observatories(self) -> Dict[str, EarthLocation]:
        """Initialize observatory locations with O4 specifications"""
        return {
            name: EarthLocation.from_geodetic(
                lon=params['lon']*u.deg,
                lat=params['lat']*u.deg,
                height=params['elevation']*u.m
            )
            for name, params in self.config['observatories'].items()
        }
    
    def start(self):
        """Start the alert monitoring"""
        if self.alert_thread and self.alert_thread.is_alive():
            logger.warning("Alert handler already running")
            return
        
        self.monitoring_active = True
        self.alert_thread = threading.Thread(target=self._run_monitoring_loop, daemon=True)
        self.alert_thread.start()
        logger.info("Alert monitoring started")
    
    def pause(self):
        """Pause alert monitoring"""
        self.monitoring_active = False
        logger.info("Alert monitoring paused")
    
    def resume(self):
        """Resume alert monitoring"""
        self.monitoring_active = True
        logger.info("Alert monitoring resumed")
    
    def shutdown(self):
        """Shutdown the alert handler"""
        self.monitoring_active = False
        if self.alert_thread and self.alert_thread.is_alive():
            self.alert_thread.join(timeout=1.0)
        logger.info("Shutting down alert handler")
    
    def _run_monitoring_loop(self):
        """Main monitoring loop running in its own thread"""
        try:
            if HAS_GCN:
                # Real GCN monitoring
                self._monitor_gcn_stream()
            else:
                # Simulation mode for development
                while self.monitoring_active:
                    # Poll for alerts in a standard thread-safe way
                    self._poll_for_alerts()
                    
                    # Sleep to avoid high CPU usage
                    time.sleep(1)
        except Exception as e:
            logger.error(f"Alert monitoring error: {str(e)}")
    
    def _monitor_gcn_stream(self):
        """Monitor GCN stream in non-async way"""
        # Define topics to monitor
        topics = [
            'gcn.classic.voevent.LVC_INITIAL',
            'gcn.classic.voevent.LVC_UPDATE',
            'gcn.classic.voevent.FERMI_GBM_FLT',
            'gcn.classic.voevent.SWIFT_BAT_ALERT'
        ]
        
        # Create consumer
        consumer = Consumer(
            client_id='nir_grb_pipeline',
            client_secret=self.config['gcn_streams']['secret'],
            domain=self.config['gcn_streams']['domain']
        )
        
        while self.monitoring_active:
            try:
                # This is a blocking call, not async anymore
                for msg in consumer.consume(topics, timeout=5):
                    if not self.monitoring_active:
                        break
                    
                    voevent = load(msg.value())
                    self._process_voevent(voevent)
            except Exception as e:
                logger.error(f"Stream error: {str(e)}")
                time.sleep(5)
    
    def _poll_for_alerts(self):
        """Poll for alerts (non-async version)"""
        # In a real application, this would connect to GCN or other alert services
        # For now, just simulate random alerts occasionally
        
        # Simulate random alert arrival (10% chance each check)
        if np.random.random() < 0.1 and self.monitoring_active:
            alert_type = np.random.choice(["GRB", "GW"])
            
            # Create a mock alert
            alert = {
                'time': Time.now().iso,
                'type': alert_type,
                'id': f"{alert_type}_{int(Time.now().mjd)}",
                'ra': np.random.uniform(0, 360),
                'dec': np.random.uniform(-90, 90),
                'confidence': np.random.uniform(0, 1)
            }
            
            # Add additional fields based on alert type
            if alert_type == "GRB":
                alert['t90'] = np.random.uniform(0.1, 100)
                alert['fluence'] = np.random.uniform(1e-7, 1e-5)
                alert['energy_band'] = "15-150 keV"
                alert['detector'] = np.random.choice(["Swift", "Fermi", "INTEGRAL"])
                
            elif alert_type == "GW":
                alert['far'] = 10**(-np.random.uniform(1, 10))  # False Alarm Rate
                alert['distance'] = np.random.uniform(40, 400)  # Mpc
                alert['error_region'] = np.random.uniform(50, 2000)  # square degrees
                alert['classification'] = np.random.choice(["BNS", "BBH", "NSBH", "MassGap", "Terrestrial"])
            
            # Process the alert
            with self.alert_lock:
                if alert_type == "GW":
                    self.gw_cache.append(alert)
                    self.metrics['gw_events'] += 1
                else:
                    self.grb_cache.append(alert)
                    self.metrics['grb_events'] += 1
                self.metrics['processed'] += 1
            
            logger.info(f"Received {alert_type} alert: {alert['id']}")
            
            # Notify UI if callback is registered
            if self.alert_callback:
                try:
                    self.alert_callback(alert)
                except Exception as e:
                    logger.error(f"Error in alert callback: {str(e)}")
            
            # Check for coincidences
            self._check_coincidences()
    
    def _process_voevent(self, voevent):
        """Process IVOA-compliant VOEvent with classification"""
        event_type = self._classify_voevent(voevent)
        logger.info(f"Processing {event_type} alert: {voevent.attrib['ivorn']}")
        
        if event_type == 'gw':
            self._handle_gw_alert(voevent)
        elif event_type == 'grb':
            self._handle_grb_alert(voevent)
        
        with self.alert_lock:
            self.metrics['processed'] += 1
    
    def _classify_voevent(self, voevent) -> str:
        """Classify VOEvent using LVC classification tree"""
        ivorn = voevent.attrib['ivorn']
        if 'lvc' in ivorn.lower():
            return 'gw'
        if any(s in ivorn.lower() for s in ['fermi', 'swift']):
            return 'grb'
        return 'unknown'
    
    def _handle_gw_alert(self, voevent):
        """Process GW alert with rapid sky localization"""
        with self.alert_lock:
            self.metrics['gw_events'] += 1
        
        # Extract skymap URL and metadata
        try:
            params = {
                'graceid': voevent.params['GraceID'].value,
                'far': float(voevent.params['FAR'].value),
                'instruments': voevent.params['Instruments'].value.split(','),
                'group': voevent.params['Group'].value
            }
            
            # Try to get classification if available
            try:
                classification = {}
                for param in voevent.iterchildren('Param'):
                    if param.attrib['name'].startswith('Classification'):
                        key = param.attrib['name'].replace('Classification', '')
                        classification[key] = float(param.value)
                params['classification'] = classification
            except Exception as e:
                logger.warning(f"Error extracting classification: {str(e)}")
                params['classification'] = {}
            
            # Get trigger time
            try:
                trigger_time = Time(voevent.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Time.TimeInstant.ISOTime.text)
                params['time'] = trigger_time
            except Exception as e:
                logger.warning(f"Error extracting trigger time: {str(e)}")
                params['time'] = Time.now()
            
            # Extract skymap URL
            skymap_url = next(
                (elem.attrib['url']
                 for elem in voevent.iterchildren('Param')
                 if elem.attrib['name'] == 'skymap_fits'), None
            )
            
            if skymap_url:
                try:
                    skymap = self._download_skymap(skymap_url)
                    self._store_skymap(params['graceid'], skymap)
                    params['skymap'] = skymap
                    
                    # Save event to cache
                    with self.alert_lock:
                        self.gw_cache.append(params)
                    
                    # Check for coincidences with GRBs
                    self._check_coincidences()
                    
                    # Trigger followup if appropriate
                    self._trigger_followup(params, skymap)
                except Exception as e:
                    logger.error(f"Error processing skymap: {str(e)}")
            else:
                logger.warning(f"No skymap URL found for {params['graceid']}")
                
        except Exception as e:
            logger.error(f"Error processing GW alert: {str(e)}")
    
    def _download_skymap(self, url: str) -> QTable:
        """Download and parse HEALPix skymap with validation"""
        try:
            import requests
            
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            with gzip.open(io.BytesIO(response.content)) as f:
                skymap = QTable.read(f)
            
            self._validate_skymap(skymap)
            return skymap
        except Exception as e:
            logger.error(f"Skymap download failed: {str(e)}")
            # Create a mock skymap for testing
            return self._create_mock_skymap()
    
    def _create_mock_skymap(self) -> QTable:
        """Generate a mock HEALPix skymap for testing"""
        nside = 64
        
        # Initialize HEALPix with ICRS frame
        from astropy.coordinates import ICRS
        hp = ah.HEALPix(nside=nside, order='nested', frame=ICRS())
        npix = hp.npix
        
        # Create a dummy skymap centered at a random point
        ra = np.random.uniform(0, 360)
        dec = np.random.uniform(-90, 90)
        
        # Calculate pixel indices using proper coordinate frame
        center_coords = SkyCoord(ra*u.deg, dec*u.deg)
        center_idx = hp.skycoord_to_healpix(center_coords)
        
        # Create probability distribution (simple Gaussian-like)
        prob = np.zeros(npix)
        all_coords = hp.healpix_to_skycoord(np.arange(npix))
        separations = center_coords.separation(all_coords).deg
        prob = np.exp(-0.5 * (separations/5)**2)
        
        # Normalize
        prob /= np.sum(prob)
        
        # Create distance information
        distmu = np.full(npix, 100.0)  # Mpc
        distsigma = np.full(npix, 20.0)  # Mpc
        
        # Create QTable
        skymap = QTable()
        skymap['PROB'] = prob
        skymap['DISTMU'] = distmu
        skymap['DISTSIGMA'] = distsigma
        
        logger.warning("Created mock skymap (for testing only)")
        return skymap
    
    def _validate_skymap(self, skymap: QTable):
        """Validate HEALPix skymap against O4 specifications"""
        required = ['PROB']
        missing = [col for col in required if col not in skymap.colnames]
        if missing:
            raise ValueError(f"Missing required skymap columns: {missing}")
        
        # Check if number of pixels is valid for HEALPix
        valid_npix = [12*4**n for n in range(0, 14)]
        if len(skymap) not in valid_npix:
            logger.warning(f"Unusual HEALPix resolution: {len(skymap)} pixels")
    
    def _store_skymap(self, graceid: str, skymap: QTable):
        """Cache skymap with LRU eviction policy"""
        if len(self.skymap_cache) >= self.config.get('skymap_cache_size', 10):
            self.skymap_cache.popitem(last=False)
        
        self.skymap_cache[graceid] = skymap
        
        # Also save to disk
        output_path = self.output_dir / f"skymap_{graceid}.fits"
        try:
            skymap.write(output_path, overwrite=True)
            logger.info(f"Saved skymap to {output_path}")
        except Exception as e:
            logger.error(f"Error saving skymap: {str(e)}")
    
    def _trigger_followup(self, params: Dict, skymap: QTable):
        """Initiate multi-messenger followup sequence"""
        logger.info(f"Initiating followup for {params['graceid']}")
        
        # Calculate observatory visibility windows
        visibility = {}
        for obs_name, location in self.locations.items():
            visibility[obs_name] = self.calculate_visibility(
                skymap, location
            )
        
        # Prioritize based on FOV and instrument capabilities
        followup_plan = self._prioritize_followup(visibility)
        
        # Trigger automated observations
        self._dispatch_observations(followup_plan, params)
    
    def calculate_visibility(self, skymap: QTable, location: EarthLocation) -> Dict:
        """Compute visibility metrics for 3D skymap"""
        # Determine HEALPix resolution
        npix = len(skymap)
        nside = ah.npix2nside(npix)
        hp = HEALPix(nside=nside, order='nested')
        
        # Create time grid (24 hours)
        time_grid = Time.now() + np.linspace(0, 24, 25)*u.hour
        
        # Calculate probability-weighted visibility
        prob_visible = np.zeros(len(time_grid))
        max_prob = 0
        max_prob_coords = None
        
        # Get coordinates for all HEALPix pixels
        all_coords = hp.healpix_to_skycoord(np.arange(npix))
        
        for i, t in enumerate(time_grid):
            # Convert to AltAz
            frame = AltAz(obstime=t, location=location)
            altaz = all_coords.transform_to(frame)
            
            # Check sun altitude for night time
            sun_altaz = get_sun(t).transform_to(frame)
            is_night = sun_altaz.alt < -12*u.deg
            
            # Check moon separation (ideally > 30 deg)
            moon_coords = get_moon(t)
            moon_sep = all_coords.separation(moon_coords)
            
            # Define visibility mask (above horizon + night + away from moon)
            visible = (altaz.alt > 30*u.deg) & is_night & (moon_sep > 30*u.deg)
            
            # Calculate probability in visible area
            visible_prob = np.sum(skymap['PROB'][visible])
            prob_visible[i] = visible_prob
            
            # Track maximum probability point for pointing
            if visible_prob > max_prob:
                max_prob = visible_prob
                
                # Find highest probability pixel that's visible
                visible_idx = np.where(visible)[0]
                if len(visible_idx) > 0:
                    max_prob_idx = visible_idx[np.argmax(skymap['PROB'][visible_idx])]
                    max_prob_coords = all_coords[max_prob_idx]
        
        # Time of best visibility
        best_idx = np.argmax(prob_visible)
        
        return {
            'time': time_grid,
            'probability': prob_visible,
            'peak_time': time_grid[best_idx],
            'peak_prob': prob_visible[best_idx],
            'best_coords': max_prob_coords
        }
    
    def _prioritize_followup(self, visibility: Dict) -> Dict:
        """Generate observation plan using facility-specific thresholds"""
        plan = {}
        threshold = self.config.get('probability_threshold', 0.1)
        
        for obs_name, vis in visibility.items():
            if vis['peak_prob'] > threshold:
                start_time = vis['peak_time'] - 0.5*u.hour
                end_time = vis['peak_time'] + 2.0*u.hour
                
                plan[obs_name] = {
                    'start': start_time.iso,
                    'end': end_time.iso,
                    'priority': float(vis['peak_prob']),
                    'coordinates': None if vis['best_coords'] is None else {
                        'ra': vis['best_coords'].ra.deg,
                        'dec': vis['best_coords'].dec.deg
                    }
                }
        
        # Sort by priority (highest first)
        return dict(sorted(plan.items(), key=lambda x: -x[1]['priority']))
    
    def _dispatch_observations(self, plan: Dict, params: Dict):
        """Send observation requests to partner facilities"""
        # Create a complete observation plan
        obs_plan = {
            'event_id': params.get('graceid', f"sim_{int(Time.now().mjd)}"),
            'event_time': params.get('time', Time.now()).iso,
            'event_type': 'GW',
            'classification': params.get('classification', {}),
            'far': params.get('far', 0),
            'observatories': plan,
            'created': datetime.utcnow().isoformat()
        }
        
        # Save the plan to disk
        output_path = self.output_dir / f"plan_{obs_plan['event_id']}.json"
        with open(output_path, 'w') as f:
            json.dump(obs_plan, f, indent=2)
        
        logger.info(f"Saved observation plan to {output_path}")
        
        # In a real implementation, this would send the plan to observatories
        # via their APIs or other interfaces
        logger.info(f"Dispatching observation plan for {obs_plan['event_id']}")
        
        # If configured, send notification email
        if 'notification_email' in self.config:
            self._send_notification(obs_plan)
    
    def _send_notification(self, plan: Dict):
        """Send notification email about observation plan"""
        # This would use an email API or SMTP in a real implementation
        email = self.config['notification_email']
        logger.info(f"Would send notification to {email} about {plan['event_id']}")
    
    def _handle_grb_alert(self, voevent):
        """Process GRB alert with joint GW correlation"""
        with self.alert_lock:
            self.metrics['grb_events'] += 1
        
        # Extract GRB parameters
        try:
            # Look for standard GRB parameters in different formats
            params = {}
            for param in voevent.iterchildren('Param'):
                params[param.attrib['name']] = param.attrib['value']
            
            # Extract position information - handle different formats
            try:
                ra = float(params.get('RA', params.get('ra', 0)))
                dec = float(params.get('Dec', params.get('dec', params.get('DEC', 0))))
                error_radius = float(params.get('Error', params.get('error', params.get('ERROR', 1.0))))
            except (ValueError, KeyError) as e:
                logger.warning(f"Could not extract position from GRB alert: {str(e)}")
                
                # Use RA/Dec from the WhereWhen section if available
                try:
                    ra = float(voevent.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Position2D.Value2.C1.text)
                    dec = float(voevent.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Position2D.Value2.C2.text)
                    error_radius = 1.0  # Default
                except Exception:
                    logger.error("Could not extract RA/Dec from WhereWhen section")
                    ra, dec, error_radius = 0.0, 0.0, 10.0
            
            # Extract time
            try:
                trigger_time = Time(voevent.WhereWhen.ObsDataLocation.ObservationLocation.AstroCoords.Time.TimeInstant.ISOTime.text)
            except Exception:
                logger.warning("Could not extract trigger time, using current time")
                trigger_time = Time.now()
            
            # Extract other metadata
            instrument = params.get('Instrument', params.get('instrument', 'Unknown'))
            energy_low = float(params.get('Energy_low', params.get('LowEnergy', 0)))
            energy_high = float(params.get('Energy_high', params.get('HighEnergy', 0)))
            duration = float(params.get('Duration', params.get('T90', 0)))
            
            # Create GRB event record
            grb_event = {
                'ra': ra,
                'dec': dec,
                'error_radius': error_radius,
                'time': trigger_time,
                'instrument': instrument,
                'energy_range': (energy_low, energy_high),
                'duration': duration,
                'coordinates': SkyCoord(ra*u.deg, dec*u.deg),
                'id': f"GRB_{trigger_time.mjd:.5f}",
                'voevent_id': voevent.attrib['ivorn']
            }
            
            # Add to GRB cache for coincidence detection
            with self.alert_lock:
                self.grb_cache.append(grb_event)
            
            # Check for coincidences with GW events
            self._check_coincidences()
            
            # Save GRB information
            self._save_grb_event(grb_event)
            
            # Check for observability
            self._check_grb_observability(grb_event)
            
            logger.info(f"Processed GRB alert: {grb_event['id']} at RA={ra:.2f}, Dec={dec:.2f}")
        except Exception as e:
            logger.error(f"Error processing GRB alert: {str(e)}")
    
    def _save_grb_event(self, grb_event: Dict):
        """Save GRB event information to disk"""
        # Create a copy without SkyCoord object (not JSON serializable)
        event_copy = grb_event.copy()
        if 'coordinates' in event_copy:
            event_copy['coordinates'] = {
                'ra': float(event_copy['coordinates'].ra.deg),
                'dec': float(event_copy['coordinates'].dec.deg)
            }
        
        if isinstance(event_copy['time'], Time):
            event_copy['time'] = event_copy['time'].iso
        
        # Save to JSON file
        output_path = self.output_dir / f"grb_{grb_event['id']}.json"
        with open(output_path, 'w') as f:
            json.dump(event_copy, f, indent=2)
        
        logger.info(f"Saved GRB event to {output_path}")
    
    def _check_grb_observability(self, grb_event: Dict):
        """Check if GRB is observable from configured observatories"""
        results = {}
        
        # Create SkyCoord for GRB position
        position = grb_event['coordinates']
        
        # Check each observatory
        for obs_name, location in self.locations.items():
            # Calculate current visibility
            current_time = Time.now()
            frame = AltAz(obstime=current_time, location=location)
            altaz = position.transform_to(frame)
            
            # Check sun altitude
            sun_altaz = get_sun(current_time).transform_to(frame)
            is_night = sun_altaz.alt < -12*u.deg
            
            # Check moon separation
            moon_coords = get_moon(current_time)
            moon_sep = position.separation(moon_coords)
            
            # Calculate next observable window
            observable_window = self._calculate_observable_window(
                position, location, current_time
            )
            
            results[obs_name] = {
                'current_alt': float(altaz.alt.deg),
                'current_az': float(altaz.az.deg),
                'is_night': bool(is_night),
                'moon_separation': float(moon_sep.deg),
                'is_observable_now': bool(altaz.alt > 30*u.deg and is_night and moon_sep > 30*u.deg),
                'next_window_start': observable_window['start'].iso if observable_window['start'] else None,
                'next_window_end': observable_window['end'].iso if observable_window['end'] else None,
                'observatory': {
                    'name': obs_name,
                    'longitude': float(location.lon.deg),
                    'latitude': float(location.lat.deg),
                    'elevation': float(location.height.to(u.m).value)
                }
            }
        
        # Save observability results
        grb_event['observability'] = results
        self._save_grb_event(grb_event)
        
        # Log results
        observable_now = any(r['is_observable_now'] for r in results.values())
        if observable_now:
            logger.info(f"GRB {grb_event['id']} is currently observable!")
        else:
            logger.info(f"GRB {grb_event['id']} is not currently observable")
        
        return results
    
    def _calculate_observable_window(self, position: SkyCoord, location: EarthLocation,
                                     start_time: Time = None) -> Dict:
        """Calculate next observable window for a target"""
        if start_time is None:
            start_time = Time.now()
        
        # Create time grid (next 24 hours)
        time_grid = start_time + np.linspace(0, 24, 145)*u.hour  # 10-minute intervals
        
        # Initialize observable flags
        is_observable = np.zeros(len(time_grid), dtype=bool)
        
        # Check each time point
        for i, t in enumerate(time_grid):
            frame = AltAz(obstime=t, location=location)
            altaz = position.transform_to(frame)
            sun_altaz = get_sun(t).transform_to(frame)
            moon_coords = get_moon(t)
            moon_sep = position.separation(moon_coords)
            
            is_observable[i] = (
                altaz.alt > 30*u.deg and
                sun_altaz.alt < -12*u.deg and
                moon_sep > 30*u.deg
            )
        
        # Find continuous intervals
        window_start = None
        window_end = None
        
        for i, observable in enumerate(is_observable):
            if observable and window_start is None:
                window_start = time_grid[i]
            elif not observable and window_start is not None:
                window_end = time_grid[i-1]
                break
        
        # If still observable at end of time grid
        if window_start is not None and window_end is None:
            window_end = time_grid[-1]
        
        return {
            'start': window_start,
            'end': window_end
        }
    
    def _check_coincidences(self):
        """Check for temporal and spatial coincidences between GW and GRB events"""
        # Set coincidence window (seconds)
        time_window = self.config.get('coincidence_window', 600)  # 10 minutes
        
        with self.alert_lock:
            # Loop through GW events
            for gw_event in self.gw_cache:
                gw_time = gw_event.get('time')
                if not gw_time:
                    continue
                
                # Loop through GRB events
                for grb_event in self.grb_cache:
                    grb_time = grb_event.get('time')
                    if not grb_time:
                        continue
                    
                    # Check temporal coincidence
                    if isinstance(gw_time, Time) and isinstance(grb_time, Time):
                        dt = abs((gw_time - grb_time).sec)
                    else:
                        # If either time is not a Time object, try to convert from string
                        try:
                            if not isinstance(gw_time, Time):
                                gw_time = Time(gw_time)
                            if not isinstance(grb_time, Time):
                                grb_time = Time(grb_time)
                            dt = abs((gw_time - grb_time).sec)
                        except:
                            dt = float('inf')  # Can't compare, set to infinity
                    
                    if dt <= time_window:
                        # Check spatial coincidence if skymap available
                        spatial_match = False
                        match_prob = 0.0
                        
                        if 'skymap' in gw_event:
                            skymap = gw_event['skymap']
                            
                            if 'coordinates' in grb_event:
                                coords = grb_event['coordinates']
                                error_radius = grb_event.get('error_radius', 1.0)
                                
                                # Evaluate the spatial coincidence
                                match_prob = self._evaluate_spatial_coincidence(
                                    skymap, coords, error_radius
                                )
                                
                                if match_prob > 0.02:  # 2% containment threshold
                                    spatial_match = True
                        
                        if spatial_match or 'skymap' not in gw_event:
                            # Record coincidence
                            coincidence = {
                                'gw_id': gw_event.get('graceid', 'unknown'),
                                'grb_id': grb_event.get('id', 'unknown'),
                                'time_difference': dt,
                                'spatial_match': spatial_match,
                                'spatial_probability': match_prob,
                                'combined_significance': match_prob * (1.0 / (1 + dt/60)),  # Simple metric
                                'detected': datetime.utcnow().isoformat()
                            }
                            
                            self.coincidences.append(coincidence)
                            self.metrics['coincident'] += 1
                            
                            # Save coincidence information
                            self._save_coincidence(coincidence)
                            
                            logger.info(f"Coincidence detected: GW {coincidence['gw_id']} - GRB {coincidence['grb_id']}")
                            logger.info(f"Time difference: {dt:.1f}s, Spatial probability: {match_prob:.4f}")
                            
                            # Generate detailed report
                            self._generate_coincidence_report(coincidence, gw_event, grb_event)
    
    def _evaluate_spatial_coincidence(self, skymap: QTable, coords: SkyCoord,
                                      error_radius: float) -> float:
        """Evaluate spatial coincidence between GW skymap and GRB position"""
        try:
            # Determine HEALPix resolution
            npix = len(skymap)
            
            # Initialize HEALPix with ICRS frame
            from astropy.coordinates import ICRS
            hp = ah.HEALPix(nside=ah.npix_to_nside(npix), order='nested', frame=ICRS())
            
            # Find the closest HEALPix pixel to the GRB position
            ipix = hp.skycoord_to_healpix(coords)
            
            # Get probability at that pixel
            prob = skymap['PROB'][ipix]
            
            # For a more accurate estimate, we should consider the GRB error circle
            # Convert error radius to radians with units
            error_rad = error_radius * u.deg.to(u.rad)
            
            # Get all pixels within the error circle using coordinate-aware method
            disc_ipix = hp.cone_search_skycoord(coords, radius=error_rad*u.rad)
            
            # Sum probability within error circle
            total_prob = np.sum(skymap['PROB'][disc_ipix])
            
            # Return probability mass contained in the error circle
            return float(total_prob)
        except Exception as e:
            logger.error(f"Error evaluating spatial coincidence: {str(e)}")
            return 0.0
    
    def _save_coincidence(self, coincidence: Dict):
        """Save coincidence information to disk"""
        output_path = self.output_dir / f"coincidence_{coincidence['gw_id']}_{coincidence['grb_id']}.json"
        
        try:
            with open(output_path, 'w') as f:
                json.dump(coincidence, f, indent=2)
            logger.info(f"Saved coincidence to {output_path}")
        except Exception as e:
            logger.error(f"Error saving coincidence: {str(e)}")
    
    def _generate_coincidence_report(self, coincidence: Dict, gw_event: Dict, grb_event: Dict):
        """Generate detailed report for coincident events"""
        report = {
            'coincidence': coincidence,
            'gw_event': self._sanitize_for_json(gw_event),
            'grb_event': self._sanitize_for_json(grb_event),
            'analysis': {
                'joint_significance': self._calculate_joint_significance(gw_event, grb_event),
                'false_alarm_rate': self._estimate_false_alarm_rate(coincidence)
            },
            'observability': self._get_observability_for_coincidence(gw_event, grb_event)
        }
        
        # Save detailed report
        output_path = self.output_dir / f"report_{coincidence['gw_id']}_{coincidence['grb_id']}.json"
        
        try:
            with open(output_path, 'w') as f:
                json.dump(report, f, indent=2)
            logger.info(f"Saved detailed report to {output_path}")
        except Exception as e:
            logger.error(f"Error saving report: {str(e)}")
    
    def _sanitize_for_json(self, event_data: Dict) -> Dict:
        """Sanitize event data for JSON serialization"""
        if not isinstance(event_data, dict):
            logger.warning(f"Expected dict, got {type(event_data)}")
            return {}
            
        # Create a deep copy to avoid modifying the original
        result = {}
        
        for key, value in event_data.items():
            # Skip skymap arrays which are large
            if key == 'skymap':
                result[key] = {'available': True, 'size': len(value) if hasattr(value, '__len__') else 'unknown'}
                continue
                
            # Handle special types
            if isinstance(value, Time):
                result[key] = value.iso
            elif isinstance(value, SkyCoord):
                result[key] = {'ra': value.ra.deg, 'dec': value.dec.deg}
            elif isinstance(value, (np.ndarray, QTable, Table)):
                # Summarize array-like objects
                result[key] = {'type': str(type(value)), 'length': len(value) if hasattr(value, '__len__') else 'unknown'}
            elif isinstance(value, dict):
                # Recursively sanitize nested dictionaries
                result[key] = self._sanitize_for_json(value)
            elif isinstance(value, (str, int, float, bool, type(None))):
                # Basic types can be kept as-is
                result[key] = value
            else:
                # Convert other types to strings
                result[key] = str(value)
                
        return result
    
    def _calculate_joint_significance(self, gw_event: Dict, grb_event: Dict) -> float:
        """Calculate joint significance for coincident events"""
        # Default to a simple formula if we don't have sophisticated data
        temporal_term = 1.0
        spatial_term = 0.05  # Default value
        
        # Get GW FAR if available
        gw_far = gw_event.get('far', 1.0)
        if gw_far < 1e-10:  # Avoid division by zero
            gw_far = 1e-10
            
        # Get GRB significance if available
        grb_snr = grb_event.get('snr', 3.0)
        if 'fluence' in grb_event:
            grb_snr = max(3.0, 10.0 * float(grb_event['fluence']) / 1e-7)
        
        # Calculate temporal term
        if 'time' in gw_event and 'time' in grb_event:
            try:
                t1 = Time(gw_event['time']) if not isinstance(gw_event['time'], Time) else gw_event['time']
                t2 = Time(grb_event['time']) if not isinstance(grb_event['time'], Time) else grb_event['time']
                dt = abs((t1 - t2).sec)
                # Temporal term peaks at dt=0, decays with larger dt
                temporal_term = np.exp(-0.5 * (dt / 60.0)**2)
            except Exception as e:
                logger.warning(f"Error calculating temporal term: {str(e)}")
        
        # Calculate spatial term if skymap available
        if 'skymap' in gw_event and 'coordinates' in grb_event:
            try:
                skymap = gw_event['skymap']
                coords = grb_event['coordinates']
                error_radius = grb_event.get('error_radius', 1.0)
                
                spatial_term = self._evaluate_spatial_coincidence(skymap, coords, error_radius)
            except Exception as e:
                logger.warning(f"Error calculating spatial term: {str(e)}")
        
        # Combine terms
        joint_significance = -np.log10(gw_far) * (grb_snr / 5.0) * temporal_term * spatial_term
        
        return float(joint_significance)
    
    def _estimate_false_alarm_rate(self, coincidence: Dict) -> float:
        """Estimate false alarm rate for the coincidence"""
        # Simplified formula based on spatial and temporal coincidence
        spatial_prob = coincidence.get('spatial_probability', 0.01)
        time_diff = coincidence.get('time_difference', 600.0)
        
        # Calculate probability of chance coincidence
        # Assuming GW events at 1 per day, GRB events at 1 per 6 hours
        gw_rate = 1.0 / (24 * 3600)  # per second
        grb_rate = 1.0 / (6 * 3600)  # per second
        
        # Time window (use the configuration value)
        time_window = self.config.get('coincidence_window', 600)  # seconds
        
        # Probability of random coincidence in time
        p_time = (2 * time_window) * gw_rate * grb_rate  # per second^2
        
        # Multiply by spatial coincidence probability
        p_coincidence = p_time * spatial_prob
        
        # Convert to false alarm rate per year
        far_per_year = p_coincidence * (365 * 24 * 3600)
        
        return float(far_per_year)
    
    def _get_observability_for_coincidence(self, gw_event: Dict, grb_event: Dict) -> Dict:
        """Calculate observability for coincident events"""
        # Default to GRB position if available, otherwise use GW skymap
        if 'coordinates' in grb_event:
            position = grb_event['coordinates']
            error_radius = grb_event.get('error_radius', 1.0)
            
            # Calculate observability based on position
            result = {}
            for obs_name, location in self.locations.items():
                window = self._calculate_observable_window(position, location)
                result[obs_name] = {
                    'position_type': 'grb',
                    'ra': position.ra.deg,
                    'dec': position.dec.deg,
                    'error_radius': error_radius,
                    'window_start': window['start'].iso if window['start'] else None,
                    'window_end': window['end'].iso if window['end'] else None
                }
            
            return result
        
        elif 'skymap' in gw_event:
            # Calculate observability based on skymap
            skymap = gw_event['skymap']
            
            result = {}
            for obs_name, location in self.locations.items():
                visibility = self.calculate_visibility(skymap, location)
                
                result[obs_name] = {
                    'position_type': 'skymap',
                    'peak_time': visibility['peak_time'].iso,
                    'peak_probability': float(visibility['peak_prob']),
                    'ra': visibility['best_coords'].ra.deg if visibility['best_coords'] else None,
                    'dec': visibility['best_coords'].dec.deg if visibility['best_coords'] else None
                }
            
            return result
        
        # Fallback if neither is available
        return {'note': 'Insufficient position information for observability calculation'}
    
    def get_metrics(self) -> Dict:
        """Get current alert metrics"""
        with self.alert_lock:
            return self.metrics.copy()
    
    def get_coincidences(self) -> List[Dict]:
        """Get list of detected coincidences"""
        with self.alert_lock:
            return self.coincidences.copy()
    
    def simulate_alerts(self, num_gw=5, num_grb=10, num_coincidences=2):
        """Simulate a set of alerts for testing"""
        try:
            logger.info(f"Generating {num_gw} GW alerts, {num_grb} GRB alerts with {num_coincidences} coincidences")
            
            alerts = []
            now = Time.now()
            
            # Generate GW alerts
            for i in range(num_gw):
                time_offset = np.random.uniform(-12, 12)  # hours
                alert_time = now + time_offset * u.hour
                
                alert = {
                    'time': alert_time,
                    'type': "GW",
                    'id': f"GW_{int(alert_time.mjd)}_{i}",
                    'graceid': f"S{int(alert_time.mjd)}a-{i}",
                    'ra': np.random.uniform(0, 360),
                    'dec': np.random.uniform(-90, 90),
                    'confidence': np.random.uniform(0, 1),
                    'far': 10**(-np.random.uniform(1, 10)),
                    'distance': np.random.uniform(40, 400),
                    'error_region': np.random.uniform(50, 2000),
                    'classification': np.random.choice(["BNS", "BBH", "NSBH", "MassGap", "Terrestrial"])
                }
                
                # Generate a simple mock skymap
                alert['skymap'] = self._create_mock_skymap()
                
                # Add to GW cache
                with self.alert_lock:
                    self.gw_cache.append(alert)
                    self.metrics['gw_events'] += 1
                    self.metrics['processed'] += 1
            
            # Generate GRB alerts (some coincident with GW alerts for multi-messenger testing)
            coincidence_indices = np.random.choice(range(num_gw), size=min(num_coincidences, num_gw), replace=False)
            
            for i in range(num_grb):
                if i < num_coincidences and i < num_gw:
                    # Make this GRB coincident with a GW
                    gw_alert = self.gw_cache[coincidence_indices[i]]
                    gw_time = Time(gw_alert['time']) if not isinstance(gw_alert['time'], Time) else gw_alert['time']
                    
                    # GRB typically happens slightly after GW for BNS mergers
                    time_delay = np.random.uniform(0, 2)  # seconds
                    alert_time = gw_time + time_delay * u.second
                    
                    # Position should be within GW error region
                    ra_offset = np.random.uniform(-5, 5)  # degrees
                    dec_offset = np.random.uniform(-5, 5)  # degrees
                    ra = (gw_alert['ra'] + ra_offset) % 360
                    dec = max(-90, min(90, gw_alert['dec'] + dec_offset))
                else:
                    # Independent GRB
                    time_offset = np.random.uniform(-12, 12)  # hours
                    alert_time = now + time_offset * u.hour
                    ra = np.random.uniform(0, 360)
                    dec = np.random.uniform(-90, 90)
                
                alert = {
                    'time': alert_time,
                    'type': "GRB",
                    'id': f"GRB_{int(alert_time.mjd)}_{i}",
                    'ra': ra,
                    'dec': dec,
                    'coordinates': SkyCoord(ra*u.deg, dec*u.deg),
                    'error_radius': np.random.uniform(0.1, 3.0),
                    'confidence': np.random.uniform(0, 1),
                    't90': np.random.uniform(0.1, 100),
                    'fluence': np.random.uniform(1e-7, 1e-5),
                    'energy_band': "15-150 keV",
                    'detector': np.random.choice(["Swift", "Fermi", "INTEGRAL"])
                }
                
                # Add to GRB cache
                with self.alert_lock:
                    self.grb_cache.append(alert)
                    self.metrics['grb_events'] += 1
                    self.metrics['processed'] += 1
            
            # Check for coincidences
            self._check_coincidences()
            
            return len(self.gw_cache) + len(self.grb_cache)
        
        except Exception as e:
            logger.error(f"Error simulating alerts: {str(e)}")
            raise
