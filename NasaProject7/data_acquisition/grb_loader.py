# data_acquisition/grb_loader.py
import os
import json
import requests
from pathlib import Path
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
import pandas as pd
import logging
from typing import Union, Dict, List
import sqlite3
from astroquery.heasarc import Heasarc

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('GRBLoader')

class GRBLoader:
    """Modern GRB loader with multiple fallback sources and local caching"""
    
    def __init__(self, catalog_path: str = None, config_path: str = None):
        self.catalog = None
        self.config = self._load_config(config_path)
        self.local_catalog_path = catalog_path or self._get_default_catalog_path()
        self.heasarc = Heasarc()
        self._init_local_cache()
        
        # Initialize catalogs
        self._load_local_catalog()
        self._init_online_sources()

    def _init_local_cache(self):
        """Initialize SQLite cache database"""
        self.cache_db = Path(self.config['data_dir']) / 'grb_cache.db'
        with sqlite3.connect(self.cache_db) as conn:
            conn.execute('''CREATE TABLE IF NOT EXISTS grbs
                         (name TEXT PRIMARY KEY,
                          data TEXT,
                          timestamp DATETIME DEFAULT CURRENT_TIMESTAMP)''')

    def _load_config(self, config_path: str) -> Dict:
        """Load updated pipeline configuration"""
        default_config = {
            'data_dir': str(Path(__file__).parent.parent / 'data'),
            'catalogs': {
                'swift': 'https://swift.gsfc.nasa.gov/archive/grb_table/current_table.html',
                'fermi': 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/',
                'heasarc_tap': 'https://heasarc.gsfc.nasa.gov/xamin/vo/tap'
            }
        }
        return default_config

    def _get_default_catalog_path(self) -> Path:
        """Get default path for GRB catalog file"""
        return Path(self.config['data_dir']) / 'grb_catalog.json'

    def _load_local_catalog(self):
        """Load GRB catalog from local JSON file"""
        try:
            if self.local_catalog_path.exists():
                with open(self.local_catalog_path) as f:
                    self.catalog = json.load(f)
                logger.info(f"Loaded local GRB catalog from {self.local_catalog_path}")
            else:
                self.catalog = {}
                logger.info("No local GRB catalog found - initializing empty catalog")
        except Exception as e:
            self.catalog = {}
            logger.error(f"Failed to load local catalog: {e} - initializing empty catalog")

    def _init_online_sources(self):
        """Initialize configuration for online data sources"""
        self.online_sources = {
            'swift': {
                'url': self.config['catalogs']['swift'],
                'active': True
            },
            'fermi': {
                'url': self.config['catalogs']['fermi'],
                'active': True
            },
            'heasarc': {
                'url': self.config['catalogs']['heasarc_tap'],
                'active': True
            }
        }
        logger.info("Initialized online GRB data sources")

    def _query_heasarc_tap(self, grb_name: str) -> Dict:
        """Query HEASARC for GRB data"""
        try:
            # Use standard query method instead of TAP
            table = self.heasarc.query_object(grb_name, mission='fermigbrst')
            if len(table) > 0:
                return self._parse_heasarc_data(table[0])
        except Exception as e:
            logger.error(f"HEASARC query failed: {e}")
        return None

    def _query_swift(self, grb_name: str) -> Dict:
        """Query Swift BAT catalog for GRB data"""
        try:
            response = requests.get(self.online_sources['swift']['url'])
            response.raise_for_status()
            # Parse HTML table response here
            # This is a simplified placeholder - actual implementation would parse the HTML
            return {
                'name': grb_name,
                'ra': 0.0,
                'dec': 0.0,
                'error_radius': 0.0,
                'time': Time.now(),
                'duration': 0.0,
                'mission': 'Swift',
                'instrument': 'BAT',
                'redshift': 0.0
            }
        except Exception as e:
            logger.error(f"Swift query failed: {e}")
        return None

    def _parse_heasarc_data(self, data: Dict) -> Dict:
        """Parse HEASARC TAP response"""
        return {
            'name': data['trigger_name'],
            'ra': float(data['ra']),
            'dec': float(data['dec']),
            'error_radius': float(data['error_radius']),
            'time': Time(data['trigger_time'], format='isot'),
            'duration': float(data['t90']),
            'mission': 'Fermi',
            'instrument': 'GBM',
            'redshift': 0.0,
            'associated_gw': data.get('associated_gw')
        }

    def load_event(self, grb_name: str) -> Dict:
        """Enhanced GRB loading with multiple fallbacks"""
        # Check local cache first
        if cached := self._get_cached_grb(grb_name):
            return cached
            
        # Try updated HEASARC TAP service
        if heasarc_data := self._query_heasarc_tap(grb_name):
            self._cache_grb(grb_name, heasarc_data)
            return heasarc_data
            
        # Fallback to Swift BAT catalog
        if swift_data := self._query_swift(grb_name):
            self._cache_grb(grb_name, swift_data)
            return swift_data
            
        raise ValueError(f"GRB {grb_name} not found in any available sources")

    def _get_cached_grb(self, grb_name: str) -> Union[Dict, None]:
        """Retrieve from SQLite cache"""
        with sqlite3.connect(self.cache_db) as conn:
            result = conn.execute('''SELECT data FROM grbs 
                                  WHERE name = ?''', (grb_name,)).fetchone()
            if result:
                return json.loads(result[0])
        return None

    def _cache_grb(self, grb_name: str, data: Dict):
        """Store in SQLite cache"""
        with sqlite3.connect(self.cache_db) as conn:
            conn.execute('''INSERT OR REPLACE INTO grbs (name, data)
                         VALUES (?, ?)''', 
                         (grb_name, json.dumps(data)))

    # Other methods remain similar but with improved error handling
    # and updated API endpoints

if __name__ == "__main__":
    loader = GRBLoader()
    try:
        grb = loader.load_event('GRB170817A')
        print(f"Loaded GRB: {grb['name']}")
        print(f"Coordinates: {grb['coord'].to_string('hmsdms')}")
        print(f"Time: {grb['time'].isot}")
    except Exception as e:
        logger.error(f"Critical error: {e}")
        logger.info("Attempting to use archival data...")
        # Additional fallback logic here
