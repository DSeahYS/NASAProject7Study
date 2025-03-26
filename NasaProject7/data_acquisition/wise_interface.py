# data_acquisition/wise_interface.py
import os
import duckdb
import numpy as np
from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from astropy.utils.data import get_readable_fileobj

class WISEInterface:
    """Interface for accessing WISE/NEOWISE NIR data"""
    
    def __init__(self, data_dir):
        """
        Initialize WISE data interface
        
        Parameters
        ----------
        data_dir : str
            Path to directory containing WISE data files
        """
        self.data_dir = Path(data_dir)
        self.conn = self._create_query_engine()
        
    def _create_query_engine(self):
        """Initialize DuckDB in-memory database with astronomical extensions"""
        conn = duckdb.connect(database=':memory:', config={'threads': '8'})
        
        # Setup astronomical extensions
        try:
            conn.execute("INSTALL astronomer; LOAD astronomer;")
            print("Loaded astronomer extension for DuckDB")
        except:
            print("Warning: astronomer extension not available, some features will be limited")
            
        # Register data directory
        if self.data_dir.exists():
            conn.execute(f"""
                CREATE VIEW wise_files AS 
                SELECT * FROM read_csv('{self.data_dir}/wise_files_index.csv', 
                                      auto_detect=true);
            """)
            
        return conn
    
    def query_by_position(self, ra, dec, radius_arcsec=300, time_start=None, time_end=None):
        """
        Query WISE data at specified sky position
        
        Parameters
        ----------
        ra : float
            Right ascension in degrees
        dec : float
            Declination in degrees
        radius_arcsec : float
            Search radius in arcseconds
        time_start : str or astropy.time.Time
            Start time for query (ISO format)
        time_end : str or astropy.time.Time
            End time for query (ISO format)
            
        Returns
        -------
        astropy.table.Table
            Table of matching WISE sources
        """
        # Convert radius to degrees
        radius_deg = radius_arcsec / 3600.0
        
        # Build query conditions
        where_clauses = [
            f"ra BETWEEN {ra-radius_deg} AND {ra+radius_deg}",
            f"dec BETWEEN {dec-radius_deg} AND {dec+radius_deg}"
        ]
        
        # Add time constraints if provided
        if time_start:
            where_clauses.append(f"mjd >= '{time_start}'")
        if time_end:
            where_clauses.append(f"mjd <= '{time_end}'")
        
        # Build final query
        where_clause = " AND ".join(where_clauses)
        
        # First query the file index to find relevant files
        try:
            files_result = self.conn.execute(f"""
                SELECT filepath 
                FROM wise_files
                WHERE {where_clause}
            """).fetchall()
            
            file_paths = [row[0] for row in files_result]
            
            # If no matching files found via database, fall back to file system search
            if not file_paths:
                file_paths = list(self.data_dir.glob("*.fits"))
                
            # Process each file and combine results
            results = []
            for filepath in file_paths:
                try:
                    with fits.open(filepath) as hdul:
                        data = Table(hdul[1].data)
                        
                        # Filter by position more precisely
                        pos_mask = self._position_filter(data, ra, dec, radius_deg)
                        filtered_data = data[pos_mask]
                        
                        # Apply time filter if needed
                        if time_start or time_end:
                            time_mask = self._time_filter(filtered_data, time_start, time_end)
                            filtered_data = filtered_data[time_mask]
                            
                        if len(filtered_data) > 0:
                            results.append(filtered_data)
                except Exception as e:
                    print(f"Error processing file {filepath}: {e}")
            
            # Combine results
            if results:
                return Table.vstack(results)
            else:
                return Table()
                
        except Exception as e:
            print(f"Database query failed: {e}")
            return Table()
    
    def _position_filter(self, data, ra, dec, radius_deg):
        """Create mask for sources within radius of position"""
        try:
            # Use precise spherical trigonometry for position matching
            ra_rad = np.radians(data['ra'])
            dec_rad = np.radians(data['dec'])
            target_ra_rad = np.radians(ra)
            target_dec_rad = np.radians(dec)
            
            # Haversine formula
            dlon = ra_rad - target_ra_rad
            dlat = dec_rad - target_dec_rad
            a = np.sin(dlat/2)**2 + np.cos(dec_rad) * np.cos(target_dec_rad) * np.sin(dlon/2)**2
            dist_rad = 2 * np.arcsin(np.sqrt(a))
            dist_deg = np.degrees(dist_rad)
            
            return dist_deg <= radius_deg
        except:
            # Fall back to simpler filter if columns don't match expected names
            ra_key = 'ra' if 'ra' in data.colnames else 'RA'
            dec_key = 'dec' if 'dec' in data.colnames else 'DEC'
            
            ra_diff = np.abs(data[ra_key] - ra)
            # Handle RA wrap-around at 0/360
            ra_diff = np.minimum(ra_diff, 360 - ra_diff)
            dec_diff = np.abs(data[dec_key] - dec)
            
            # Simple box filter
            return (ra_diff <= radius_deg) & (dec_diff <= radius_deg)
    
    def _time_filter(self, data, time_start, time_end):
        """Create mask for sources within time range"""
        time_key = 'mjd' if 'mjd' in data.colnames else 'MJD'
        
        mask = np.ones(len(data), dtype=bool)
        
        if time_start:
            if isinstance(time_start, str):
                from astropy.time import Time
                time_start = Time(time_start).mjd
            mask &= data[time_key] >= time_start
            
        if time_end:
            if isinstance(time_end, str):
                from astropy.time import Time
                time_end = Time(time_end).mjd
            mask &= data[time_key] <= time_end
            
        return mask
