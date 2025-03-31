from astroquery.gaia import Gaia
from astropy.io import fits
import numpy as np
import json
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import curve_fit
import emcee
import corner

class GaiaValidationEngine:
    """GAIA DR3 validation engine for astronomical data"""
    
    def __init__(self):
        self.gaia_table = None
        self.cross_matched = None
        
    def validate_fits_file(self, fits_path):
        """
        Validate a FITS file against GAIA DR3 catalog
        
        Args:
            fits_path (str): Path to FITS file to validate
            
        Returns:
            dict: Validation report with metrics and statistics
            
        Raises:
            FileNotFoundError: If FITS file doesn't exist
            ValueError: If FITS file is invalid
        """
        if not os.path.exists(fits_path):
            raise FileNotFoundError(f"FITS file not found: {fits_path}")
            
        try:
            # Test if file is valid FITS
            with fits.open(fits_path) as test:
                if len(test) < 1:
                    raise ValueError("Invalid FITS file - no HDUs found")
        except Exception as e:
            raise ValueError(f"Error validating FITS file: {str(e)}")
        # Load FITS data
        with fits.open(fits_path) as hdul:
            data = hdul[1].data
            header = hdul[0].header
            
        # Get coordinates from header
        ra = header.get('RA', np.mean(data['RA']))
        dec = header.get('DEC', np.mean(data['DEC']))
        radius = 0.2 * u.deg  # Search radius
        
        # Query GAIA DR3
        coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg))
        self.gaia_table = Gaia.query_object_async(
            coordinate=coord, 
            radius=radius,
            table='gaiadr3.gaia_source'
        )
        
        # Cross-match with input catalog
        self._cross_match(data)
        
        # Perform validations
        report = {
            'astrometry': self._validate_astrometry(),
            'photometry': self._validate_photometry(),
            'matching': self._get_matching_stats(),
            'visualizations': self._generate_plots()
        }
        
        return report
    
    def _cross_match(self, data):
        """Cross-match input catalog with GAIA sources"""
        # Simple positional cross-match within 1 arcsec
        input_coords = SkyCoord(ra=data['RA']*u.deg, dec=data['DEC']*u.deg)
        gaia_coords = SkyCoord(ra=self.gaia_table['ra'], 
                              dec=self.gaia_table['dec'])
        
        idx, sep, _ = input_coords.match_to_catalog_sky(gaia_coords)
        match_mask = sep < 1 * u.arcsec
        
        self.cross_matched = {
            'input_idx': np.where(match_mask)[0],
            'gaia_idx': idx[match_mask],
            'separations': sep[match_mask]
        }
    
    def _validate_astrometry(self):
        """Perform astrometric validation"""
        if not self.cross_matched or len(self.cross_matched['input_idx']) == 0:
            return {'n_matches': 0, 'status': 'No matches found'}
            
        # Calculate positional offsets
        offsets = self.cross_matched['separations'].to(u.arcsec)
        
        # Calculate statistics
        return {
            'total_rms': np.sqrt(np.mean(offsets**2)).value,
            'offset': np.mean(offsets).value,
            'std_dev': np.std(offsets).value,
            'n_matches': len(offsets),
            'status': 'OK'
        }
    
    def _validate_photometry(self):
        """Perform photometric validation"""
        if not self.cross_matched or len(self.cross_matched['input_idx']) == 0:
            return {'n_matches': 0, 'status': 'No matches found'}
            
        # Get magnitudes (simplified - would need actual band mapping)
        input_mags = self.cross_matched['input_idx']['MAG']  # Assuming MAG column
        gaia_mags = self.gaia_table['phot_g_mean_mag'][self.cross_matched['gaia_idx']]
        
        # Calculate zero point offset
        mag_diff = input_mags - gaia_mags
        zp_offset = np.median(mag_diff)
        
        return {
            'zp_offset': zp_offset,
            'mag_offset': np.mean(mag_diff),
            'std_dev': np.std(mag_diff),
            'n_matches': len(mag_diff),
            'status': 'OK'
        }
    
    def _get_matching_stats(self):
        """Get cross-matching statistics"""
        if not self.cross_matched:
            return {'n_matches': 0, 'completeness': 0}
            
        n_input = len(self.cross_matched['input_idx'])
        n_gaia = len(self.gaia_table)
        
        return {
            'n_matches': n_input,
            'completeness': n_input / n_gaia if n_gaia > 0 else 0
        }
    
    def _generate_plots(self):
        """Generate validation plots"""
        plot_paths = []
        
        # Astrometry plot
        if self.cross_matched and len(self.cross_matched['input_idx']) > 0:
            plt.figure()
            plt.hist(self.cross_matched['separations'].arcsec, bins=20)
            plt.xlabel('Offset (arcsec)')
            plt.ylabel('Count')
            plt.title('Positional Offsets')
            plot_paths.append('astrometry_plot.png')
            plt.savefig(plot_paths[-1])
            plt.close()
            
        # Photometry plot
        if self.cross_matched and len(self.cross_matched['input_idx']) > 0:
            plt.figure()
            input_mags = self.cross_matched['input_idx']['MAG']
            gaia_mags = self.gaia_table['phot_g_mean_mag'][self.cross_matched['gaia_idx']]
            plt.scatter(gaia_mags, input_mags - gaia_mags)
            plt.xlabel('GAIA G mag')
            plt.ylabel('Δ mag (input - GAIA)')
            plt.title('Photometric Comparison')
            plot_paths.append('photometry_plot.png')
            plt.savefig(plot_paths[-1])
            plt.close()
            
        return plot_paths
