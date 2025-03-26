import numpy as np
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import median_filter, gaussian_filter

class NIRBackground:
    """Class for NIR-specific background estimation and subtraction"""
    
    def __init__(self, method='median', mask_threshold=3.0, filter_size=5):
        """
        Initialize NIR background handler
        
        Parameters
        ----------
        method : str
            Background estimation method ('median', 'mesh', or 'polynomial')
        mask_threshold : float
            Threshold for source masking in sigma units
        filter_size : int
            Size of the filter for mesh method
        """
        self.method = method
        self.mask_threshold = mask_threshold
        self.filter_size = filter_size
    
    def create_source_mask(self, image_data):
        """Create a mask of detected sources in the image"""
        mean, median, std = sigma_clipped_stats(image_data)
        threshold = median + (self.mask_threshold * std)
        return image_data > threshold
    
    def estimate_background(self, image_data, source_mask=None):
        """
        Estimate background level in NIR image
        
        Parameters
        ----------
        image_data : numpy.ndarray
            2D image array
        source_mask : numpy.ndarray, optional
            Boolean mask of source pixels to exclude
            
        Returns
        -------
        numpy.ndarray
            Estimated background map
        """
        # Create source mask if not provided
        if source_mask is None:
            source_mask = self.create_source_mask(image_data)
        
        # Make a copy of the data for masking
        masked_data = image_data.copy()
        masked_data[source_mask] = np.nan
        
        if self.method == 'median':
            # Simple median background
            mean, median, std = sigma_clipped_stats(image_data, mask=source_mask)
            background = np.ones_like(image_data) * median
            
        elif self.method == 'mesh':
            # Mesh-based background estimation
            # Fill masked values with local median for filtering
            temp_data = masked_data.copy()
            temp_data[np.isnan(temp_data)] = np.nanmedian(temp_data)
            
            # Apply median filter
            background = median_filter(temp_data, size=self.filter_size)
            
            # Smooth the background
            background = gaussian_filter(background, sigma=self.filter_size/2)
            
        elif self.method == 'polynomial':
            # Polynomial background fitting
            y, x = np.indices(image_data.shape)
            y_flat = y[~source_mask].flatten()
            x_flat = x[~source_mask].flatten()
            z_flat = masked_data[~source_mask].flatten()
            
            # Remove NaN values
            valid = ~np.isnan(z_flat)
            y_fit = y_flat[valid]
            x_fit = x_flat[valid]
            z_fit = z_flat[valid]
            
            # Fit polynomial (degree 2)
            from sklearn.preprocessing import PolynomialFeatures
            from sklearn.linear_model import LinearRegression
            
            # Create polynomial features
            poly = PolynomialFeatures(degree=2)
            X = poly.fit_transform(np.column_stack((x_fit, y_fit)))
            
            # Fit linear regression model
            model = LinearRegression()
            model.fit(X, z_fit)
            
            # Predict background for all pixels
            X_all = poly.transform(np.column_stack((x.flatten(), y.flatten())))
            background = model.predict(X_all).reshape(image_data.shape)
            
        else:
            raise ValueError(f"Unknown background method: {self.method}")
        
        return background
    
    def subtract_background(self, image_data, return_background=False):
        """
        Subtract estimated background from image
        
        Parameters
        ----------
        image_data : numpy.ndarray
            2D image array
        return_background : bool
            If True, return both the subtracted image and the background
            
        Returns
        -------
        numpy.ndarray or tuple
            Background-subtracted image, or tuple of (subtracted_image, background)
        """
        background = self.estimate_background(image_data)
        subtracted = image_data - background
        
        if return_background:
            return subtracted, background
        return subtracted
