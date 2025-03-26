# examples/gw170817_analysis.py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
import sys

# Add project root to path
project_root = Path(__file__).parents[1]
sys.path.append(str(project_root))

from data_acquisition.grb_loader import GRBLoader
from data_acquisition.wise_interface import WISEInterface
from analysis.models.kilonova import KilonovaModel
from visualization.lightcurve_plot import plot_multiband_lightcurve

def main():
    print("NASA Project 7: GW170817 Kilonova Analysis")
    print("===========================================")
    
    # Initialize data paths
    data_dir = Path("C:/Users/Dave/PycharmProjects/NasaProject7/data")
    
    # Ensure the directories exist
    (data_dir / "processed").mkdir(exist_ok=True)
    
    # 1. Load GRB data
    grb_loader = GRBLoader()
    grb = grb_loader.load_event('GRB170817A')
    print(f"Loaded GRB: {grb['name']} at position {grb['coord'].to_string('hmsdms')}")
    
    # 2. Query WISE data for the GRB position
    wise = WISEInterface(data_dir)
    ra, dec = grb['coord'].ra.degree, grb['coord'].dec.degree
    
    # Normally we'd query real data, but for this example we'll create synthetic observations
    print(f"Would query WISE data at RA={ra:.4f}, Dec={dec:.4f}")
    print("Using synthetic NIR observations for demonstration")
    
    # Create synthetic observations for GW170817
    # Based on published light curve data from Villar et al. 2017
    nir_data = {
        'time': np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 7.0]),  # days
        'J': np.array([17.2, 17.5, 17.9, 18.3, 19.0, 19.5, 20.8]),
        'H': np.array([17.5, 17.7, 18.0, 18.1, 18.5, 19.0, 20.2]),
        'K': np.array([18.0, 18.1, 18.3, 18.4, 18.7, 19.1, 19.8]),
        'J_err': np.array([0.1, 0.1, 0.15, 0.15, 0.2, 0.25, 0.4]),
        'H_err': np.array([0.1, 0.1, 0.15, 0.15, 0.2, 0.25, 0.4]),
        'K_err': np.array([0.15, 0.15, 0.2, 0.2, 0.25, 0.3, 0.5])
    }
    
    # 3. Initialize kilonova model
    print("Initializing kilonova models")
    
    # Create default model
    default_model = KilonovaModel(
        ejecta_mass=0.05, 
        velocity=0.2,
        opacity=10.0,
        lanthanide_fraction=0.01,
        viewing_angle=30.0,
        model_type='default',
        data_dir=str(data_dir / "reference")
    )
    
    # Create AT2017gfo-specific model
    at2017gfo_model = KilonovaModel(
        ejecta_mass=0.05,
        velocity=0.2,
        opacity=10.0,
        lanthanide_fraction=0.01,
        viewing_angle=30.0,
        model_type='AT2017gfo',
        data_dir=str(data_dir / "reference")
    )
    
    # Create BH-NS model for comparison
    bhns_model = KilonovaModel(
        ejecta_mass=0.1,
        velocity=0.15,
        opacity=15.0,
        lanthanide_fraction=0.1,
        viewing_angle=30.0,
        model_type='BH-NS',
        data_dir=str(data_dir / "reference")
    )
    
    # 4. Fit default model to data
    print("Fitting kilonova model to observations")
    default_model.fit(nir_data, params_to_fit=['ejecta_mass', 'velocity', 'lanthanide_fraction'])
    
    print(f"Best-fit parameters:")
    for param, value in default_model.best_fit_params.items():
        print(f"  {param}: {value:.4f}")
    print(f"Chi-squared: {default_model.chi2:.2f}")
    
    # 5. Generate light curves
    print("Generating light curves")
    time_grid = np.linspace(0.1, 10.0, 100)
    
    # Calculate light curves for all models
    default_lc = default_model.light_curve(time_grid, bands=['J', 'H', 'K'])
    at2017gfo_lc = at2017gfo_model.light_curve(time_grid, bands=['J', 'H', 'K'])
    bhns_lc = bhns_model.light_curve(time_grid, bands=['J', 'H', 'K'])
    
    # 6. Visualize results
    print("Creating visualizations")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Plot J-band
    ax = axes[0]
    ax.errorbar(nir_data['time'], nir_data['J'], yerr=nir_data['J_err'], 
                fmt='o', color='blue', label='GW170817 J-band')
    ax.plot(time_grid, default_lc['J'], 'b-', label='Best-fit model')
    ax.plot(time_grid, at2017gfo_lc['J'], 'b--', label='AT2017gfo template')
    ax.plot(time_grid, bhns_lc['J'], 'b:', label='BH-NS model')
    ax.set_xlabel('Time since merger (days)')
    ax.set_ylabel('J-band magnitude')
    ax.set_ylim(22, 16)  # Inverted y-axis
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title('J-band')
    
    # Plot H-band
    ax = axes[1]
    ax.errorbar(nir_data['time'], nir_data['H'], yerr=nir_data['H_err'], 
                fmt='o', color='green', label='GW170817 H-band')
    ax.plot(time_grid, default_lc['H'], 'g-', label='Best-fit model')
    ax.plot(time_grid, at2017gfo_lc['H'], 'g--', label='AT2017gfo template')
    ax.plot(time_grid, bhns_lc['H'], 'g:', label='BH-NS model')
    ax.set_xlabel('Time since merger (days)')
    ax.set_ylabel('H-band magnitude')
    ax.set_ylim(22, 16)  # Inverted y-axis
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title('H-band')
    
    # Plot K-band
    ax = axes[2]
    ax.errorbar(nir_data['time'], nir_data['K'], yerr=nir_data['K_err'], 
                fmt='o', color='red', label='GW170817 K-band')
    ax.plot(time_grid, default_lc['K'], 'r-', label='Best-fit model')
    ax.plot(time_grid, at2017gfo_lc['K'], 'r--', label='AT2017gfo template')
    ax.plot(time_grid, bhns_lc['K'], 'r:', label='BH-NS model')
    ax.set_xlabel('Time since merger (days)')
    ax.set_ylabel('K-band magnitude')
    ax.set_ylim(22, 16)  # Inverted y-axis
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.set_title('K-band')
    
    plt.tight_layout()
    plt.savefig(data_dir / "processed" / "gw170817_nir_analysis.png", dpi=300)
    
    print(f"Analysis complete. Results saved to {data_dir / 'processed'}")

if __name__ == "__main__":
    main()
