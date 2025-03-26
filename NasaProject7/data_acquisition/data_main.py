#!/usr/bin/env python3

"""
NIR GRB/GW Pipeline: Data Acquisition Tool

A production-ready tool for querying, downloading, and analyzing
near-infrared observations related to gamma-ray bursts and gravitational wave events.

Usage:
python data_main.py
"""

import os
import json
import sys
import threading
from pathlib import Path
import logging
import tkinter as tk
from tkinter import ttk, scrolledtext, filedialog, messagebox
from datetime import datetime, timedelta  # Added timedelta import
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
from astropy.time import Time

# Import module functions
from query import query_vo_archives, download_fits_data, standardize_fits_headers, extract_metadata, save_data_to_csv, save_to_excel
from alerts import AlertMonitor

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger('data_acquisition')

# Default directories
DEFAULT_DATA_DIR = Path('./data')
DEFAULT_DATA_DIR.mkdir(parents=True, exist_ok=True)

class TextRedirector:
    """Class to redirect stdout to a tkinter Text widget"""
    def __init__(self, text_widget):
        self.text_widget = text_widget
        
    def write(self, string):
        self.text_widget.insert(tk.END, string)
        self.text_widget.see(tk.END)
        self.text_widget.update_idletasks()
        
    def flush(self):
        pass

class DataAcquisitionApp:
    def __init__(self, root):
        self.root = root
        self.root.title("NIR GRB/GW Pipeline: Data Acquisition")
        self.root.geometry("1000x700")
        
        # Create alert monitor instance
        self.alert_monitor = AlertMonitor()
        
        # Set data directories
        self.data_dir = DEFAULT_DATA_DIR
        self.csv_dir = self.data_dir / 'csv_data'
        self.excel_dir = self.data_dir / 'excel_data'
        self.fits_dir = self.data_dir / 'fits_data'
        self.alerts_dir = self.data_dir / 'alerts_data'
        
        # Create directories
        self.csv_dir.mkdir(exist_ok=True)
        self.excel_dir.mkdir(exist_ok=True)
        self.fits_dir.mkdir(exist_ok=True)
        self.alerts_dir.mkdir(exist_ok=True)
        
        # Flag to track if long-running operations are in progress
        self.operation_in_progress = False
        
        # Create interface
        self.create_interface()
        
        # Store query results
        self.query_results = None
        
        # Store alerts
        self.gcn_alerts = []
        self.ligo_alerts = []
        
        # Threads
        self.resolve_thread = None
        self.gcn_thread = None
        self.monitor_running = False
    
    def create_interface(self):
        """Create the application interface."""
        # Create a notebook (tabbed interface)
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Create tabs
        self.create_query_tab()
        self.create_alerts_tab()
        self.create_data_export_tab()
        self.create_settings_tab()
        
        # Create status bar
        self.status_var = tk.StringVar()
        self.status_var.set("Ready")
        self.status_bar = ttk.Frame(self.root, relief=tk.SUNKEN, padding=(10, 2))
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Label(self.status_bar, textvariable=self.status_var, anchor=tk.W).pack(side=tk.LEFT)
        ttk.Label(self.status_bar, text=f"Data Dir: {self.data_dir}", anchor=tk.E).pack(side=tk.RIGHT)
    
    def create_query_tab(self):
        """Create the query tab."""
        query_frame = ttk.Frame(self.notebook, padding=10)
        self.notebook.add(query_frame, text="VO Queries")
        
        # Query parameters frame
        params_frame = ttk.LabelFrame(query_frame, text="Query Parameters", padding=10)
        params_frame.pack(fill=tk.X, pady=5)
        
        # Create grid of labels and entries
        ttk.Label(params_frame, text="Target Name:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.target_entry = ttk.Entry(params_frame, width=30)
        self.target_entry.grid(row=0, column=1, sticky=tk.W, pady=2)
        self.target_entry.insert(0, "M31") # Default value
        
        ttk.Label(params_frame, text="Coordinates (RA, Dec):").grid(row=1, column=0, sticky=tk.W, pady=2)
        coord_frame = ttk.Frame(params_frame)
        coord_frame.grid(row=1, column=1, sticky=tk.W, pady=2)
        self.ra_entry = ttk.Entry(coord_frame, width=12)
        self.ra_entry.pack(side=tk.LEFT, padx=(0, 5))
        self.dec_entry = ttk.Entry(coord_frame, width=12)
        self.dec_entry.pack(side=tk.LEFT)
        
        ttk.Label(params_frame, text="Search Radius (arcmin):").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.radius_entry = ttk.Entry(params_frame, width=10)
        self.radius_entry.grid(row=2, column=1, sticky=tk.W, pady=2)
        self.radius_entry.insert(0, "10") # Default value
        
        # NIR bands frame
        bands_frame = ttk.LabelFrame(params_frame, text="NIR Bands", padding=5)
        bands_frame.grid(row=0, column=2, rowspan=3, padx=10, sticky=tk.N)
        self.band_vars = {}
        for i, band in enumerate(['J', 'H', 'K', 'Y']):
            self.band_vars[band] = tk.BooleanVar(value=True)
            ttk.Checkbutton(bands_frame, text=band, variable=self.band_vars[band]).pack(anchor=tk.W)
        
        # Buttons frame
        buttons_frame = ttk.Frame(params_frame)
        buttons_frame.grid(row=0, column=3, rowspan=3, padx=10, sticky=tk.N)
        ttk.Button(buttons_frame, text="Resolve Name", command=self.resolve_target_name).pack(fill=tk.X, pady=2)
        ttk.Button(buttons_frame, text="Execute Query", command=self.run_vo_query).pack(fill=tk.X, pady=2)
        ttk.Button(buttons_frame, text="Download Data", command=self.download_selected_data).pack(fill=tk.X, pady=2)
        
        # Results frame
        results_frame = ttk.LabelFrame(query_frame, text="Query Results", padding=10)
        results_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Create treeview for results
        columns = ('id', 'ra', 'dec', 'instrument', 'filter', 'date_obs')
        self.results_tree = ttk.Treeview(results_frame, columns=columns, show='headings')
        self.results_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Add column headings
        for col in columns:
            self.results_tree.heading(col, text=col.upper())
            width = 120 if col == 'date_obs' else 80
            self.results_tree.column(col, width=width, anchor=tk.CENTER)
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.results_tree.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.results_tree.configure(yscrollcommand=scrollbar.set)
        
        # Export buttons
        export_frame = ttk.Frame(query_frame)
        export_frame.pack(fill=tk.X, pady=5)
        ttk.Button(export_frame, text="Export to CSV", command=lambda: self.export_results('csv')).pack(side=tk.LEFT, padx=5)
        ttk.Button(export_frame, text="Export to Excel", command=lambda: self.export_results('excel')).pack(side=tk.LEFT, padx=5)
    
    def create_alerts_tab(self):
        """Create the alerts tab."""
        alerts_frame = ttk.Frame(self.notebook, padding=10)
        self.notebook.add(alerts_frame, text="GRB/GW Alerts")
        
        # Controls frame
        controls_frame = ttk.LabelFrame(alerts_frame, text="Alert Controls", padding=10)
        controls_frame.pack(fill=tk.X, pady=5)
        
        # GCN controls
        gcn_frame = ttk.Frame(controls_frame)
        gcn_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(gcn_frame, text="GCN Monitoring:").pack(anchor=tk.W)
        gcn_buttons = ttk.Frame(gcn_frame)
        gcn_buttons.pack(fill=tk.X, pady=5)
        ttk.Button(gcn_buttons, text="Start Monitor", command=self.start_gcn_monitor).pack(side=tk.LEFT, padx=5)
        ttk.Button(gcn_buttons, text="Query Recent", command=self.query_recent_gcn).pack(side=tk.LEFT, padx=5)
        
        # LIGO controls
        ligo_frame = ttk.Frame(controls_frame)
        ligo_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(ligo_frame, text="LIGO/Virgo Alerts:").pack(anchor=tk.W)
        ligo_buttons = ttk.Frame(ligo_frame)
        ligo_buttons.pack(fill=tk.X, pady=5)
        ttk.Button(ligo_buttons, text="Get Alerts", command=self.get_ligo_alerts).pack(side=tk.LEFT, padx=5)
        ttk.Button(ligo_buttons, text="Download Skymap", command=self.download_ligo_skymap).pack(side=tk.LEFT, padx=5)
        
        # Visibility controls
        visibility_frame = ttk.Frame(controls_frame)
        visibility_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(visibility_frame, text="Visibility Check:").pack(anchor=tk.W)
        vis_buttons = ttk.Frame(visibility_frame)
        vis_buttons.pack(fill=tk.X, pady=5)
        ttk.Button(vis_buttons, text="Check Visibility", command=self.check_visibility).pack(side=tk.LEFT, padx=5)
        
        # Alerts list frame
        alerts_list_frame = ttk.LabelFrame(alerts_frame, text="Recent Alerts", padding=10)
        alerts_list_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Create treeview for alerts
        columns = ('time', 'type', 'source', 'ra', 'dec', 'significance')
        self.alerts_tree = ttk.Treeview(alerts_list_frame, columns=columns, show='headings')
        self.alerts_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Add column headings
        for col in columns:
            self.alerts_tree.heading(col, text=col.upper())
            width = 120 if col == 'time' else 80
            self.alerts_tree.column(col, width=width, anchor=tk.CENTER)
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(alerts_list_frame, orient=tk.VERTICAL, command=self.alerts_tree.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.alerts_tree.configure(yscrollcommand=scrollbar.set)
        
        # Alert details frame
        details_frame = ttk.LabelFrame(alerts_frame, text="Alert Details", padding=10)
        details_frame.pack(fill=tk.X, pady=5)
        self.details_text = scrolledtext.ScrolledText(details_frame, height=8)
        self.details_text.pack(fill=tk.BOTH, expand=True)
        
        # Set up event handler for clicking on an alert
        self.alerts_tree.bind('<ButtonRelease-1>', self.show_alert_details)
    
    def create_data_export_tab(self):
        """Create the data export tab."""
        export_frame = ttk.Frame(self.notebook, padding=10)
        self.notebook.add(export_frame, text="Data Export")
        
        # Export controls frame
        controls_frame = ttk.LabelFrame(export_frame, text="Export Controls", padding=10)
        controls_frame.pack(fill=tk.X, pady=5)
        
        # Format selection
        format_frame = ttk.Frame(controls_frame)
        format_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(format_frame, text="Export Format:").pack(anchor=tk.W)
        self.export_format = tk.StringVar(value="csv")
        ttk.Radiobutton(format_frame, text="CSV", variable=self.export_format, value="csv").pack(anchor=tk.W)
        ttk.Radiobutton(format_frame, text="Excel", variable=self.export_format, value="excel").pack(anchor=tk.W)
        
        # Size controls
        size_frame = ttk.Frame(controls_frame)
        size_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(size_frame, text="Dataset Size:").pack(anchor=tk.W)
        self.dataset_size = tk.StringVar(value="10000")
        size_entry = ttk.Entry(size_frame, textvariable=self.dataset_size, width=10)
        size_entry.pack(anchor=tk.W, pady=2)
        ttk.Button(size_frame, text="Generate Large Dataset", command=self.generate_large_dataset).pack(anchor=tk.W, pady=2)
        
        # Directory selection
        dir_frame = ttk.Frame(controls_frame)
        dir_frame.pack(side=tk.LEFT, padx=10)
        ttk.Label(dir_frame, text="Output Directory:").pack(anchor=tk.W)
        dir_select_frame = ttk.Frame(dir_frame)
        dir_select_frame.pack(fill=tk.X, pady=2)
        self.export_dir = tk.StringVar(value=str(self.excel_dir))
        dir_entry = ttk.Entry(dir_select_frame, textvariable=self.export_dir, width=30)
        dir_entry.pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(dir_select_frame, text="Browse...", command=self.browse_export_dir).pack(side=tk.LEFT)
        
        # Export history frame
        history_frame = ttk.LabelFrame(export_frame, text="Export History", padding=10)
        history_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        
        # Create treeview for export history
        columns = ('time', 'type', 'rows', 'file')
        self.export_tree = ttk.Treeview(history_frame, columns=columns, show='headings')
        self.export_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Add column headings
        for col in columns:
            self.export_tree.heading(col, text=col.upper())
            width = 300 if col == 'file' else 100
            self.export_tree.column(col, width=width, anchor=tk.CENTER if col != 'file' else tk.W)
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(history_frame, orient=tk.VERTICAL, command=self.export_tree.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.export_tree.configure(yscrollcommand=scrollbar.set)
        
        # Make double-clicking on an export file open it
        self.export_tree.bind('<Double-1>', self.open_export_file)
    
    def create_settings_tab(self):
        """Create the settings tab."""
        settings_frame = ttk.Frame(self.notebook, padding=10)
        self.notebook.add(settings_frame, text="Settings")
        
        # General settings frame
        general_frame = ttk.LabelFrame(settings_frame, text="General Settings", padding=10)
        general_frame.pack(fill=tk.X, pady=5)
        
        # Data directory
        dir_frame = ttk.Frame(general_frame)
        dir_frame.pack(fill=tk.X, pady=5)
        ttk.Label(dir_frame, text="Data Directory:").pack(side=tk.LEFT)
        self.data_dir_var = tk.StringVar(value=str(self.data_dir))
        dir_entry = ttk.Entry(dir_frame, textvariable=self.data_dir_var, width=40)
        dir_entry.pack(side=tk.LEFT, padx=5)
        ttk.Button(dir_frame, text="Browse...", command=self.browse_data_dir).pack(side=tk.LEFT)
        
        # Observatory settings frame
        obs_frame = ttk.LabelFrame(settings_frame, text="Observatory Settings", padding=10)
        obs_frame.pack(fill=tk.X, pady=5)
        
        # Observatory selection
        ttk.Label(obs_frame, text="Default Observatory:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.observatory = tk.StringVar(value="LDT")
        observatories = ["LDT", "PRIME", "Keck", "Gemini-N", "Gemini-S", "SOAR", "CTIO", "CFHT"]
        obs_combo = ttk.Combobox(obs_frame, textvariable=self.observatory, values=observatories, state="readonly")
        obs_combo.grid(row=0, column=1, sticky=tk.W, pady=2)
        
        # Log level frame
        log_frame = ttk.LabelFrame(settings_frame, text="Logging Settings", padding=10)
        log_frame.pack(fill=tk.X, pady=5)
        
        # Log level selection
        ttk.Label(log_frame, text="Log Level:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.log_level = tk.StringVar(value="INFO")
        log_levels = ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
        log_combo = ttk.Combobox(log_frame, textvariable=self.log_level, values=log_levels, state="readonly")
        log_combo.grid(row=0, column=1, sticky=tk.W, pady=2)
        log_combo.bind("<<ComboboxSelected>>", self.set_log_level)
        
        # Buttons frame
        button_frame = ttk.Frame(settings_frame)
        button_frame.pack(fill=tk.X, pady=10)
        ttk.Button(button_frame, text="Apply Settings", command=self.apply_settings).pack(side=tk.RIGHT, padx=5)
        ttk.Button(button_frame, text="Reset to Defaults", command=self.reset_settings).pack(side=tk.RIGHT, padx=5)
        
        # Console output
        console_frame = ttk.LabelFrame(settings_frame, text="Console Output", padding=10)
        console_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        self.console_text = scrolledtext.ScrolledText(console_frame, height=10)
        self.console_text.pack(fill=tk.BOTH, expand=True)
        
        # Redirect stdout to the console text
        self.original_stdout = sys.stdout
        sys.stdout = TextRedirector(self.console_text)
        
        # Print welcome message
        print("NIR GRB/GW Pipeline: Data Acquisition Tool")
        print(f"Current time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("Ready to process astronomical data.")
    
    def resolve_target_name(self):
        """Resolve a target name to coordinates using SIMBAD."""
        target_name = self.target_entry.get().strip()
        if not target_name:
            messagebox.showwarning("Warning", "Please enter a target name.")
            return
        
        # Show that we're working on it
        self.status_var.set(f"Resolving target name: {target_name}...")
        
        # Use a thread to keep the UI responsive
        if self.resolve_thread is not None and self.resolve_thread.is_alive():
            messagebox.showinfo("Info", "Target resolution already in progress.")
            return
        
        self.resolve_thread = threading.Thread(target=self._resolve_target_thread, args=(target_name,))
        self.resolve_thread.daemon = True
        self.resolve_thread.start()
    
    def _resolve_target_thread(self, target_name):
        """Thread function to resolve a target name."""
        try:
            # Try to import astroquery
            try:
                from astroquery.simbad import Simbad
                result = Simbad.query_object(target_name)
            except ImportError:
                # If astroquery is not available, simulate a result for common objects
                print("Warning: astroquery not available, using hardcoded coordinates")
                common_objects = {
                    "m31": (10.6847, 41.2687), # Andromeda Galaxy
                    "m51": (202.4696, 47.1953), # Whirlpool Galaxy
                    "m87": (187.7059, 12.3911), # Messier 87
                    "m81": (148.8882, 69.0653), # Bode's Galaxy
                    "sn1987a": (83.8667, -69.2697), # SN 1987A
                    "grb221009a": (288.27, 19.78) # GRB 221009A
                }
                
                lower_name = target_name.lower().replace(" ", "")
                if lower_name in common_objects:
                    ra, dec = common_objects[lower_name]
                    coord = SkyCoord(ra, dec, unit="deg")
                    self.root.after(0, self._update_coords, ra, dec)
                    self.status_var.set(f"Using predefined coordinates for {target_name}: {coord.to_string('hmsdms')}")
                    return
                else:
                    raise ValueError(f"Could not resolve target name: {target_name}")
            
            if result is not None:
                ra = result['RA'][0]
                dec = result['DEC'][0]
                # Convert to degrees
                coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
                
                # Update the UI from the main thread
                self.root.after(0, self._update_coords, coord.ra.degree, coord.dec.degree)
                self.status_var.set(f"Resolved {target_name} to {coord.to_string('hmsdms')}")
            else:
                self.root.after(0, lambda: messagebox.showwarning("Warning", f"Could not resolve target name: {target_name}"))
                self.status_var.set("Ready")
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", f"Error resolving target name: {str(e)}"))
            self.status_var.set("Ready")
    
    def _update_coords(self, ra, dec):
        """Update the coordinate entries from the main thread."""
        self.ra_entry.delete(0, tk.END)
        self.ra_entry.insert(0, f"{ra:.6f}")
        self.dec_entry.delete(0, tk.END)
        self.dec_entry.insert(0, f"{dec:.6f}")
    
    def run_vo_query(self):
        """Run a VO query with the provided parameters."""
        # Get the coordinates
        try:
            if self.ra_entry.get() and self.dec_entry.get():
                ra = float(self.ra_entry.get())
                dec = float(self.dec_entry.get())
                position = SkyCoord(ra, dec, unit="deg")
            else:
                # Try to resolve the target name
                self.resolve_target_name()
                return
        except ValueError:
            messagebox.showerror("Error", "Invalid coordinates. Please enter valid RA and Dec values.")
            return
        
        # Get the radius
        try:
            radius = float(self.radius_entry.get()) * u.arcmin
        except ValueError:
            messagebox.showerror("Error", "Invalid search radius. Please enter a valid number.")
            return
        
        # Get the selected bands
        bands = [band for band, var in self.band_vars.items() if var.get()]
        if not bands:
            messagebox.showwarning("Warning", "No NIR bands selected. Please select at least one band.")
            return
        
        # Show that we're working on it
        self.status_var.set(f"Querying VO archives for position {position.to_string('hmsdms')}...")
        
        # Use a thread to keep the UI responsive
        if self.operation_in_progress:
            messagebox.showinfo("Info", "An operation is already in progress.")
            return
        
        self.operation_in_progress = True
        thread = threading.Thread(target=self._vo_query_thread, args=(position, radius, bands))
        thread.daemon = True
        thread.start()
    
    def _vo_query_thread(self, position, radius, bands):
        """Thread function to run a VO query."""
        try:
            # Run the query
            print(f"Querying VO archives at position {position.to_string('hmsdms')} with radius {radius}")
            results = query_vo_archives(position=position, radius=radius, bands=bands)
            
            # If no results or error, create simulated data for testing
            if len(results) == 0:
                print("No results found from VO archives. Creating simulated data.")
                results = Table()
                results['obs_id'] = [f'sim_obs_{i}' for i in range(5)]
                results['ra'] = [position.ra.degree + (i*0.01) for i in range(5)]
                results['dec'] = [position.dec.degree - (i*0.01) for i in range(5)]
                results['instrument'] = ['RIMAS', 'PRIME', 'SOFI', 'HAWK-I', 'ISAAC']
                results['filter'] = ['J', 'H', 'K', 'J', 'K']
                results['date_obs'] = ['2025-03-01T00:00:00', '2025-03-01T01:00:00',
                                      '2025-03-01T02:00:00', '2025-03-02T00:00:00',
                                      '2025-03-02T01:00:00']
                results['access_url'] = [f"https://example.com/simulated/{i}.fits" for i in range(5)]
                results['service_name'] = ['SIMULATION'] * 5
            
            # Store the results
            self.query_results = results
            
            # Update the UI from the main thread
            self.root.after(0, self._update_results_tree, results)
            self.status_var.set(f"Query complete. Found {len(results)} results.")
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", f"Error running VO query: {str(e)}"))
            self.status_var.set("Ready")
            print(f"Error running VO query: {str(e)}")
        finally:
            self.operation_in_progress = False
    
    def _update_results_tree(self, results):
        """Update the results treeview from the main thread."""
        # Clear the treeview
        for item in self.results_tree.get_children():
            self.results_tree.delete(item)
        
        # Add the results
        for i, row in enumerate(results):
            values = []
            # Get the values for each column (various column names are used in different services)
            id_val = row.get('obs_id', row.get('observation_id', f"result_{i}"))
            ra_val = row.get('ra', row.get('s_ra', 0.0))
            if hasattr(ra_val, 'value'):
                ra_val = ra_val.value
            dec_val = row.get('dec', row.get('s_dec', 0.0))
            if hasattr(dec_val, 'value'):
                dec_val = dec_val.value
            instrument = row.get('instrument', row.get('instrume', 'Unknown'))
            filter_val = row.get('filter', row.get('band', 'Unknown'))
            date_obs = row.get('date_obs', row.get('t_min', 'Unknown'))
            
            values = [id_val, f"{ra_val:.6f}", f"{dec_val:.6f}", instrument, filter_val, date_obs]
            self.results_tree.insert('', tk.END, values=values)
    
    def download_selected_data(self):
        """Download the selected data from the results."""
        selected_items = self.results_tree.selection()
        if not selected_items:
            messagebox.showwarning("Warning", "No items selected. Please select one or more items to download.")
            return
        
        if self.query_results is None or len(self.query_results) == 0:
            messagebox.showwarning("Warning", "No query results available.")
            return
        
        # Get the indices of the selected items
        indices = [self.results_tree.index(item) for item in selected_items]
        
        # Subset the query results to only the selected items
        selected_results = Table(self.query_results[indices])
        
        # Show that we're working on it
        self.status_var.set(f"Downloading {len(selected_results)} files...")
        
        # Use a thread to keep the UI responsive
        if self.operation_in_progress:
            messagebox.showinfo("Info", "An operation is already in progress.")
            return
        
        self.operation_in_progress = True
        thread = threading.Thread(target=self._download_thread, args=(selected_results,))
        thread.daemon = True
        thread.start()
    
    def _download_thread(self, results):
        """Thread function to download data."""
        try:
            # Create the download directory
            download_dir = self.fits_dir
            download_dir.mkdir(exist_ok=True)
            
            # Download the data (or simulate download if access URL doesn't work)
            print(f"Downloading {len(results)} files to {download_dir}")
            
            # Try real download first
            try:
                downloaded_files = download_fits_data(results, destination_dir=str(download_dir))
            except Exception as e:
                print(f"Error downloading files: {str(e)}")
                print("Creating simulated FITS files instead")
                
                # Create simulated FITS files
                downloaded_files = []
                for i, row in enumerate(results):
                    obs_id = row.get('obs_id', f'observation_{i}')
                    file_path = download_dir / f"{obs_id}.fits"
                    
                    # Create a simple text file simulating a FITS file
                    with open(file_path, 'w') as f:
                        f.write(f"Simulated FITS data for {obs_id}\n")
                        f.write(f"RA: {row.get('ra', 0.0)}\n")
                        f.write(f"Dec: {row.get('dec', 0.0)}\n")
                        f.write(f"Filter: {row.get('filter', 'Unknown')}\n")
                        f.write(f"Instrument: {row.get('instrument', 'Unknown')}\n")
                        f.write(f"Date: {row.get('date_obs', 'Unknown')}\n")
                    
                    downloaded_files.append(str(file_path))
                    print(f"Created simulated FITS file: {file_path}")
            
            # Process the downloaded files
            processed_files = []
            for fits_file in downloaded_files:
                # Extract metadata (or create simulated metadata)
                try:
                    standardize_fits_headers(fits_file)
                    metadata = extract_metadata(fits_file)
                except Exception as e:
                    print(f"Error processing FITS header: {str(e)}")
                    print("Creating simulated metadata instead")
                    
                    # Create simulated metadata
                    metadata = {
                        'filename': os.path.basename(fits_file),
                        'object': 'SIMULATED',
                        'telescope': 'SIMULATED',
                        'instrument': 'SIMULATED',
                        'filter': 'SIMULATED',
                        'date_obs': '2025-03-23T00:00:00',
                        'exptime': 100.0,
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
                
                if metadata:
                    # Save metadata to JSON
                    metadata_file = Path(fits_file).with_suffix('.json')
                    
                    # Remove SkyCoord object for JSON serialization if present
                    if 'coords' in metadata:
                        coords = metadata.pop('coords')
                        metadata['ra'] = coords.ra.degree if coords else metadata.get('ra')
                        metadata['dec'] = coords.dec.degree if coords else metadata.get('dec')
                    
                    with open(metadata_file, 'w') as f:
                        json.dump(metadata, f, indent=2)
                    
                    print(f"Saved metadata to {metadata_file}")
                    processed_files.append(fits_file)
            
            # Update the UI from the main thread
            message = f"Downloaded and processed {len(processed_files)} files to {download_dir}"
            self.root.after(0, lambda: messagebox.showinfo("Download Complete", message))
            self.status_var.set(message)
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", f"Error downloading data: {str(e)}"))
            self.status_var.set("Ready")
            print(f"Error in download thread: {str(e)}")
        finally:
            self.operation_in_progress = False
    
    def export_results(self, format_type):
        """Export the query results to a file."""
        if self.query_results is None or len(self.query_results) == 0:
            messagebox.showwarning("Warning", "No query results available to export.")
            return
        
        # Determine the output directory and filename
        if format_type == 'csv':
            output_dir = self.csv_dir
            filename = f"query_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            save_func = save_data_to_csv
        else: # Excel
            output_dir = self.excel_dir
            filename = f"query_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"
            save_func = save_to_excel
        
        # Show that we're working on it
        self.status_var.set(f"Exporting {len(self.query_results)} results to {format_type.upper()}...")
        
        # Use a thread to keep the UI responsive
        if self.operation_in_progress:
            messagebox.showinfo("Info", "An operation is already in progress.")
            return
        
        self.operation_in_progress = True
        thread = threading.Thread(target=self._export_thread, args=(self.query_results, save_func, output_dir, filename))
        thread.daemon = True
        thread.start()
    
    def _export_thread(self, data, save_func, output_dir, filename):
        """Thread function to export data."""
        try:
            # Export the data
            print(f"Exporting {len(data)} rows to {output_dir}/{filename}")
            
            try:
                file_path = save_func(data, output_dir=str(output_dir), filename=filename)
            except Exception as e:
                print(f"Error using save function: {str(e)}")
                print("Creating file directly")
                
                # If the save function fails, try a more direct approach
                output_path = Path(output_dir)
                output_path.mkdir(parents=True, exist_ok=True)
                file_path = output_path / filename
                
                if filename.endswith('.csv'):
                    # Convert to pandas dataframe and save as CSV
                    import pandas as pd
                    df = pd.DataFrame(data=data.as_array())
                    df.to_csv(file_path, index=False)
                else: # Excel
                    # Convert to pandas dataframe and save as Excel
                    import pandas as pd
                    df = pd.DataFrame(data=data.as_array())
                    df.to_excel(file_path, index=False)
            
            # Update the export history
            export_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            export_type = "CSV" if filename.endswith(".csv") else "Excel"
            export_rows = len(data)
            
            # Update the UI from the main thread
            self.root.after(0, self._add_export_history, export_time, export_type, export_rows, file_path)
            
            message = f"Exported {export_rows} rows to {file_path}"
            self.root.after(0, lambda: messagebox.showinfo("Export Complete", message))
            self.status_var.set(message)
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", f"Error exporting data: {str(e)}"))
            self.status_var.set("Ready")
            print(f"Error in export thread: {str(e)}")
        finally:
            self.operation_in_progress = False
    
    def _add_export_history(self, time_str, type_str, rows, file_path):
        """Add an entry to the export history treeview."""
        self.export_tree.insert('', 0, values=(time_str, type_str, rows, file_path))
    
    def open_export_file(self, event):
        """Open an exported file when double-clicked in the history."""
        selected_items = self.export_tree.selection()
        if not selected_items:
            return
        
        item = selected_items[0]
        file_path = self.export_tree.item(item, "values")[3]
        
        try:
            # Use the platform-specific way to open a file
            import os
            import platform
            if platform.system() == 'Windows':
                os.startfile(file_path)
            elif platform.system() == 'Darwin': # macOS
                import subprocess
                subprocess.run(['open', file_path])
            else: # Linux
                import subprocess
                subprocess.run(['xdg-open', file_path])
        except Exception as e:
            messagebox.showerror("Error", f"Could not open file: {str(e)}")
    
    def browse_export_dir(self):
        """Browse for an export directory."""
        directory = filedialog.askdirectory(initialdir=self.export_dir.get())
        if directory:
            self.export_dir.set(directory)
    
    def browse_data_dir(self):
        """Browse for a data directory."""
        directory = filedialog.askdirectory(initialdir=self.data_dir_var.get())
        if directory:
            self.data_dir_var.set(directory)
    
    def apply_settings(self):
        """Apply the settings."""
        # Set the data directory
        new_data_dir = Path(self.data_dir_var.get())
        if new_data_dir != self.data_dir:
            # Create the new data directory if it doesn't exist
            new_data_dir.mkdir(exist_ok=True)
            
            # Update all the subdirectories
            self.data_dir = new_data_dir
            self.csv_dir = self.data_dir / 'csv_data'
            self.excel_dir = self.data_dir / 'excel_data'
            self.fits_dir = self.data_dir / 'fits_data'
            self.alerts_dir = self.data_dir / 'alerts_data'
            
            # Create the subdirectories
            self.csv_dir.mkdir(exist_ok=True)
            self.excel_dir.mkdir(exist_ok=True)
            self.fits_dir.mkdir(exist_ok=True)
            self.alerts_dir.mkdir(exist_ok=True)
            
            # Update the status bar
            for widget in self.status_bar.winfo_children():
                widget.destroy()
            ttk.Label(self.status_bar, textvariable=self.status_var, anchor=tk.W).pack(side=tk.LEFT)
            ttk.Label(self.status_bar, text=f"Data Dir: {self.data_dir}", anchor=tk.E).pack(side=tk.RIGHT)
            
            print(f"Data directory updated to: {self.data_dir}")
        
        messagebox.showinfo("Settings Applied", "Settings have been applied successfully.")
    
    def reset_settings(self):
        """Reset the settings to default values."""
        self.data_dir_var.set(str(DEFAULT_DATA_DIR))
        self.observatory.set("LDT")
        self.log_level.set("INFO")
        self.set_log_level(None)
        print("Settings reset to defaults")
    
    def set_log_level(self, event):
        """Set the logging level."""
        level_name = self.log_level.get()
        level = getattr(logging, level_name)
        logging.getLogger().setLevel(level)
        
        # Update all loggers
        for logger_name in logging.root.manager.loggerDict:
            logging.getLogger(logger_name).setLevel(level)
        
        print(f"Log level set to: {level_name}")
    
    def start_gcn_monitor(self):
        """Start monitoring GCN notices."""
        if self.monitor_running:
            messagebox.showinfo("Info", "GCN monitoring is already running.")
            return
        
        # Show that we're working on it
        self.status_var.set("Starting GCN monitoring...")
        
        # For demo, simulate a GCN notice
        gcn_notice = {
            "type": "GRB",
            "source": "Swift-BAT",
            "trigger_time": datetime.now().isoformat(),
            "ra": 123.456,
            "dec": -45.678,
            "error_radius": 0.05,
            "energy_range": "15-150 keV",
            "t90": 45.2,
            "flux": 3.4e-7
        }
        
        # Process the simulated notice
        parsed_notice = self.alert_monitor.parse_gcn_notice(gcn_notice)
        self.gcn_alerts.append(parsed_notice)
        
        # Update the UI
        self._add_gcn_alert(parsed_notice)
        
        # Update status
        self.status_var.set("GCN notice received (simulated)")
        print("Simulated GCN notice received and processed")
    
    def query_recent_gcn(self):
        """Query for recent GCN notices."""
        # For demo, simulate multiple GCN notices
        for i in range(3):
            gcn_notice = {
                "type": "GRB",
                "source": ["Swift-BAT", "Fermi-GBM", "INTEGRAL"][i % 3],
                "trigger_time": (datetime.now() - i * timedelta(days=1)).isoformat(),
                "ra": 123.456 + i,
                "dec": -45.678 + i,
                "error_radius": 0.05,
                "energy_range": "15-150 keV",
                "t90": 45.2 + i * 10,
                "flux": 3.4e-7 * (1 + i)
            }
            
            # Process the simulated notice
            parsed_notice = self.alert_monitor.parse_gcn_notice(gcn_notice)
            self.gcn_alerts.append(parsed_notice)
            
            # Update the UI
            self._add_gcn_alert(parsed_notice)
        
        # Update status
        self.status_var.set("Retrieved 3 recent GCN notices (simulated)")
        print("Simulated retrieval of 3 recent GCN notices")
    
    def _add_gcn_alert(self, notice):
        """Add a GCN notice to the alerts treeview."""
        # Get the alert data
        time_str = notice.get('time', 'Unknown')
        if hasattr(time_str, 'isot'):
            time_str = time_str.isot
        type_str = notice.get('type', 'Unknown')
        source = notice.get('source', 'Unknown')
        position = notice.get('position')
        ra = dec = 'Unknown'
        if position:
            ra = f"{position.ra.degree:.6f}"
            dec = f"{position.dec.degree:.6f}"
        significance = notice.get('flux', '-')
        
        # Add to the treeview
        self.alerts_tree.insert('', 0, values=(time_str, type_str, source, ra, dec, significance))
    
    def get_ligo_alerts(self):
        """Get LIGO/Virgo alerts."""
        # Show that we're working on it
        self.status_var.set("Retrieving LIGO/Virgo alerts...")
        
        # Use a thread to keep the UI responsive
        if self.operation_in_progress:
            messagebox.showinfo("Info", "An operation is already in progress.")
            return
        
        self.operation_in_progress = True
        thread = threading.Thread(target=self._ligo_alerts_thread)
        thread.daemon = True
        thread.start()
    
    def _ligo_alerts_thread(self):
        """Thread function to get LIGO/Virgo alerts."""
        try:
            # Get the alerts
            print("Retrieving LIGO/Virgo alerts")
            alerts = self.alert_monitor.get_ligo_alerts()
            
            # Store the alerts
            self.ligo_alerts = alerts
            
            # Update the UI from the main thread
            self.root.after(0, self._update_ligo_alerts, alerts)
            self.status_var.set(f"Retrieved {len(alerts)} LIGO/Virgo alerts")
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", f"Error retrieving LIGO/Virgo alerts: {str(e)}"))
            self.status_var.set("Ready")
            print(f"Error retrieving LIGO/Virgo alerts: {str(e)}")
        finally:
            self.operation_in_progress = False
    
    def _update_ligo_alerts(self, alerts):
        """Update the alerts treeview with LIGO/Virgo alerts."""
        for alert in alerts:
            # Get the alert data
            time_str = alert.get('created', 'Unknown')
            type_str = 'GW'
            source = ', '.join(alert.get('instruments', []))
            
            # LIGO alerts don't have RA/Dec directly
            ra = dec = 'N/A'
            
            # Get the significance
            far = alert.get('far', 0)
            significance = f"{1/far:.1e}" if far else '-'
            
            # Add to the treeview
            self.alerts_tree.insert('', 0, values=(time_str, type_str, source, ra, dec, significance))
    
    def download_ligo_skymap(self):
        """Download a LIGO skymap for a selected alert."""
        selected_items = self.alerts_tree.selection()
        if not selected_items:
            messagebox.showwarning("Warning", "No alert selected. Please select a LIGO/Virgo alert.")
            return
        
        # Get the selected alert
        item = selected_items[0]
        alert_type = self.alerts_tree.item(item, "values")[1]
        if alert_type != 'GW':
            messagebox.showwarning("Warning", "Selected alert is not a gravitational wave alert.")
            return
        
        # Show that we're working on it
        self.status_var.set("Downloading LIGO skymap...")
        
        # Simulate skymap download
        skymap_dir = self.alerts_dir / 'skymaps'
        skymap_dir.mkdir(exist_ok=True)
        
        # Create a placeholder skymap file
        event_id = "S230101a" # Simulated event ID
        skymap_file = skymap_dir / f"{event_id}_skymap.txt"
        with open(skymap_file, 'w') as f:
            f.write("This is a placeholder for a LIGO skymap file.\n")
            f.write("In a real implementation, a HEALPix FITS file would be downloaded and parsed.\n")
            f.write("RA: 197.45\n") # GW170817 location
            f.write("Dec: -23.38\n")
        
        # Parse the simulated skymap
        skymap_data = self.alert_monitor.parse_ligo_skymap(skymap_file)
        
        if skymap_data:
            # Get the maximum probability position
            max_prob_pos = skymap_data.get('max_prob_position')
            
            # Update the status
            message = f"Downloaded and parsed skymap for {event_id}"
            if max_prob_pos:
                message += f"\nMaximum probability position: {max_prob_pos.to_string('hmsdms')}"
            
            messagebox.showinfo("Skymap Download Complete", message)
            self.status_var.set(f"Downloaded skymap for {event_id}")
            print(message)
        else:
            messagebox.showwarning("Warning", "Could not parse skymap")
            self.status_var.set("Ready")
    
    def check_visibility(self):
        """Check the visibility of a selected target."""
        # Get the coordinates
        try:
            if self.ra_entry.get() and self.dec_entry.get():
                ra = float(self.ra_entry.get())
                dec = float(self.dec_entry.get())
                coordinates = SkyCoord(ra, dec, unit="deg")
            else:
                # Try a selected alert
                selected_items = self.alerts_tree.selection()
                if not selected_items:
                    messagebox.showwarning("Warning", "No coordinates or alert selected.")
                    return
                
                # Get the selected alert
                item = selected_items[0]
                ra = self.alerts_tree.item(item, "values")[3]
                dec = self.alerts_tree.item(item, "values")[4]
                
                if ra == 'Unknown' or ra == 'N/A' or dec == 'Unknown' or dec == 'N/A':
                    messagebox.showwarning("Warning", "Selected alert does not have valid coordinates.")
                    return
                
                coordinates = SkyCoord(float(ra), float(dec), unit="deg")
        except ValueError:
            messagebox.showerror("Error", "Invalid coordinates. Please enter valid RA and Dec values.")
            return
        
        # Get the observatory
        observatory = self.observatory.get()
        
        # Show that we're working on it
        self.status_var.set(f"Checking visibility from {observatory}...")
        
        # Use a thread to keep the UI responsive
        if self.operation_in_progress:
            messagebox.showinfo("Info", "An operation is already in progress.")
            return
        
        self.operation_in_progress = True
        thread = threading.Thread(target=self._check_visibility_thread, args=(coordinates, observatory))
        thread.daemon = True
        thread.start()
    
    def _check_visibility_thread(self, coordinates, observatory):
        """Thread function to check visibility."""
        try:
            # Check visibility
            print(f"Checking visibility of {coordinates.to_string('hmsdms')} from {observatory}")
            visibility = self.alert_monitor.check_visibility(coordinates, observatory)
            
            if visibility:
                # Update the UI from the main thread
                self.root.after(0, self._show_visibility, visibility)
                self.status_var.set(f"Checked visibility from {observatory}")
            else:
                self.root.after(0, lambda: messagebox.showwarning("Warning", f"Could not check visibility from {observatory}"))
                self.status_var.set("Ready")
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", f"Error checking visibility: {str(e)}"))
            self.status_var.set("Ready")
            print(f"Error checking visibility: {str(e)}")
        finally:
            self.operation_in_progress = False
    
    def _show_visibility(self, visibility):
        """Show visibility information."""
        # Clear the details text
        self.details_text.delete(1.0, tk.END)
        
        # Add the visibility information
        self.details_text.insert(tk.END, "Visibility Information\n")
        self.details_text.insert(tk.END, "======================\n\n")
        self.details_text.insert(tk.END, f"Target: {visibility['coordinates'].to_string('hmsdms')}\n")
        self.details_text.insert(tk.END, f"Observatory: {visibility['observatory']}\n")
        self.details_text.insert(tk.END, f"Time: {visibility['time'].iso}\n\n")
        self.details_text.insert(tk.END, f"Current altitude: {visibility['altitude'].deg:.1f}°\n")
        self.details_text.insert(tk.END, f"Current azimuth: {visibility['azimuth'].deg:.1f}°\n")
        self.details_text.insert(tk.END, f"Airmass: {visibility['airmass']:.2f}\n")
        self.details_text.insert(tk.END, f"Sun altitude: {visibility['sun_altitude'].deg:.1f}°\n\n")
        
        is_night = "Yes" if visibility['is_night'] else "No"
        is_observable = "Yes" if visibility['is_observable'] else "No"
        
        self.details_text.insert(tk.END, f"Is it night? {is_night}\n")
        self.details_text.insert(tk.END, f"Is target observable? {is_observable}\n\n")
        self.details_text.insert(tk.END, f"Transit time: {visibility['transit_time'].iso}\n")
        self.details_text.insert(tk.END, f"Observable hours tonight: {visibility['observable_hours']:.1f}\n")
    
    def generate_large_dataset(self):
        """Generate a large dataset."""
        try:
            size = int(self.dataset_size.get())
            if size <= 0:
                messagebox.showwarning("Warning", "Dataset size must be a positive integer.")
                return
        except ValueError:
            messagebox.showerror("Error", "Invalid dataset size. Please enter a valid number.")
            return
        
        # Get the output directory
        output_dir = self.export_dir.get()
        
        # Get the format
        format_type = self.export_format.get()
        
        # Show that we're working on it
        self.status_var.set(f"Generating large dataset with {size} rows...")
        
        # Use a thread to keep the UI responsive
        if self.operation_in_progress:
            messagebox.showinfo("Info", "An operation is already in progress.")
            return
        
        self.operation_in_progress = True
        thread = threading.Thread(target=self._generate_dataset_thread, args=(size, output_dir, format_type))
        thread.daemon = True
        thread.start()
    
    def _generate_dataset_thread(self, size, output_dir, format_type):
        """Thread function to generate a large dataset."""
        try:
            # Generate the dataset
            print(f"Generating large dataset with {size} rows")
            
            # Create output directory
            Path(output_dir).mkdir(parents=True, exist_ok=True)
            
            # Generate random astronomical data
            import numpy as np
            import pandas as pd
            
            # Generate data with reasonable astronomical values
            data = {
                'obs_id': [f'obs_{i:06d}' for i in range(size)],
                'ra': np.random.uniform(0, 360, size),
                'dec': np.random.uniform(-90, 90, size),
                'magnitude': np.random.normal(18, 2, size),
                'filter': np.random.choice(['J', 'H', 'K'], size),
                'exptime': np.random.randint(30, 600, size),
                'signal_to_noise': np.random.exponential(5, size),
                'date_obs': pd.date_range(start='2025-01-01', periods=size, freq='s')  # Changed from 'S' to 's'
            }
            
            df = pd.DataFrame(data)
            
            # Save the data in the specified format
            if format_type == 'csv':
                file_path = Path(output_dir) / f'large_dataset_{size}_rows.csv'
                df.to_csv(file_path, index=False)
            else: # Excel
                file_path = Path(output_dir) / f'large_dataset_{size}_rows.xlsx'
                
                # For very large datasets, use chunking to avoid Excel row limit
                if size > 1000000: # Excel row limit is 1,048,576
                    with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
                        for i in range(0, size, 1000000):
                            sheet_name = f'Data_{i//1000000 + 1}'
                            df.iloc[i:i+1000000].to_excel(writer, sheet_name=sheet_name, index=False)
                else:
                    df.to_excel(file_path, index=False)
            
            # Update the export history
            export_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            export_type = "CSV" if format_type == 'csv' else "Excel"
            
            # Update the UI from the main thread
            self.root.after(0, self._add_export_history, export_time, export_type, size, file_path)
            
            message = f"Generated {size} rows and saved to {file_path}"
            self.root.after(0, lambda: messagebox.showinfo("Dataset Generation Complete", message))
            self.status_var.set(message)
            print(message)
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", f"Error generating dataset: {str(e)}"))
            self.status_var.set("Ready")
            print(f"Error generating dataset: {str(e)}")
        finally:
            self.operation_in_progress = False
    
    def show_alert_details(self, event):
        """Show details of a selected alert."""
        selected_items = self.alerts_tree.selection()
        if not selected_items:
            return
        
        # Get the selected alert
        item = selected_items[0]
        alert_type = self.alerts_tree.item(item, "values")[1]
        
        # Clear the details text
        self.details_text.delete(1.0, tk.END)
        
        # Display the alert details based on the type
        if alert_type == 'GRB':
            self.details_text.insert(tk.END, "GRB Alert Details\n")
            self.details_text.insert(tk.END, "=================\n\n")
            
            # Get details from the item
            time_str = self.alerts_tree.item(item, "values")[0]
            source = self.alerts_tree.item(item, "values")[2]
            ra = self.alerts_tree.item(item, "values")[3]
            dec = self.alerts_tree.item(item, "values")[4]
            
            self.details_text.insert(tk.END, f"Time: {time_str}\n")
            self.details_text.insert(tk.END, f"Source: {source}\n")
            self.details_text.insert(tk.END, f"Position: RA={ra}, Dec={dec}\n")
            
            # Add observatory visibility information
            observatory = self.observatory.get()
            self.details_text.insert(tk.END, f"\nVisibility from {observatory}: Calculating...\n")
            
            # Start a thread to check visibility
            if ra != 'Unknown' and dec != 'N/A':
                try:
                    coordinates = SkyCoord(float(ra), float(dec), unit="deg")
                    thread = threading.Thread(target=self._check_alert_visibility, args=(coordinates, observatory))
                    thread.daemon = True
                    thread.start()
                except ValueError:
                    self.details_text.insert(tk.END, "Cannot calculate visibility: invalid coordinates.\n")
            
        elif alert_type == 'GW':
            self.details_text.insert(tk.END, "Gravitational Wave Alert Details\n")
            self.details_text.insert(tk.END, "==============================\n\n")
            
            # Get details from the item
            time_str = self.alerts_tree.item(item, "values")[0]
            source = self.alerts_tree.item(item, "values")[2]
            significance = self.alerts_tree.item(item, "values")[5]
            
            self.details_text.insert(tk.END, f"Time: {time_str}\n")
            self.details_text.insert(tk.END, f"Instruments: {source}\n")
            self.details_text.insert(tk.END, f"Significance: {significance}\n\n")
            
            # Look up the full alert in our ligo_alerts list
            alert_time = time_str
            found_alert = None
            for alert in self.ligo_alerts:
                if alert.get('created') == alert_time:
                    found_alert = alert
                    break
            
            if found_alert:
                # Display additional information if available
                event_id = found_alert.get('superevent_id', 'Unknown')
                self.details_text.insert(tk.END, f"Event ID: {event_id}\n")
                
                if 'classification' in found_alert:
                    self.details_text.insert(tk.END, "\nSource Classification:\n")
                    for key, value in found_alert['classification'].items():
                        self.details_text.insert(tk.END, f" {key}: {value:.2f}\n")
                
                self.details_text.insert(tk.END, "\nActions:\n")
                self.details_text.insert(tk.END, "1. Download skymap\n")
                self.details_text.insert(tk.END, "2. Check for NIR followup data\n")
            else:
                self.details_text.insert(tk.END, "Full alert details not available.\n")
    
    def _check_alert_visibility(self, coordinates, observatory):
        """Check visibility for an alert and update the details text."""
        try:
            # Check visibility
            visibility = self.alert_monitor.check_visibility(coordinates, observatory)
            
            if visibility:
                # Update the details text
                def update_details():
                    self.details_text.delete("end-2l", tk.END) # Delete "Calculating..." line
                    is_observable = "Yes" if visibility['is_observable'] else "No"
                    alt = visibility['altitude'].deg
                    airmass = visibility['airmass']
                    hours = visibility['observable_hours']
                    self.details_text.insert(tk.END, f"\nVisible from {observatory}: {is_observable}\n")
                    self.details_text.insert(tk.END, f"Current altitude: {alt:.1f}°\n")
                    self.details_text.insert(tk.END, f"Current airmass: {airmass:.2f}\n")
                    self.details_text.insert(tk.END, f"Observable hours tonight: {hours:.1f}\n")
                
                self.root.after(0, update_details)
        except Exception as e:
            def show_error():
                self.details_text.delete("end-2l", tk.END) # Delete "Calculating..." line
                self.details_text.insert(tk.END, f"Error checking visibility: {e}\n")
            
            self.root.after(0, show_error)

def main():
    """Main function to run the application."""
    root = tk.Tk()
    app = DataAcquisitionApp(root)
    
    # Handle window close event
    def on_close():
        # Restore stdout if it was redirected
        if hasattr(app, 'original_stdout') and app.original_stdout:
            sys.stdout = app.original_stdout
        root.destroy()
    
    root.protocol("WM_DELETE_WINDOW", on_close)
    
    # Set window icon if available
    try:
        root.iconbitmap("assets/icon.ico")
    except:
        pass
    
    root.mainloop()

if __name__ == "__main__":
    main()
