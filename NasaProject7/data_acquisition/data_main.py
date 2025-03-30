import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import threading
import logging
import json
import os
import sys
import random
import time
from astropy.time import Time
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

# Import project modules
from data_acquisition.alerts import MultiMessengerAlertHandler
from data_acquisition.query import HybridDataAcquisition as DataAcquisition

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

class NIRCounterpartPipeline:
    """Main application for NIR GRB/GW Counterpart Pipeline"""
    
    def __init__(self, root):
        """Initialize the application"""
        self.root = root
        self.root.title("NIR GRB/GW Counterpart Pipeline")
        self.root.geometry("900x700")
        
        # Initialize status variables
        self.status_var = tk.StringVar()
        self.status_var.set("Ready")
        
        # Initialize data
        self.alert_list = []
        self.current_results = []
        self.alert_handler = None
        self.data_acquisition = DataAcquisition()
        
        # Create interface
        self._create_interface()
        
        # Create status bar
        status_frame = ttk.Frame(self.root)
        status_frame.pack(side=tk.BOTTOM, fill=tk.X)
        status_label = ttk.Label(status_frame, textvariable=self.status_var)
        status_label.pack(side=tk.LEFT, padx=5, pady=2)
        
        # Initialize alert handler
        self._initialize_alert_handler()
    
    def _create_interface(self):
        """Build the Tkinter GUI interface"""
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True)
        
        # Create tabs
        self._create_query_tab()
        self._create_alerts_tab()
        self._create_analysis_tab()  # Keep this as a simple visualizer
        self._create_settings_tab()
    
    def _create_query_tab(self):
        """Create the query tab for data acquisition"""
        query_frame = ttk.Frame(self.notebook)
        self.notebook.add(query_frame, text="Query")
        
        # Left panel for query parameters
        left_frame = ttk.LabelFrame(query_frame, text="Query Parameters")
        left_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=False, padx=5, pady=5)
        
        # Object ID/Coordinates
        ttk.Label(left_frame, text="Object ID or Coordinates:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.coord_var = tk.StringVar(value="GRB121128A")
        ttk.Entry(left_frame, textvariable=self.coord_var, width=30).grid(row=0, column=1, padx=5, pady=5)
        
        # Radius
        ttk.Label(left_frame, text="Search Radius (deg):").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.radius_var = tk.StringVar(value="0.2")
        ttk.Entry(left_frame, textvariable=self.radius_var, width=10).grid(row=1, column=1, sticky=tk.W, padx=5, pady=5)
        
        # Instruments frame
        inst_frame = ttk.LabelFrame(left_frame, text="Instruments")
        inst_frame.grid(row=2, column=0, columnspan=2, sticky=tk.W+tk.E, padx=5, pady=5)
        
        # Instrument checkboxes
        self.instruments = {}
        instruments_list = ["UKIRT", "VISTA", "WFCAM", "VIRCAM", "HAWK-I", "SOFI"]
        for i, inst in enumerate(instruments_list):
            self.instruments[inst] = tk.BooleanVar(value=True)
            ttk.Checkbutton(inst_frame, text=inst, variable=self.instruments[inst]).grid(
                row=i//2, column=i%2, sticky=tk.W, padx=5, pady=2
            )
        
        # Bands frame
        bands_frame = ttk.LabelFrame(left_frame, text="NIR Bands")
        bands_frame.grid(row=3, column=0, columnspan=2, sticky=tk.W+tk.E, padx=5, pady=5)
        
        # Band checkboxes
        self.bands = {}
        bands_list = ["J", "H", "K", "Y", "z", "other"]
        for i, band in enumerate(bands_list):
            self.bands[band] = tk.BooleanVar(value=True)
            ttk.Checkbutton(bands_frame, text=band, variable=self.bands[band]).grid(
                row=i//3, column=i%3, sticky=tk.W, padx=5, pady=2
            )
        
        # Query buttons
        button_frame = ttk.Frame(left_frame)
        button_frame.grid(row=4, column=0, columnspan=2, pady=10)
        
        ttk.Button(button_frame, text="Execute Query", command=self._execute_query).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear Results", command=self._clear_results).pack(side=tk.LEFT, padx=5)
        
        # Right panel for results
        right_frame = ttk.LabelFrame(query_frame, text="Query Results")
        right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Results text area
        self.results_text = scrolledtext.ScrolledText(right_frame, wrap=tk.WORD, width=60, height=30)
        self.results_text.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Download button
        ttk.Button(right_frame, text="Download Selected", command=self._download_selected).pack(pady=5)
    
    def _create_alerts_tab(self):
        """Create the alerts tab for monitoring GRB and GW events"""
        alerts_frame = ttk.Frame(self.notebook)
        self.notebook.add(alerts_frame, text="Alerts")
        
        # Control panel
        control_frame = ttk.LabelFrame(alerts_frame, text="Alert Controls")
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        
        # Control buttons
        self.start_button = ttk.Button(control_frame, text="Start Monitoring", command=self._start_alert_monitor)
        self.start_button.pack(side=tk.LEFT, padx=5, pady=5)
        
        self.pause_button = ttk.Button(control_frame, text="Pause", command=self._pause_alerts)
        self.pause_button.pack(side=tk.LEFT, padx=5, pady=5)
        
        self.resume_button = ttk.Button(control_frame, text="Resume", command=self._resume_alerts)
        self.resume_button.pack(side=tk.LEFT, padx=5, pady=5)
        
        # Simulation controls
        sim_frame = ttk.LabelFrame(control_frame, text="Simulation")
        sim_frame.pack(side=tk.RIGHT, padx=5, pady=5)
        
        ttk.Button(sim_frame, text="Generate Test Alerts", command=self._generate_test_alerts).pack(side=tk.LEFT, padx=5, pady=5)
        
        # Alerts display
        alerts_display = ttk.LabelFrame(alerts_frame, text="Recent Alerts")
        alerts_display.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Create treeview for alerts
        self.alerts_tree = ttk.Treeview(alerts_display, columns=("time", "type", "id", "ra", "dec"), show="headings")
        self.alerts_tree.heading("time", text="Time")
        self.alerts_tree.heading("type", text="Type")
        self.alerts_tree.heading("id", text="Alert ID")
        self.alerts_tree.heading("ra", text="RA")
        self.alerts_tree.heading("dec", text="Dec")
        
        self.alerts_tree.column("time", width=150)
        self.alerts_tree.column("type", width=50)
        self.alerts_tree.column("id", width=150)
        self.alerts_tree.column("ra", width=80)
        self.alerts_tree.column("dec", width=80)
        
        self.alerts_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Add scrollbar
        scrollbar = ttk.Scrollbar(alerts_display, orient=tk.VERTICAL, command=self.alerts_tree.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.alerts_tree.configure(yscrollcommand=scrollbar.set)
        
        # Alert details
        details_frame = ttk.LabelFrame(alerts_frame, text="Alert Details")
        details_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.alert_details = scrolledtext.ScrolledText(details_frame, wrap=tk.WORD, height=10)
        self.alert_details.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Bind selection event
        self.alerts_tree.bind("<<TreeviewSelect>>", self._show_alert_details)
    
    def _create_analysis_tab(self):
        """Create a simple analysis tab for visualization"""
        analysis_frame = ttk.Frame(self.notebook)
        self.notebook.add(analysis_frame, text="Analysis")
        
        # Control panel
        control_frame = ttk.LabelFrame(analysis_frame, text="Visualization Controls")
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=5, pady=5)
        
        # Plot type selector
        ttk.Label(control_frame, text="Plot Type:").pack(side=tk.LEFT, padx=5, pady=5)
        self.plot_type = tk.StringVar(value="Light Curve")
        plot_combo = ttk.Combobox(control_frame, textvariable=self.plot_type, 
                                 values=["Light Curve", "Finding Chart", "Color-Magnitude"])
        plot_combo.pack(side=tk.LEFT, padx=5, pady=5)
        plot_combo.bind("<<ComboboxSelected>>", self._update_plot)
        
        # Plot button
        ttk.Button(control_frame, text="Generate Plot", command=self._generate_plot).pack(side=tk.LEFT, padx=5, pady=5)
        
        # Canvas for matplotlib
        self.figure = plt.Figure(figsize=(6, 4), dpi=100)
        self.plot_canvas = FigureCanvasTkAgg(self.figure, analysis_frame)
        self.plot_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
    
    def _create_settings_tab(self):
        """Create the settings tab"""
        settings_frame = ttk.Frame(self.notebook)
        self.notebook.add(settings_frame, text="Settings")
        
        # General settings
        general_frame = ttk.LabelFrame(settings_frame, text="General Settings")
        general_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Data directory
        ttk.Label(general_frame, text="Data Directory:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.data_dir_var = tk.StringVar(value="./data")
        ttk.Entry(general_frame, textvariable=self.data_dir_var, width=40).grid(row=0, column=1, padx=5, pady=5)
        ttk.Button(general_frame, text="Browse...", command=self._browse_data_dir).grid(row=0, column=2, padx=5, pady=5)
        
        # Log level
        ttk.Label(general_frame, text="Log Level:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.log_level = tk.StringVar(value="INFO")
        ttk.Combobox(general_frame, textvariable=self.log_level, values=["DEBUG", "INFO", "WARNING", "ERROR"]).grid(
            row=1, column=1, sticky=tk.W, padx=5, pady=5
        )
        
        # Alert settings
        alert_frame = ttk.LabelFrame(settings_frame, text="Alert Settings")
        alert_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Alert sources
        ttk.Label(alert_frame, text="GRB Alert Source:").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.grb_source = tk.StringVar(value="GCN")
        ttk.Combobox(alert_frame, textvariable=self.grb_source, values=["GCN", "SWIFT", "Fermi", "Simulation"]).grid(
            row=0, column=1, sticky=tk.W, padx=5, pady=5
        )
        
        ttk.Label(alert_frame, text="GW Alert Source:").grid(row=1, column=0, sticky=tk.W, padx=5, pady=5)
        self.gw_source = tk.StringVar(value="LIGO/Virgo")
        ttk.Combobox(alert_frame, textvariable=self.gw_source, values=["LIGO/Virgo", "Simulation"]).grid(
            row=1, column=1, sticky=tk.W, padx=5, pady=5
        )
        
        # Save and load buttons
        button_frame = ttk.Frame(settings_frame)
        button_frame.pack(pady=10)
        
        ttk.Button(button_frame, text="Save Settings", command=self._save_settings).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Load Settings", command=self._load_settings).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Reset to Defaults", command=self._reset_settings).pack(side=tk.LEFT, padx=5)
    
    # --- Functionality methods ---
    
    def _initialize_alert_handler(self):
        """Initialize the alert handler"""
        try:
            self.alert_handler = MultiMessengerAlertHandler()
            logger.info("Alert handler initialized")
        except Exception as e:
            self.status_var.set(f"Error initializing alert handler: {str(e)}")
    
    def _execute_query(self):
        """Execute a data query based on form values"""
        try:
            # Get query parameters
            coord = self.coord_var.get().strip()
            radius = float(self.radius_var.get())
            
            # Get selected instruments and bands
            selected_instruments = [inst for inst, var in self.instruments.items() if var.get()]
            selected_bands = [band for band, var in self.bands.items() if var.get()]
            
            if not coord:
                messagebox.showerror("Input Error", "Please enter an object ID or coordinates")
                return
            
            # Clear previous results
            self.results_text.delete(1.0, tk.END)
            self.status_var.set(f"Querying data for {coord}")
            
            # Run query in background thread
            threading.Thread(
                target=self._execute_query_thread,
                args=(coord, radius, selected_instruments, selected_bands),
                daemon=True
            ).start()
            
        except ValueError as e:
            messagebox.showerror("Input Error", f"Invalid input: {str(e)}")
        except Exception as e:
            messagebox.showerror("Query Error", f"Error executing query: {str(e)}")
    
    def _execute_query_thread(self, coord, radius, instruments, bands):
        """Background thread for query execution"""
        try:
            # Check network connectivity
            results = self.data_acquisition.query_archives(coord, radius, instruments, bands)
            
            # Store results for later use
            self.current_results = results
            
            # Update UI in main thread
            def _update_ui():
                self._display_results(results)
                self.status_var.set(f"Found {len(results)} observations")
            
            self.root.after(0, _update_ui)
            
        except Exception as e:
            def _show_error(error_msg):
                self.status_var.set(f"Query error: {error_msg}")
                messagebox.showerror("Query Error", f"Error: {error_msg}")
            error_msg = str(e)
            self.root.after(0, lambda msg=error_msg: _show_error(msg))
            error_msg = str(e)
            self.root.after(0, lambda msg=error_msg: _show_error(msg))
    
    def _display_results(self, results):
        """Display query results in the text area"""
        if not results:
            self.results_text.insert(tk.END, "No results found.\n")
            return
        
        for i, item in enumerate(results):
            self.results_text.insert(tk.END, f"Result #{i+1}:\n")
            self.results_text.insert(tk.END, f"  Instrument: {item.get('instrument', 'Unknown')}\n")
            self.results_text.insert(tk.END, f"  Band: {item.get('band', 'Unknown')}\n")
            self.results_text.insert(tk.END, f"  Date: {item.get('date', 'Unknown')}\n")
            self.results_text.insert(tk.END, f"  Exposure: {item.get('exposure', 'Unknown')} seconds\n")
            self.results_text.insert(tk.END, f"  URL: {item.get('url', 'N/A')}\n")
            self.results_text.insert(tk.END, "-" * 50 + "\n")
    
    def _clear_results(self):
        """Clear the results display"""
        self.results_text.delete(1.0, tk.END)
        self.current_results = []
        self.status_var.set("Results cleared")
    
    def _download_selected(self):
        """Download selected data products"""
        if not self.current_results:
            messagebox.showinfo("No Data", "No query results to download")
            return
        
        # In a real implementation, you would:
        # 1. Show a selection dialog for which results to download
        # 2. Download the selected files in a background thread
        # 3. Update the UI with progress
        
        # For now, just show a message
        self.status_var.set("Downloading data...")
        
        # Simulate download in a thread
        threading.Thread(target=self._simulate_download, daemon=True).start()
    
    def _simulate_download(self):
        """Simulate downloading data"""
        # Simulate a download process
        time.sleep(2)
        
        # Update UI in main thread
        self.root.after(0, lambda: self.status_var.set("Download complete"))
    
    def _start_alert_monitor(self):
        """Start alert monitoring with proper thread handling"""
        try:
            # Use proper threading instead of raw asyncio
            self.alert_thread = threading.Thread(
                target=self._run_alert_monitor_thread,
                daemon=True
            )
            self.alert_thread.start()
            self.status_var.set("Alert monitoring started")
        except Exception as e:
            def _set_error_status(error_msg):
                self.status_var.set(f"Error starting alert monitor: {error_msg}")
            error_msg = str(e)
            self.root.after(0, lambda msg=error_msg: _set_error_status(msg))
    
    def _run_alert_monitor_thread(self):
        """Run alert monitor in a separate thread"""
        try:
            # Start the alert handler
            if hasattr(self.alert_handler, 'start'):
                self.alert_handler.start()
            else:
                # Fallback to manual alert generation
                while True:
                    if random.random() < 0.2:  # 20% chance to generate an alert
                        alert_type = random.choice(["GRB", "GW"])
                        alert_data = {
                            'time': Time.now().iso,
                            'type': alert_type,
                            'id': f"{alert_type}_{int(Time.now().mjd)}",
                            'ra': random.uniform(0, 360),
                            'dec': random.uniform(-90, 90),
                            'confidence': random.uniform(0, 1)
                        }
                        # Update UI on main thread
                        self.root.after(0, lambda a=alert_data: self._add_alert_to_ui(a))
                    time.sleep(3)
                
        except Exception as e:
            def _set_error_status(error_msg):
                self.status_var.set(f"Alert monitor error: {error_msg}")
            error_msg = str(e)
            self.root.after(0, lambda msg=error_msg: _set_error_status(msg))
    
    def _pause_alerts(self):
        """Pause alert monitoring"""
        try:
            if hasattr(self.alert_handler, 'pause'):
                self.alert_handler.pause()
                self.status_var.set("Alert monitoring paused")
            else:
                # Fallback if pause method is missing
                self.status_var.set("Pausing alert monitoring")
        except Exception as e:
            def _set_error_status(error_msg):
                self.status_var.set(f"Error pausing alerts: {error_msg}")
            error_msg = str(e)
            self.root.after(0, lambda msg=error_msg: _set_error_status(msg))
    
    def _resume_alerts(self):
        """Resume alert monitoring"""
        try:
            if hasattr(self.alert_handler, 'resume'):
                self.alert_handler.resume()
                self.status_var.set("Alert monitoring resumed")
            else:
                # Fallback if resume method is missing
                self.status_var.set("Resuming alert monitoring")
        except Exception as e:
            def _set_error_status(error_msg):
                self.status_var.set(f"Error resuming alerts: {error_msg}")
            error_msg = str(e)
            self.root.after(0, lambda msg=error_msg: _set_error_status(msg))
    
    def _generate_test_alerts(self):
        """Generate test alerts for development"""
        try:
            self.status_var.set("Generating test alerts...")
            
            # Set up parameters
            num_gw = 5
            num_grb = 10
            num_coincidences = 2
            
            # Start a thread to generate the alerts
            threading.Thread(
                target=self._generate_alerts_thread,
                args=(num_gw, num_grb, num_coincidences),
                daemon=True
            ).start()
            
        except Exception as e:
            def _set_error_status(error_msg):
                self.status_var.set(f"Simulation failed: {error_msg}")
            error_msg = str(e)
            self.root.after(0, lambda msg=error_msg: _set_error_status(msg))
    
    def _generate_alerts_thread(self, num_gw, num_grb, num_coincidences):
        """Generate simulated alerts in background thread"""
        try:
            # Log the simulation start
            logger.info(f"Generating {num_gw} GW alerts, {num_grb} GRB alerts with {num_coincidences} coincidences")
            
            # Try using the alert handler's simulate method if available
            if hasattr(self.alert_handler, 'simulate_alerts'):
                self.alert_handler.simulate_alerts(num_gw, num_grb, num_coincidences)
                return
                
            # Fallback to generating our own alerts
            for i in range(num_gw + num_grb):
                alert_type = "GW" if i < num_gw else "GRB"
                
                alert_data = {
                    'time': Time.now().iso,
                    'type': alert_type,
                    'id': f"{alert_type}_{int(Time.now().mjd)}_{i}",
                    'ra': random.uniform(0, 360),
                    'dec': random.uniform(-90, 90),
                    'confidence': random.uniform(0, 1)
                }
                
                # Add to UI
                self.root.after(0, lambda a=alert_data: self._add_alert_to_ui(a))
                
                # Brief delay to space out alerts
                time.sleep(0.5)
            
            # Update status when done
            self.root.after(0, lambda: self.status_var.set(f"Generated {num_gw} GW and {num_grb} test alerts"))
            
        except Exception as e:
            def _set_error_status(error_msg):
                self.status_var.set(f"Simulation failed: {error_msg}")
            error_msg = str(e)
            self.root.after(0, lambda msg=error_msg: _set_error_status(msg))
    
    def _add_alert_to_ui(self, alert):
        """Add an alert to the treeview"""
        # Add to the list of alerts
        self.alert_list.append(alert)
        
        # Format coordinates for display
        ra_str = f"{alert['ra']:.2f}"
        dec_str = f"{alert['dec']:.2f}"
        
        # Add to the treeview
        self.alerts_tree.insert(
            "", tk.END, 
            values=(alert['time'], alert['type'], alert['id'], ra_str, dec_str)
        )
    
    def _show_alert_details(self, event):
        """Show details for the selected alert"""
        selection = self.alerts_tree.selection()
        if not selection:
            return
        
        item = self.alerts_tree.item(selection[0])
        values = item['values']
        if not values:
            return
        
        # Find the corresponding alert data
        alert_id = values[2]  # Index 2 should be the ID
        alert_data = None
        for alert in self.alert_list:
            if alert['id'] == alert_id:
                alert_data = alert
                break
        
        if not alert_data:
            return
        
        # Display the details
        self.alert_details.delete(1.0, tk.END)
        self.alert_details.insert(tk.END, f"Alert ID: {alert_id}\n")
        self.alert_details.insert(tk.END, f"Type: {alert_data['type']}\n")
        self.alert_details.insert(tk.END, f"Time: {alert_data['time']}\n")
        self.alert_details.insert(tk.END, f"RA: {alert_data['ra']:.4f}\n")
        self.alert_details.insert(tk.END, f"Dec: {alert_data['dec']:.4f}\n")
        
        if 'confidence' in alert_data:
            self.alert_details.insert(tk.END, f"Confidence: {alert_data['confidence']:.2f}\n")
        
        # Add additional details if available
        for key, value in alert_data.items():
            if key not in ['id', 'type', 'time', 'ra', 'dec', 'confidence']:
                self.alert_details.insert(tk.END, f"{key}: {value}\n")
    
    def _generate_plot(self):
        """Generate plot based on current selection"""
        plot_type = self.plot_type.get()
        
        # Clear the figure
        self.figure.clear()
        
        # Add a subplot
        ax = self.figure.add_subplot(111)
        
        if plot_type == "Light Curve":
            self._plot_light_curve(ax)
        elif plot_type == "Finding Chart":
            self._plot_finding_chart(ax)
        elif plot_type == "Color-Magnitude":
            self._plot_color_magnitude(ax)
        
        # Update the canvas
        self.plot_canvas.draw()
    
    def _update_plot(self, event):
        """Update plot when selection changes"""
        self._generate_plot()
    
    def _plot_light_curve(self, ax):
        """Plot a sample light curve"""
        # Sample data - in a real app, this would come from your data
        times = np.linspace(0, 10, 20)
        flux = 10 * np.exp(-0.3 * times) + 0.5 * np.random.randn(len(times))
        
        # Plot
        ax.errorbar(times, flux, yerr=0.2, fmt='o-', color='blue', ecolor='gray', capsize=5)
        ax.set_xlabel('Time (days)')
        ax.set_ylabel('Flux (μJy)')
        ax.set_title('GRB Afterglow Light Curve (J-band)')
        ax.grid(True, linestyle='--', alpha=0.7)
    
    def _plot_finding_chart(self, ax):
        """Plot a sample finding chart"""
        # Sample data - in a real app, this would be an image
        np.random.seed(42)  # For reproducible random noise
        image = np.random.normal(loc=100, scale=15, size=(100, 100))
        
        # Add a fake source
        x, y = 50, 50
        for i in range(-3, 4):
            for j in range(-3, 4):
                if 0 <= x+i < 100 and 0 <= y+j < 100:
                    r = np.sqrt(i*i + j*j)
                    if r < 3:
                        image[y+j, x+i] += 100 * np.exp(-r*r/2)
        
        # Plot
        im = ax.imshow(image, cmap='viridis', origin='lower')
        ax.set_title('NIR Finding Chart')
        ax.set_xlabel('Pixel X')
        ax.set_ylabel('Pixel Y')
        
        # Add a circle around the source
        circle = plt.Circle((x, y), 5, color='r', fill=False, linestyle='--')
        ax.add_patch(circle)
        
        # Add colorbar
        self.figure.colorbar(im, ax=ax, label='Counts')
    
    def _plot_color_magnitude(self, ax):
        """Plot a sample color-magnitude diagram"""
        # Sample data - in a real app, this would come from your catalog
        np.random.seed(0)
        
        # Generate some random stars
        n_stars = 200
        j_mag = np.random.uniform(16, 22, n_stars)
        k_mag = j_mag - 0.8 + 0.4 * np.random.randn(n_stars)
        
        # Add a special object (our GRB)
        grb_j, grb_k = 18.5, 17.2
        
        # Plot
        ax.scatter(j_mag - k_mag, j_mag, alpha=0.6, s=10, color='gray', label='Field Stars')
        ax.scatter(grb_j - grb_k, grb_j, color='red', s=100, marker='*', label='GRB Host')
        
        ax.set_xlabel('J-K Color')
        ax.set_ylabel('J Magnitude')
        ax.set_title('Color-Magnitude Diagram')
        ax.legend()
        
        # Invert y-axis (astronomical convention)
        ax.invert_yaxis()
        
        # Add grid
        ax.grid(True, linestyle='--', alpha=0.3)
    
    def _browse_data_dir(self):
        """Browse for data directory"""
        directory = filedialog.askdirectory()
        if directory:
            self.data_dir_var.set(directory)
    
    def _save_settings(self):
        """Save settings to a JSON file"""
        settings = {
            'data_dir': self.data_dir_var.get(),
            'log_level': self.log_level.get(),
            'grb_source': self.grb_source.get(),
            'gw_source': self.gw_source.get()
        }
        
        try:
            os.makedirs('config', exist_ok=True)
            with open('config/settings.json', 'w') as f:
                json.dump(settings, f, indent=2)
            self.status_var.set("Settings saved")
        except Exception as e:
            messagebox.showerror("Save Error", f"Could not save settings: {e}")
    
    def _load_settings(self):
        """Load settings from a JSON file"""
        try:
            with open('config/settings.json', 'r') as f:
                settings = json.load(f)
            
            self.data_dir_var.set(settings.get('data_dir', './data'))
            self.log_level.set(settings.get('log_level', 'INFO'))
            self.grb_source.set(settings.get('grb_source', 'GCN'))
            self.gw_source.set(settings.get('gw_source', 'LIGO/Virgo'))
            
            self.status_var.set("Settings loaded")
        except FileNotFoundError:
            messagebox.showinfo("Settings", "No saved settings found")
        except Exception as e:
            messagebox.showerror("Load Error", f"Could not load settings: {e}")
    
    def _reset_settings(self):
        """Reset settings to defaults"""
        self.data_dir_var.set('./data')
        self.log_level.set('INFO')
        self.grb_source.set('GCN')
        self.gw_source.set('LIGO/Virgo')
        self.status_var.set("Settings reset to defaults")

# Run the application
if __name__ == "__main__":
    root = tk.Tk()
    app = NIRCounterpartPipeline(root)
    root.mainloop()
