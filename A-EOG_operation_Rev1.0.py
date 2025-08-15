import pygame
import sys
import subprocess
import os
import numpy as np
import pandas as pd
import math
from scipy import linalg, signal
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import threading
import json
from typing import List, Tuple

class AEOG_GUI:
    def __init__(self, root): # called whenever creating an instance of this class, self referring to the instance being created, root being a parameter representing another object instance, see in main()
        self.root = root
        self.root.title("A-EOG System Operation Interface")
        self.root.geometry("1200x800")
        
        # Initialize variables
        self.screen_distance_cm = tk.DoubleVar(value=70.0)
        self.screen_width_cm = tk.DoubleVar(value=52.0)
        self.screen_height_cm = tk.DoubleVar(value=30.0)
        self.screen_width_pix = tk.DoubleVar()
        self.screen_height_pix = tk.DoubleVar()
        self.cal_spot_duration_s = tk.IntVar(value=1) # use self.cal_spot_duration_s.get() to extract the integer value
        self.data_acquisition_duration_s = tk.IntVar(value=0) # in seconds, 0 for indefinite acquisition duration # use self.data_acquisition_duration_s.get() to extract the integer value
        self.baudrate_kbps = tk.StringVar(value="1000") # in kbits/s
        self.file_path = tk.StringVar(value=os.getcwd()) # default files path = current working directory
        self.display_filtered_angles = tk.BooleanVar(value=True)
        self.first_abs_stamp = tk.DoubleVar(value=0)  # Will be set when loading data
        self.last_abs_stamp = tk.DoubleVar(value=1000) # Will be set when loading data

        # Calibration data storage
        self.is_calibrating = False
        self.is_acquiring = False
        pygame.init()
        self.setup_ui()

        # Try to load existing filter data
        if hasattr(self, 'filter_designer'):
            self.filter_designer.load_lpf_filter_data(self.file_path.get())
            self.filter_designer.load_notch_filter_data(self.file_path.get())
            self.filter_designer.load_hpf_filter_data(self.file_path.get())

    def setup_ui(self):
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill='both', expand=True, padx=10, pady=10)
        
        # Configuration, Calibration and Acquisition tab
        config_frame = ttk.Frame(notebook)
        notebook.add(config_frame, text="Configuration, Calibration and Acquisition")
        self.setup_config_and_acqu_tab(config_frame)
        
        # Calibration results tab
        results_frame = ttk.Frame(notebook)
        notebook.add(results_frame, text="Calibration Results")
        self.setup_calibration_results_tab(results_frame)

        # Data acquisition results tab
        data_acquisition_results_frame = ttk.Frame(notebook)
        notebook.add(data_acquisition_results_frame, text="Data Acquisition Results")
        self.setup_data_acquisition_results_tab(data_acquisition_results_frame)
        
        # Data filter design tab
        filter_frame = ttk.Frame(notebook)
        notebook.add(filter_frame, text="Angular responses data filtering")
        self.filter_designer = FiltersDesigner(filter_frame, self)

        # Signals FFT analysis tab
        Debug_FFT(notebook)

        # Signals time domain analysis tab
        Debug_temporal(notebook)

    def setup_config_and_acqu_tab(self, parent):
        main_frame = ttk.Frame(parent)
        main_frame.pack(fill='both', expand=True, padx=10, pady=10)

        # --- PATCH START: Adjust grid weights and status frame width ---

        # Main content frame (left) and status frame (right)
        content_frame = ttk.Frame(main_frame)
        content_frame.grid(row=0, column=0, sticky='nsew')
        status_frame = ttk.LabelFrame(main_frame, text="Status", padding=6)
        status_frame.grid(row=0, column=1, sticky='ns', padx=(10, 0), pady=0)

        # Make left column (content) expand, right column (status) fixed width
        main_frame.columnconfigure(0, weight=3)  # Increase weight for content
        main_frame.columnconfigure(1, weight=0, minsize=280)  # Set minsize for status frame

        main_frame.rowconfigure(0, weight=1)
        content_frame.rowconfigure(0, weight=1)

        # --- Data Acquisition and Storage Configuration ---
        data_config_frame = ttk.LabelFrame(content_frame, text="Data Acquisition and Storage Configuration", padding=6)
        data_config_frame.pack(fill='x', pady=(0, 4))

        ttk.Label(data_config_frame, text="Storage Path:").grid(row=0, column=0, sticky='w', padx=(0, 5), pady=2)
        path_entry = ttk.Entry(data_config_frame, textvariable=self.file_path, width=40)
        path_entry.grid(row=0, column=1, padx=(0, 3), sticky='ew', pady=2)
        ttk.Button(data_config_frame, text="Browse", command=self.browse_folder).grid(row=0, column=2, padx=(3, 0), pady=2)
        ttk.Label(data_config_frame, text="Data Rate (kbit/second):").grid(row=1, column=0, sticky='w', padx=(0, 5), pady=2)
        ttk.Combobox(data_config_frame, textvariable=self.baudrate_kbps, values=["1000", "3090"], 
                     state="readonly", width=10).grid(row=1, column=1, sticky='w', pady=2)
        data_config_frame.columnconfigure(1, weight=1)

        # --- Electrodes configuration ---
        electrodes_frame = ttk.LabelFrame(content_frame, text="Electrodes Configuration", padding=6)
        electrodes_frame.pack(fill='x', pady=(0, 4))
        self.setup_electrodes_ui(electrodes_frame)

        # --- Calibration frame ---
        calibration_frame = ttk.LabelFrame(content_frame, text="Calibration", padding=6)
        calibration_frame.pack(fill='x', pady=(0, 4))

        # Compact: all screen params in one row
        ttk.Label(calibration_frame, text="Screen Distance (cm):").grid(row=0, column=0, sticky='w', padx=(0, 3), pady=2)
        ttk.Entry(calibration_frame, textvariable=self.screen_distance_cm, width=8).grid(row=0, column=1, padx=(0, 8), pady=2)
        ttk.Label(calibration_frame, text="Screen Width (cm):").grid(row=0, column=2, sticky='w', padx=(0, 3), pady=2)
        ttk.Entry(calibration_frame, textvariable=self.screen_width_cm, width=8).grid(row=0, column=3, padx=(0, 8), pady=2)
        ttk.Label(calibration_frame, text="Screen Height (cm):").grid(row=0, column=4, sticky='w', padx=(0, 3), pady=2)
        ttk.Entry(calibration_frame, textvariable=self.screen_height_cm, width=8).grid(row=0, column=5, padx=(0, 8), pady=2)

        # Calibration spot duration and button in one row
        ttk.Label(calibration_frame, text="Calibration Spot Duration (s):").grid(
            row=1, column=0, sticky='w', padx=(0, 3), pady=2
        )
        ttk.Entry(calibration_frame, textvariable=self.cal_spot_duration_s, width=8).grid(
            row=1, column=1, padx=(0, 8), pady=2, sticky='w'
        )
        self.calibrate_btn = ttk.Button(
            calibration_frame, text="Start Calibration",
            command=self.start_calibration, style='Accent.TButton'
        )
        self.calibrate_btn.grid(row=1, column=2, padx=(8, 0), pady=2, sticky='w')

        # --- A-EOG Data Acquisition frame ---
        acquisition_frame = ttk.LabelFrame(content_frame, text="A-EOG Data Acquisition", padding=6)
        acquisition_frame.pack(fill='x', pady=(0, 4))

        ttk.Label(acquisition_frame, text="A-EOG Data Acquisition Duration (s, 0=indefinite):").grid(row=0, column=0, sticky='w', padx=(0, 5), pady=2)
        ttk.Entry(acquisition_frame, textvariable=self.data_acquisition_duration_s, width=8).grid(row=0, column=1, padx=(0, 8), pady=2)
        self.acquire_btn = ttk.Button(acquisition_frame, text="Start A-EOG Data Acquisition",
                                      command=self.start_data_acquisition, style='Accent.TButton')
        self.acquire_btn.grid(row=0, column=2, padx=(8, 0), pady=2, sticky='w')

        # --- Status frame on the right ---
        self.status_text = tk.Text(status_frame, height=18, wrap='word', state='disabled', width=36)
        status_scrollbar = ttk.Scrollbar(status_frame, orient='vertical', command=self.status_text.yview)
        self.status_text.configure(yscrollcommand=status_scrollbar.set)
        self.status_text.pack(side='left', fill='both', expand=True, padx=(0, 0), pady=(0, 0))
        status_scrollbar.pack(side='right', fill='y', padx=(0, 0), pady=(0, 0))

        # Make sure the status frame expands vertically
        status_frame.rowconfigure(0, weight=1)
        status_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(0, weight=1)
        content_frame.rowconfigure(0, weight=1)

        # --- PATCH END ---
        
    def setup_electrodes_ui(self, parent):
        # Initialize electrode selection variables
        if not hasattr(self, 'electrodes_selection_status'):
            self.electrodes_selection_status = {}
            for electrode in ['L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6']:
                self.electrodes_selection_status[electrode] = tk.BooleanVar(value=True)
        
        container_frame = ttk.Frame(parent)
        container_frame.pack(fill='x', pady=10)
        
        # Try to use PNG electrode layout, fallback to simple grid
        png_path = "A-EOG_electrodes_configuration.png"
        if self.setup_image_electrode_interface(container_frame, png_path):
            instruction_text = "Click on electrodes labels in the diagram to select/deselect them"
        else:
            self.setup_simple_electrode_grid(container_frame)
            instruction_text = "Click on electrode buttons to select/deselect active electrodes"
        
        # --- PATCH START: Horizontal row for instructions and button ---
        instruction_row = ttk.Frame(container_frame)
        instruction_row.pack(fill='x', pady=(10, 5))
        
        # Add instruction text (left)
        instruction_label = ttk.Label(instruction_row, text=instruction_text)
        instruction_label.pack(side='left', padx=(0, 10))
        
        # Add selection summary (left, after instruction)
        self.electrode_summary = ttk.Label(instruction_row, text="Selected electrodes: All")
        self.electrode_summary.pack(side='left', padx=(0, 10))
        
        # Add Re-calculate Calibration button (right)
        ttk.Button(
            instruction_row,
            text="Re-calculate Calibration",
            command=self.recalculate_calibration
        ).pack(side='right', padx=(10, 0))
        # --- PATCH END ---
        
        # Update initial display
        self.update_electrode_display()

    def setup_simple_electrode_grid(self, parent): # in case the png file illustrating the electrodes around the eyes is not found in the same folder as the python file...
        """Create a simple grid-based electrode selection interface"""
        # Create main frame for electrode layout
        electrode_main_frame = ttk.Frame(parent)
        electrode_main_frame.pack(pady=10)
        
        # Title
        ttk.Label(electrode_main_frame, text="Select Active Electrodes", font=('Arial', 12, 'bold')).pack(pady=(0, 10))
        
        # Create two columns for left and right electrodes
        columns_frame = ttk.Frame(electrode_main_frame)
        columns_frame.pack()
        
        # Left electrodes frame
        left_frame = ttk.LabelFrame(columns_frame, text="Left Eye Electrodes", padding=10)
        left_frame.pack(side='left', padx=(0, 20))
        
        # Right electrodes frame  
        right_frame = ttk.LabelFrame(columns_frame, text="Right Eye Electrodes", padding=10)
        right_frame.pack(side='right')
        
        # Store electrode widgets for later reference
        self.electrode_widgets = {}
        
        # Left electrodes (L1-L5)
        left_electrodes = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6'] # L6 is the electrode placed at the right eye cantus as measured/used for the left eye
        for i, electrode in enumerate(left_electrodes):
            var = self.electrodes_selection_status[electrode]
            btn = ttk.Checkbutton(left_frame, text=f"{electrode}", variable=var, 
                                command=lambda e=electrode: self.on_electrode_change(e))
            btn.grid(row=i, column=0, sticky='w', pady=2, padx=5)
            self.electrode_widgets[electrode] = {'button': btn}
        
        # Right electrodes (R1-R5)  
        right_electrodes = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6'] # R6 is the electrode placed at the left eye cantus as measured/used for the right eye
        for i, electrode in enumerate(right_electrodes):
            var = self.electrodes_selection_status[electrode]
            btn = ttk.Checkbutton(right_frame, text=f"{electrode}", variable=var,
                                command=lambda e=electrode: self.on_electrode_change(e))
            btn.grid(row=i, column=0, sticky='w', pady=2, padx=5)
            self.electrode_widgets[electrode] = {'button': btn}

    def on_electrode_change(self, electrode):
        """Handle electrode selection change"""
        self.update_electrode_display()

    def setup_image_electrode_interface(self, parent, image_path):
        """Setup image-based electrode selection with clickable regions"""
        try:
            from PIL import Image, ImageTk
            import os
            
            if not os.path.exists(image_path):
                print(f"Electrode image not found: {image_path}")
                return False
                
            # Load and display the electrode image
            img = Image.open(image_path)
            
            # Resize image to fit nicely in the interface (maintain aspect ratio)
            max_width, max_height = 600, 400
            img_width, img_height = img.size
            
            # Calculate scaling factor
            scale_width = max_width / img_width
            scale_height = max_height / img_height  
            scale = min(scale_width, scale_height)
            
            new_width = int(img_width * scale)
            new_height = int(img_height * scale)
            
            img = img.resize((new_width, new_height), Image.Resampling.LANCZOS)
            self.electrode_photo = ImageTk.PhotoImage(img)
            
            # Store scaling factor for coordinate adjustment
            self.image_scale = scale
            self.image_width = new_width
            self.image_height = new_height
            
            # Create canvas frame
            canvas_frame = ttk.Frame(parent)
            canvas_frame.pack(pady=10)
            
            # Create canvas with proper size
            electrode_canvas = tk.Canvas(canvas_frame, width=new_width, height=new_height, bg='white')
            electrode_canvas.pack()
            electrode_canvas.create_image(new_width//2, new_height//2, image=self.electrode_photo)
            
            # Get electrode regions (adjusted for scaling)
            electrode_regions = self.get_electrode_click_regions_scaled()
            
            # Create clickable regions
            self.electrode_widgets = {}
            for electrode, region_coords in electrode_regions.items():
                if len(region_coords) == 4:  # Rectangle: (x1, y1, x2, y2)
                    x1, y1, x2, y2 = region_coords
                    # Create invisible clickable rectangle
                    rect = electrode_canvas.create_rectangle(x1, y1, x2, y2, 
                                                        outline='', fill='', width=0)
                    # Create visible indicator rectangle (initially hidden)
                    indicator = electrode_canvas.create_rectangle(x1, y1, x2, y2,
                                                                outline='lime', fill='', width=3,
                                                                state='hidden')
                elif len(region_coords) >= 6:  # Polygon
                    rect = electrode_canvas.create_polygon(region_coords,
                                                        outline='', fill='', width=0)
                    indicator = electrode_canvas.create_polygon(region_coords,
                                                            outline='lime', fill='', width=3,
                                                            state='hidden')
                
                # Make clickable
                electrode_canvas.tag_bind(rect, '<Button-1>', 
                                        lambda e, elec=electrode: self.toggle_electrode(elec))
                
                # Store references
                self.electrode_widgets[electrode] = {
                    'rect': rect,
                    'indicator': indicator,
                    'canvas': electrode_canvas,
                    'coords': region_coords
                }
            
            print(f"Successfully loaded electrode image: {image_path}")
            print("Click on electrode positions in the image to select/deselect them")
            return True
            
        except ImportError:
            print("PIL (Pillow) not available. Install with: pip install pillow")
            return False
        except Exception as e:
            print(f"Error loading electrode image: {e}")
            return False

    def get_electrode_click_regions_scaled(self):
        """Get electrode click regions scaled to current image size"""
        # Get base regions
        base_regions = self.get_electrode_click_regions()
        
        # Scale coordinates if image was resized
        if hasattr(self, 'image_scale'):
            scaled_regions = {}
            for electrode, coords in base_regions.items():
                scaled_coords = tuple(int(coord * self.image_scale) for coord in coords)
                scaled_regions[electrode] = scaled_coords
            return scaled_regions
        else:
            return base_regions

    def toggle_electrode(self, electrode):
        """Toggle electrode selection state"""
        current_state = self.electrodes_selection_status[electrode].get()
        self.electrodes_selection_status[electrode].set(not current_state)
        self.update_electrode_display()

    def get_electrode_click_regions(self):
        """
        Defines clickable regions for each electrode based on the electrodes layout defined in the PNG file.
        Returns dictionary with electrode names as keys and click regions as values.
        Regions are rectangles (x1, y1, x2, y2) but could be polygons (x1, y1, x2, y2, x3, y3, ...)
        """
        electrode_regions = {
            # Left eye electrodes
            'L1': (365, 3, 394, 32),    # Rectangle: (x1, y1, x2, y2)
            'L2': (442, 3, 471, 32),    
            'L3': (556, 86, 585, 115),   
            'L4': (442, 210, 471, 239),   
            'L5': (365, 210, 394, 239),   
            'L6': (8, 131, 37, 160),    # Right eye cantus electrode as measured/used for the left eye
            
            # Right eye electrodes  
            'R1': (197, 5, 226, 34),    
            'R2': (120, 5, 149, 34),   
            'R3': (8, 80, 37, 109),   
            'R4': (120, 210, 149, 239),   
            'R5': (197, 210, 226, 239),   
            'R6': (556, 137, 585, 166)    # Left eye cantus electrode as measured/used for the right eye
        }
        
        return electrode_regions

    def update_electrode_display(self):
        """Update the visual display of electrode selections"""
        selected_electrodes = []
        
        for electrode, var in self.electrodes_selection_status.items():
            if var.get():
                selected_electrodes.append(electrode)
                
                # Show selection indicator if using image interface
                if (hasattr(self, 'electrode_widgets') and 
                    electrode in self.electrode_widgets and 
                    'indicator' in self.electrode_widgets[electrode]):
                    canvas = self.electrode_widgets[electrode]['canvas']
                    indicator = self.electrode_widgets[electrode]['indicator']
                    canvas.itemconfig(indicator, state='normal')
            else:
                # Hide selection indicator if using image interface  
                if (hasattr(self, 'electrode_widgets') and 
                    electrode in self.electrode_widgets and 
                    'indicator' in self.electrode_widgets[electrode]):
                    canvas = self.electrode_widgets[electrode]['canvas']
                    indicator = self.electrode_widgets[electrode]['indicator']
                    canvas.itemconfig(indicator, state='hidden')
        
        # Update summary
        if selected_electrodes:
            summary_text = f"Selected electrodes: {', '.join(sorted(selected_electrodes))}"
        else:
            summary_text = "Selected electrodes: None"
        
        self.electrode_summary.config(text=summary_text)

    def get_selected_electrodes(self):
        """Return list of selected electrodes"""
        return [electrode for electrode, var in self.electrodes_selection_status.items() if var.get()]
        
    def setup_calibration_results_tab(self, parent):
        self.cal_fig = Figure(figsize=(12, 5))
        self.cal_ax_left = self.cal_fig.add_subplot(121)
        self.cal_ax_right = self.cal_fig.add_subplot(122)
        self.setup_plots('Calibration', self.cal_ax_left, self.cal_ax_right)
        
        canvas = FigureCanvasTkAgg(self.cal_fig, parent)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True, padx=10, pady=10)
        
        summary_frame = ttk.LabelFrame(parent, text="Calibration Summary", padding=10)
        summary_frame.pack(fill='x', padx=10, pady=(0, 10))
        
        self.cal_summary_text = tk.Text(summary_frame, height=6, wrap='word', state='disabled')
        summary_scrollbar = ttk.Scrollbar(summary_frame, orient='vertical', command=self.cal_summary_text.yview)
        self.cal_summary_text.configure(yscrollcommand=summary_scrollbar.set)
        self.cal_summary_text.pack(side='left', fill='both', expand=True)
        summary_scrollbar.pack(side='right', fill='y')
        
    def setup_data_acquisition_results_tab(self, parent):
        self.acq_fig = Figure(figsize=(12, 5), constrained_layout=True)
        self.acq_ax_left = self.acq_fig.add_subplot(121)
        self.acq_ax_right = self.acq_fig.add_subplot(122)
        self.setup_plots('Data Acquisition', self.acq_ax_left, self.acq_ax_right)

        self.acq_canvas = FigureCanvasTkAgg(self.acq_fig, parent)
        self.acq_canvas.draw()
        self.acq_canvas.get_tk_widget().pack(fill='both', expand=True, padx=10, pady=10)

        # --- Time Interval Selection Frame (with filter button to the right) ---
        self.time_interval_frame = ttk.LabelFrame(parent, text="Select Two Time Intervals (ms)", padding=10)
        self.time_interval_frame.pack(fill='x', padx=10, pady=5)

        # Create a frame for time bars and filter button in a row
        time_and_filter_row = ttk.Frame(self.time_interval_frame)
        time_and_filter_row.pack(fill='x')

        # Frame for the two time bars (left side)
        time_bars_frame = ttk.Frame(time_and_filter_row)
        time_bars_frame.pack(side='left', fill='x', expand=True)

        # Default relative time min and max values (will be updated when data is loaded)
        self.time_rel_min = tk.DoubleVar(value=0)
        self.time_rel_max = tk.DoubleVar(value=1000)

        # Relative time Interval 1 variables
        self.interval1_rel_start = tk.DoubleVar(value=self.time_rel_min)
        self.interval1_rel_end = tk.DoubleVar(value=self.time_rel_max)
        # Relative time Interval 2 variables
        self.interval2_rel_start = tk.DoubleVar(value=self.time_rel_min)
        self.interval2_rel_end = tk.DoubleVar(value=self.time_rel_max)

        def create_rel_time_bar_row(label, start_var, end_var):
            row = ttk.Frame(time_bars_frame)
            row.pack(fill='x', pady=2)
            ttk.Label(row, text=label).pack(side='left', padx=(0, 5))
            start_entry = ttk.Entry(row, textvariable=start_var, width=8)
            start_entry.pack(side='left')
            # Start scale (min cursor, leftmost bar)
            start_scale = tk.Scale(row, from_=self.time_rel_min.get(), to=self.time_rel_max.get(),
                                orient='horizontal', variable=start_var, showvalue=0, length=300,
                                command=lambda v, sv=start_var, ev=end_var: self._sync_time_entry(sv, ev, True))
            start_scale.pack(side='left', padx=2)
            # End scale (max cursor, rightmost bar)
            end_scale = tk.Scale(row, from_=self.time_rel_min.get(), to=self.time_rel_max.get(),
                                orient='horizontal', variable=end_var, showvalue=0, length=300,
                                command=lambda v, sv=start_var, ev=end_var: self._sync_time_entry(sv, ev, False))
            end_scale.pack(side='left', padx=2)
            end_entry = ttk.Entry(row, textvariable=end_var, width=8)
            end_entry.pack(side='left')
            return (start_scale, end_scale, start_entry, end_entry)

        # Create both relative time interval bars
        self.interval1_rel_widgets = create_rel_time_bar_row("Interval 1:", self.interval1_rel_start, self.interval1_rel_end)
        self.interval2_rel_widgets = create_rel_time_bar_row("Interval 2:", self.interval2_rel_start, self.interval2_rel_end)

        # Filter button frame (right side)
        filter_btn_frame = ttk.Frame(time_and_filter_row)
        filter_btn_frame.pack(side='right', padx=10, pady=5, anchor='e')
        ttk.Checkbutton(
            filter_btn_frame,
            text="Display Filtered Angles",
            variable=self.display_filtered_angles,
            command=self.update_acquisition_plots
        ).pack(side='right')
        ttk.Button(
            filter_btn_frame,
            text="Update Plots",
            command=self.update_acquisition_plots
        ).pack(side='right', padx=(10, 0))

        # --- End Time Interval Frame ---

        summary_frame = ttk.LabelFrame(parent, text="Data Acquisition Summary", padding=10)
        summary_frame.pack(fill='x', padx=10, pady=(0, 10))

        self.acq_summary_text = tk.Text(summary_frame, height=6, wrap='word', state='disabled')
        summary_scrollbar = ttk.Scrollbar(summary_frame, orient='vertical', command=self.acq_summary_text.yview)
        self.acq_summary_text.configure(yscrollcommand=summary_scrollbar.set)
        self.acq_summary_text.pack(side='left', fill='both', expand=True)
        summary_scrollbar.pack(side='right', fill='y')

    def _sync_time_entry(self, start_var, end_var, is_start):
        """Ensure start <= end and update text fields when sliders move."""
        if is_start:
            if start_var.get() > end_var.get():
                start_var.set(end_var.get())
        else:
            if end_var.get() < start_var.get():
                end_var.set(start_var.get())

    def update_rel_time_interval_selectors(self, t_rel_min, t_rel_max):
        """Call this after loading new data to update the time bar ranges."""
        self.time_rel_min.set(t_rel_min)
        self.time_rel_max.set(t_rel_max)
        # Update all scales' range
        for scale, _, _, _ in [self.interval1_rel_widgets, self.interval2_rel_widgets]:
            scale.config(from_=t_rel_min, to=t_rel_max)
        for _, scale, _, _ in [self.interval1_rel_widgets, self.interval2_rel_widgets]:
            scale.config(from_=t_rel_min, to=t_rel_max)
        # Reset values to the new range
        self.interval1_rel_start.set(t_rel_min)
        self.interval1_rel_end.set(t_rel_max)
        self.interval2_rel_start.set(t_rel_min)
        self.interval2_rel_end.set(t_rel_max)

    def setup_plots(self, context_suffix, ax_left, ax_right): # context_suffix is appended to the graphs titles
        for ax, title in [(ax_left, 'Left Eye'), (ax_right, 'Right Eye')]:
            ax.clear()
            ax.set_title(f'{title} {context_suffix}')
            ax.set_xlabel('Azimuth (degrees)')
            ax.set_ylabel('Elevation (degrees)')
            ax.grid(True, alpha=0.3)
            ax.set_aspect('equal')
        
    def log_status(self, message):
        self.status_text.config(state='normal')
        self.status_text.insert(tk.END, f"{message}\n")
        self.status_text.see(tk.END)
        self.status_text.config(state='disabled')
        self.root.update_idletasks()
        
    def browse_folder(self):
        # --- Clear all output graphs and filter info ---
        # Calibration Results tab
        if hasattr(self, 'cal_ax_left') and hasattr(self, 'cal_ax_right') and hasattr(self, 'cal_fig'):
            self.cal_ax_left.clear()
            self.cal_ax_right.clear()
            self.cal_fig.canvas.draw()
            if hasattr(self, 'cal_summary_text'):
                self.cal_summary_text.config(state='normal')
                self.cal_summary_text.delete(1.0, tk.END)
                self.cal_summary_text.config(state='disabled')

        # Data Acquisition Results tab
        if hasattr(self, 'acq_ax_left') and hasattr(self, 'acq_ax_right') and hasattr(self, 'acq_fig'):
            self.acq_ax_left.clear()
            self.acq_ax_right.clear()
            self.acq_fig.canvas.draw()
            if hasattr(self, 'acq_summary_text'):
                self.acq_summary_text.config(state='normal')
                self.acq_summary_text.delete(1.0, tk.END)
                self.acq_summary_text.config(state='disabled')

        # Angular responses data filtering tab (filter response and info)
        if hasattr(self, 'filter_designer'):
            # Clear filter info
            if hasattr(self.filter_designer, 'lpf_info'):
                self.filter_designer.lpf_info.config(state='normal')
                self.filter_designer.lpf_info.delete(1.0, tk.END)
                self.filter_designer.lpf_info.config(state='disabled')
            # Clear filter response plot
            if hasattr(self.filter_designer, 'lpf_response_frame'):
                for widget in self.filter_designer.lpf_response_frame.winfo_children():
                    widget.destroy()
                self.filter_designer.lpf_response_frame.pack_forget()

        if self.is_calibrating or self.is_acquiring:
            messagebox.showwarning("Warning", "Calibration or A-EOG Data Acquisition already in progress!")
            return
        folder = filedialog.askdirectory(initialdir=self.file_path.get())
        if folder:
            self.file_path.set(folder)
            # Auto-load calibration data if available
            cal_file = os.path.join(folder, "calibration_data.json")
            if os.path.exists(cal_file):
                self.load_calibration_from_path(cal_file)
            else:
                self.log_status("No calibration data file found in the folder.")
            
            # Try to load filters data from new path
            filter_file = os.path.join(folder, "lpf_filter_data.json")
            if os.path.exists(filter_file):
                self.filter_designer.load_lpf_filter_data(folder)
            else:
                self.log_status("No Low-Pass Filter data file found in the folder.")
            
            filter_file = os.path.join(folder, "notch_filter_data.json")
            if os.path.exists(filter_file):
                self.filter_designer.load_notch_filter_data(folder)
            else:
                self.log_status("No Notch Filter data file found in the folder.")

            filter_file = os.path.join(folder, "hpf_filter_data.json")
            if os.path.exists(filter_file):
                self.filter_designer.load_hpf_filter_data(folder)
            else:
                self.log_status("No High-Pass Filter data file found in the folder.")

            temp_first_abs_stamp, temp_last_abs_stamp = self.get_first_and_last_timestamp(os.path.join(folder, "Processed_angles.csv"))
            self.first_abs_stamp.set(temp_first_abs_stamp)
            self.last_abs_stamp.set(temp_last_abs_stamp)
#            self.update_rel_time_interval_selectors(self.first_abs_stamp.get(), self.last_abs_stamp.get())
            self.update_rel_time_interval_selectors(self.first_abs_stamp.get()-self.first_abs_stamp.get(), self.last_abs_stamp.get()-self.first_abs_stamp.get()) # Convert to relative time

            self.update_acquisition_plots() # Update data output graphs

    def get_first_and_last_timestamp(self, csv_path):
        # Read only the first row (after header)
        first_row = pd.read_csv(csv_path, nrows=1)
        first_timestamp = first_row.iloc[0, 0]

        # Read only the last row efficiently
        # Use iterator to avoid loading the whole file
        last_timestamp = None
        for chunk in pd.read_csv(csv_path, usecols=[0], chunksize=100000):
            last_timestamp = chunk.iloc[-1, 0]

        return first_timestamp, last_timestamp

    def start_data_acquisition(self):
        if self.is_calibrating or self.is_acquiring:
            messagebox.showwarning("Warning", "Calibration or A-EOG Data Acquisition already in progress!")
            return

            # Validate and create directory
        if not os.path.exists(self.file_path.get()):
            try:
                os.makedirs(self.file_path.get())
            except Exception as e:
                messagebox.showerror("Error", f"Cannot create directory: {e}")
                return
                
        # Check executables
        current_dir = os.path.dirname(os.path.realpath(__file__))
        executables = {
            "acq": os.path.join(current_dir, "A-EOG_acqu_generic.exe"),
            "csv": os.path.join(current_dir, "bin2csv.exe")
        }
        
        for name, path in executables.items():
            if not os.path.exists(path):
                messagebox.showerror("Error", f"{name.upper()} executable not found: {path}")
                return
                
        self.is_acquiring = True
        self.acquire_btn.config(state='disabled', text='Acquiring...')
        self.log_status("Starting A-EOG Data Acquisition...")
        
        threading.Thread(target=self.run_AEOG_acquisition, daemon=True).start()
    
    def run_AEOG_acquisition(self):
        """Acquisition method with proper error handling."""
        try:
            spots = self.create_calibration_spots()
            screen = pygame.display.get_surface() or pygame.display.set_mode((0, 0), pygame.FULLSCREEN)
            pygame.display.set_caption("A-EOG Data Acquisition")
        
            colors = {'blue': (0, 0, 255), 'red': (255, 0, 0), 'green': (0, 255, 0), 'black': (0, 0, 0)}
            current_dir = os.path.dirname(os.path.realpath(__file__))
            
            # Generate filenames
            acqu_bin_file = f"Acquisition_data.bin"
            acqu_csv_file = f"Acquisition_data.csv"
            acqu_processed_angles_file = f"Processed_angles.csv"

            # Display calibration spots
            screen.fill(colors['black'])
            for idx, spot in enumerate(spots, 1):
                x, y = spot['position']
                diameter = spot['diameter']
                pygame.draw.circle(screen, colors['green'], (int(x), int(y)), diameter // 2)
            pygame.display.flip()
            
            self.wait_for_spacebar()
            
            # Change spots to red during acquisition
            for idx, spot in enumerate(spots, 1):
                x, y = spot['position']
                diameter = spot['diameter']
                pygame.draw.circle(screen, colors['red'], (int(x), int(y)), diameter // 2)
            pygame.display.flip()

            # Acquire data
            if self.acquire_data(current_dir, acqu_bin_file, self.data_acquisition_duration_s.get()):
                # Change spots color to blue after acquisition
                for idx, spot in enumerate(spots, 1):
                    x, y = spot['position']
                    diameter = spot['diameter']
                    pygame.draw.circle(screen, colors['blue'], (int(x), int(y)), diameter // 2)
                pygame.display.flip()

                # Process the acquired data
                if self.process_acqu_data(os.path.join(self.file_path.get(), acqu_csv_file),os.path.join(self.file_path.get(), acqu_processed_angles_file)):

                    self.log_status("Acquisition and processing completed successfully!")
                    temp_first_abs_stamp, temp_last_abs_stamp = self.get_first_and_last_timestamp(os.path.join(self.file_path.get(), acqu_processed_angles_file))
                    self.first_abs_stamp.set(temp_first_abs_stamp)
                    self.last_abs_stamp.set(temp_last_abs_stamp)
                    self.update_rel_time_interval_selectors(self.first_abs_stamp.get()-self.first_abs_stamp.get(), self.last_abs_stamp.get()-self.first_abs_stamp.get()) # Convert to relative time
                    self.update_acquisition_plots()
                else:
                    self.log_status("Acquisition completed but processing failed!")
            else:
                self.log_status("Data acquisition failed!")

            pygame.quit()
            
        except Exception as e:
            self.log_status(f"Acquisition failed: {str(e)}")
            messagebox.showerror("Error", f"Acquisition failed: {str(e)}")
        finally:
            self.is_acquiring = False
            self.root.after(0, lambda: self.acquire_btn.config(state='normal', text='Start A-EOG Data Acquisition'))
    
    def start_calibration(self):
        if self.is_calibrating or self.is_acquiring:
            messagebox.showwarning("Warning", "Calibration or A-EOG Data Acquisition already in progress!")
            return
            
        # Validate and create directory
        if not os.path.exists(self.file_path.get()):
            try:
                os.makedirs(self.file_path.get())
            except Exception as e:
                messagebox.showerror("Error", f"Cannot create directory: {e}")
                return
                
        # Check executables
        current_dir = os.path.dirname(os.path.realpath(__file__))
        executables = {
            "acq": os.path.join(current_dir, "A-EOG_acqu_generic.exe"),
            "csv": os.path.join(current_dir, "bin2csv.exe")
        }
        
        for name, path in executables.items():
            if not os.path.exists(path):
                messagebox.showerror("Error", f"{name.upper()} executable not found: {path}")
                return
                
        self.is_calibrating = True
        self.calibrate_btn.config(state='disabled', text='Calibrating...')
        self.log_status("Starting calibration...")
        
        threading.Thread(target=self.run_calibration, daemon=True).start()

    def run_calibration(self):
        try:
            self.initialize_calibration_data()
            spots = self.create_calibration_spots()
            self.run_calibration_sequence(spots)
            self.process_calibration_data()
            self.root.after(0, self.update_cal_plots_and_outcome('Calibration'))
            self.log_status("Calibration completed successfully!")
            self.save_calibration()
            
        except Exception as e:
            self.log_status(f"Calibration failed: {str(e)}")
            messagebox.showerror("Error", f"Calibration failed: {str(e)}")
        finally:
            self.is_calibrating = False
            self.root.after(0, lambda: self.calibrate_btn.config(state='normal', text='Start Calibration'))
            
    def initialize_calibration_data(self):
        '''Initialize all calibration data arrays ot the class to empty lists'''
        arrays = ['LAz_calibration_angles', 'RAz_calibration_angles', 'LEl_calibration_angles', 
                 'REl_calibration_angles', 'LAz_voltages_combinations', 'LEl_voltages_combinations',
                 'RAz_voltages_combinations', 'REl_voltages_combinations', 'all_voltages_left', 'all_voltages_right']
        for array in arrays:
            setattr(self, array, [])    # Initialize all calibration data arrays of the class to empty lists
        
    def create_calibration_spots(self):
        screen = pygame.display.set_mode((0, 0), pygame.FULLSCREEN)
        width, height = screen.get_size()
        self.screen_width_pix.set(width)
        self.screen_height_pix.set(height)
        # Create 5x5 grid and select 10-point calibration pattern
        rows, cols, diameter = 5, 5, 20
        cal_spots = [[None for _ in range(cols)] for _ in range(rows)]

        # Use .get() to access DoubleVar values
        width_pix = self.screen_width_pix.get()
        height_pix = self.screen_height_pix.get()

        for r in range(rows):
            for c in range(cols):
                x = width_pix / 2 + (width_pix - diameter) * (c - (cols - 1) / 2) / (cols - 1)
                y = height_pix / 2 + (height_pix - diameter) * (r - (rows - 1) / 2) / (rows - 1)
                cal_spots[r][c] = {'position': (x, y), 'diameter': diameter}

        # 10-point calibration sequence
        center = [cal_spots[2][2]]
        periphery = [cal_spots[0][0], cal_spots[2][4], cal_spots[4][0], cal_spots[0][2],
                    cal_spots[4][4], cal_spots[2][0], cal_spots[0][4], cal_spots[4][2]]

        return center + periphery + center
        
    def run_calibration_sequence(self, spots):
        screen = pygame.display.get_surface() or pygame.display.set_mode((0, 0), pygame.FULLSCREEN)
        pygame.display.set_caption("EOG Calibration")
        
        colors = {'blue': (0, 0, 255), 'red': (255, 0, 0), 'green': (0, 255, 0), 'black': (0, 0, 0)}
        current_dir = os.path.dirname(os.path.realpath(__file__))
        
        for idx, spot in enumerate(spots, 1):
            self.log_status(f"Point {idx}/{len(spots)}")
            
            x, y = spot['position']
            diameter = spot['diameter']
            
            # Generate filenames
            cal_bin_file = f"cal_seq{idx}_x{int(x)}_y{int(y)}.bin"
            cal_csv_file = f"cal_seq{idx}_x{int(x)}_y{int(y)}.csv"
            
            # Show green spot and wait
            screen.fill(colors['black'])
            pygame.draw.circle(screen, colors['green'], (int(x), int(y)), diameter // 2)
            pygame.display.flip()
            self.wait_for_spacebar()
            
            # Show red spot during acquisition
            pygame.draw.circle(screen, colors['red'], (int(x), int(y)), diameter // 2)
            pygame.display.flip()
            
            if self.acquire_data(current_dir, cal_bin_file,self.cal_spot_duration_s.get()):
                # Turn spot color to blue after successful calibration spot data acquisition
                pygame.draw.circle(screen, colors['blue'], (int(x), int(y)), diameter // 2)
                pygame.display.flip()

                self.process_spot_data(idx, x, y, cal_csv_file)
            
        pygame.quit()

    def recalculate_calibration(self):
        try:
            # Check if calibration data exists
            if not hasattr(self, 'all_voltages_left') or not self.all_voltages_left:
                messagebox.showwarning("Warning", "No calibration data available! Please run calibration first or load existing calibration.")
                return
            
            self.log_status("Re-calculating calibration with current electrode selection...")
            
            # Get current electrode selection
            selected_electrodes = self.get_selected_electrodes()
            left_labels = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6']
            right_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']
            
            left_indices = [i for i, lbl in enumerate(left_labels) if lbl in selected_electrodes]
            right_indices = [i for i, lbl in enumerate(right_labels) if lbl in selected_electrodes]
            
            # Re-build voltage combinations with current electrode selection
            self.LAz_voltages_combinations = []
            self.LEl_voltages_combinations = []
            self.RAz_voltages_combinations = []
            self.REl_voltages_combinations = []
            
            for i in range(len(self.all_voltages_left)):
                # Filter voltages based on current selection
                left_selected = [self.all_voltages_left[i][j] for j in left_indices]
                right_selected = [self.all_voltages_right[i][j] for j in right_indices]
                
                # Add bias term
                left_coef = left_selected + [1]
                right_coef = right_selected + [1]
                
                self.LAz_voltages_combinations.append(left_coef)
                self.LEl_voltages_combinations.append(left_coef)
                self.RAz_voltages_combinations.append(right_coef)
                self.REl_voltages_combinations.append(right_coef)
            
            # Convert to numpy arrays and re-process
            for array in ['LAz_voltages_combinations', 'LEl_voltages_combinations', 
                        'RAz_voltages_combinations', 'REl_voltages_combinations']:
                setattr(self, array, np.array(getattr(self, array)))
            
            # Re-calculate calibration parameters
            self.process_calibration_data()
            # Update calibration plots and summary and save calibration data
            self.update_cal_plots_and_outcome('Calibration')
            self.log_status("Calibration re-calculated successfully!")
            self.save_calibration()
            
            # re-process acquisition data to raw angles
            self.process_acqu_data(os.path.join(self.file_path.get(), "Acquisition_data.csv"),os.path.join(self.file_path.get(), "Processed_angles.csv"))
            # generate filtered angles file if filters are designed and active
            if not hasattr(self, 'filter_designer') or self.filter_designer.lpf_filter_coeffs is None:
                self.log_status("No filtered data file produced since no filter data was found.")
                return
            self.apply_filters_to_data()
            self.update_acquisition_plots()

        except Exception as e:
            self.log_status(f"Re-calculation failed: {str(e)}")
            messagebox.showerror("Error", f"Re-calculation failed: {str(e)}")

    def wait_for_spacebar(self):
        while True:
            for event in pygame.event.get():
                if event.type == pygame.QUIT or (event.type == pygame.KEYDOWN and event.key == pygame.K_ESCAPE):
                    pygame.quit()
                    sys.exit()
                elif event.type == pygame.KEYDOWN and event.key == pygame.K_SPACE:
                    return

    def acquire_data(self, current_dir, bin_file, duration_s):
        try:
            # Run acquisition
            acq_cmd = [os.path.join(current_dir, "A-EOG_acqu_generic.exe"),
                      self.baudrate_kbps.get(), str(duration_s), self.file_path.get(), bin_file]
            if duration_s > 0:
                result = subprocess.run(acq_cmd, capture_output=True, text=True, timeout=duration_s*3)
                if result.returncode != 0:
                    self.log_status(f"Acquisition failed: {result.returncode}")
                    return False
            else:
                process = subprocess.Popen(acq_cmd)
                self.wait_for_spacebar()
                process.terminate()
                process.wait()

            # Convert to CSV
            csv_cmd = [os.path.join(current_dir, "bin2csv.exe"), 
                      os.path.join(self.file_path.get(), bin_file)]
            if duration_s > 0:
                result = subprocess.run(csv_cmd, capture_output=True, text=True, timeout=duration_s*3)
            else:
                result = subprocess.run(csv_cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                self.log_status(f"CSV conversion failed: {result.returncode}")
                return False
                
            return True
            
        except subprocess.TimeoutExpired:
            self.log_status("Process timeout")
            return False
        except Exception as e:
            self.log_status(f"Acquisition error: {str(e)}")
            return False
            
    def process_spot_data(self, idx, x, y, csv_file):
        try:
            data = pd.read_csv(os.path.join(self.file_path.get(), csv_file))
            
            # Calculate averages for all channels
            voltages_left = [data[f'Ch{i} (uV)'].mean() for i in range(6)]  # voltages for left electrodes averaged for the calibration spot
            voltages_right = [data[f'Ch{i} (uV)'].mean() for i in range(6, 12)] # voltages for right electrodes averaged for the calibration spot
            
            self.all_voltages_left.append(voltages_left)    #append calibration spot voltages to the previous ones
            self.all_voltages_right.append(voltages_right)  
            
            # Use .get() for all tk.Variable values
            width_pix = self.screen_width_pix.get()
            height_pix = self.screen_height_pix.get()
            width_cm = self.screen_width_cm.get()
            height_cm = self.screen_height_cm.get()
            distance_cm = self.screen_distance_cm.get()

            az_rad = math.atan((x - width_pix / 2) * width_cm / width_pix / distance_cm)
            el_rad = math.atan((height_pix / 2 - y) * height_cm / height_pix / distance_cm)
            az_deg, el_deg = math.degrees(az_rad), math.degrees(el_rad)
            
            # Store angles (choice/approximation made so far: same for both eyes, corresponding to the gaze angle)
            for angle_list in [self.LAz_calibration_angles, self.RAz_calibration_angles]:
                angle_list.append(az_deg)
            for angle_list in [self.LEl_calibration_angles, self.REl_calibration_angles]:
                angle_list.append(el_deg)
            
            # Store voltage combinations
            # Identify selected electrodes
            selected_electrodes = self.get_selected_electrodes()
            left_labels = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6']
            right_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']

            left_indices = [i for i, lbl in enumerate(left_labels) if lbl in selected_electrodes]
            right_indices = [i for i, lbl in enumerate(right_labels) if lbl in selected_electrodes]

            # Select/pick voltages to be taken into account for calibration
            left_selected = [voltages_left[i] for i in left_indices]
            right_selected = [voltages_right[i] for i in right_indices]

            # Append bias term
            left_coef = left_selected + [1]
            right_coef = right_selected + [1]
            
            self.LAz_voltages_combinations.append(left_coef)    #append one row to the calibration voltages matrix for this additional calibration spot
            self.LEl_voltages_combinations.append(left_coef)
            self.RAz_voltages_combinations.append(right_coef)
            self.REl_voltages_combinations.append(right_coef)
            
        except Exception as e:
            self.log_status(f"Error processing spot {idx}: {str(e)}")
            
    def process_calibration_data(self):
        try:
            # Convert all relevant arrays to numpy arrays
            arrays = [
                'LAz_voltages_combinations', 'LEl_voltages_combinations',
                'RAz_voltages_combinations', 'REl_voltages_combinations',
                'LAz_calibration_angles', 'LEl_calibration_angles',
                'RAz_calibration_angles', 'REl_calibration_angles'
            ]
            for array in arrays:
                setattr(self, array, np.array(getattr(self, array)))

            # Process each calibration component
            params = ['LAz', 'LEl', 'RAz', 'REl']
            for param in params:
                voltages = getattr(self, f'{param}_voltages_combinations')
                angles = getattr(self, f'{param}_calibration_angles')

                # Ensure array shapes are valid
                voltages = np.atleast_2d(voltages)
                angles = np.atleast_1d(angles)

                # Log matrix shape
                n_samples, n_features = voltages.shape
                self.log_status(f"Solving least squares for {param}: {n_samples} equations, {n_features} unknowns")

                # Check matrix rank
                rank = np.linalg.matrix_rank(voltages)
                if rank < n_features:
                    self.log_status(f"⚠️ {param}: Rank deficient (rank={rank} < {n_features}). Solution may be non-unique.")

                # Check condition number
                cond_number = np.linalg.cond(voltages)
                self.log_status(f"{param} condition number: {cond_number:.2e}")
                if cond_number > 1e8:
                    self.log_status(f"⚠️ {param}: Ill-conditioned matrix — solution may be unstable.")

                # Solve least squares
                calibration_params, residuals, rank_used, s = linalg.lstsq(voltages, angles)
                setattr(self, f'{param}_calibration_parameters', calibration_params)

                # Evaluate fit
                verified = voltages @ calibration_params
                mse = np.mean((verified - angles) ** 2)
                setattr(self, f'{param}_mse', mse)

                # Log residuals safely
                try:
                    if isinstance(residuals, (list, np.ndarray)) and len(residuals) > 0:
                        self.log_status(f"{param} residual sum of squares: {residuals[0]:.6f}")
                    else:
                        self.log_status(f"{param} residuals not returned — possibly rank-deficient or underdetermined.")
                except Exception as e:
                    self.log_status(f"{param} residual logging failed: {str(e)}")

            # Log summary of MSEs
            mse_values = [float(getattr(self, f'{param}_mse')) for param in params]
            mse_strs = [f"{p}:{v:.6f}" for p, v in zip(params, mse_values)]
            self.log_status("MSE summary: " + ", ".join(mse_strs))

        except Exception as e:
            import traceback
            tb = traceback.format_exc()
            self.log_status(f"Error processing calibration:\n{tb}")
            raise
            
    def update_cal_plots_and_outcome(self, context_update):
        """Update calibration plots with the latest data, and the calibration summary text."""
        try:
            self.setup_plots(context_update, self.cal_ax_left, self.cal_ax_right)
            
            if hasattr(self, 'LAz_calibration_angles') and len(self.LAz_calibration_angles) > 0:
                # Plot data for both eyes
                eye_data = [
                    (self.cal_ax_left, 'LAz', 'LEl'),
                    (self.cal_ax_right, 'RAz', 'REl')
                ]
                
                for ax, az_param, el_param in eye_data:
                    az_angles = getattr(self, f'{az_param}_calibration_angles')
                    el_angles = getattr(self, f'{el_param}_calibration_angles')
                    az_voltages = getattr(self, f'{az_param}_voltages_combinations')
                    el_voltages = getattr(self, f'{el_param}_voltages_combinations')
                    az_params = getattr(self, f'{az_param}_calibration_parameters')
                    el_params = getattr(self, f'{el_param}_calibration_parameters')
                    
                    # Target orientation angles
                    ax.scatter(az_angles, el_angles, c='blue', s=50, label='Target', alpha=0.7)
                    
                    # Verified orientation angles
                    az_verif = az_voltages @ az_params
                    el_verif = el_voltages @ el_params
                    ax.scatter(az_verif, el_verif, c='red', s=30, marker='x', label='Verified')
                    
                    # Connect target to verified
                    for i in range(len(az_angles)):
                        ax.plot([az_angles[i], az_verif[i]], [el_angles[i], el_verif[i]], 
                               'gray', linewidth=0.5, alpha=0.5)
                    
                    ax.legend()
                
            self.cal_fig.tight_layout()
            self.cal_fig.canvas.draw()
            self.update_cal_summary()   # Update calibration summary text as well...
            
        except Exception as e:
            self.log_status(f"Error updating plots: {str(e)}")
            
    def update_cal_summary(self):
        try:
            if hasattr(self, 'LAz_mse'):
                mse_data = {param: getattr(self, f'{param}_mse') for param in ['LAz', 'LEl', 'RAz', 'REl']}
                az_range = (np.min(self.LAz_calibration_angles), np.max(self.LAz_calibration_angles))
                el_range = (np.min(self.LEl_calibration_angles), np.max(self.LEl_calibration_angles))

                # Always show selected electrodes + Bias as last label
                def param_lines(param, base_labels):
                    arr = getattr(self, f"{param}_calibration_parameters")
                    if param.startswith('L'):
                        selected = [lbl for lbl in base_labels if self.electrodes_selection_status[lbl].get()]
                    else:
                        selected = [lbl for lbl in base_labels if self.electrodes_selection_status[lbl].get()]
                    labels = selected + ['Bias']
                    # Always use the last label as "Bias" if there are more parameters than selected electrodes
                    if len(arr) == len(selected) + 1:
                        used_labels = selected + ['Bias']
                    else:
                        used_labels = labels[:len(arr)]
                    return chr(10).join(
                        f'- {lbl}: {arr[i]:.6f}' for i, lbl in enumerate(used_labels)
                    )

                summary = f"""Calibration Summary:
Points: {len(self.LAz_calibration_angles)}
MSE:
{chr(10).join(f'- {param}: {mse:.6f}' for param, mse in mse_data.items())}

LAz Calibration Parameters:
{param_lines('LAz', ['L1', 'L2', 'L3', 'L4', 'L5', 'L6'])}
LEl Calibration Parameters:
{param_lines('LEl', ['L1', 'L2', 'L3', 'L4', 'L5', 'L6'])}
RAz Calibration Parameters:
{param_lines('RAz', ['R1', 'R2', 'R3', 'R4', 'R5', 'R6'])}
REl Calibration Parameters:
{param_lines('REl', ['R1', 'R2', 'R3', 'R4', 'R5', 'R6'])}

Range:
- Azimuth: {az_range[0]:.2f}° to {az_range[1]:.2f}°
- Elevation: {el_range[0]:.2f}° to {el_range[1]:.2f}°
"""

                self.cal_summary_text.config(state='normal')
                self.cal_summary_text.delete(1.0, tk.END)
                self.cal_summary_text.insert(1.0, summary)
                self.cal_summary_text.config(state='disabled')

        except Exception as e:
            self.log_status(f"Error updating summary: {str(e)}")
            
    def save_calibration(self):
        try:
            if not hasattr(self, 'LAz_calibration_parameters'):
                messagebox.showwarning("Warning", "No calibration data to save!")
                return
                
            # Prepare data for saving
            arrays_to_save = ['calibration_angles', 'voltages_combinations', 'calibration_parameters']
            params = ['LAz', 'LEl', 'RAz', 'REl']
            
            calibration_data = {}
            
            for param in params:
                for array_type in arrays_to_save:
                    key = f"{param}_{array_type}"
                    if hasattr(self, key):
                        calibration_data[key] = getattr(self, key).tolist()
                        
                # Add MSE
                if hasattr(self, f"{param}_mse"):
                    calibration_data[f"{param}_mse"] = getattr(self, f"{param}_mse")
            
            # Add additional data
            calibration_data.update({
                "all_voltages_left": self.all_voltages_left,
                "all_voltages_right": self.all_voltages_right,
                "screen_parameters": {
                    "distance_cm": self.screen_distance_cm.get(),
                    "width_cm": self.screen_width_cm.get(),
                    "height_cm": self.screen_height_cm.get(),
                    "width_pix": self.screen_width_pix.get(),
                    "height_pix": self.screen_height_pix.get()
                },
                "selected_electrodes": self.get_selected_electrodes()
            })
            # Add electrode configuration
            calibration_data["electrode_configuration"] = {
                electrode: var.get() for electrode, var in self.electrodes_selection_status.items()
            }
            
            # Save files
            base_path = self.file_path.get()
            np.save(os.path.join(base_path, "calibration_data.npy"), calibration_data)
            
            with open(os.path.join(base_path, "calibration_data.json"), 'w') as f:
                json.dump(calibration_data, f, indent=2)
                
            self.log_status("Calibration saved successfully!")
            
        except Exception as e:
            self.log_status(f"Save error: {str(e)}")
            messagebox.showerror("Error", f"Failed to save: {str(e)}")
            
    def load_calibration_from_path(self, file_path):
        """Load calibration from specific file path"""
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
                    
            # Restore arrays
            params = ['LAz', 'LEl', 'RAz', 'REl']
            array_types = ['calibration_angles', 'voltages_combinations', 'calibration_parameters']
            
            for param in params:
                for array_type in array_types:
                    key = f"{param}_{array_type}"
                    if key in data:
                        setattr(self, key, np.array(data[key]))
                        
                # Restore MSE
                mse_key = f"{param}_mse"
                if mse_key in data:
                    setattr(self, mse_key, data[mse_key])
            
            # Restore additional data
            if "all_voltages_left" in data:
                self.all_voltages_left = data["all_voltages_left"]
                self.all_voltages_right = data["all_voltages_right"]
                
            if "screen_parameters" in data:
                screen_params = data["screen_parameters"]
                self.screen_width_cm.set(screen_params["width_cm"])
                self.screen_height_cm.set(screen_params["height_cm"])
                self.screen_distance_cm.set(screen_params["distance_cm"])
                self.screen_width_pix.set(screen_params["width_pix"])
                self.screen_height_pix.set(screen_params["height_pix"])
            
            # Restore electrodes configuration
            if "electrode_configuration" in data:
                electrode_config = data["electrode_configuration"]
                for electrode, state in electrode_config.items():
                    if electrode in self.electrodes_selection_status:
                        self.electrodes_selection_status[electrode].set(state)
                self.update_electrode_display()

            self.update_cal_plots_and_outcome('Calibration')
            self.log_status("Calibration loaded automatically!")
            
        except Exception as e:
            self.log_status(f"Load error: {str(e)}")
            messagebox.showerror("Error", f"Failed to load: {str(e)}")

    def process_acqu_data(self, csv_file: str, processed_angles_file: str) -> bool:
        """
        Process acquired A-EOG data using calibration parameters.
        
        Args:
            csv_file: Name of the CSV file to process
            processed_angles_file: Name of the output file for processed angles
            
        Returns:
            bool: True if processing succeeded, False otherwise
        """
        if not hasattr(self, 'LAz_calibration_parameters'):
            self.log_status("No calibration data available for processing!")
            return False
        
        try:
            # Get selected indices
            selected_left = self.get_selected_electrodes()
            left_labels = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6']
            right_labels = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6']

            left_indices = [i for i, lbl in enumerate(left_labels) if lbl in selected_left]
            right_indices = [i for i, lbl in enumerate(right_labels) if lbl in selected_left]

            cal_params = {
                'LAz_params': self.LAz_calibration_parameters.tolist(),
                'LEl_params': self.LEl_calibration_parameters.tolist(), 
                'RAz_params': self.RAz_calibration_parameters.tolist(),
                'REl_params': self.REl_calibration_parameters.tolist(),
                'left_indices': left_indices,
                'right_indices': right_indices
            }
            processor = EOGDataProcessor(cal_params)
            
            # Build file paths
            input_path = os.path.join(self.file_path.get(), csv_file)
            output_path = os.path.join(self.file_path.get(), processed_angles_file)
            
            # Check if input file exists
            if not os.path.exists(input_path):
                self.log_status(f"Input file not found: {input_path}")
                return False
            
            # Process the data
            self.log_status(f"Processing data: {csv_file} using current calibration parameters")
            processor.process_csv_chunked(input_path, output_path)
            self.log_status(f"Data processing complete: {processed_angles_file}")
            
            return True
            
        except Exception as e:
            self.log_status(f"Processing failed: {str(e)}")
            return False
    
    def update_acquisition_plots(self):
        """Update the data acquisition results tab with processed data, showing only data within selected intervals."""
        try:
            if self.display_filtered_angles.get():
                # Check if filtered file exists
                angles_file_path = os.path.join(self.file_path.get(), "Filtered_angles.csv")
                if not os.path.exists(angles_file_path):
                    self.display_filtered_angles.set(False)  # Reset checkbox
                    messagebox.showwarning("Warning", "No filtered data available. Please design and apply a filter first.")
                    angles_file_path = os.path.join(self.file_path.get(), "Processed_angles.csv")
            else:
                angles_file_path = os.path.join(self.file_path.get(), "Processed_angles.csv")
            # Read processed data
            data = pd.read_csv(angles_file_path)

            # Get interval selections
            t1_abs_start = self.interval1_rel_start.get() + self.first_abs_stamp.get()  # Convert to absolute time
            t1_abs_end = self.interval1_rel_end.get() + self.first_abs_stamp.get()  # Convert to absolute time
            t2_abs_start = self.interval2_rel_start.get() + self.first_abs_stamp.get()  # Convert to absolute time
            t2_abs_end = self.interval2_rel_end.get() + self.first_abs_stamp.get()  # Convert to absolute time

            # Ensure intervals are ordered
            t1_abs_min, t1_abs_max = min(t1_abs_start, t1_abs_end), max(t1_abs_start, t1_abs_end)
            t2_abs_min, t2_abs_max = min(t2_abs_start, t2_abs_end), max(t2_abs_start, t2_abs_end)

            # Build absolute time mask for both intervals, avoiding duplication if they overlap
            abs_mask1 = (data['timestamp'] >= t1_abs_min) & (data['timestamp'] <= t1_abs_max)
            abs_mask2 = (data['timestamp'] >= t2_abs_min) & (data['timestamp'] <= t2_abs_max)
            combined_abs_mask = abs_mask1 | abs_mask2

            # Remove duplicates if intervals overlap (keep order)
#            filtered_data = data.loc[combined_abs_mask].drop_duplicates(subset='timestamp', keep='first')
            filtered_data = data.loc[combined_abs_mask]

            # Clear previous plots but reuse the existing axes
            self.acq_ax_left.clear()
            self.acq_ax_right.clear()
            self.setup_plots('Data Acquisition', self.acq_ax_left, self.acq_ax_right)

            if not filtered_data.empty:
                # Create a color gradient based on progression (index)
                colors_gradient = np.linspace(0, 1, len(filtered_data))

                # Plot left eye with color gradient
                self.acq_ax_left.scatter(
                    filtered_data['LAz_angle'], filtered_data['LEl_angle'],
                    c=colors_gradient, cmap='viridis', s=1, alpha=0.8, label='Gaze points'
                )
                # Plot right eye with color gradient
                self.acq_ax_right.scatter(
                    filtered_data['RAz_angle'], filtered_data['REl_angle'],
                    c=colors_gradient, cmap='viridis', s=1, alpha=0.8, label='Gaze points'
                )

                # Add legends
                self.acq_ax_left.legend()
                self.acq_ax_right.legend()

                # Remove and recreate colorbar if needed
                if hasattr(self, '_colorbar') and self._colorbar:
                    try:
                        self._colorbar.remove()
                        self._colorbar = None
                    except Exception as e:
                        self.log_status(f"Warning: failed to remove previous colorbar: {e}")

                norm = mcolors.Normalize(vmin=0, vmax=len(filtered_data))
                sm = cm.ScalarMappable(cmap='viridis', norm=norm)
                sm.set_array([])
                self._colorbar = self.acq_fig.colorbar(sm, ax=[self.acq_ax_left, self.acq_ax_right],
                                                    shrink=0.8, pad=0.1, label='Sample Progression')

                # Compute global axis limits
                all_az = np.concatenate([filtered_data['LAz_angle'], filtered_data['RAz_angle']])
                all_el = np.concatenate([filtered_data['LEl_angle'], filtered_data['REl_angle']])
                az_min, az_max = np.min(all_az), np.max(all_az)
                el_min, el_max = np.min(all_el), np.max(all_el)

                # Add a small margin
                az_margin = (az_max - az_min) * 0.05
                el_margin = (el_max - el_min) * 0.05

                self.acq_ax_left.set_xlim(az_min - az_margin, az_max + az_margin)
                self.acq_ax_right.set_xlim(az_min - az_margin, az_max + az_margin)
                self.acq_ax_left.set_ylim(el_min - el_margin, el_max + el_margin)
                self.acq_ax_right.set_ylim(el_min - el_margin, el_max + el_margin)

            self.acq_fig.canvas.draw()
            self.log_status(f"Acquisition plots updated")

            # Update summary
            self.update_acquisition_summary(filtered_data)

        except Exception as e:
            self.log_status(f"Error updating acquisition plots: {str(e)}")

    def update_acquisition_summary(self, data: pd.DataFrame):
        """Update the acquisition summary with statistics."""
        try:
            # CORRECTED: Handle potential timestamp format issues
            if data['timestamp'].dtype == 'object':
                # If timestamps are strings, try to convert
                timestamps = pd.to_numeric(data['timestamp'], errors='coerce')
            else:
                timestamps = data['timestamp']
            
            duration = (timestamps.iloc[-1] - timestamps.iloc[0]) / 1000.0  # Convert ms to seconds if needed
            
            summary = f"""Data Acquisition Summary:

Total samples: {len(data):,}
Duration: {duration:.2f} seconds
Sampling rate: {len(data)/duration:.1f} Hz

Left Eye Statistics:
- Azimuth range: {data['LAz_angle'].min():.2f}° to {data['LAz_angle'].max():.2f}°
- Elevation range: {data['LEl_angle'].min():.2f}° to {data['LEl_angle'].max():.2f}°
- Azimuth std: {data['LAz_angle'].std():.3f}°
- Elevation std: {data['LEl_angle'].std():.3f}°

Right Eye Statistics:
- Azimuth range: {data['RAz_angle'].min():.2f}° to {data['RAz_angle'].max():.2f}°
- Elevation range: {data['REl_angle'].min():.2f}° to {data['REl_angle'].max():.2f}°
- Azimuth std: {data['RAz_angle'].std():.3f}°
- Elevation std: {data['REl_angle'].std():.3f}°
"""
            self.acq_summary_text.config(state='normal')
            self.acq_summary_text.delete(1.0, tk.END)
            self.acq_summary_text.insert(1.0, summary)
            self.acq_summary_text.config(state='disabled')
            
        except Exception as e:
            self.log_status(f"Error updating acquisition summary: {str(e)}")

    def apply_filters_to_data(self):
        try:
            fd = self.filter_designer
            # Check if at least one filter is both designed and active
            no_filter_active = (
                (fd.lpf_filter_coeffs is None or not fd.apply_lpf.get()) and
                (fd.notch_filter_coeffs is None or not fd.apply_notch.get()) and
                (fd.hpf_filter_coeffs is None or not fd.apply_hpf.get())
            )
            if not hasattr(self, 'filter_designer') or no_filter_active:
                messagebox.showwarning("Warning", "Please design and/or activate a filter first!")
                return
                    
            # Get the latest processed data file
            processed_file = os.path.join(self.file_path.get(), "Processed_angles.csv")
            if not os.path.exists(processed_file):
                messagebox.showwarning("Warning", "No processed angle data found!")
                return
                
            # Load data
            data = pd.read_csv(processed_file)
            
            # Apply filter with group delay compensation
            filtered_data = self.filter_designer.apply_filters_compensated(data)
            
            # Save filtered data
            filtered_file = os.path.join(self.file_path.get(), "Filtered_angles.csv")
            filtered_data.to_csv(filtered_file, index=False)
            
            self.log_status(f"Filter applied successfully! Saved to: Filtered_angles.csv")
                
        except Exception as e:
            self.log_status(f"Filter application failed: {str(e)}")
            messagebox.showerror("Error", f"Filter application failed: {str(e)}")

class EOGDataProcessor:
    def __init__(self, calibration_params: dict):
        """Initialize with calibration parameters from your EOG calibration."""
        self.cal_params = calibration_params
        
    def process_csv_chunked(self, input_file: str, output_file: str, 
                          chunk_size: int = 100000) -> None:
        """Process large CSV file in chunks for memory efficiency."""

        # Extract calibration parameters
        laz_params = np.array(self.cal_params['LAz_params'])
        lel_params = np.array(self.cal_params['LEl_params'])  
        raz_params = np.array(self.cal_params['RAz_params'])
        rel_params = np.array(self.cal_params['REl_params'])
        
        # Process file in chunks
        first_chunk = True
        total_rows = 0
        
        for chunk in pd.read_csv(input_file, chunksize=chunk_size): # names=input_cols, chunksize=chunk_size, skiprows=1):
            # Extract voltage data (first 5 voltages for left eye, next 5 for right eye)
            # Define indices based on electrode selection
            left_indices = self.cal_params.get('left_indices', list(range(5)))
            right_indices = self.cal_params.get('right_indices', list(range(5)))

            # Adjust to correct column offsets (Ch0–Ch5 → left: Ch1–Ch5 → indices 1–5, Ch6–Ch11 → right: 7–11)
            left_cols = [1 + i for i in left_indices]     # Ch1–Ch5
            right_cols = [7 + i for i in right_indices]   # Ch7–Ch11

            left_voltages = chunk.iloc[:, left_cols].values
            right_voltages = chunk.iloc[:, right_cols].values

            # Add ones for bias
            left_voltages_with_bias = np.column_stack([left_voltages, np.ones(len(chunk))])
            right_voltages_with_bias = np.column_stack([right_voltages, np.ones(len(chunk))])
            
            # Compute linear combinations (vectorized)
            laz_angles = left_voltages_with_bias @ laz_params
            lel_angles = left_voltages_with_bias @ lel_params
            raz_angles = right_voltages_with_bias @ raz_params
            rel_angles = right_voltages_with_bias @ rel_params
            
            # Create output dataframe
            output_chunk = pd.DataFrame({
                'timestamp': chunk.iloc[:, 0],
                'LAz_angle': laz_angles,
                'LEl_angle': lel_angles, 
                'RAz_angle': raz_angles,
                'REl_angle': rel_angles
            })
            
            # Write to output file
            output_chunk.to_csv(output_file, mode='w' if first_chunk else 'a',
                              header=first_chunk, index=False)
            
            first_chunk = False
            total_rows += len(chunk)
            
            # Progress update
            if total_rows % 500000 == 0:
                print(f"Processed {total_rows:,} rows...")
        
        print(f"Processing complete! Total rows: {total_rows:,}")
        pass

    def process_csv_memory_efficient(self, input_file: str, output_file: str) -> None:
        """Alternative approach using NumPy for maximum speed with reasonable memory usage."""
        print("Loading data...")
        
        # Load full data (we'll subselect columns manually)
        data = np.loadtxt(input_file, delimiter=',', skiprows=1)
        
        print(f"Loaded {len(data):,} rows. Processing...")
        
        # Extract components
        timestamps = data[:, 0]

        # Get selected channel indices (relative to 0-based offsets for file)
        left_indices = self.cal_params.get('left_indices', list(range(5)))
        right_indices = self.cal_params.get('right_indices', list(range(5)))

        # Voltage columns: left = Ch1–Ch5 → cols 1–5, right = Ch7–Ch11 → cols 7–11
        left_cols = [1 + i for i in left_indices]
        right_cols = [7 + i for i in right_indices]

        # Extract only selected channels
        left_voltages = data[:, left_cols]
        right_voltages = data[:, right_cols]
        
        # Add bias column
        left_with_bias = np.column_stack([left_voltages, np.ones(len(data))])
        right_with_bias = np.column_stack([right_voltages, np.ones(len(data))])
        
        # Extract parameters
        laz_params = np.array(self.cal_params['LAz_params'])
        lel_params = np.array(self.cal_params['LEl_params'])
        raz_params = np.array(self.cal_params['RAz_params'])
        rel_params = np.array(self.cal_params['REl_params'])
        
        # Compute all transformations at once (very fast)
        laz_angles = left_with_bias @ laz_params
        lel_angles = left_with_bias @ lel_params
        raz_angles = right_with_bias @ raz_params
        rel_angles = right_with_bias @ rel_params
        
        print("Saving results...")
        
        # Combine results
        output_data = np.column_stack([timestamps, laz_angles, lel_angles, raz_angles, rel_angles])
        
        # Save with headers
        header = 'timestamp,LAz_angle,LEl_angle,RAz_angle,REl_angle'
        np.savetxt(output_file, output_data, delimiter=',', header=header, comments='')
        
        print(f"Processing complete! Output saved to {output_file}")
        pass

class FiltersDesigner:
    def __init__(self, parent_frame, main_app):
        self.parent = parent_frame
        self.main_app = main_app  # Store reference to main AEOG_GUI instance
        self.lpf_filter_coeffs = None
        self.notch_filter_coeffs = None
        self.hpf_filter_coeffs = None
        self.lpf_passband_freq = tk.DoubleVar(value=30.0)
        self.stopband_freq = tk.DoubleVar(value=50.0)
        self.passband_ripple = tk.DoubleVar(value=0.1)
        self.stopband_atten = tk.DoubleVar(value=60.0)
        self.lpf_group_delay = 0
        self.notch_min_freq = tk.DoubleVar(value=48.0)
        self.notch_max_freq = tk.DoubleVar(value=52.0)
        self.notch_transition_bw = tk.DoubleVar(value=0.8)
        self.notch_atten = tk.DoubleVar(value=50.0)
        self.notch_passband_ripple = tk.DoubleVar(value=0.1)
        self.notch_group_delay = 0
        self.hpf_group_delay = 0
        self.apply_lpf = tk.BooleanVar(value=True)
        self.apply_notch = tk.BooleanVar(value=True)
        self.apply_hpf = tk.BooleanVar(value=True)
        self.sampling_rate = tk.DoubleVar(value=1000.0) # Shared sampling rate variable for all filters
        self.setup_filters_ui()

    def setup_filters_ui(self):

        # --- Parent frame for all three filters ---
        filters_row = ttk.Frame(self.parent)
        filters_row.pack(fill='both', expand=True)

        for i in range(3):
            filters_row.columnconfigure(i, weight=1)
        filters_row.rowconfigure(0, weight=1)

        # --- Low-Pass FIR Filter Design (left third) ---
        lpf_frame = ttk.LabelFrame(filters_row, text="Low-Pass FIR Filter Design (Parks-McClellan)", padding=6)
        lpf_frame.grid(row=0, column=0, sticky='nsew', padx=(0, 4), pady=2)
        lpf_frame.columnconfigure(0, weight=1)

        # Sampling Rate (shared)
        lpf_sr_row = ttk.Frame(lpf_frame)
        lpf_sr_row.pack(fill='x', pady=(2, 0))
        ttk.Label(lpf_sr_row, text="Sampling Rate (Hz):").pack(side='left', padx=(0, 2))
        lpf_sr_entry = ttk.Entry(lpf_sr_row, textvariable=self.sampling_rate, width=10)
        lpf_sr_entry.pack(side='left')

        params_frame = ttk.Frame(lpf_frame)
        params_frame.pack(fill='x', pady=2)
        labels_vars = [
            ("Passband Cutoff (Hz):", self.lpf_passband_freq),
            ("Stopband Cutoff (Hz):", self.stopband_freq),
            ("Passband Ripple (dB):", self.passband_ripple),
            ("Stopband Attenuation (dB):", self.stopband_atten)
        ]
        for i, (label, var) in enumerate(labels_vars):
            ttk.Label(params_frame, text=label).grid(row=i, column=0, sticky='w', padx=2, pady=1)
            ttk.Entry(params_frame, textvariable=var, width=10).grid(row=i, column=1, padx=2, pady=1)
        filter_design_button_frame = ttk.Frame(lpf_frame)
        filter_design_button_frame.pack(fill='x', pady=2)
        ttk.Button(filter_design_button_frame, text="Design and Save Low-Pass Filter",
                command=self.design_and_save_lpf_filter).pack(side='left', padx=2)
        ttk.Checkbutton(
            lpf_frame,
            text="Apply Low-Pass Filter",
            variable=self.apply_lpf,
            command=self.on_apply_lpf_changed
        ).pack(anchor='w')
        self.lpf_info = tk.Text(lpf_frame, height=6, wrap='word', state='disabled')
        self.lpf_info.pack(fill='x', pady=2)
        self.lpf_response_frame = ttk.LabelFrame(lpf_frame, text="Filter Response", padding=2)
        # Don't pack initially

        # --- Notch Filter (center third) ---
        notch_frame = ttk.LabelFrame(filters_row, text="Notch Filter", padding=6)
        notch_frame.grid(row=0, column=1, sticky='nsew', padx=2, pady=2)
        notch_frame.columnconfigure(0, weight=1)

        # Sampling Rate (shared)
        notch_sr_row = ttk.Frame(notch_frame)
        notch_sr_row.pack(fill='x', pady=(2, 0))
        ttk.Label(notch_sr_row, text="Sampling Rate (Hz):").pack(side='left', padx=(0, 2))
        notch_sr_entry = ttk.Entry(notch_sr_row, textvariable=self.sampling_rate, width=10)
        notch_sr_entry.pack(side='left')

        notch_params_frame = ttk.Frame(notch_frame)
        notch_params_frame.pack(fill='x', pady=2)
        notch_labels_vars = [
            ("Min Cutoff (Hz):", self.notch_min_freq),
            ("Max Cutoff (Hz):", self.notch_max_freq),
            ("Transition BW (Hz):", self.notch_transition_bw),
            ("Attenuation (dB):", self.notch_atten),
            ("Passband Ripple (dB):", self.notch_passband_ripple)
        ]
        for i, (label, var) in enumerate(notch_labels_vars):
            ttk.Label(notch_params_frame, text=label).grid(row=i, column=0, sticky='w', padx=2, pady=1)
            ttk.Entry(notch_params_frame, textvariable=var, width=10).grid(row=i, column=1, padx=2, pady=1)
        notch_button_frame = ttk.Frame(notch_frame)
        notch_button_frame.pack(fill='x', pady=2)
        ttk.Button(notch_button_frame, text="Design and Save Notch Filter",
                    command=self.design_and_save_notch_filter).pack(side='left', padx=2)
        ttk.Checkbutton(
            notch_frame,
            text="Apply Notch Filter",
            variable=self.apply_notch,
            command=self.on_apply_notch_changed
        ).pack(anchor='w')
        self.notch_info = tk.Text(notch_frame, height=6, wrap='word', state='disabled')
        self.notch_info.pack(fill='x', pady=2)
        # --- Add a dedicated response frame for the notch filter ---
        self.notch_response_frame = ttk.LabelFrame(notch_frame, text="Notch Filter Response", padding=2)
        # Don't pack initially

        # --- High-Pass FIR Filter (right third) ---
        hpf_frame = ttk.LabelFrame(filters_row, text="High-Pass FIR Filter", padding=6)
        hpf_frame.grid(row=0, column=2, sticky='nsew', padx=(4, 0), pady=2)
        hpf_frame.columnconfigure(0, weight=1)

        # Sampling Rate (shared)
        hpf_sr_row = ttk.Frame(hpf_frame)
        hpf_sr_row.pack(fill='x', pady=(2, 0))
        ttk.Label(hpf_sr_row, text="Sampling Rate (Hz):").pack(side='left', padx=(0, 2))
        hpf_sr_entry = ttk.Entry(hpf_sr_row, textvariable=self.sampling_rate, width=10)
        hpf_sr_entry.pack(side='left')

        hpf_params_frame = ttk.Frame(hpf_frame)
        hpf_params_frame.pack(fill='x', pady=2)
        self.hpf_passband_freq = tk.DoubleVar(value=1.0)
        self.hpf_stopband_freq = tk.DoubleVar(value=0.5)
        self.hpf_passband_ripple = tk.DoubleVar(value=0.1)
        self.hpf_stopband_atten = tk.DoubleVar(value=60.0)
        hpf_labels_vars = [
            ("Passband Cutoff (Hz):", self.hpf_passband_freq),
            ("Stopband Cutoff (Hz):", self.hpf_stopband_freq),
            ("Passband Ripple (dB):", self.hpf_passband_ripple),
            ("Stopband Attenuation (dB):", self.hpf_stopband_atten)
        ]
        for i, (label, var) in enumerate(hpf_labels_vars):
            ttk.Label(hpf_params_frame, text=label).grid(row=i, column=0, sticky='w', padx=2, pady=1)
            ttk.Entry(hpf_params_frame, textvariable=var, width=10).grid(row=i, column=1, padx=2, pady=1)
        hpf_button_frame = ttk.Frame(hpf_frame)
        hpf_button_frame.pack(fill='x', pady=2)
        ttk.Button(hpf_button_frame, text="Design and Save High-Pass Filter",
                   command=self.design_and_save_hpf_filter).pack(side='left', padx=2)
        ttk.Checkbutton(
            hpf_frame,
            text="Apply High-Pass Filter",
            variable=self.apply_hpf,
            command=self.on_apply_hpf_changed
        ).pack(anchor='w')

        self.hpf_info = tk.Text(hpf_frame, height=6, wrap='word', state='disabled')
        self.hpf_info.pack(fill='x', pady=2)
        self.hpf_response_frame = ttk.LabelFrame(hpf_frame, text="High-Pass Filter Response", padding=2)
        # Don't pack initially

    def design_and_save_lpf_filter(self):
        """Design LPF filter, save it and apply it to data if enabled."""
        try:
            # Get parameters
            fs = self.sampling_rate.get()
            fp = self.lpf_passband_freq.get()
            fa = self.stopband_freq.get() 
            rp = self.passband_ripple.get()
            rs = self.stopband_atten.get()
            
            # Validate inputs
            if fp >= fa:
                messagebox.showerror("Error", "Passband frequency must be < stopband frequency")
                return
            if fp >= fs/2 or fa >= fs/2:
                messagebox.showerror("Error", "Cutoff frequencies must be < Nyquist frequency")
                return
                
            # Convert to normalized frequencies
            wp = fp / (fs/2)
            wa = fa / (fs/2)
            
            # Convert dB to linear scale
            delta_p = (10**(rp/20) - 1) / (10**(rp/20) + 1)  # Passband ripple
            delta_s = 10**(-rs/20)                            # Stopband ripple
            
            # Estimate filter order using Kaiser formula
            df = (fa - fp) / fs  # Transition width
            A = -20 * np.log10(delta_s)
            
            if A > 50:
                beta = 0.1102 * (A - 8.7)
            elif A >= 21:
                beta = 0.5842 * (A - 21)**0.4 + 0.07886 * (A - 21)
            else:
                beta = 0
                
            N_est = int(np.ceil((A - 8) / (2.285 * 2 * np.pi * df))) + 1
            
            # Ensure odd order for linear phase Type I filter
            if N_est % 2 == 0:
                N_est += 1
                
            # Limit to maximum 10001 taps
            N_est = min(N_est, 10001)
            
            # Design using Parks-McClellan (Remez exchange)
            # Try different orders around estimate
            for N in range(max(3, N_est-10), min(10002, N_est+20), 2):  # Only odd orders
                try:
                    bands = [0, wp, wa, 1]
                    desired = [1, 0]
                    weight = [1/delta_p, 1/delta_s]  # Weight inversely proportional to allowed error
                    
                    h = signal.remez(N, bands, desired, weight=weight, fs=2)
                    
                    # Check if specifications are met
                    w, H = signal.freqz(h, worN=8192, fs=fs)
                    H_mag = np.abs(H)
                    H_db = 20 * np.log10(H_mag + 1e-12)
                    
                    # Find passband and stopband regions
                    pb_mask = w <= fp
                    sb_mask = w >= fa
                    
                    if np.any(pb_mask) and np.any(sb_mask):
                        pb_ripple = np.max(H_db[pb_mask]) - np.min(H_db[pb_mask])
                        sb_atten = -np.max(H_db[sb_mask])
                        
                        if pb_ripple <= rp and sb_atten >= rs:
                            self.lpf_filter_coeffs = h
                            self.lpf_group_delay = (len(h) - 1) // 2  # Group delay in samples
                            break
                            
                except Exception:
                    continue
            else:
                # If no solution found, use best effort
                N = min(N_est, 1001)
                if N % 2 == 0:
                    N += 1
                bands = [0, wp, wa, 1]
                desired = [1, 0]  
                weight = [1, 1]
                self.lpf_filter_coeffs = signal.remez(N, bands, desired, weight=weight, fs=2)
                self.lpf_group_delay = (len(self.lpf_filter_coeffs) - 1) // 2
            
            # Update info display
            self.update_lpf_filter_info()
            
            # After successful filter design, show response and save
            self.update_lpf_filter_info()
            self.show_lpf_response_inline()
            self.save_lpf_filter_data()

            # Auto-apply to data if data is available and if LPF filter is enabled
            processed_file = os.path.join(self.main_app.file_path.get(), "Processed_angles.csv")
            if os.path.exists(processed_file) & self.apply_lpf.get():
                self.main_app.apply_filters_to_data()
                # Update display if filtered display is enabled
                if self.main_app.display_filtered_angles.get():
                    self.main_app.update_acquisition_plots()
            
        except Exception as e:
            messagebox.showerror("Error", f"Low-Pass Filter design failed: {str(e)}")
            
    def update_lpf_filter_info(self):
        if self.lpf_filter_coeffs is None:
            return
            
        # Analyze filter performance
        fs = self.sampling_rate.get()
        w, H = signal.freqz(self.lpf_filter_coeffs, worN=8192, fs=fs)
        H_mag = np.abs(H)
        H_db = 20 * np.log10(H_mag + 1e-12)
        
        fp = self.lpf_passband_freq.get()
        fa = self.stopband_freq.get()
        
        pb_mask = w <= fp
        sb_mask = w >= fa
        
        if np.any(pb_mask) and np.any(sb_mask):
            pb_ripple = np.max(H_db[pb_mask]) - np.min(H_db[pb_mask])
            sb_atten = -np.max(H_db[sb_mask])
        else:
            pb_ripple = sb_atten = 0
            
        info_text = f"""Filter Design Results:
Taps: {len(self.lpf_filter_coeffs)}
Group Delay: {self.lpf_group_delay} samples ({self.lpf_group_delay/fs*1000:.2f} ms)
Achieved Passband Ripple: {pb_ripple:.3f} dB
Achieved Stopband Attenuation: {sb_atten:.1f} dB"""
        
        self.lpf_info.config(state='normal')
        self.lpf_info.delete(1.0, tk.END)
        self.lpf_info.insert(1.0, info_text)
        self.lpf_info.config(state='disabled')
        
    def show_lpf_response_inline(self):
        """Show LPF filter response in the GUI tab."""
        if self.lpf_filter_coeffs is None:
            return

        self.lpf_response_frame.pack(fill='both', expand=True, pady=5)

        for widget in self.lpf_response_frame.winfo_children():
            widget.destroy()

        fig = Figure(figsize=(4, 2.5))

        # Frequency response
        ax1 = fig.add_subplot(121)
        fs = self.sampling_rate.get()
        w, H = signal.freqz(self.lpf_filter_coeffs, worN=8192, fs=fs)
        ax1.plot(w, 20*np.log10(np.abs(H)))
        ax1.axvline(self.lpf_passband_freq.get(), color='r', linestyle='--', alpha=0.7)
        ax1.axvline(self.stopband_freq.get(), color='r', linestyle='--', alpha=0.7)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Magnitude (dB)')
        ax1.set_title('Magnitude Response')
        ax1.grid(True, alpha=0.3)

        # Impulse response
        ax2 = fig.add_subplot(122)
        ax2.stem(range(len(self.lpf_filter_coeffs)), self.lpf_filter_coeffs, basefmt=' ')
        ax2.set_xlabel('Sample')
        ax2.set_ylabel('Amplitude')
        ax2.set_title('Impulse Response')
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, self.lpf_response_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)

    def save_lpf_filter_data(self):
        """Save filter parameters and design to file."""
        try:
            # Get parent AEOG_GUI instance to access file path
            main_app = self.main_app
            if hasattr(main_app, 'file_path'):
                file_path = main_app.file_path.get()
            else:
                file_path = os.getcwd()
                
            filter_file = os.path.join(file_path, "lpf_filter_data.json")
            
            filter_data = {
                'sampling_rate': self.sampling_rate.get(),
                'lpf_passband_freq': self.lpf_passband_freq.get(),
                'stopband_freq': self.stopband_freq.get(),
                'passband_ripple': self.passband_ripple.get(),
                'stopband_atten': self.stopband_atten.get(),
                'apply_lpf': self.apply_lpf.get()  # <-- Save checkbox state
            }
            
            # Add filter coefficients if filter exists
            if self.lpf_filter_coeffs is not None:
                filter_data['lpf_filter_coeffs'] = self.lpf_filter_coeffs.tolist()
                filter_data['lpf_group_delay'] = self.lpf_group_delay
            else:
                filter_data['lpf_filter_coeffs'] = None
                filter_data['lpf_group_delay'] = 0
                
            with open(filter_file, 'w') as f:
                json.dump(filter_data, f, indent=2)
                
        except Exception as e:
            print(f"Error saving filter data: {e}")

    def load_lpf_filter_data(self, file_path):
        """Load Low-Pass filter parameters from file and update GUI."""
        try:
            filter_file = os.path.join(file_path, "lpf_filter_data.json")
            if not os.path.exists(filter_file):
                return False
                
            with open(filter_file, 'r') as f:
                filter_data = json.load(f)
                
            # Restore parameters
            self.sampling_rate.set(filter_data['sampling_rate'])
            self.lpf_passband_freq.set(filter_data['lpf_passband_freq'])
            self.stopband_freq.set(filter_data['stopband_freq'])
            self.passband_ripple.set(filter_data['passband_ripple'])
            self.stopband_atten.set(filter_data['stopband_atten'])
            # Restore checkbox state (default to True if not present for backward compatibility)
            self.apply_lpf.set(filter_data.get('apply_lpf', True))
            # Restore filter coefficients if they exist
            if filter_data.get('lpf_filter_coeffs') is not None:
                self.lpf_filter_coeffs = np.array(filter_data['lpf_filter_coeffs'])
                self.lpf_group_delay = filter_data['lpf_group_delay']
            else:
                self.lpf_filter_coeffs = None
                self.lpf_group_delay = 0
                
            # Update display
            self.update_lpf_filter_info()
            if self.lpf_filter_coeffs is not None:
                self.show_lpf_response_inline()
            
            self.main_app.log_status("Low-Pass Filter data loaded automatically!")

            return True
            
        except Exception as e:
            print(f"Error loading Low-Pass Filter data: {e}")
            return False

    def on_apply_lpf_changed(self):
        """Callback for when the apply LPF checkbox is toggled."""
        self.save_lpf_filter_data()  # Save state when toggled
        # Auto-apply to data if data is available
        processed_file = os.path.join(self.main_app.file_path.get(), "Processed_angles.csv")
        if os.path.exists(processed_file):
            self.main_app.apply_filters_to_data()
            # Update display if filtered display is enabled
            if self.main_app.display_filtered_angles.get():
                self.main_app.update_acquisition_plots()

    def design_and_save_notch_filter(self):
        """Design a linear phase FIR notch filter and show its response, iteratively increasing order if needed."""
        try:
            fs = self.sampling_rate.get()
            f1 = self.notch_min_freq.get()
            f2 = self.notch_max_freq.get()
            trans_bw = self.notch_transition_bw.get()
            attn = self.notch_atten.get()
            rp = self.notch_passband_ripple.get()

            # Validate
            if f1 <= 0 or f2 <= f1 or f2 >= fs/2:
                messagebox.showerror("Error", "Min/Max cutoff must be >0, min < max, and max < Nyquist.")
                return
            if trans_bw <= 0:
                messagebox.showerror("Error", "Transition bandwidth must be > 0.")
                return

            bands = [
                0,
                max(0, f1 - trans_bw),
                f1,
                f2,
                min(fs/2, f2 + trans_bw),
                fs/2
            ]
            desired = [1, 0, 1]
            delta_p = (10**(rp/20) - 1) / (10**(rp/20) + 1)
            weight = [1/delta_p, 10**(attn/20), 1/delta_p]
            bands_norm = [b / (fs/2) for b in bands]

            # Start with a reasonable order and increase if needed
            N_min = 301
            N_max = 10001
            N_step = 50
            for N in range(N_min, N_max+1, N_step):
                if N % 2 == 0:
                    N += 1
                try:
                    coeffs = signal.remez(N, bands_norm, desired, weight=weight, fs=2)
                    w, H = signal.freqz(coeffs, worN=8192, fs=fs)
                    H_db = 20 * np.log10(np.abs(H) + 1e-12)
                    pb1_mask = w <= max(0, f1 - trans_bw)
                    pb2_mask = w >= min(fs/2, f2 + trans_bw)
                    sb_mask = (w >= f1) & (w <= f2)
                    min_attn_sb = -np.max(H_db[sb_mask]) if np.any(sb_mask) else float('nan')
                    pb1_ripple = np.max(H_db[pb1_mask]) - np.min(H_db[pb1_mask]) if np.any(pb1_mask) else float('nan')
                    pb2_ripple = np.max(H_db[pb2_mask]) - np.min(H_db[pb2_mask]) if np.any(pb2_mask) else float('nan')
                    # Check if specs are met
                    if (min_attn_sb >= attn) and (pb1_ripple <= rp) and (pb2_ripple <= rp):
                        self.notch_filter_coeffs = coeffs
                        self.notch_group_delay = (len(coeffs) - 1) // 2
                        break
                except Exception:
                    continue
            else:
                # If not met, use last attempt
                self.notch_filter_coeffs = coeffs
                self.notch_group_delay = (len(coeffs) - 1) // 2

            self.update_notch_filter_info()
            self.show_notch_response_inline()
            self.save_notch_filter_data()

            # Auto-apply to data if data is available and if notch filter is enabled
            processed_file = os.path.join(self.main_app.file_path.get(), "Processed_angles.csv")
            if os.path.exists(processed_file) and self.apply_notch.get():
                self.main_app.apply_filters_to_data()
                if self.main_app.display_filtered_angles.get():
                    self.main_app.update_acquisition_plots()

        except Exception as e:
            messagebox.showerror("Error", f"Notch filter design failed: {str(e)}")

    def update_notch_filter_info(self):
        if not hasattr(self, 'notch_filter_coeffs') or self.notch_filter_coeffs is None:
            return
        fs = self.sampling_rate.get()
        w, H = signal.freqz(self.notch_filter_coeffs, worN=8192, fs=fs)
        H_db = 20 * np.log10(np.abs(H) + 1e-12)
        f1 = self.notch_min_freq.get()
        f2 = self.notch_max_freq.get()
        trans_bw = self.notch_transition_bw.get()
        # Passbands: below notch and above notch
        pb1_mask = w <= max(0, f1 - trans_bw)
        pb2_mask = w >= min(fs/2, f2 + trans_bw)
        # Stopband: notch region
        sb_mask = (w >= f1) & (w <= f2)
        # Minimum attenuation in stopband (dB)
        min_attn_sb = -np.max(H_db[sb_mask]) if np.any(sb_mask) else float('nan')
        # Attenuation at notch center
        f0 = (f1 + f2) / 2
        idx = np.argmin(np.abs(w - f0))
        notch_attn = -H_db[idx]
        # Passband ripple calculations
        pb1_ripple = np.max(H_db[pb1_mask]) - np.min(H_db[pb1_mask]) if np.any(pb1_mask) else float('nan')
        pb2_ripple = np.max(H_db[pb2_mask]) - np.min(H_db[pb2_mask]) if np.any(pb2_mask) else float('nan')
        info_text = f"""Notch Filter Design Results:
    Taps: {len(self.notch_filter_coeffs)}
    Group Delay: {self.notch_group_delay} samples ({self.notch_group_delay/fs*1000:.2f} ms)
    Attenuation at {f0:.1f} Hz: {notch_attn:.1f} dB
    Minimum attenuation in stopband: {min_attn_sb:.2f} dB
    Max ripple below / above notch: {pb1_ripple:.3f} / {pb2_ripple:.3f} dB"""
        self.notch_info.config(state='normal')
        self.notch_info.delete(1.0, tk.END)
        self.notch_info.insert(1.0, info_text)
        self.notch_info.config(state='disabled')

    def show_notch_response_inline(self):
        """Show notch filter response in the central GUI frame."""
        if not hasattr(self, 'notch_filter_coeffs') or self.notch_filter_coeffs is None:
            return

        self.notch_response_frame.pack(fill='both', expand=True, pady=5)

        for widget in self.notch_response_frame.winfo_children():
            widget.destroy()

        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        fig = Figure(figsize=(4, 2.5))  # Compact

        # Frequency response
        ax1 = fig.add_subplot(121)
        fs = self.sampling_rate.get()
        w, H = signal.freqz(self.notch_filter_coeffs, worN=8192, fs=fs)
        ax1.plot(w, 20*np.log10(np.abs(H)))
#        ax1.axvline(self.notch_min_freq.get(), color='r', linestyle='--', alpha=0.7)
#        ax1.axvline(self.notch_max_freq.get(), color='r', linestyle='--', alpha=0.7)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Magnitude (dB)')
        ax1.set_title('Magnitude Response')
        ax1.grid(True, alpha=0.3)

        # Impulse response
        ax2 = fig.add_subplot(122)
        ax2.stem(range(len(self.notch_filter_coeffs)), self.notch_filter_coeffs, basefmt=' ')
        ax2.set_xlabel('Sample')
        ax2.set_ylabel('Amplitude')
        ax2.set_title('Impulse Response')
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, self.notch_response_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)

    def save_notch_filter_data(self):
        """Save notch filter parameters and design to file."""
        try:
            # Get parent AEOG_GUI instance to access file path
            main_app = self.main_app
            if hasattr(main_app, 'file_path'):
                file_path = main_app.file_path.get()
            else:
                file_path = os.getcwd()

            filter_file = os.path.join(file_path, "notch_filter_data.json")

            filter_data = {
                'sampling_rate': self.sampling_rate.get(),
                'notch_min_freq': self.notch_min_freq.get(),
                'notch_max_freq': self.notch_max_freq.get(),
                'notch_transition_bw': self.notch_transition_bw.get(),
                'notch_atten': self.notch_atten.get(),
                'notch_passband_ripple': self.notch_passband_ripple.get(),
                'apply_notch': self.apply_notch.get()  # <-- Save checkbox state
            }

            # Add filter coefficients if filter exists
            if self.notch_filter_coeffs is not None:
                filter_data['notch_filter_coeffs'] = self.notch_filter_coeffs.tolist()
                filter_data['notch_group_delay'] = self.notch_group_delay
            else:
                filter_data['notch_filter_coeffs'] = None
                filter_data['notch_group_delay'] = 0

            with open(filter_file, 'w') as f:
                json.dump(filter_data, f, indent=2)

        except Exception as e:
            print(f"Error saving notch filter data: {e}")

    def load_notch_filter_data(self, file_path):
        """Load Notch filter parameters from file and update GUI."""
        try:
            filter_file = os.path.join(file_path, "notch_filter_data.json")
            if not os.path.exists(filter_file):
                return False

            with open(filter_file, 'r') as f:
                filter_data = json.load(f)

            # Restore parameters
            self.sampling_rate.set(filter_data['sampling_rate'])
            self.notch_min_freq.set(filter_data['notch_min_freq'])
            self.notch_max_freq.set(filter_data['notch_max_freq'])
            self.notch_transition_bw.set(filter_data['notch_transition_bw'])
            self.notch_atten.set(filter_data['notch_atten'])
            self.notch_passband_ripple.set(filter_data['notch_passband_ripple'])
            # Restore checkbox state (default to True if not present for backward compatibility)
            self.apply_notch.set(filter_data.get('apply_notch', True))
            self.parent.update()  # <-- Force UI update if needed

            # Restore filter coefficients if they exist
            if filter_data.get('notch_filter_coeffs') is not None:
                self.notch_filter_coeffs = np.array(filter_data['notch_filter_coeffs'])
                self.notch_group_delay = filter_data.get('notch_group_delay', 0)
            else:
                self.notch_filter_coeffs = None
                self.notch_group_delay = 0

            # Update display
            self.update_notch_filter_info()
            if self.notch_filter_coeffs is not None:
                self.show_notch_response_inline()

            self.main_app.log_status("Notch Filter data loaded automatically!")

            return True

        except Exception as e:
            print(f"Error loading Notch Filter data: {e}")
            return False

    def on_apply_notch_changed(self):
        """Callback for when the apply notch filter checkbox is toggled."""
        self.save_notch_filter_data()  # Save state when toggled
        # Auto-apply to data if data is available
        processed_file = os.path.join(self.main_app.file_path.get(), "Processed_angles.csv")
        if os.path.exists(processed_file):
            self.main_app.apply_filters_to_data()
            # Update display if filtered display is enabled
            if self.main_app.display_filtered_angles.get():
                self.main_app.update_acquisition_plots()

    def design_and_save_hpf_filter(self):
        """Design a high-pass FIR filter and show its response."""
        try:
            fs = self.sampling_rate.get()
            fp = self.hpf_passband_freq.get()
            fa = self.hpf_stopband_freq.get()
            rp = self.hpf_passband_ripple.get()
            rs = self.hpf_stopband_atten.get()

            # Validate
            if fa >= fp:
                messagebox.showerror("Error", "Stopband frequency must be < passband frequency")
                return
            if fp >= fs/2 or fa >= fs/2:
                messagebox.showerror("Error", "Cutoff frequencies must be < Nyquist frequency")
                return

            # Convert to normalized frequencies
            wp = fp / (fs/2)
            wa = fa / (fs/2)

            # Convert dB to linear scale
            delta_p = (10**(rp/20) - 1) / (10**(rp/20) + 1)  # Passband ripple
            delta_s = 10**(-rs/20)                            # Stopband ripple

            # Estimate filter order using Kaiser formula
            df = (fp - fa) / fs  # Transition width
            A = -20 * np.log10(delta_s)

            if A > 50:
                beta = 0.1102 * (A - 8.7)
            elif A >= 21:
                beta = 0.5842 * (A - 21)**0.4 + 0.07886 * (A - 21)
            else:
                beta = 0

            N_est = int(np.ceil((A - 8) / (2.285 * 2 * np.pi * df))) + 1

            # Ensure odd order for linear phase Type I filter
            if N_est % 2 == 0:
                N_est += 1

            # Limit to maximum 50001 taps
            N_est = min(N_est, 50001)

            # Design using Parks-McClellan (Remez exchange)
            for N in range(max(3, N_est-10), min(1002, N_est+20), 2):  # Only odd orders
                try:
                    bands = [0, wa, wp, 1]
                    desired = [0, 1]
                    weight = [1/delta_s, 1/delta_p]  # Weight inversely proportional to allowed error

                    h = signal.remez(N, bands, desired, weight=weight, fs=2)

                    # Check if specifications are met
                    w, H = signal.freqz(h, worN=8192, fs=fs)
                    H_mag = np.abs(H)
                    H_db = 20 * np.log10(H_mag + 1e-12)

                    # Find passband and stopband regions
                    pb_mask = w >= fp
                    sb_mask = w <= fa

                    if np.any(pb_mask) and np.any(sb_mask):
                        pb_ripple = np.max(H_db[pb_mask]) - np.min(H_db[pb_mask])
                        sb_atten = -np.max(H_db[sb_mask])

                        if pb_ripple <= rp and sb_atten >= rs:
                            self.hpf_filter_coeffs = h
                            self.hpf_group_delay = (len(h) - 1) // 2
                            break

                except Exception:
                    continue
            else:
                # If no solution found, use best effort
                N = min(N_est, 50001)
                if N % 2 == 0:
                    N += 1
                bands = [0, wa, wp, 1]
                desired = [0, 1]
                weight = [1, 1]
                self.hpf_filter_coeffs = signal.remez(N, bands, desired, weight=weight, fs=2)
                self.hpf_group_delay = (len(self.hpf_filter_coeffs) - 1) // 2

            # Update info display
            self.update_hpf_filter_info()
            self.show_hpf_response_inline()
            self.save_hpf_filter_data()

            # Auto-apply to data if data is available and if HPF filter is enabled
            processed_file = os.path.join(self.main_app.file_path.get(), "Processed_angles.csv")
            if os.path.exists(processed_file) & self.apply_hpf.get():
                self.main_app.apply_filters_to_data()
                # Update display if filtered display is enabled
                if self.main_app.display_filtered_angles.get():
                    self.main_app.update_acquisition_plots()

        except Exception as e:
            messagebox.showerror("Error", f"High-Pass filter design failed: {str(e)}")

    def update_hpf_filter_info(self):
        if not hasattr(self, 'hpf_filter_coeffs') or self.hpf_filter_coeffs is None:
            return
        fs = self.sampling_rate.get()
        w, H = signal.freqz(self.hpf_filter_coeffs, worN=8192, fs=fs)
        H_db = 20 * np.log10(np.abs(H) + 1e-12)
        fp = self.hpf_passband_freq.get()
        fa = self.hpf_stopband_freq.get()
        # Passband: above passband freq
        pb_mask = w >= fp
        # Stopband: below stopband freq
        sb_mask = w <= fa
        # Passband ripple
        pb_ripple = np.max(H_db[pb_mask]) - np.min(H_db[pb_mask]) if np.any(pb_mask) else float('nan')
        # Stopband attenuation
        sb_atten = -np.max(H_db[sb_mask]) if np.any(sb_mask) else float('nan')
        info_text = f"""High-Pass Filter Design Results:
Taps: {len(self.hpf_filter_coeffs)}
Group Delay: {self.hpf_group_delay} samples ({self.hpf_group_delay/fs*1000:.2f} ms)
Achieved Passband Ripple: {pb_ripple:.3f} dB
Achieved Stopband Attenuation: {sb_atten:.1f} dB"""
        self.hpf_info.config(state='normal')
        self.hpf_info.delete(1.0, tk.END)
        self.hpf_info.insert(1.0, info_text)
        self.hpf_info.config(state='disabled')

    def show_hpf_response_inline(self):
        """Show high-pass filter response in the right GUI frame."""
        if not hasattr(self, 'hpf_filter_coeffs') or self.hpf_filter_coeffs is None:
            return

        self.hpf_response_frame.pack(fill='both', expand=True, pady=5)

        for widget in self.hpf_response_frame.winfo_children():
            widget.destroy()

        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        fig = Figure(figsize=(4, 2.5))  # Compact

        # Frequency response
        ax1 = fig.add_subplot(121)
        fs = self.sampling_rate.get()
        w, H = signal.freqz(self.hpf_filter_coeffs, worN=8192, fs=fs)
        ax1.plot(w, 20*np.log10(np.abs(H)))
        ax1.axvline(self.hpf_passband_freq.get(), color='r', linestyle='--', alpha=0.7)
        ax1.axvline(self.hpf_stopband_freq.get(), color='r', linestyle='--', alpha=0.7)
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Magnitude (dB)')
        ax1.set_title('Magnitude Response')
        ax1.grid(True, alpha=0.3)

        # Impulse response
        ax2 = fig.add_subplot(122)
        ax2.stem(range(len(self.hpf_filter_coeffs)), self.hpf_filter_coeffs, basefmt=' ')
        ax2.set_xlabel('Sample')
        ax2.set_ylabel('Amplitude')
        ax2.set_title('Impulse Response')
        ax2.grid(True, alpha=0.3)

        fig.tight_layout()

        canvas = FigureCanvasTkAgg(fig, self.hpf_response_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill='both', expand=True)

    def save_hpf_filter_data(self):
        """Save high-pass filter parameters and design to file."""
        try:
            # Get parent AEOG_GUI instance to access file path
            main_app = self.main_app
            if hasattr(main_app, 'file_path'):
                file_path = main_app.file_path.get()
            else:
                file_path = os.getcwd()

            filter_file = os.path.join(file_path, "hpf_filter_data.json")

            filter_data = {
                'sampling_rate': self.sampling_rate.get(),
                'hpf_passband_freq': self.hpf_passband_freq.get(),
                'hpf_stopband_freq': self.hpf_stopband_freq.get(),
                'hpf_passband_ripple': self.hpf_passband_ripple.get(),
                'hpf_stopband_atten': self.hpf_stopband_atten.get(),
                'apply_hpf': self.apply_hpf.get()  # <-- Save checkbox state
            }

            # Add filter coefficients if filter exists
            if self.hpf_filter_coeffs is not None:
                filter_data['hpf_filter_coeffs'] = self.hpf_filter_coeffs.tolist()
                filter_data['hpf_group_delay'] = self.hpf_group_delay
            else:
                filter_data['hpf_filter_coeffs'] = None
                filter_data['hpf_group_delay'] = 0

            with open(filter_file, 'w') as f:
                json.dump(filter_data, f, indent=2)

        except Exception as e:
            print(f"Error saving high-pass filter data: {e}")

    def load_hpf_filter_data(self, file_path):
        """Load High-Pass filter parameters from file and update GUI."""
        try:
            filter_file = os.path.join(file_path, "hpf_filter_data.json")
            if not os.path.exists(filter_file):
                return False

            with open(filter_file, 'r') as f:
                filter_data = json.load(f)

            # Restore parameters
            self.sampling_rate.set(filter_data['sampling_rate'])
            self.hpf_passband_freq.set(filter_data['hpf_passband_freq'])
            self.hpf_stopband_freq.set(filter_data['hpf_stopband_freq'])
            self.hpf_passband_ripple.set(filter_data['hpf_passband_ripple'])
            self.hpf_stopband_atten.set(filter_data['hpf_stopband_atten'])
            # Restore checkbox state (default to True if not present for backward compatibility)
            self.apply_hpf.set(filter_data.get('apply_hpf', True))
            # Restore filter coefficients if they exist
            if filter_data.get('hpf_filter_coeffs') is not None:
                self.hpf_filter_coeffs = np.array(filter_data['hpf_filter_coeffs'])
                self.hpf_group_delay = filter_data.get('hpf_group_delay', 0)
            else:
                self.hpf_filter_coeffs = None
                self.hpf_group_delay = 0

            # Update display
            self.update_hpf_filter_info()
            if self.hpf_filter_coeffs is not None:
                self.show_hpf_response_inline()

            self.main_app.log_status("High-Pass Filter data loaded automatically!")

            return True

        except Exception as e:
            print(f"Error loading High-Pass Filter data: {e}")
            return False

    def on_apply_hpf_changed(self):
        """Callback for when the apply HPF checkbox is toggled."""
        self.save_hpf_filter_data()  # Save state when toggled
        # Auto-apply to data if data is available
        processed_file = os.path.join(self.main_app.file_path.get(), "Processed_angles.csv")
        if os.path.exists(processed_file):
            self.main_app.apply_filters_to_data()
            # Update display if filtered display is enabled
            if self.main_app.display_filtered_angles.get():
                self.main_app.update_acquisition_plots()

    def apply_filters_compensated(self, data, column_names=['LAz_angle', 'LEl_angle', 'RAz_angle', 'REl_angle']):
        """Apply all active filters with group delay compensation to angle data."""
        filtered_data = data.copy()

        for col in column_names:
            if col in filtered_data.columns:
                signal_data = filtered_data[col].values

                # Apply High-Pass Filter if designed and active
                if self.hpf_filter_coeffs is not None and self.apply_hpf.get():
                    signal_data = signal.filtfilt(self.hpf_filter_coeffs, 1, signal_data)

                # Apply Notch Filter if designed and active
                if self.notch_filter_coeffs is not None and self.apply_notch.get():
                    signal_data = signal.filtfilt(self.notch_filter_coeffs, 1, signal_data)

                # Apply Low-Pass Filter if designed and active
                if self.lpf_filter_coeffs is not None and self.apply_lpf.get():
                    signal_data = signal.filtfilt(self.lpf_filter_coeffs, 1, signal_data)

                filtered_data[col] = signal_data

        return filtered_data

class Debug_FFT:
    def __init__(self, parent_notebook):
        self.parent_notebook = parent_notebook
        self.frame = ttk.Frame(parent_notebook)
        parent_notebook.add(self.frame, text="Debug: FFT Explorer")

        # Variables
        self.csv_path = tk.StringVar()
        self.sampling_rate = tk.DoubleVar(value=1000.0)
        self.fft_samples = tk.IntVar(value=512)
        self.timepoint = tk.DoubleVar(value=0)
        self.time_min = 0
        self.time_max = 0
        self.time_offset = 0
        self.data = None

        # UI
        self.setup_ui()

    def setup_ui(self):
        # File selection
        file_frame = ttk.LabelFrame(self.frame, text="CSV File Selection", padding=6)
        file_frame.pack(fill='x', padx=8, pady=4)
        ttk.Label(file_frame, text="CSV File:").pack(side='left', padx=(0, 4))
        file_entry = ttk.Entry(file_frame, textvariable=self.csv_path, width=40, state='readonly')
        file_entry.pack(side='left', padx=(0, 4))
        ttk.Button(file_frame, text="Browse", command=self.browse_csv).pack(side='left', padx=(0, 4))

        # Sampling rate and FFT size
        param_frame = ttk.Frame(self.frame)
        param_frame.pack(fill='x', padx=8, pady=4)
        ttk.Label(param_frame, text="Sampling Rate (Hz):").pack(side='left', padx=(0, 2))
        self.sampling_rate_entry = ttk.Entry(param_frame, width=8)
        self.sampling_rate_entry.pack(side='left', padx=(0, 8))
        self.sampling_rate_entry.insert(0, f"{self.sampling_rate.get():.0f}")
        self.sampling_rate_entry.bind('<Return>', lambda e: self.on_sampling_rate_entry_validate())
        self.sampling_rate_entry.bind('<FocusOut>', lambda e: self.on_sampling_rate_entry_validate())
        self.sampling_rate_entry.bind('<FocusIn>', lambda e: setattr(self, '_last_focus_widget', self.sampling_rate_entry))

        ttk.Label(param_frame, text="FFT Samples:").pack(side='left', padx=(0, 2))
        self.fft_samples_entry = ttk.Entry(param_frame, width=8)
        self.fft_samples_entry.pack(side='left', padx=(0, 8))
        self.fft_samples_entry.insert(0, f"{self.fft_samples.get()}")
        self.fft_samples_entry.bind('<Return>', lambda e: self.on_fft_samples_entry_validate())
        self.fft_samples_entry.bind('<FocusOut>', lambda e: self.on_fft_samples_entry_validate())
        self.fft_samples_entry.bind('<FocusIn>', lambda e: setattr(self, '_last_focus_widget', self.fft_samples_entry))

        # Timepoint slider and entry
        self.time_slider_frame = ttk.Frame(self.frame)
        self.time_slider_frame.pack(fill='x', padx=8, pady=4)
        slider_row = ttk.Frame(self.time_slider_frame)
        slider_row.pack(fill='x')

        ttk.Label(slider_row, text="Timepoint:").pack(side='left', padx=(0, 2))
        self.timepoint_entry = ttk.Entry(slider_row, width=10)
        self.timepoint_entry.pack(side='left', padx=(0, 8))
        self.timepoint_entry.bind('<Return>', lambda e: self.on_timepoint_entry_validate())
        self.timepoint_entry.bind('<FocusOut>', lambda e: self.on_timepoint_entry_validate())
        self.timepoint_entry.bind('<FocusIn>', lambda e: setattr(self, '_last_focus_widget', self.timepoint_entry))

        self.time_slider = tk.Scale(
            slider_row, from_=0, to=1000, orient='horizontal',
            variable=self.timepoint, showvalue=1, resolution=1, length=600
        )
        self.time_slider.pack(side='left', fill='x', expand=True)
        self.time_slider.bind('<ButtonRelease-1>', lambda e: self.on_timepoint_slider_validate())
        self.time_slider_label = ttk.Label(self.time_slider_frame, text="Timepoint: 0 ms")
        self.time_slider_label.pack()

        # FFT Figure (unchanged)
        from matplotlib.gridspec import GridSpec
        self.fig = Figure(figsize=(12, 8))
        gs = GridSpec(6, 2, width_ratios=[1, 1], figure=self.fig, wspace=0.25)
        self.axes = []
        for i in range(6):
            ax_left = self.fig.add_subplot(gs[i, 0])
            ax_left.set_ylabel(f"L{i+1}\n(uV/Hz)", rotation=0, labelpad=28, va='center', ha='right', fontsize=10)
            ax_left.yaxis.set_label_coords(-0.13, 0.5)
            if i != 5:
                ax_left.set_xticklabels([])
            else:
                ax_left.set_xlabel("Frequency (Hz)")
            self.axes.append((ax_left, None))
            ax_right = self.fig.add_subplot(gs[i, 1])
            ax_right.set_ylabel(f"R{i+1}\n(uV/Hz)", rotation=0, labelpad=28, va='center', ha='right', fontsize=10)
            ax_right.yaxis.set_label_coords(-0.13, 0.5)
            if i != 5:
                ax_right.set_xticklabels([])
            else:
                ax_right.set_xlabel("Frequency (Hz)")
            self.axes[i] = (ax_left, ax_right)
        self.fig.subplots_adjust(left=0.13, right=0.98, top=0.97, bottom=0.06, wspace=0.25, hspace=0.35)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.get_tk_widget().pack(fill='both', expand=True, padx=8, pady=8)

        # --- Global click-outside validation for entries ---
        self._last_focus_widget = None
        self.frame.bind_all('<Button-1>', self.global_click_validate, add='+')

    def browse_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.csv_path.set(path)
            self.load_csv(path)

    def load_csv(self, path):
        try:
            self.data = pd.read_csv(path)
            filename = os.path.basename(path)
            if filename.endswith("angles.csv"):
                self._fft_mode = "angles"
                if self.data.shape[1] < 5:
                    messagebox.showerror("Error", "Angles CSV must have at least 5 columns (timestamp + 4 angles).")
                    self.data = None
                    return
            else:
                self._fft_mode = "raw"
                if self.data.shape[1] < 13:
                    messagebox.showerror("Error", "CSV must have at least 13 columns (timestamp + 12 channels).")
                    self.data = None
                    return

            self.time_min = float(self.data.iloc[0, 0])
            self.time_max = float(self.data.iloc[-1, 0])
            self.time_offset = self.time_min

            duration = self.time_max - self.time_min
            self.time_slider.config(from_=0, to=duration)
            self.timepoint.set(0)
            self.time_slider.set(0)
            self.time_slider_label.config(text=f"Timepoint: 0 ms")
            self.timepoint_entry.delete(0, tk.END)
            self.timepoint_entry.insert(0, "0")
            self.update_fft()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load CSV: {e}")
            self.data = None

    def on_sampling_rate_entry_validate(self):
        entry_val = self.sampling_rate_entry.get()
        try:
            val = float(entry_val)
            if val <= 0:
                raise ValueError
        except Exception:
            val = self.sampling_rate.get()
            self.sampling_rate_entry.delete(0, tk.END)
            self.sampling_rate_entry.insert(0, f"{val:.0f}")
            return
        self.sampling_rate.set(val)
        self.sampling_rate_entry.delete(0, tk.END)
        self.sampling_rate_entry.insert(0, f"{val:.0f}")
        self.update_fft()

    def on_fft_samples_entry_validate(self):
        entry_val = self.fft_samples_entry.get()
        try:
            val = int(float(entry_val))
            if val <= 0:
                raise ValueError
        except Exception:
            val = self.fft_samples.get()
            self.fft_samples_entry.delete(0, tk.END)
            self.fft_samples_entry.insert(0, f"{val}")
            return
        self.fft_samples.set(val)
        self.fft_samples_entry.delete(0, tk.END)
        self.fft_samples_entry.insert(0, f"{val}")
        self.update_fft()

    def on_timepoint_entry_validate(self):
        entry_val = self.timepoint_entry.get()
        try:
            val = float(entry_val)
        except Exception:
            val = self.timepoint.get()
            self.timepoint_entry.delete(0, tk.END)
            self.timepoint_entry.insert(0, f"{val:.0f}")
            self.time_slider_label.config(text=f"Timepoint: {val:.0f} ms")
            return
        duration = self.time_max - self.time_min
        val = max(0, min(val, duration))
        self.timepoint.set(val)
        self.time_slider.set(val)
        self.timepoint_entry.delete(0, tk.END)
        self.timepoint_entry.insert(0, f"{val:.0f}")
        self.time_slider_label.config(text=f"Timepoint: {val:.0f} ms")
        self.update_fft()

    def on_timepoint_slider_validate(self):
        val = self.timepoint.get()
        self.timepoint_entry.delete(0, tk.END)
        self.timepoint_entry.insert(0, f"{val:.0f}")
        self.time_slider_label.config(text=f"Timepoint: {val:.0f} ms")
        self.update_fft()

    def update_fft(self):
        if self.data is None:
            for ax_left, ax_right in self.axes:
                ax_left.clear()
                ax_right.clear()
            self.canvas.draw()
            return

        try:
            fs = float(self.sampling_rate_entry.get())
            if fs <= 0:
                raise ValueError
        except Exception:
            fs = self.sampling_rate.get()
        try:
            N = int(float(self.fft_samples_entry.get()))
            if N <= 0:
                raise ValueError
        except Exception:
            N = self.fft_samples.get()

        t0 = float(self.timepoint.get()) + getattr(self, "time_offset", 0)
        idx_center = np.argmin(np.abs(self.data.iloc[:, 0] - t0))
        idx_start = max(0, idx_center - N // 2)
        idx_end = min(len(self.data), idx_start + N)
        if idx_end - idx_start < N:
            idx_start = max(0, idx_end - N)
        segment = self.data.iloc[idx_start:idx_end, :]
        if self._fft_mode == "raw":
            segment = segment.iloc[:, 1:13].values
            if segment.shape[0] < N:
                for ax_left, ax_right in self.axes:
                    ax_left.clear()
                    ax_right.clear()
                self.canvas.draw()
                return

            freqs = np.fft.rfftfreq(N, d=1/fs)
            fft_vals = np.fft.rfft(segment * np.hanning(N)[:, None], axis=0)
            fft_mag = np.abs(fft_vals) / N * 2

            min_plot_freq = 5.0
            freq_mask = freqs >= min_plot_freq
            if not np.any(freq_mask):
                freq_mask = np.ones_like(freqs, dtype=bool)

            for i, (ax_left, ax_right) in enumerate(self.axes):
                ax_left.clear()
                ax_right.clear()
                # Plot only frequencies >= 1 Hz
                plot_mask = freqs >= 1.0
                ax_left.plot(freqs[plot_mask], fft_mag[plot_mask, i])
                ax_right.plot(freqs[plot_mask], fft_mag[plot_mask, i+6])
                ax_left.set_xlim(1.0, fs/2)
                ax_right.set_xlim(1.0, fs/2)
                # Compute max for y-limits only for freqs >= 5 Hz
                max_mask = freqs >= 5.0
                left_max = np.max(fft_mag[max_mask, i]) if np.any(max_mask) and np.any(fft_mag[max_mask, i] > 0) else np.max(fft_mag[:, i])
                right_max = np.max(fft_mag[max_mask, i+6]) if np.any(max_mask) and np.any(fft_mag[max_mask, i+6] > 0) else np.max(fft_mag[:, i+6])
                ax_left.set_ylim(0, left_max * 1.05 if left_max > 0 else 1)
                ax_right.set_ylim(0, right_max * 1.05 if right_max > 0 else 1)
                ax_left.set_yticks(np.linspace(0, ax_left.get_ylim()[1], 3))
                ax_right.set_yticks(np.linspace(0, ax_right.get_ylim()[1], 3))
                if i == 5:
                    ax_left.set_xlabel("Frequency (Hz)")
                    ax_right.set_xlabel("Frequency (Hz)")
                else:
                    ax_left.set_xticklabels([])
                    ax_right.set_xticklabels([])
                ax_left.set_ylabel(f"L{i+1}\n(uV/Hz)", rotation=0, labelpad=28, va='center', ha='right', fontsize=10)
                ax_left.yaxis.set_label_coords(-0.13, 0.5)
                ax_right.set_ylabel(f"R{i+1}\n(uV/Hz)", rotation=0, labelpad=28, va='center', ha='right', fontsize=10)
                ax_right.yaxis.set_label_coords(-0.13, 0.5)

        else:  # angles mode
            # Only use first 4 columns after timestamp
            segment = segment.iloc[:, 1:5].values
            if segment.shape[0] < N:
                for ax_left, ax_right in self.axes:
                    ax_left.clear()
                    ax_right.clear()
                self.canvas.draw()
                return

            freqs = np.fft.rfftfreq(N, d=1/fs)
            fft_vals = np.fft.rfft(segment * np.hanning(N)[:, None], axis=0)
            fft_mag = np.abs(fft_vals) / N * 2

            min_plot_freq = 0.1
            freq_mask = freqs >= min_plot_freq
            if not np.any(freq_mask):
                freq_mask = np.ones_like(freqs, dtype=bool)

            # Show only 4 graphs (use axes 0-3 left, clear right)
            angle_labels = ["LAz (°/Hz)", "LEl (°/Hz)", "RAz (°/Hz)", "REl (°/Hz)"]
            for i in range(4):
                ax_left, ax_right = self.axes[i]
                ax_left.clear()
                ax_right.clear()
                # Plot only frequencies >= 1 Hz
                plot_mask = freqs >= 1.0
                ax_left.plot(freqs[plot_mask], fft_mag[plot_mask, i])
                ax_left.set_xlim(1.0, fs/2)
                # Compute max for y-limits only for freqs >= 10 Hz
                max_mask = freqs >= 10.0
                if np.any(max_mask) and np.any(fft_mag[max_mask, i] > 0):
                    left_max = np.max(fft_mag[max_mask, i])
                else:
                    left_max = np.max(fft_mag[:, i])
                ax_left.set_ylim(0, left_max * 1.05 if left_max > 0 else 1)
                ax_left.set_yticks(np.linspace(0, ax_left.get_ylim()[1], 3))
                if i == 3:
                    ax_left.set_xlabel("Frequency (Hz)")
                else:
                    ax_left.set_xticklabels([])
                ax_left.set_ylabel(angle_labels[i], rotation=0, labelpad=28, va='center', ha='right', fontsize=10)
                ax_left.yaxis.set_label_coords(-0.13, 0.5)
                ax_left.set_title("")
                ax_right.set_xticklabels([])
                ax_right.set_yticklabels([])
                ax_right.set_title("")
                ax_right.set_ylabel("")
                ax_right.set_xlabel("")

            # Clear remaining axes
            for i in range(4, 6):
                ax_left, ax_right = self.axes[i]
                ax_left.clear()
                ax_right.clear()
                ax_left.set_xticklabels([])
                ax_left.set_yticklabels([])
                ax_left.set_title("")
                ax_left.set_ylabel("")
                ax_left.set_xlabel("")
                ax_right.set_xticklabels([])
                ax_right.set_yticklabels([])
                ax_right.set_title("")
                ax_right.set_ylabel("")
                ax_right.set_xlabel("")

        self.fig.subplots_adjust(left=0.13, right=0.98, top=0.97, bottom=0.06, wspace=0.25, hspace=0.35)
        self.canvas.draw()

    def global_click_validate(self, event=None):
        # Use after_idle to ensure focus_get() is updated after the click
        self.frame.after_idle(self._validate_if_entry_lost_focus)

    def _validate_if_entry_lost_focus(self):
        entries = []
        # List all entry widgets for this class
        for name in ['sampling_rate_entry', 'fft_samples_entry', 'timepoint_entry', 'start_entry', 'end_entry']:
            if hasattr(self, name):
                entries.append(getattr(self, name))
        focus = self.frame.focus_get()
        # If focus is not on any entry, validate the last one if needed
        if hasattr(self, '_last_focus_widget') and self._last_focus_widget in entries and focus not in entries:
            last = self._last_focus_widget
            if last == getattr(self, 'sampling_rate_entry', None):
                self.on_sampling_rate_entry_validate()
            elif last == getattr(self, 'fft_samples_entry', None):
                self.on_fft_samples_entry_validate()
            elif last == getattr(self, 'timepoint_entry', None):
                self.on_timepoint_entry_validate()
            elif last == getattr(self, 'start_entry', None):
                self.on_time_entry_validate('start')
            elif last == getattr(self, 'end_entry', None):
                self.on_time_entry_validate('end')
            self._last_focus_widget = None  # Reset after validation

class Debug_temporal:
    def __init__(self, parent_notebook):
        self.parent_notebook = parent_notebook
        self.frame = ttk.Frame(parent_notebook)
        parent_notebook.add(self.frame, text="Debug: Temporal Explorer")

        # Variables
        self.csv_path = tk.StringVar()
        self.sampling_rate = tk.DoubleVar(value=1000.0)
        self.time_min = 0
        self.time_max = 1000
        self.time_start = tk.DoubleVar(value=0)
        self.time_end = tk.DoubleVar(value=1000)
        self.data = None
        self.time_offset = 0

        self.setup_ui()

    def setup_ui(self):
        # File selection
        file_frame = ttk.LabelFrame(self.frame, text="CSV File Selection", padding=6)
        file_frame.pack(fill='x', padx=8, pady=4)
        ttk.Label(file_frame, text="CSV File:").pack(side='left', padx=(0, 4))
        file_entry = ttk.Entry(file_frame, textvariable=self.csv_path, width=40, state='readonly')
        file_entry.pack(side='left', padx=(0, 4))
        ttk.Button(file_frame, text="Browse", command=self.browse_csv).pack(side='left', padx=(0, 4))

        # Sampling rate
        param_frame = ttk.Frame(self.frame)
        param_frame.pack(fill='x', padx=8, pady=4)
        ttk.Label(param_frame, text="Sampling Rate (Hz):").pack(side='left', padx=(0, 2))
        self.sampling_rate_entry = ttk.Entry(param_frame, width=8)
        self.sampling_rate_entry.pack(side='left', padx=(0, 8))
        self.sampling_rate_entry.insert(0, f"{self.sampling_rate.get():.0f}")
        self.sampling_rate_entry.bind('<Return>', lambda e: self.on_sampling_rate_entry_validate())
        self.sampling_rate_entry.bind('<FocusOut>', lambda e: self.on_sampling_rate_entry_validate())
        self.sampling_rate_entry.bind('<FocusIn>', lambda e: setattr(self, '_last_focus_widget', self.sampling_rate_entry))

        # Time interval selection with double slider and entry fields
        self.time_slider_frame = ttk.LabelFrame(self.frame, text="Select Time Interval", padding=6)
        self.time_slider_frame.pack(fill='x', padx=8, pady=4)
        self.time_slider_label = ttk.Label(self.time_slider_frame, text="Interval: 0 ms - 1000 ms")
        self.time_slider_label.pack()

        slider_row = ttk.Frame(self.time_slider_frame)
        slider_row.pack(fill='x')

        # Start time entry
        ttk.Label(slider_row, text="Start:").pack(side='left', padx=(0, 2))
        self.start_entry = ttk.Entry(slider_row, width=8)
        self.start_entry.pack(side='left', padx=(0, 8))
        self.start_entry.insert(0, f"{self.time_start.get():.0f}")
        self.start_entry.bind('<Return>', lambda e: self.on_time_entry_validate('start'))
        self.start_entry.bind('<FocusOut>', lambda e: self.on_time_entry_validate('start'))
        self.start_entry.bind('<FocusIn>', lambda e: setattr(self, '_last_focus_widget', self.start_entry))

        # Start slider
        self.start_scale = tk.Scale(
            slider_row, from_=self.time_min, to=self.time_max, orient='horizontal',
            variable=self.time_start, showvalue=0, resolution=1, length=250
        )
        self.start_scale.pack(side='left', fill='x', expand=True, padx=(0, 5))
        self.start_scale.bind('<ButtonRelease-1>', lambda e: self.on_slider_validate('start'))

        # End time entry
        ttk.Label(slider_row, text="End:").pack(side='left', padx=(0, 2))
        self.end_entry = ttk.Entry(slider_row, width=8)
        self.end_entry.pack(side='left', padx=(0, 8))
        self.end_entry.insert(0, f"{self.time_end.get():.0f}")
        self.end_entry.bind('<Return>', lambda e: self.on_time_entry_validate('end'))
        self.end_entry.bind('<FocusOut>', lambda e: self.on_time_entry_validate('end'))
        self.end_entry.bind('<FocusIn>', lambda e: setattr(self, '_last_focus_widget', self.end_entry))

        # End slider
        self.end_scale = tk.Scale(
            slider_row, from_=self.time_min, to=self.time_max, orient='horizontal',
            variable=self.time_end, showvalue=0, resolution=1, length=250
        )
        self.end_scale.pack(side='left', fill='x', expand=True, padx=(0, 0))
        self.end_scale.bind('<ButtonRelease-1>', lambda e: self.on_slider_validate('end'))

        # Matplotlib Figure: 6 rows, 4 columns (label, left, label, right)
        from matplotlib.gridspec import GridSpec
        self.fig = Figure(figsize=(12, 8))
        gs = GridSpec(6, 4, width_ratios=[0.10, 1, 0.10, 1], figure=self.fig, wspace=0.25)
        self.axes = []
        self.label_axes_left = []
        self.label_axes_right = []

        for i in range(6):
            ax_label_left = self.fig.add_subplot(gs[i, 0])
            ax_label_left.axis('off')
            ax_label_left.text(
                0.98, 0.5, f"L{i+1}", va='center', ha='right', fontsize=11, fontweight='bold'
            )
            ax_label_left.text(
                0.98, 0.15, "uV", va='center', ha='right', fontsize=9, rotation=0, color='gray'
            )
            self.label_axes_left.append(ax_label_left)
            ax_left = self.fig.add_subplot(gs[i, 1])
            if i == 5:
                ax_left.set_xlabel("Time (ms)")
            else:
                ax_left.set_xticklabels([])
            self.axes.append((ax_left, None))
            ax_label_right = self.fig.add_subplot(gs[i, 2])
            ax_label_right.axis('off')
            ax_label_right.text(
                0.02, 0.5, f"R{i+1}", va='center', ha='left', fontsize=11, fontweight='bold'
            )
            ax_label_right.text(
                0.02, 0.15, "uV", va='center', ha='left', fontsize=9, rotation=0, color='gray'
            )
            self.label_axes_right.append(ax_label_right)
            ax_right = self.fig.add_subplot(gs[i, 3])
            if i == 5:
                ax_right.set_xlabel("Time (ms)")
            else:
                ax_right.set_xticklabels([])
            self.axes[i] = (ax_left, ax_right)

        self.fig.subplots_adjust(left=0.08, right=0.98, top=0.97, bottom=0.06, wspace=0.35, hspace=0.25)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame)
        self.canvas.get_tk_widget().pack(fill='both', expand=True, padx=8, pady=8)

        # --- Global click-outside validation for entries ---
        self._last_focus_widget = None
        self.frame.bind_all('<Button-1>', self.global_click_validate, add='+')

    def browse_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.csv_path.set(path)
            self.load_csv(path)

    def load_csv(self, path):
        try:
            self.data = pd.read_csv(path)
            if self.data.shape[1] < 13:
                messagebox.showerror("Error", "CSV must have at least 13 columns (timestamp + 12 channels).")
                self.data = None
                return
            self.time_min = float(self.data.iloc[0, 0])
            self.time_max = float(self.data.iloc[-1, 0])
            self.time_offset = self.time_min  # Set offset to first timestamp

            # Set slider ranges to RELATIVE values
            duration = self.time_max - self.time_min
            self.start_scale.config(from_=0, to=duration)
            self.end_scale.config(from_=0, to=duration)
            # Set values to min and max (relative)
            self.time_start.set(0)
            self.time_end.set(duration)
            self.start_scale.set(0)
            self.end_scale.set(duration)
            # PATCH: Update entry fields to match slider values
            self.start_entry.delete(0, tk.END)
            self.start_entry.insert(0, "0")
            self.end_entry.delete(0, tk.END)
            self.end_entry.insert(0, f"{int(duration)}")
            # Update label and plot
            self.time_slider_label.config(
                text=f"Interval: {int(self.time_start.get())} ms - {int(self.time_end.get())} ms"
            )
            self.update_temporal()
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load CSV: {e}")
            self.data = None

    def update_time_slider(self):
        # Update slider ranges and values
        self.start_scale.config(from_=self.time_min, to=self.time_max)
        self.end_scale.config(from_=self.time_min, to=self.time_max)
        self.start_scale.set(self.time_start.get())
        self.end_scale.set(self.time_end.get())
        self.on_time_slider_change()

    def on_time_slider_change(self):
        # Ensure start <= end
        if self.time_start.get() > self.time_end.get():
            self.time_start.set(self.time_end.get())
        if self.time_end.get() < self.time_start.get():
            self.time_end.set(self.time_start.get())
        self.time_slider_label.config(
            text=f"Interval: {int(self.time_start.get())} ms - {int(self.time_end.get())} ms"
        )
        self.update_temporal()

    def on_sampling_rate_entry_validate(self):
        entry_val = self.sampling_rate_entry.get()
        try:
            val = float(entry_val)
            if val <= 0:
                raise ValueError
        except Exception:
            val = self.sampling_rate.get()
            self.sampling_rate_entry.delete(0, tk.END)
            self.sampling_rate_entry.insert(0, f"{val:.0f}")
            return
        self.sampling_rate.set(val)
        self.sampling_rate_entry.delete(0, tk.END)
        self.sampling_rate_entry.insert(0, f"{val:.0f}")
        self.update_temporal()

    def on_time_entry_validate(self, which):
        try:
            start = float(self.start_entry.get())
            end = float(self.end_entry.get())
        except Exception:
            # If invalid, revert to last valid values
            start = self.time_start.get()
            end = self.time_end.get()
            self.start_entry.delete(0, tk.END)
            self.start_entry.insert(0, f"{start:.0f}")
            self.end_entry.delete(0, tk.END)
            self.end_entry.insert(0, f"{end:.0f}")
            return

        duration = self.time_max - self.time_min
        start = max(0, min(start, duration))
        end = max(0, min(end, duration))

        if which == 'start':
            if start > end:
                start = end
            self.time_start.set(start)
            self.start_scale.set(start)
            self.start_entry.delete(0, tk.END)
            self.start_entry.insert(0, f"{start:.0f}")
        elif which == 'end':
            if end < start:
                end = start
            self.time_end.set(end)
            self.end_scale.set(end)
            self.end_entry.delete(0, tk.END)
            self.end_entry.insert(0, f"{end:.0f}")

        self.time_slider_label.config(
            text=f"Interval: {int(self.time_start.get())} ms - {int(self.time_end.get())} ms"
        )
        self.update_temporal()

    def on_slider_validate(self, which):
        start = self.time_start.get()
        end = self.time_end.get()
        if which == 'start':
            if start > end:
                self.time_start.set(end)
                self.start_scale.set(end)
                start = end  # update local variable
            # PATCH: update start_entry field to match slider
            self.start_entry.delete(0, tk.END)
            self.start_entry.insert(0, f"{start:.0f}")
        elif which == 'end':
            if end < start:
                self.time_end.set(start)
                self.end_scale.set(start)
                end = start  # update local variable
            # PATCH: update end_entry field to match slider
            self.end_entry.delete(0, tk.END)
            self.end_entry.insert(0, f"{end:.0f}")
        self.time_slider_label.config(
            text=f"Interval: {int(self.time_start.get())} ms - {int(self.time_end.get())} ms"
        )
        self.update_temporal()

    def update_temporal(self):
        if self.data is None:
            for ax_left, ax_right in self.axes:
                ax_left.clear()
                ax_right.clear()
            self.canvas.draw()
            return

        # Convert relative slider values to absolute timestamps
        t_start = self.time_start.get() + self.time_offset
        t_end = self.time_end.get() + self.time_offset
        if t_start == t_end:
            t_start -= 1
            t_end += 1
        mask = (self.data.iloc[:, 0] >= t_start) & (self.data.iloc[:, 0] <= t_end)
        segment = self.data.loc[mask]
        # Display times as relative to offset
        times = segment.iloc[:, 0].values - self.time_offset

        for i, (ax_left, ax_right) in enumerate(self.axes):
            ax_left.clear()
            ax_right.clear()
            if len(times) > 0:
                ax_left.plot(times, segment.iloc[:, 1 + i], color='tab:blue')
                ax_right.plot(times, segment.iloc[:, 7 + i], color='tab:orange')
            ax_left.set_xlim(self.time_start.get(), self.time_end.get())
            ax_right.set_xlim(self.time_start.get(), self.time_end.get())
            if i == 5:
                ax_left.set_xlabel("Time (ms)")
                ax_right.set_xlabel("Time (ms)")
            else:
                ax_left.set_xticklabels([])
                ax_right.set_xticklabels([])
            ax_left.set_ylabel("")
            ax_right.set_ylabel("")
            ax_left.set_title("")
            ax_right.set_title("")

        self.fig.subplots_adjust(left=0.11, right=0.98, top=0.97, bottom=0.06, wspace=0.15, hspace=0.25)
        self.canvas.draw()

    def global_click_validate(self, event=None):
        # Use after_idle to ensure focus_get() is updated after the click
        self.frame.after_idle(self._validate_if_entry_lost_focus)

    def _validate_if_entry_lost_focus(self):
        entries = [self.sampling_rate_entry, self.start_entry, self.end_entry]
        focus = self.frame.focus_get()
        if hasattr(self, '_last_focus_widget') and self._last_focus_widget in entries and focus not in entries:
            last = self._last_focus_widget
            if last == self.sampling_rate_entry:
                self.on_sampling_rate_entry_validate()
            elif last == self.start_entry:
                self.on_time_entry_validate('start')
            elif last == self.end_entry:
                self.on_time_entry_validate('end')
            self._last_focus_widget = None  # Reset after validation

def main():
    root = tk.Tk()
    
    # Set style
    style = ttk.Style()
    if "clam" in style.theme_names():
        style.theme_use("clam")
    
    app = AEOG_GUI(root)
    
    def on_closing():
        if app.is_calibrating or app.is_acquiring:
            if messagebox.askokcancel("Quit", "Calibration or A-EOG Data Acquisition in progress. Quit anyway?"):
                pygame.quit()
                root.destroy()
        else:
            pygame.quit()
            root.destroy()
    
    root.protocol("WM_DELETE_WINDOW", on_closing)
    
    try:
        root.mainloop()
    except KeyboardInterrupt:
        pygame.quit()
        root.destroy()


if __name__ == "__main__":
    main()