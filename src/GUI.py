import json
import os
import tkinter as tk
import tkinter.filedialog as fd
import tkinter.ttk as ttk
from PIL import ImageTk, Image
import sv_ttk
import copy
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter
from collections import defaultdict


from src.frames import *
from src.options import * 
from src.data_reduction import reduce_data

import tkinter as tk
from tkinter import ttk
from collections import defaultdict

def filter_by_common_frame_shape(bias_list, flat_list, shifted_flat_list, science_list, complamp_list):
	# Group frames by their data.shape
	shape_dict = defaultdict(lambda: {'bias': [], 'flat': [], 'science': [], 'complamp': [], 'shifted_flat': []})

	for frm in bias_list:
		shape_dict[frm.data.shape]['bias'].append(frm)
	for frm in flat_list:
		shape_dict[frm.data.shape]['flat'].append(frm)
	for frm in science_list:
		shape_dict[frm.data.shape]['science'].append(frm)
	for frm in complamp_list:
		shape_dict[frm.data.shape]['complamp'].append(frm)
	for frm in shifted_flat_list:
		shape_dict[frm.data.shape]['shifted_flat'].append(frm)

	# Filter valid groups (those with at least one frame in each of the 4 main lists)
	valid_shapes = {shape: frames for shape, frames in shape_dict.items()
					if all(frames[key] for key in ['bias', 'flat', 'science', 'complamp'])}

	if not valid_shapes:
		print("No common valid frame shape group found.")
		return  # or raise error

	# If only one valid shape, select it automatically
	if len(valid_shapes) == 1:
		selected_shape = next(iter(valid_shapes))
	else:
		selected_shape = _prompt_user_for_shape_selection(valid_shapes)

	# Replace original lists with filtered ones
	selected_group = valid_shapes[selected_shape]
	bias_list[:] = selected_group['bias']
	flat_list[:] = selected_group['flat']
	science_list[:] = selected_group['science']
	complamp_list[:] = selected_group['complamp']
	shifted_flat_list[:] = selected_group['shifted_flat']

def _prompt_user_for_shape_selection(valid_shapes):
	selected = {"shape": None}

	def on_confirm():
		selection = shape_var.get()
		selected["shape"] = eval(selection)
		popup.destroy()

	popup = tk.Tk()
	popup.title("Select Frame Shape Group")
	tk.Label(popup, text="Multiple valid frame shape groups found. Please select one:").pack(padx=10, pady=10)

	shape_var = tk.StringVar()
	dropdown = ttk.Combobox(popup, textvariable=shape_var, state="readonly", width=60)

	options = []
	for shape, frames in valid_shapes.items():
		count_summary = f"{len(frames['bias'])} bias, {len(frames['flat'])} flat, {len(frames['science'])} science, {len(frames['complamp'])} complamp"
		options.append(f"{shape} → {count_summary}")
	dropdown['values'] = options
	dropdown.pack(padx=10, pady=5)
	dropdown.current(0)

	tk.Button(popup, text="OK", command=on_confirm).pack(pady=10)
	popup.mainloop()

	return eval(shape_var.get().split(' → ')[0])  # e.g., "(1024, 1024)"



class ConfigWindow(tk.Toplevel):
	def __init__(self, master):
		self.frame_config = master.frame_config
		self.reduction_options = master.reduction_options
		super().__init__(master)
		self.title("MIDIR Configuration")
		self.geometry("600x640+550+200")
		self.resizable(False, False)
		try:
			imgicon = ImageTk.PhotoImage(Image.open("src/MIDIR_icon.png"))
			self.iconphoto(False, imgicon)
			self._icon_img_ref = imgicon
		except Exception as e:
			print(e)

		self.protocol("WM_DELETE_WINDOW", self.on_close)

		notebook = ttk.Notebook(self)
		notebook.pack(fill='both', expand=True, padx=10, pady=10)

		# File Selection Tab
		file_selection_tab = ttk.Frame(notebook)
		notebook.add(file_selection_tab, text='File Selection')

		# Auto-Detect button
		autodetect_button = ttk.Button(file_selection_tab, text="Auto-Detect", command=self.auto_detect)
		#autodetect_button["state"] = "disabled"
		autodetect_button.pack(fill='x', pady=(10, 20), padx=20)

		# Visual separator labeled "OR"
		separator_frame = ttk.Frame(file_selection_tab)
		separator_frame.pack(fill='x', padx=20, pady=(0, 10))
		ttk.Separator(separator_frame, orient='horizontal').pack(fill='x', side='left', expand=True, padx=(0, 5))
		ttk.Label(separator_frame, text="OR").pack(side='left')
		ttk.Separator(separator_frame, orient='horizontal').pack(fill='x', side='left', expand=True, padx=(5, 0))


		# Buttons for selecting different types of frames
		self.file_buttons = {}
		self.file_count_labels = {}
		frame_types = ["Bias Frames", "Flat Frames", "Shifted Flat Frames (Optional)", "Arc Frames", "Science Frames"]

		for frame_type in frame_types:
			frame = ttk.Frame(file_selection_tab)
			frame.pack(fill='x', padx=20, pady=3)

			button = ttk.Button(frame, text=f"Select {frame_type}", 
								command=lambda ft=frame_type: self.select_files(ft))
			button.pack(side='left', fill='x', expand=True)

			count_label = ttk.Label(frame, text="0 files", foreground="red")
			count_label.pack(side='right', padx=10)

			self.file_buttons[frame_type] = button
			self.file_count_labels[frame_type] = count_label

		separator_frame2 = ttk.Frame(file_selection_tab)
		separator_frame2.pack(fill='x', padx=20, pady=(10, 0))
		ttk.Separator(separator_frame2, orient='horizontal').pack(fill='x', side='left', expand=True, padx=(0, 0))

		# Header/Image HDU Index row
		hdu_frame = ttk.Frame(file_selection_tab)
		hdu_frame.pack(fill='x', padx=20, pady=(15, 5))

		ttk.Label(hdu_frame, text="Header HDU Index:").pack(side='left')
		self.header_hdu_var = tk.IntVar(value=getattr(self.reduction_options, "header_hdu_ind", 0))
		ttk.Entry(hdu_frame, textvariable=self.header_hdu_var, width=5).pack(side='left', padx=(5, 20))
		self.header_hdu_var.trace_add("write", lambda *_, v=self.header_hdu_var: setattr(self.reduction_options, "header_hdu_ind", v.get()))

		ttk.Label(hdu_frame, text="Image HDU Index:").pack(side='left')
		self.image_hdu_var = tk.IntVar(value=getattr(self.reduction_options, "image_hdu_ind", 0))
		ttk.Entry(hdu_frame, textvariable=self.image_hdu_var, width=5).pack(side='left', padx=(5, 0))
		self.image_hdu_var.trace_add("write", lambda *_, v=self.image_hdu_var: setattr(self.reduction_options, "image_hdu_ind", v.get()))

		ttk.Label(file_selection_tab, text="Comparison Lamp Matching").pack(anchor='w', padx=20, pady=(15, 5))

		# Frame for grouping
		match_frame = ttk.Frame(file_selection_tab)
		match_frame.pack(fill='x', padx=20)
		# --- Comparison Matching Settings ---

		# Mode dropdown
		mode_dict = {"count": "Select N Comparison Frames", "time": "Select by Time Window (minutes)"}
		match_modes = ["Select N Comparison Frames", "Select by Time Window (minutes)"]
		self.match_mode_var = tk.StringVar(value=mode_dict[self.reduction_options.comparison_match_mode])
		ttk.Label(match_frame, text="Match Mode:").grid(row=0, column=0, sticky='w', pady=2)
		match_mode_dropdown = ttk.Combobox(match_frame, textvariable=self.match_mode_var, values=match_modes, state="readonly")
		match_mode_dropdown.grid(row=0, column=1, sticky='ew', padx=5, pady=2)
		match_frame.columnconfigure(1, weight=1)

		# Position dropdown
		match_positions = ["BEFORE", "AFTER", "AROUND"]
		self.match_position_var = tk.StringVar(value=self.reduction_options.comparison_match_position.upper())
		ttk.Label(match_frame, text="Relative Position:").grid(row=1, column=0, sticky='w', pady=2)
		match_position_dropdown = ttk.Combobox(match_frame, textvariable=self.match_position_var, values=match_positions, state="readonly")
		match_position_dropdown.grid(row=1, column=1, sticky='ew', padx=5, pady=2)

		# Value input (N or Minutes)
		self.match_value_var = tk.IntVar(value=self.reduction_options.comparison_match_value)
		ttk.Label(match_frame, text="N / Minutes:").grid(row=2, column=0, sticky='w', pady=2)
		ttk.Entry(match_frame, textvariable=self.match_value_var).grid(row=2, column=1, sticky='ew', padx=5, pady=2)

		# Bind to reduction_options
		def update_comparison_match_settings(*args):
			self.reduction_options.comparison_match_mode = (
				"count" if self.match_mode_var.get() == "Select N Comparison Frames" else "time"
			)
			self.reduction_options.comparison_match_value = self.match_value_var.get()
			self.reduction_options.comparison_match_position = self.match_position_var.get().lower()

		self.match_mode_var.trace_add("write", update_comparison_match_settings)
		self.match_position_var.trace_add("write", update_comparison_match_settings)
		self.match_value_var.trace_add("write", update_comparison_match_settings)

		# --- General Settings Tab (GRID layout) ---
		general_settings_tab = ttk.Frame(notebook)
		notebook.add(general_settings_tab, text='General Settings')
		# Wrap upper content in a frame
		content_frame = ttk.Frame(general_settings_tab)
		content_frame.pack(fill='both', expand=True)

		# Output Directory Path
		ttk.Label(content_frame, text="Output Directory Path").pack(anchor='w', padx=20, pady=(20, 5))
		self.output_dir_var = tk.StringVar(value=self.reduction_options.output_dict)
		output_dir_entry = ttk.Entry(content_frame, textvariable=self.output_dir_var)
		output_dir_entry.pack(fill='x', padx=20)
		self.output_dir_var.trace_add("write", lambda *args: setattr(self.reduction_options, "output_dict", self.output_dir_var.get()))

		# Output Table Path
		ttk.Label(content_frame, text="Output Table Path").pack(anchor='w', padx=20, pady=(10, 5))
		self.output_table_var = tk.StringVar(value=self.reduction_options.outputfile_path)
		output_table_entry = ttk.Entry(content_frame, textvariable=self.output_table_var)
		output_table_entry.pack(fill='x', padx=20)
		self.output_table_var.trace_add("write", lambda *args: setattr(self.reduction_options, "outputfile_path", self.output_table_var.get()))

		# Rotate Images By Dropdown
		ttk.Label(content_frame, text="Rotate Images by").pack(anchor='w', padx=20, pady=(10, 5))
		rotation_options = ["0°", "90°", "180°", "270°"]
		rotation_degrees_map = {"0°": 0, "90°": 90, "180°": 180, "270°": 270}
		inverse_rotation_map = {v: k for k, v in rotation_degrees_map.items()}

		# Use string label, but store int value
		self.rotation_var = tk.StringVar(value=inverse_rotation_map.get(getattr(self.reduction_options, "rotationangle", 0), "0°"))
		rotation_dropdown = ttk.Combobox(content_frame, textvariable=self.rotation_var, values=rotation_options, state="readonly")
		rotation_dropdown.pack(fill='x', padx=20)

		def update_rotation(*_):
			selected_label = self.rotation_var.get()
			angle = rotation_degrees_map.get(selected_label, 0)
			setattr(self.reduction_options, "rotationangle", angle)

		self.rotation_var.trace_add("write", update_rotation)

		# Skyflux Separation
		ttk.Label(content_frame, text="Skyflux Separation (Pixels)").pack(anchor='w', padx=20, pady=(10, 5))
		self.skyfluxsep_var = tk.IntVar(value=getattr(self.reduction_options, "skyfluxsep", 10))
		ttk.Entry(content_frame, textvariable=self.skyfluxsep_var).pack(fill='x', padx=20)
		self.skyfluxsep_var.trace_add("write", lambda *_, var=self.skyfluxsep_var: setattr(self.reduction_options, "skyfluxsep", var.get()))

		# Telescope Location Dropdown
		ttk.Label(content_frame, text="Telescope Location").pack(anchor='w', padx=20, pady=(10, 5))
		location_options = ["La Silla Observatory", "Cerro Pachon", "Roque de los Muchachos"]
		self.location_var = tk.StringVar(value=getattr(self.reduction_options, "telescope_location", "Cerro Pachon"))
		location_dropdown = ttk.Combobox(content_frame, textvariable=self.location_var, values=location_options, state="readonly")
		location_dropdown.pack(fill='x', padx=20)
		self.location_var.trace_add("write", lambda *_, var=self.location_var: setattr(self.reduction_options, "telescope_location", var.get()))

		# Inline Crop Frames Manually Checkbox with Entry Fields
		crop_frame = ttk.Frame(content_frame)
		crop_frame.pack(fill='x', padx=20, pady=(10, 0))

		self.cropframes_var = tk.BooleanVar(value=getattr(self.reduction_options, "manual_crop", False))

		def update_crop_status():
			self.reduction_options.manual_crop = self.cropframes_var.get()
			state = "normal" if self.cropframes_var.get() else "disabled"
			for entry in (self.xlo_entry, self.xhi_entry, self.ylo_entry, self.yhi_entry):
				entry.configure(state=state)

		ttk.Checkbutton(
			crop_frame,
			text="Crop Frames Manually",
			variable=self.cropframes_var,
			command=update_crop_status
		).pack(side='left', padx=(0, 10))

		# X range
		ttk.Label(crop_frame, text="X").pack(side='left')
		self.xlo_var = tk.IntVar(value=getattr(self.reduction_options, "x_lo", 0))
		self.xlo_entry = ttk.Entry(crop_frame, textvariable=self.xlo_var, width=5)
		self.xlo_entry.pack(side='left', padx=(5, 0))
		ttk.Label(crop_frame, text="-").pack(side='left', padx=(3, 3))
		self.xhi_var = tk.IntVar(value=getattr(self.reduction_options, "x_hi", 0))
		self.xhi_entry = ttk.Entry(crop_frame, textvariable=self.xhi_var, width=5)
		self.xhi_entry.pack(side='left', padx=(0, 10))

		# Y range
		ttk.Label(crop_frame, text="Y").pack(side='left')
		self.ylo_var = tk.IntVar(value=getattr(self.reduction_options, "y_lo", 0))
		self.ylo_entry = ttk.Entry(crop_frame, textvariable=self.ylo_var, width=5)
		self.ylo_entry.pack(side='left', padx=(5, 0))
		ttk.Label(crop_frame, text="-").pack(side='left', padx=(3, 3))
		self.yhi_var = tk.IntVar(value=getattr(self.reduction_options, "y_hi", 0))
		self.yhi_entry = ttk.Entry(crop_frame, textvariable=self.yhi_var, width=5)
		self.yhi_entry.pack(side='left')

		# Bind vars to reduction_options
		self.xlo_var.trace_add("write", lambda *_, v=self.xlo_var: setattr(self.reduction_options, "x_lo", v.get()))
		self.xhi_var.trace_add("write", lambda *_, v=self.xhi_var: setattr(self.reduction_options, "x_hi", v.get()))
		self.ylo_var.trace_add("write", lambda *_, v=self.ylo_var: setattr(self.reduction_options, "y_lo", v.get()))
		self.yhi_var.trace_add("write", lambda *_, v=self.yhi_var: setattr(self.reduction_options, "y_hi", v.get()))

		# Set initial enabled/disabled state
		update_crop_status()
				
		# Checkbox: Cosmic Rejection
		self.cosmic_reject_var = tk.BooleanVar(value=getattr(self.reduction_options, "cosmicrejection", False))
		ttk.Checkbutton(
			content_frame, text="Attempt Automatic Cosmic Rejection",
			variable=self.cosmic_reject_var,
			command=lambda: setattr(self.reduction_options, "cosmicrejection", self.cosmic_reject_var.get())
		).pack(anchor='w', padx=20, pady=(10, 0))

		# Checkbox: Boxcut Extraction
		self.boxcut_var = tk.BooleanVar(value=getattr(self.reduction_options, "use_boxcut", True))
		ttk.Checkbutton(
			content_frame, text="Use Boxcut Extraction",
			variable=self.boxcut_var,
			command=lambda: setattr(self.reduction_options, "use_boxcut", self.boxcut_var.get())
		).pack(anchor='w', padx=20, pady=(0, 0))

		# Show Plots checkbox
		self.show_plots_var = tk.BooleanVar(value=getattr(self.reduction_options, "debugimages", False))
		ttk.Checkbutton(
			content_frame, text="Show Plots", variable=self.show_plots_var,
			command=lambda: setattr(self.reduction_options, "debugimages", self.show_plots_var.get())
		).pack(anchor='w', padx=20, pady=(0, 10))



		# Buttons frame that sticks to the bottom
		preset_button_frame = ttk.Frame(general_settings_tab)
		preset_button_frame.pack(side="bottom", fill="x", padx=20, pady=20)
		preset_button_frame.columnconfigure((0, 1), weight=1)


		def load_preset():
			preset_path = fd.askopenfilename(
				title="Load Preset",
				initialdir=os.path.join(os.getcwd(), "presets"),
				filetypes=[("JSON files", "*.json")]
			)
			if preset_path:
				self.reduction_options.load_from_json(json_path=preset_path)
				self.refresh_ui_fields()

		def save_preset():
			preset_path = fd.asksaveasfilename(
				title="Save Preset",
				initialdir=os.path.join(os.getcwd(), "presets"),
				defaultextension=".json",
				filetypes=[("JSON files", "*.json")]
			)
			if preset_path:
				self.reduction_options.save_to_json(json_path=preset_path)

		ttk.Button(preset_button_frame, text="Load Preset", command=load_preset).pack(side="left", expand=True, fill="x", padx=(0, 5))
		ttk.Button(preset_button_frame, text="Save Preset", command=save_preset).pack(side="left", expand=True, fill="x", padx=(5, 0))


		wavelength_reduction_tab = ttk.Frame(notebook)
		notebook.add(wavelength_reduction_tab, text='Wavelength Reduction Settings')

		# --- Wavelength Reduction Options ---

		# Linelist Dropdown Selector
		ttk.Label(wavelength_reduction_tab, text="Select Linelist").pack(anchor='w', padx=20, pady=(20, 5))

		linelist_folder = os.path.join(os.getcwd(), "linelists")
		linelist_files = ["linelists/"+f for f in os.listdir(linelist_folder) if os.path.isfile(os.path.join(linelist_folder, f))]

		self.linelist_var = tk.StringVar(value=self.reduction_options.linelist)
		linelist_dropdown = ttk.Combobox(
			wavelength_reduction_tab,
			textvariable=self.linelist_var,
			values=linelist_files,
			state="readonly"
		)
		linelist_dropdown.pack(fill='x', padx=20)
		self.linelist_var.trace_add("write", lambda *_, var=self.linelist_var: setattr(self.reduction_options, "linelist", var.get()))

		ttk.Label(wavelength_reduction_tab, text="Amount of MC Samples (per Core)").pack(anchor='w', padx=20, pady=(20, 5))
		self.mc_sample_var = tk.IntVar(value=self.reduction_options.sampleamt)
		ttk.Entry(wavelength_reduction_tab, textvariable=self.mc_sample_var).pack(fill='x', padx=20)
		self.mc_sample_var.trace_add("write", lambda *_, var=self.mc_sample_var: setattr(self.reduction_options, "sampleamt", var.get()))

		ttk.Label(wavelength_reduction_tab, text="Acceptance Rate Modifier").pack(anchor='w', padx=20, pady=(10, 5))
		self.accept_param_var = tk.DoubleVar(value=self.reduction_options.accept_param)
		ttk.Entry(wavelength_reduction_tab, textvariable=self.accept_param_var).pack(fill='x', padx=20)
		self.accept_param_var.trace_add("write", lambda *_, var=self.accept_param_var: setattr(self.reduction_options, "accept_param", var.get()))

		# Spacer and disclaimer
		ttk.Label(wavelength_reduction_tab, text="").pack(pady=5)
		ttk.Label(wavelength_reduction_tab, text="Modify these parameters only if you know what you are doing!", foreground='red').pack(pady=(0, 10))

		# Grid for MCMC polynomial parameters
		param_frame = ttk.Frame(wavelength_reduction_tab)
		param_frame.pack(fill='x', padx=20)

		# Column headers
		ttk.Label(param_frame, text="", anchor="center").grid(row=0, column=0, padx=5, pady=5)  # Empty top-left
		headers = ["offset", "linear", "quadratic", "cubic"]
		for col, text in enumerate(headers, start=1):
			ttk.Label(param_frame, text=text, anchor="center").grid(row=0, column=col, padx=5, pady=5, sticky="ew")

		# Row labels
		row_labels = ["Zero Point", "Stepsize", "Bounds"]
		for row, text in enumerate(row_labels, start=1):
			ttk.Label(param_frame, text=text, anchor="e").grid(row=row, column=0, padx=5, pady=2, sticky="e")

		# Row 1: Zero points
		self.zero_point_vars = []
		zero_attrs = ["offset_zero", "extent_zero", "quad_zero", "cube_zero"]
		for col, attr in enumerate(zero_attrs, start=1):
			if attr != "offset_zero":
				var = tk.DoubleVar(value=getattr(self.reduction_options, attr))
			else:
				var = tk.StringVar(value=getattr(self.reduction_options, attr))
			entry = ttk.Entry(param_frame, textvariable=var)
			entry.grid(row=1, column=col, padx=5, pady=2, sticky="ew")
			var.trace_add("write", lambda *_, a=attr, v=var: setattr(self.reduction_options, a, v.get()))
			self.zero_point_vars.append(var)


		# Row 2: Stepsizes
		self.stepsize_vars = []
		stepsize_attrs = ["offset_stepsize", "linear_stepsize", "quad_stepsize", "cube_stepsize"]
		for col, attr in enumerate(stepsize_attrs, start=1):
			var = tk.DoubleVar(value=getattr(self.reduction_options, attr))
			entry = ttk.Entry(param_frame, textvariable=var)
			entry.grid(row=2, column=col, padx=5, pady=2, sticky="ew")
			var.trace_add("write", lambda *_, a=attr, v=var: setattr(self.reduction_options, a, v.get()))
			self.stepsize_vars.append(var)

		# Row 3: Bounds (tuples as string input)
		self.bounds_vars = []
		bounds_attrs = ["c_cov", "s_cov", "q_cov", "cub_cov"]
		for col, attr in enumerate(bounds_attrs, start=1):
			var = tk.StringVar(value=str(getattr(self.reduction_options, attr)))
			entry = ttk.Entry(param_frame, textvariable=var)
			entry.grid(row=3, column=col, padx=5, pady=2, sticky="ew")

			def update_bounds(attr=attr, var=var):
				try:
					val = eval(var.get(), {"__builtins__": None}, {})
					setattr(self.reduction_options, attr, val)
				except:
					pass

			var.trace_add("write", lambda *_, a=attr, v=var: update_bounds(a, v))
			self.bounds_vars.append(var)


		# Make all parameter columns expand equally
		for col in range(1, 5):
			param_frame.grid_columnconfigure(col, weight=1)

		# Arc Lamp Maximum Filter Window Size
		ttk.Label(wavelength_reduction_tab, text="Arc Lamp Maximum Filter Window Size").pack(anchor='w', padx=20, pady=(20, 5))
		self.lampfilterwindow_var = tk.IntVar(value=getattr(self.reduction_options, "lampfilterwindow", 15))
		ttk.Entry(wavelength_reduction_tab, textvariable=self.lampfilterwindow_var).pack(fill='x', padx=20)
		self.lampfilterwindow_var.trace_add("write", lambda *_, var=self.lampfilterwindow_var: setattr(self.reduction_options, "lampfilterwindow", var.get()))


		# Populate count indicators on window load
		for ft in self.file_buttons.keys():
			existing = getattr(self.frame_config, {
				"Bias Frames": "biases",
				"Flat Frames": "flats",
				"Shifted Flat Frames (Optional)": "shiftedflats",
				"Arc Frames": "comparisonframes",
				"Science Frames": "scienceframes"
			}[ft])
			try:
				self.update_file_count_label(ft, len(existing.filepath_list()))
			except AttributeError:
				self.update_file_count_label(ft, 0)

	def refresh_ui_fields(self):
		# Match Mode
		match_mode_display = "Select N Comparison Frames" if self.reduction_options.comparison_match_mode == "count" else "Select by Time Window (minutes)"
		self.match_mode_var.set(match_mode_display)
		self.match_position_var.set(self.reduction_options.comparison_match_position.upper())
		self.match_value_var.set(self.reduction_options.comparison_match_value)

		# General Settings
		self.output_dir_var.set(self.reduction_options.output_dict)
		self.output_table_var.set(self.reduction_options.outputfile_path)
		self.rotation_var.set({0: "0°", 90: "90°", 180: "180°", 270: "270°"}.get(self.reduction_options.rotationangle, "0°"))
		self.skyfluxsep_var.set(self.reduction_options.skyfluxsep)
		self.location_var.set(self.reduction_options.telescope_location)

		self.cosmic_reject_var.set(self.reduction_options.cosmicrejection)
		self.boxcut_var.set(self.reduction_options.use_boxcut)
		self.show_plots_var.set(self.reduction_options.debugimages)

		self.cropframes_var.set(self.reduction_options.manual_crop)
		self.xlo_var.set(self.reduction_options.x_lo)
		self.xhi_var.set(self.reduction_options.x_hi)
		self.ylo_var.set(self.reduction_options.y_lo)
		self.yhi_var.set(self.reduction_options.y_hi)

		self.header_hdu_var.set(self.reduction_options.header_hdu_ind)
		self.image_hdu_var.set(self.reduction_options.image_hdu_ind)

		# Wavelength Reduction
		self.linelist_var.set(self.reduction_options.linelist)
		self.mc_sample_var.set(self.reduction_options.sampleamt)
		self.accept_param_var.set(self.reduction_options.accept_param)

		# Polynomial Zero Points
		for attr, var in zip(["offset_zero", "extent_zero", "quad_zero", "cube_zero"], self.zero_point_vars):
			var.set(getattr(self.reduction_options, attr))

		# Stepsizes
		for attr, var in zip(["offset_stepsize", "linear_stepsize", "quad_stepsize", "cube_stepsize"], self.stepsize_vars):
			var.set(getattr(self.reduction_options, attr))

		# Bounds (as string)
		for attr, var in zip(["c_cov", "s_cov", "q_cov", "cub_cov"], self.bounds_vars):
			var.set(str(getattr(self.reduction_options, attr)))

		self.lampfilterwindow_var.set(self.reduction_options.lampfilterwindow)

	def user_review_frames(self, bias_list, flat_list, shifted_flat_list, science_list, complamp_list):
		review_window = tk.Toplevel()
		review_window.title("Review Detected Frames")
		review_window.geometry("900x600")
		review_window.resizable(False, False)

		frame_dict = {
			"Bias": bias_list,
			"Flat": flat_list,
			"Shifted Flat": shifted_flat_list,
			"Science": science_list,
			"Arc": complamp_list
		}

		listboxes = {}
		file_to_frame = {}

		# Create tabbed view
		notebook = ttk.Notebook(review_window)
		notebook.pack(fill="both", expand=True)

		for category, framelist in frame_dict.items():
			tab = ttk.Frame(notebook)
			notebook.add(tab, text=category)

			# Map filenames to frame objects
			frame_map = {os.path.basename(f.filepath or "Created from array"): f for f in framelist}
			file_to_frame[category] = frame_map

			# Scrollable Listbox
			scrollbar = ttk.Scrollbar(tab)
			listbox = tk.Listbox(tab, selectmode=tk.SINGLE, width=100, height=30)
			listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
			scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
			listbox.config(yscrollcommand=scrollbar.set)
			scrollbar.config(command=listbox.yview)
			listboxes[category] = listbox

			# Populate
			for filename in frame_map:
				listbox.insert(tk.END, filename)

			# Capture current category and listbox in closure
			def bind_right_click(lb, cat):
				menu = tk.Menu(lb, tearoff=0)

				def remove_selected():
					try:
						index = lb.curselection()[0]
						filename = lb.get(index)
						lb.delete(index)

						# Remove from underlying list
						frame = file_to_frame[cat].get(filename)
						if frame in frame_dict[cat]:
							frame_dict[cat].remove(frame)
					except IndexError:
						pass

				menu.add_command(label="Remove Frame", command=remove_selected)

				def show_menu(event):
					lb.selection_clear(0, tk.END)
					index = lb.nearest(event.y)
					lb.selection_set(index)
					lb.activate(index)
					menu.tk_popup(event.x_root, event.y_root)

				lb.bind("<Button-3>", show_menu)

			bind_right_click(listbox, category)

		review_window.grab_set()
		review_window.wait_window()

		return frame_dict["Bias"], frame_dict["Flat"], frame_dict["Shifted Flat"], frame_dict["Science"], frame_dict["Arc"]

	def auto_detect(self):
		folder_path = fd.askdirectory(title=f"Select Folder Containing Raw Files")
		
		find_shifted_flats = False
		all_frames = []

		for filepath in os.listdir(folder_path):
			if not ".fits" in filepath:
				continue
			frame = Frame(os.path.join(folder_path, filepath), self.reduction_options)

			### Eliminating Alignment frames (Only necessary for SOAR)
			if "GRATING" in frame.header and "TELESCOP" in frame.header:
				if frame.header["TELESCOP"] == "SOAR 4.1m":
					find_shifted_flats = True
				if frame.header["GRATING"] == "NO_GRATING" and frame.header["TELESCOP"] == "SOAR 4.1m":
					# Alignment Frames should be the only without a grating
					continue

			all_frames.append(frame)
			print(filepath)

		bias_list = []
		flat_list = []
		shifted_flat_list = []
		science_list = []
		complamp_list = []

		# xs = []
		# ys = []
		# ts = []
		for frm in all_frames:
			if "OBJECT" in frm.header:
				if "test" in frm.header["OBJECT"].lower():
					continue

				if "zero" in frm.header["OBJECT"].lower() or "bias" in frm.header["OBJECT"].lower():
					bias_list.append(frm)
					continue

				if "quartz" in frm.header["OBJECT"].lower() or "flat" in frm.header["OBJECT"].lower() or "halogen" in frm.header["OBJECT"].lower():
					flat_list.append(frm)
					continue


			frm.determine_frametype()

			if frm.type == "Arc":
				complamp_list.append(frm)
			elif frm.type == "Science":
				science_list.append(frm)
			# xs.append(x)
			# ys.append(y)
			# ts.append(t)

		colordict = {"Arc": "Blue", "Science": "Red", "Indeterminate": "Grey"}
		# plt.scatter(xs, ys, c=[colordict[t] for t in ts])
		# plt.show()
		# Group frames by grating angle
		if find_shifted_flats:
			angle_groups = defaultdict(list)
			for frm in flat_list:
				angle = frm.header["GRT_TARG"]
				angle_groups[angle].append(frm)

			if len(angle_groups) != 2:
				raise ValueError("Expected exactly 2 distinct grating angles in flat_list")

			# Sort the angles and separate the frames
			sorted_angles = sorted(angle_groups.keys())
			lower_angle, higher_angle = sorted_angles[0], sorted_angles[1]

			shifted_flat_list = angle_groups[higher_angle]
			flat_list = angle_groups[lower_angle]

		filter_by_common_frame_shape(bias_list,flat_list,shifted_flat_list,science_list,complamp_list)

		bias_list, flat_list, shifted_flat_list, science_list, complamp_list = \
			self.user_review_frames(bias_list, flat_list, shifted_flat_list, science_list, complamp_list)


		self.frame_config.biases = BiasList(bias_list, self.reduction_options)
		self.frame_config.flats = FlatList(flat_list, self.reduction_options)
		self.frame_config.shiftedflats = FlatList(shifted_flat_list, self.reduction_options)
		self.frame_config.comparisonframes = ComplampList(complamp_list, self.reduction_options)
		self.frame_config.scienceframes = ScienceList(science_list, self.reduction_options)

		for ft in self.file_buttons.keys():
			existing = getattr(self.frame_config, {
				"Bias Frames": "biases",
				"Flat Frames": "flats",
				"Shifted Flat Frames (Optional)": "shiftedflats",
				"Arc Frames": "comparisonframes",
				"Science Frames": "scienceframes"
			}[ft])
			try:
				self.update_file_count_label(ft, len(existing.filepath_list()))
			except AttributeError:
				self.update_file_count_label(ft, 0)

		self.frame_config.save_to_json()

	def update_file_count_label(self, frame_type, count):
		label = self.file_count_labels[frame_type]
		label.config(text=f"{count} file{'s' if count != 1 else ''}")

		if count > 0:
			label.config(foreground="green")
		elif "Optional" in frame_type:
			label.config(foreground="orange")
		else:
			label.config(foreground="red")

	def select_files(self, frame_type):
		filenames = fd.askopenfilenames(title=f"Select {frame_type}")
		count = len(filenames)

		if frame_type == "Bias Frames":
			self.frame_config.biases = BiasList(filenames, self.reduction_options)

		elif frame_type == "Flat Frames":
			self.frame_config.flats = FlatList(filenames, self.reduction_options)

		elif frame_type == "Shifted Flat Frames (Optional)":
			self.frame_config.shiftedflats = FlatList(filenames, self.reduction_options)

		elif frame_type == "Arc Frames":
			self.frame_config.comparisonframes = ComplampList(filenames, self.reduction_options)

		elif frame_type == "Science Frames":
			self.frame_config.scienceframes = ScienceList(filenames, self.reduction_options)

		print(f"Selected {len(filenames)} {frame_type}:")
		if len(filenames) < 8:
			for filename in filenames:
				print(filename)
		else:
			for filename in filenames[:4]:
				print(filename)
			print("...")
			for filename in filenames[-4:]:
				print(filename)

		self.update_file_count_label(frame_type, count)
		self.frame_config.save_to_json()


	def on_close(self):
		self.reduction_options.save_to_json()
		self.frame_config.save_to_json()
		self.destroy()


class ProgressWindow(tk.Toplevel):
	def __init__(self, master):
		super().__init__(master)
		self.title("Reduction Progress")
		self.geometry("400x100+200+500")
		self.resizable(False, False)

		content_frame = tk.Frame(self)
		content_frame.pack(expand=True)

		# Label and progress bar for overall progress
		self.overall_label = ttk.Label(content_frame, text="Starting reduction...")
		self.overall_label.pack(pady=(10, 0))

		self.overall_bar = ttk.Progressbar(content_frame, orient="horizontal", length=350, mode="determinate")
		self.overall_bar.pack(pady=(0, 10))

		# Label and progress bar for current frame's step
		self.current_label = ttk.Label(content_frame, text="")
		self.current_label.pack()

		self.current_bar = ttk.Progressbar(content_frame, orient="horizontal", length=350, mode="determinate")
		self.current_bar.pack(pady=(0, 10))

		self.protocol("WM_DELETE_WINDOW", self.disable_event)  # Disable manual closing

	def disable_event(self):
		pass  # Prevent user from closing the progress window manually

	def update_overall(self, text, fraction):
		self.overall_label.config(text=text)
		self.overall_bar["value"] = fraction * 100
		self.update_idletasks()

	def update_current(self, text, fraction):
		self.current_label.config(text=text)
		self.current_bar["value"] = fraction * 100
		self.update_idletasks()



class DataReductionGUI(tk.Tk):
	def __init__(self, *args, **kwargs):
		super().__init__()
		self.title("Multi-Instrument Data Input Reducer")
		try:
			imgicon = ImageTk.PhotoImage(Image.open("src/MIDIR_icon.png"))
			self.iconphoto(False, imgicon)
			self._icon_img_ref = imgicon
		except Exception as e:
			print(e)

		self.protocol("WM_DELETE_WINDOW", self.on_close)

		self.geometry("300x150+200+200")
		self.resizable(False, False)

		self.reduction_options = ReductionOptions()
		try:
			self.frame_config = FrameConfiguration(self.reduction_options)
		except AttributeError:
			print("ERROR: Could not load previous Frame config!\nThis can safely be ignored if you switched presets.")
			self.frame_config = FrameConfiguration(self.reduction_options, use_previous=False)

		config_window = ConfigWindow(self)
		config_window.focus()

		# Centered content frame
		content_frame = tk.Frame(self)
		content_frame.pack(expand=True)

		# Welcome label
		welcome_label = ttk.Label(content_frame, text="Welcome to MIDIR!", font=("Segoe UI", 12))
		welcome_label.pack(pady=(0, 10))

		# Start button
		start_button = ttk.Button(
			content_frame, 
			text="Start Reduction", 
			command=self.start_reduction
		)
		start_button.pack()

	def on_close(self):
		self.reduction_options.save_to_json()
		self.frame_config.save_to_json()
		self.destroy()

	def start_reduction(self):
		print(self.reduction_options)
		progress_window = ProgressWindow(self)
		self.update()  # Ensures the window shows immediately
		reduce_data(self.reduction_options, copy.deepcopy(self.frame_config), progress_window)
		


if __name__ == "__main__":
	root = DataReductionGUI()
	sv_ttk.set_theme("light")
	root.mainloop()
