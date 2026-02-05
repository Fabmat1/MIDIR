import json
import os
from astropy.coordinates import EarthLocation


class ReductionOptions:
	# Define defaults as class attribute for single source of truth
	DEFAULTS = {
		"output_dict": "./output",
		"outputfile_path": "./reduced_spectra.csv",
		"linelist": "./linelists/SOAR_930_M1_FeAr.txt",
		"telescope_location": "Cerro Pachon",
		"header_hdu_ind": 0,
		"image_hdu_ind": 0,
		"comparison_match_mode": "time",
		"comparison_match_value": 20,
		"comparison_match_position": "after",
		"debugimages": True,
		"offset_zero": "Auto",
		"extent_zero": 1725,
		"quad_zero": -7e-6,
		"cube_zero": 0,
		"offset_stepsize": 0.1,
		"linear_stepsize": 0.0001,
		"quad_stepsize": 1.e-7,
		"cube_stepsize": 1.e-10,
		"c_cov": 150,
		"s_cov": 0.1,
		"q_cov": 5e-5,
		"cub_cov": 3e-8,
		"sampleamt": 2500000,
		"accept_param": 1.05,
		"skyfluxsep": 100,
		"cosmicrejection": True,
		"lampfilterwindow": 50,
		"use_boxcut": False,
		"multitrace": False,
		"low_quality_mode": False,
		"rotationangle": 0,
		"manual_crop": False,
		"x_lo": 0,
		"x_hi": 0,
		"y_lo": 0,
		"y_hi": 0,
	}

	# Fields that accept multiple types
	MULTI_TYPE_FIELDS = {
		"offset_zero": (str, int, float),
	}

	def __init__(self):
		# Always initialize with defaults first
		self._set_defaults()
		
		# Then try to load from previous config
		if os.path.exists(".previousoptions.json"):
			print("Reloading Options...")
			try:
				self.load_from_json()
			except (json.JSONDecodeError, IOError) as e:
				print(f"Warning: Could not load previous options ({e}), using defaults.")
		else:
			print("No previous configuration found, using a default one.")

	def _set_defaults(self):
		"""Set all attributes to their default values"""
		for key, value in self.DEFAULTS.items():
			setattr(self, key, value)

	def _is_valid_type(self, key, value, default_value):
		"""Check if value has a valid type for the given key"""
		if default_value is None or value is None:
			return True
		
		# Check if this field accepts multiple types
		if key in self.MULTI_TYPE_FIELDS:
			return isinstance(value, self.MULTI_TYPE_FIELDS[key])
		
		expected_type = type(default_value)
		
		# Special handling for numeric types (int/float interchangeable)
		if expected_type in (int, float) and isinstance(value, (int, float)):
			return True
		
		return isinstance(value, expected_type)

	def load_from_json(self, json_path=".previousoptions.json"):
		try:
			with open(json_path, "r") as jsonfile:
				config = json.load(jsonfile)
		except (json.JSONDecodeError, IOError) as e:
			print(f"Warning: Could not parse JSON file ({e}), keeping current values.")
			return

		# Load each value with fallback to default if key is missing or value is invalid
		for key, default_value in self.DEFAULTS.items():
			try:
				value = config.get(key, default_value)
				# Type checking: ensure loaded value matches expected type
				if not self._is_valid_type(key, value, default_value):
					if key in self.MULTI_TYPE_FIELDS:
						expected = " or ".join(t.__name__ for t in self.MULTI_TYPE_FIELDS[key])
					else:
						expected = type(default_value).__name__
					print(f"Warning: Invalid type for '{key}', expected {expected}, using default.")
					value = default_value
				setattr(self, key, value)
			except Exception as e:
				print(f"Warning: Error loading '{key}' ({e}), using default.")
				setattr(self, key, default_value)

	def save_to_json(self, json_path=".previousoptions.json"):
		json_dict = {}
		for key in self.DEFAULTS.keys():
			try:
				json_dict[key] = getattr(self, key, self.DEFAULTS[key])
			except AttributeError:
				json_dict[key] = self.DEFAULTS[key]
		
		try:
			with open(json_path, "w") as jsonfile:
				json.dump(json_dict, jsonfile, indent=2)
		except IOError as e:
			print(f"Warning: Could not save options to {json_path}: {e}")

	def populate_earth_location(self):
		try:
			self.earth_location = EarthLocation.of_site(self.telescope_location)
		except Exception as e:
			print(f"Warning: Could not resolve telescope location '{self.telescope_location}': {e}")
			# Fall back to Cerro Pachon
			self.earth_location = EarthLocation.of_site("Cerro Pachon")

	def __str__(self):
		lines = ["#### REDUCTION OPTIONS SELECTED: ####"]
		for key in self.DEFAULTS.keys():
			value = getattr(self, key, self.DEFAULTS[key])
			lines.append(f" {key}: {value}")
		return "\n".join(lines)