import numpy as np
from pathlib import Path
import io
import gzip
import zipfile
import subprocess
import shutil
import os
import json
import copy
import cv2
from scipy.ndimage import sobel, median_filter, gaussian_filter, minimum_filter, maximum_filter
from astropy.io import fits
from astropy.time import Time, TimeDelta
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt


def is_compressed_by_extension(filepath):
	compressed_exts = {'.gz', '.zip', '.bz2', '.xz', '.7z', '.rar', '.z'}
	extension = Path(filepath).suffix.lower()

	if extension in compressed_exts:
		return extension
	else:
		return False

def merge_dicts(*dicts):
	merged_dict = {}
	count_dict = {}

	for d in dicts:
		for key, value in d.items():
			if isinstance(value, (int, float)):
				merged_dict[key] = merged_dict.get(key, 0) + value
				count_dict[key] = count_dict.get(key, 0) + 1
			elif key not in merged_dict:
				merged_dict[key] = value

	for key in merged_dict:
		if isinstance(merged_dict[key], (int, float)):
			merged_dict[key] /= count_dict[key]

	return merged_dict


def get_header_value(header, keys):
	for key in keys:
		if key in header:
			return header[key]
	return None  # or raise an error if mandatory


class Frame:
	def __init__(self, fits_filepath_or_array, reduction_options):
		if isinstance(fits_filepath_or_array, np.ndarray):
			self.reduction_options = reduction_options
			self.filepath = None

			self.header = {}
			self.data = fits_filepath_or_array

			self.obstime = None
			self.type = None
			self.star_info = None

		else:
			fits_filepath = fits_filepath_or_array
			self.reduction_options = reduction_options
			self.filepath = fits_filepath

			hdul = self._load_fits()

			try:
				self.header = dict(hdul[self.reduction_options.header_hdu_ind].header)
				self.data = hdul[self.reduction_options.image_hdu_ind].data.astype(np.float64)
			except IndexError:
				print("ERROR: Wrong HDU indices! You need to reload your frames!")
				self.header = None
				self.data = None

			try:
				self.obstime = Time(self.header["DATE-OBS"], format='isot', scale='utc')
				self.obstime += TimeDelta(self.header["EXPTIME"], format='sec')/2
			except (KeyError, TypeError):
				self.obstime = None
			self.type = None
			self.star_info = None

	def rotate(self):
		angle = self.reduction_options.rotationangle
		if angle == 0:
			return
		elif angle in (90, 180, 270):
			k = angle // 90  # number of 90° rotations
			self.data = np.rot90(self.data, k=k)
		else:
			raise ValueError("Rotation angle must be 0, 90, 180, or 270 degrees")

	def reverse_rotate(self):
		angle = self.reduction_options.rotationangle
		if angle == 0:
			return
		elif angle in (90, 180, 270):
			k = angle // 90  # number of 90° rotations
			self.data = np.rot90(self.data, k=-k)  # reverse direction
		else:
			raise ValueError("Rotation angle must be 0, 90, 180, or 270 degrees")

	def crop(self, crop_bounds):
		self.data = self.data[crop_bounds[1][0]:crop_bounds[1][1], crop_bounds[0][0]:crop_bounds[0][1]]


	def detect_spectral_area(self):
		# Edge detection
		minumum_truncation = 5

		image_data = self.data.astype(np.float64)[minumum_truncation:-minumum_truncation,
					 minumum_truncation:-minumum_truncation]

		x_img = sobel(image_data, axis=0, mode="nearest")
		y_img = sobel(image_data, axis=1, mode="nearest")

		edge_detection = np.sqrt(x_img ** 2 + y_img ** 2)
		edge_detection *= 1 / np.max(edge_detection)
		edge_detection[edge_detection > 0.075] = 1
		edge_detection[edge_detection < 1] = 0

		edge_detection = (255 * edge_detection / edge_detection.max()).astype(np.uint8)

		lines = cv2.HoughLinesP(edge_detection, 1, np.pi / 180, 50, None, 500, 0)

		x_intercepts = []
		y_intercepts = []
		# Loop through the detected lines
		for line in lines:
			x1, y1, x2, y2 = line[0]
			cv2.line(edge_detection, (x1, y1), (x2, y2), (255,), 1)

			if y2 == y1:
				y_intercepts.append(y1)
				continue

			if x2 == x1:
				x_intercepts.append(x1)
				continue

			m = (y2 - y1) / (x2 - x1)
			y_intercept = y1 - m * x1
			x_intercept = -y_intercept / m if m != 0 else None

			if x_intercept is not None:
				if 0 < x_intercept < edge_detection.shape[0]:
					y_intercepts.append(x_intercept + minumum_truncation)

			if 0 < y_intercept < edge_detection.shape[1]:
				x_intercepts.append(y_intercept + minumum_truncation)

		x_intercepts = np.array(x_intercepts)
		y_intercepts = np.array(y_intercepts)

		u_x = x_intercepts[x_intercepts > edge_detection.shape[1] / 2]
		l_x = x_intercepts[x_intercepts < edge_detection.shape[1] / 2]

		u_y = y_intercepts[y_intercepts > edge_detection.shape[0] / 2]
		l_y = y_intercepts[y_intercepts < edge_detection.shape[0] / 2]

		if len(u_x) == 0:
			u_x = image_data.shape[1] - 5
		else:
			u_x = np.min(u_x) - 5
		if len(l_x) == 0:
			l_x = 5
		else:
			l_x = np.max(l_x) + 5
		if len(u_y) == 0:
			u_y = image_data.shape[0] - 5
		else:
			u_y = np.min(u_y) - 5
		if len(l_y) == 0:
			l_y = 5
		else:
			l_y = np.max(l_y) + 5

		return [(l_x, u_x), (l_y, u_y)]


	def _load_fits(self):
		try:
			suffix = is_compressed_by_extension(self.filepath)

			if suffix == '.gz':
				with gzip.open(self.filepath, 'rb') as f:
					return fits.open(io.BytesIO(f.read()))

			elif suffix == '.z':
				if not shutil.which('uncompress'):
					raise RuntimeError("'uncompress' command not found for .Z file.")
				process = subprocess.Popen(['uncompress', '-c', self.filepath], stdout=subprocess.PIPE)
				decompressed_data, _ = process.communicate()
				return fits.open(io.BytesIO(decompressed_data))

			elif suffix == '.zip':
				with zipfile.ZipFile(self.filepath, 'r') as zf:
					fits_files = [name for name in zf.namelist() if name.lower().endswith('.fits')]
					if not fits_files:
						raise FileNotFoundError("No .fits file found inside the ZIP archive.")
					with zf.open(fits_files[0]) as f:
						return fits.open(io.BytesIO(f.read()))

			else:
				return fits.open(self.filepath)

		except FileNotFoundError:
			print(f"ERROR: An input file does not exist: {self.filepath}")
			raise  

		except IOError:
			print(f"ERROR: Could not read .fits file! Is the file compressed? File: {self.filepath}")
			raise

	def determine_frametype(self):
		current_frm = copy.copy(self)
		current_frm.rotate()
		if not self.reduction_options.manual_crop:
			x_lo = int(current_frm.data.shape[1]*0.15)
			x_hi = int(current_frm.data.shape[1]*0.85)
			y_lo = int(current_frm.data.shape[0]*0.15)
			y_hi = int(current_frm.data.shape[0]*0.85)
			current_frm.crop(([x_lo, x_hi],[y_lo, y_hi]))
		else:
			crop = [[self.reduction_options.x_lo, current_frm.data.shape[1]-self.reduction_options.x_hi], [self.reduction_options.y_hi, current_frm.data.shape[0]-self.reduction_options.y_lo]]
			current_frm.crop(crop)

		along_x = median_filter(np.sum(current_frm.data, axis=0).astype(np.float32), 15)
		along_x -= np.min(along_x)
		along_x /= np.max(along_x)

		along_y = median_filter(np.sum(current_frm.data, axis=1).astype(np.float32), 15)
		if "TELESCOP" in self.header:
			if self.header["TELESCOP"] == "NOT":
				along_y -= minimum_filter(along_y, 25)

		along_y -= np.min(along_y)
		along_y /= np.max(along_y)

		x_peak_median_diff = np.nanmax(along_x)-np.nanmedian(along_x)
		y_peak_median_diff = np.nanmax(along_y)-np.nanmedian(along_y)

		if x_peak_median_diff < .5 and y_peak_median_diff < .5:
			self.type = "Indeterminate"
		if y_peak_median_diff < 0.8 and x_peak_median_diff > 0.8:
			self.type = "Arc"
		else:
			self.type = "Science"


		
		# if x_peak_median_diff > .8:
		# 	plt.imshow(current_frm.data, cmap="Greys_r")
		# 	plt.show()
		# 	plt.plot(along_x)
		# 	plt.plot(along_y)
		# 	plt.show()

		# if y_peak_median_diff < 0.8 and x_peak_median_diff > 0.8:
		# 	return x_peak_median_diff, y_peak_median_diff, "Arc"
		# else:
		# 	return x_peak_median_diff, y_peak_median_diff, "Science"
			
		
		

	def get_star_info(self):
		rename_dict = {
			"Name": "name",
			"RA_ICRS": "ra",
			"DE_ICRS": "dec",
			"GaiaEDR3": "source_id",
			"SpClass": "SPEC_CLASS",
			"BP-RP": "bp_rp",
			"GGAIA": "gmag",
			"Gmag": "gmag",
			"pmRA": "pmra",
			"pmRAGAIA": "pmra",
			"pmDE": "pmdec",
			"pmDEGAIA": "pmdec",
			"e_pmRAGAIA": "pmra_error",
			"e_pmRA": "pmra_error",
			"e_pmDEGAIA": "pmdec_error",
			"e_pmDE": "pmdec_error",
			"Plx": "parallax",
			"e_Plx": "parallax_error",
		}

		ra_keys = ["OBJRA", "RA"]
		dec_keys = ["OBJDEC", "DEC"]

		ra = get_header_value(self.header, ra_keys)
		dec = get_header_value(self.header, dec_keys)

		sky_coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
		# Define the VizieR catalog ID
		catalog_id = "J/A+A/662/A40"

		# Query the VizieR catalog
		vizier = Vizier(columns=['all'], row_limit=1)
		sinfo = vizier.query_region(sky_coord, radius=30 * u.arcsec, catalog=catalog_id)

		if len(sinfo) == 0:
			# Define the coordinates
			coord = SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg))
			star_found = False
			tries_to_find_star = 0

			while not star_found:
				tries_to_find_star += 1
				try:
					# Query Gaia DR3
					width = u.Quantity(10*tries_to_find_star, u.arcsecond)
					result = Gaia.query_object(coordinate=coord, radius=width)
					star = result[0]
					sinfo = {}
					sinfo["name"] = f"Gaia DR3 {star['source_id']}"
					sinfo["source_id"] = star['source_id']
					sinfo["ra"] = star['ra']
					sinfo["dec"] = star['dec']
			
					# Fill in other values if they exist, otherwise set to "N/A"
					sinfo["file"] = self.filepath
					sinfo["SPEC_CLASS"] = "unknown"
					sinfo["bp_rp"] = star['bp_rp'] if 'bp_rp' in star.columns else "N/A"
					sinfo["gmag"] = star['phot_g_mean_mag'] if 'phot_g_mean_mag' in star.columns else "N/A"
					sinfo["nspec"] = 1
					sinfo["pmra"] = star['pmra'] if 'pmra' in star.columns else "N/A"
					sinfo["pmra_error"] = star['pmra_error'] if 'pmra_error' in star.columns else "N/A"
					sinfo["pmdec"] = star['pmdec'] if 'pmdec' in star.columns else "N/A"
					sinfo["pmdec_error"] = star['pmdec_error'] if 'pmdec_error' in star.columns else "N/A"
					sinfo["parallax"] = star['parallax'] if 'parallax' in star.columns else "N/A"
					sinfo["parallax_error"] = star['parallax_error'] if 'parallax_error' in star.columns else "N/A"
					star_found = True
				except (KeyError, IndexError) as e:
					print(f"Star not found in Gaia catalogues after {tries_to_find_star} iterations, trying with bigger radius...")

			if tries_to_find_star != 1:
				found_coordinates = SkyCoord(sinfo["ra"], sinfo["dec"], unit=(u.deg, u.deg))
				print(f"Star from file {self.filepath} was identified as {sinfo['name']} after {tries_to_find_star} iterations.\nDistance from header coordinates is {coord.separation(found_coordinates).arcsecond}\"")

		elif len(sinfo) > 0:
			if len(sinfo) == 1:
				sinfo = sinfo[0].to_pandas().to_dict(orient='records')[0]
				sinfo["name"] = "-"
			else:
				sinfo = sinfo[1].to_pandas().to_dict(orient='records')[0]
				sinfo["name"] = "-"
				sinfo["bp_rp"] = sinfo["BPGAIA"] - sinfo["RPGAIA"]
		else:
			sinfo = sinfo[0].to_pandas().to_dict(orient='records')[0]
			sinfo["name"] = "-"

		for a, b in rename_dict.items():
			if a in sinfo.keys():
				sinfo[b] = sinfo.pop(a)

		if "SPEC_CLASS" in sinfo.keys():
			if sinfo["SPEC_CLASS"] == "":
				sinfo["SPEC_CLASS"] = "unknown"
		else:
			sinfo["SPEC_CLASS"] = "unknown"
		sinfo["nspec"] = 1
		sinfo["file"] = self.filepath.split("/")[-1]

		self.star_info = sinfo

	def apply_corrections(self, master_flat, frame_config):
		master_bias = frame_config.biases.master_frame.data

		self.data[self.data < master_bias] = 0
		self.data[self.data >= master_bias] = (self.data - master_bias)[self.data >= master_bias]


		if 0.0 in master_flat.data:
			print("WARNING: Found invalid value in master flat frame, continuing anyways...")
			master_flat.data[master_flat.data == 0] = 1.0

		if 0.0 in frame_config.flats.master_frame.data:
			print("WARNING: Found invalid value in master flat frame, continuing anyways...")
			frame_config.flats.master_frame.data[frame_config.flats.master_frame.data == 0] = 1.0

		self.data = np.floor_divide(self.data.astype(np.float64), master_flat.data)
		self.data = np.floor_divide(self.data.astype(np.float64), frame_config.flats.master_frame.data)
		self.data *= 65535 / self.data.max()
		self.data = self.data.astype(np.uint16)



class FrameList:
	def __init__(self, list_of_filepaths, reduction_options):
		if len(list_of_filepaths) != 0:
			if isinstance(list_of_filepaths[0], str):
				self.reduction_options = reduction_options
				self.frames = []
				for filepath in list_of_filepaths:
					try:
						self.frames.append(Frame(filepath, reduction_options))
					except FileNotFoundError:
						print(f"WARNING: Could not load file {filepath}, skipping it...")
				self.master_frame = None
			elif isinstance(list_of_filepaths[0], Frame):
				self.reduction_options = reduction_options
				self.frames = list_of_filepaths
				self.master_frame = None

		else:
			self.reduction_options = reduction_options
			self.frames = []
			self.master_frame = None

	def filepath_list(self):
		return [frame.filepath for frame in self.frames]

	def create_master_image(self, master_bias=None, master_continuum=None, return_header=False):
		# Stack all image data into a 3D array: shape (N_frames, H, W)
		image_stack = np.stack([frame.data.astype(np.uint32) for frame in self.frames])
		
		# Collect all headers
		headers = [dict(frame.header) for frame in self.frames]

		if master_bias is not None:
			bias_data = master_bias.data.astype(np.uint32)
			# Subtract bias where image > bias, else zero
			image_stack = np.where(image_stack > bias_data, image_stack - bias_data, 0)

		# Compute mean across all frames
		master = np.mean(image_stack, axis=0)

		if master_continuum is not None:
			master = master.astype(np.float64)
			if 0.0 in master_continuum.data:
				print("WARNING: Found invalid value in master flat frame, continuing anyways...")
				master_continuum.data[master_continuum.data == 0] = 1.0
			master /= master_continuum.data
			master /= master.max()
			master *= 65535

		master = master.astype(np.uint16)
		master[np.isnan(master)] = 0.0

		self.master_frame = Frame(master, self.reduction_options)

		if not return_header:
			return master, None
		else:
			master_header = merge_dicts(*headers)
			return master, master_header


	def __iter__(self):
		return iter(self.frames)

class ScienceList(FrameList):
	pass



class ComplampList(FrameList):
	def __init__(self, list_of_frames_or_filepaths, reduction_options):
		self.reduction_options = reduction_options
		if len(list_of_frames_or_filepaths) != 0:
			if not isinstance(list_of_frames_or_filepaths[0], Frame):
				self.frames = []
				for filepath in list_of_frames_or_filepaths:
					try:
						self.frames.append(Frame(filepath, reduction_options))
					except FileNotFoundError:
						print(f"WARNING: Could not load file {filepath}, skipping it...")
					self.master_frame = None
			else:
				self.frames = list_of_frames_or_filepaths
				self.master_frame = None
		else:
			self.frames = []
			self.master_frame = None
		self.compparams = None

	def _get_complamp_key(self):
		"""Return a unique key for the current set of comp lamps, order-independent."""
		complamp_names = sorted(frame.filepath for frame in self.frames)
		return "|".join(complamp_names)

	def get_compparams(self, traceparams):
		"""Try to load polynomial parameters for current comp lamp set. Set self.compparams."""

		traceparams = [str(t) for t in traceparams]

		paramfile = ".previoussolutions.json"
		key = self._get_complamp_key() + "|".join(traceparams)

		if os.path.isfile(paramfile):
			with open(paramfile, "r") as f:
				try:
					solutions = json.load(f)
				except json.JSONDecodeError:
					solutions = {}
		else:
			solutions = {}

		self.compparams = solutions.get(key, None)
		return self.compparams

	def save_compparams(self, traceparams):
		"""Store the current self.compparams under the current set of comp lamps."""

		traceparams = [str(t) for t in traceparams]

		if self.compparams is None:
			raise ValueError("compparams is None — nothing to save.")

		paramfile = ".previoussolutions.json"
		key = self._get_complamp_key() + "|".join(traceparams)

		# Load or create solutions
		if os.path.isfile(paramfile):
			with open(paramfile, "r") as f:
				try:
					solutions = json.load(f)
				except json.JSONDecodeError:
					solutions = {}
		else:
			solutions = {}

		solutions[key] = self.compparams

		with open(paramfile, "w") as f:
			json.dump(solutions, f, indent=4)

	def create_master_bias():
		pass


class BiasList(FrameList):
	pass


class FlatList(FrameList):
	def create_master_flat(self, second_image_list=None, master_bias=None, bounds=None):
		def apply_bias_correction(images, bias):
			bias_data = bias.data.astype(np.float64)
			return np.where(images >= bias_data, images - bias_data, 0)

		# Stack all primary images
		image_stack = np.stack([img.data.astype(np.float64) for img in self.frames])
		if master_bias is not None:
			image_stack = apply_bias_correction(image_stack, master_bias)
		
		master = np.sum(image_stack, axis=0)

		# Handle second image list if present
		if second_image_list is not None:
			image_stack2 = np.stack([img.data.astype(np.float64) for img in second_image_list.frames])
			if master_bias is not None:
				image_stack2 = apply_bias_correction(image_stack2, master_bias)
			master2 = np.sum(image_stack2, axis=0)
		else:
			master2 = None

		# Normalize master
		master /= master.max()

		if second_image_list is not None:
			master2 /= master2.max()

			if bounds is not None:
				# Correct for Littrow ghost
				center_diff = np.median(master) - np.median(master2)
				master2 += center_diff

			master = np.minimum(master, master2)

		# Median smoothing
		smooth_master = median_filter(master, 25)
		smooth_master[np.isnan(smooth_master)] = 1.0

		self.master_frame = Frame(smooth_master, self.reduction_options)

		if bounds is not None:
			master /= smooth_master
			smooth_master /= smooth_master.max()
			return master, smooth_master
		else:
			master /= master.max()
			return master



class FrameConfiguration:
	def __init__(self, reduction_options, use_previous=True):
		self.reduction_options = reduction_options
		if os.path.exists(".previousframeconfig.json") and use_previous:
			print("Reloading Frame Selection...")
			self.load_from_json()

		else:
			print("First time running, creating basic configuration...")

			self.flats = FrameList([], reduction_options)
			self.shiftedflats = FrameList([], reduction_options)
			self.biases = FrameList([], reduction_options)
			self.scienceframes = FrameList([], reduction_options)
			self.comparisonframes = FrameList([], reduction_options)

			basic_json_dict = {
				"flats": [],
				"shiftedflats": [],
				"biases": [],
				"scienceframes": [],
				"comparisonframes": [],
			}

			with open(".previousframeconfig.json", "w") as jsonfile:
				json.dump(basic_json_dict, jsonfile)

	def load_from_json(self):
		with open("./.previousframeconfig.json", "r") as jsonfile:
			previous_config = json.load(jsonfile)

		self.flats = FlatList(previous_config["flats"], self.reduction_options)
		self.shiftedflats = FlatList(previous_config["shiftedflats"], self.reduction_options)
		self.biases = BiasList(previous_config["biases"], self.reduction_options)
		self.scienceframes = ScienceList(previous_config["scienceframes"], self.reduction_options)
		self.comparisonframes = ComplampList(previous_config["comparisonframes"], self.reduction_options)


	def save_to_json(self):
		json_dict = {
			"flats": self.flats.filepath_list(),
			"shiftedflats": self.shiftedflats.filepath_list(),
			"biases": self.biases.filepath_list(),
			"scienceframes": self.scienceframes.filepath_list(),
			"comparisonframes": self.comparisonframes.filepath_list()
		}	
		with open(".previousframeconfig.json", "w") as jsonfile:
			json.dump(json_dict, jsonfile)


	def __str__(self):
		returnstr = ""
		for framename, filenames in zip(["Flat Frames", "Shifted Flat Frames", "Bias Frames", "Science Frames", "Arc Lamp Frames"], [self.flats, self.shiftedflats, self.biases, self.scienceframes, self.comparisonframes]):
			returnstr += f"\nSelected {framename}:\n"
			filenames = filenames.filepath_list()
			if len(filenames) < 8:
				for filename in filenames:
					returnstr += filename + "\n"
			else:
				for filename in filenames[:4]:
					returnstr += filename + "\n"
				returnstr += "..." + "\n"
				for filename in filenames[-4:]:
					returnstr += filename + "\n"
		return returnstr
