import json
import os
from astropy.coordinates import EarthLocation

class ReductionOptions:
	def __init__(self):
		if os.path.exists(".previousoptions.json"):
			print("Reloading Options...")
			self.load_from_json()
		else:
			print("No previous configuration found, using a default one.")
			self.output_dict = "./output"
			self.outputfile_path = "./reduced_spectra.csv"
			self.linelist = "./linelists/SOAR_930_M1_FeAr.txt"
			self.telescope_location = "Cerro Pachon"
			self.header_hdu_ind = 0
			self.image_hdu_ind = 0
			self.comparison_match_mode = "time"
			self.comparison_match_value = 20
			self.comparison_match_position = "after"
			self.debugimages = True
			self.offset_zero = "Auto"
			self.extent_zero = 1725
			self.quad_zero = -7e-6
			self.cube_zero = 0
			self.offset_stepsize =  0.1
			self.linear_stepsize =  0.0001
			self.quad_stepsize = 1.e-7
			self.cube_stepsize = 1.e-10
			self.c_cov = 150
			self.s_cov = 0.1
			self.q_cov = 5e-5
			self.cub_cov = 3e-8
			self.sampleamt = 2500000
			self.accept_param = 1.05
			self.skyfluxsep = 100
			self.cosmicrejection = True
			self.lampfilterwindow = 50
			self.use_boxcut = False
			self.multitrace = False
			self.rotationangle = 0
			self.manual_crop = False
			self.x_lo = 0
			self.x_hi = 0
			self.y_lo = 0
			self.y_hi = 0

	def load_from_json(self, json_path=".previousoptions.json"):
		with open(json_path, "r") as jsonfile:
			config = json.load(jsonfile)

		self.output_dict = config["output_dict"]
		self.outputfile_path = config["outputfile_path"]
		self.linelist = config["linelist"]
		self.telescope_location = config["telescope_location"]
		self.header_hdu_ind = config["header_hdu_ind"]
		self.image_hdu_ind = config["image_hdu_ind"]
		self.comparison_match_mode = config["comparison_match_mode"]
		self.comparison_match_value = config["comparison_match_value"]
		self.comparison_match_position = config["comparison_match_position"]
		self.debugimages = config["debugimages"]
		self.offset_zero = config["offset_zero"]
		self.extent_zero = config["extent_zero"]
		self.quad_zero = config["quad_zero"]
		self.cube_zero = config["cube_zero"]
		self.offset_stepsize = config["offset_stepsize"]
		self.linear_stepsize = config["linear_stepsize"]
		self.quad_stepsize = config["quad_stepsize"]
		self.cube_stepsize = config["cube_stepsize"]
		self.c_cov = config["c_cov"]
		self.s_cov = config["s_cov"]
		self.q_cov = config["q_cov"]
		self.cub_cov = config["cub_cov"]
		self.sampleamt = config["sampleamt"]
		self.accept_param = config["accept_param"]
		self.skyfluxsep = config["skyfluxsep"]
		self.cosmicrejection = config["cosmicrejection"]
		self.lampfilterwindow = config["lampfilterwindow"]
		self.use_boxcut = config["use_boxcut"]
		self.multitrace = config["multitrace"]
		self.rotationangle = config["rotationangle"]
		self.manual_crop = config["manual_crop"]
		self.x_lo = config["x_lo"]
		self.x_hi = config["x_hi"]
		self.y_lo = config["y_lo"]
		self.y_hi = config["y_hi"]


	def save_to_json(self, json_path=".previousoptions.json"):
		json_dict = {
			"output_dict": self.output_dict,
			"outputfile_path": self.outputfile_path,
			"linelist": self.linelist,
			"telescope_location": self.telescope_location,
			"header_hdu_ind": self.header_hdu_ind,
			"image_hdu_ind": self.image_hdu_ind,
			"comparison_match_mode": self.comparison_match_mode,
			"comparison_match_value": self.comparison_match_value,
			"comparison_match_position": self.comparison_match_position,
			"debugimages": self.debugimages,
			"offset_zero": self.offset_zero,
			"extent_zero": self.extent_zero,
			"quad_zero": self.quad_zero,
			"cube_zero": self.cube_zero,
			"offset_stepsize": self.offset_stepsize,
			"linear_stepsize": self.linear_stepsize,
			"quad_stepsize": self.quad_stepsize,
			"cube_stepsize": self.cube_stepsize,
			"c_cov": self.c_cov,
			"s_cov": self.s_cov,
			"q_cov": self.q_cov,
			"cub_cov": self.cub_cov,
			"sampleamt": self.sampleamt,
			"accept_param": self.accept_param,
			"skyfluxsep": self.skyfluxsep,
			"cosmicrejection": self.cosmicrejection,
			"lampfilterwindow": self.lampfilterwindow,
			"use_boxcut": self.use_boxcut,
			"multitrace": self.multitrace,
			"rotationangle": self.rotationangle,
			"manual_crop": self.manual_crop,
			"x_lo": self.x_lo,
			"x_hi": self.x_hi,
			"y_lo": self.y_lo,
			"y_hi": self.y_hi,

		}	
		with open(json_path, "w") as jsonfile:
			json.dump(json_dict, jsonfile)

	def populate_earth_location(self):
		self.earth_location = EarthLocation.of_site(self.telescope_location)

	def __str__(self):
		return f"""#### REDUCTION OPTIONS SELECTED: ####
 output_dict: {self.output_dict}
 outputfile_path: {self.outputfile_path}
 linelist: {self.linelist}
 telescope_location: {self.telescope_location},
 header_hdu_ind: {self.header_hdu_ind} 
 image_hdu_ind: {self.image_hdu_ind}
 comparison_match_mode: {self.comparison_match_mode},
 comparison_match_value: {self.comparison_match_value},
 comparison_match_position: {self.comparison_match_position},
 debugimages: {self.debugimages}
 offset_zero: {self.offset_zero}
 extent_zero: {self.extent_zero}
 quad_zero: {self.quad_zero} 
 cube_zero: {self.cube_zero} 
 offset_stepsize: {self.offset_stepsize}
 linear_stepsize: {self.linear_stepsize}
 quad_stepsize: {self.quad_stepsize}
 cube_stepsize: {self.cube_stepsize}
 c_cov: {self.c_cov}
 s_cov: {self.s_cov}
 q_cov: {self.q_cov}
 cub_cov: {self.cub_cov}
 sampleamt: {self.sampleamt}
 accept_param: {self.accept_param}
 skyfluxsep: {self.skyfluxsep}
 cosmicrejection: {self.cosmicrejection}
 lampfilterwindow: {self.lampfilterwindow}
 use_boxcut: {self.use_boxcut}
 multitrace: {self.multitrace}
 rotationangle: {self.rotationangle}
 manual_crop: {self.manual_crop}
 x_lo: {self.x_lo},
 x_hi: {self.x_hi},
 y_lo: {self.y_lo},
 y_hi: {self.y_hi},
"""

