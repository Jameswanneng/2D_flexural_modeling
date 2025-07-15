#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Key Features:
# - Interpolates sediment thickness along a geological profile
# - Builds various mountain belt geometries (triangle, rectangle, trapezoid, etc.)
# - Computes flexural deflection due to tectonic and sedimentary loads
# - Visualization of modeling results


#    (c) China University of Geosciences (Beijing)
#            ALL RIGHTS RESERVED
#        Authors: Neng Wan, Xueyan Li


import math
import sys
import os
import numpy as np
from scipy import interpolate
import matplotlib
from matplotlib import pyplot as plt

# Configure Matplotlib for publication-quality output
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'Arial'  # 'Helvetica'


def sediment_thickness_interpolation_along_profile(well_location_list, original_sedth, displacement, displacement2, load_width, load_width2,
                                                   section_length=2000, interp_method='linear', extra_thickness=0, extra_thickness2 = 0,
												   erosion_recovery_list=None, taper_distance=800,
                                                   include_wells_involved_in_thrust_belt=0, load_shape='triangle'):
	"""
	   Interpolate sediment thickness along the modeling section.

	   Parameters:
	       well_location_list (list): Well positions in km.
	       original_sedth (np.array): Sediment thickness at each well.
	       displacement (float): Displacement of main mountain belt.
	       displacement2 (float): Displacement of secondary mountain belt.
	       load_width (float): Width of primary mountain load.
	       load_width2 (float): Width of secondary mountain load.
	       section_length (int): Total length of the modeled section in km.
	       interp_method (str): Interpolation method ('linear', 'nearest', etc.).
	       extra_thickness (float): Extra sediment added near load front.
	       extra_thickness2 (float): Additional sediment for secondary load.
	       erosion_recovery_list (list): List of sediment recovery values for wells.
	       taper_distance (float): Distance to taper off sediment thickness.
	       include_wells_involved_in_thrust_belt (bool): Whether to zero out certain regions.
	       load_shape (str): Shape of the mountain load (e.g., triangle).

	   Returns:
	       tuple: interpolated sediment thickness array and updated well positions.
	   """
	if load_width + displacement < 0:
		print("Caution!! The absolute value of load width should be greater than displacement. Check it out")
		sys.exit()

	sedth = original_sedth.copy()

	if erosion_recovery_list:
		# we do erosion recovery for wells adjacent to thrust front.
		for i in range(0, len(erosion_recovery_list), 1):
			sedth[i] = sedth[i] + erosion_recovery_list[i]

	# Insert a point at mountain belt load front
	if displacement2 > 0:
		well_location_list = np.append(well_location_list, displacement2 + well_location_list[-1])
		sedth = np.append(sedth, sedth[-1] + extra_thickness2)

	well_location_list = np.append(well_location_list, load_width2 + well_location_list[-1])
	sedth = np.append(sedth, 0)

	# Shift all well positions based on displacement and load width
	well_location_list = np.array(well_location_list) + displacement + load_width

	# Add a control point at taper distance
	well_location_list = np.append(well_location_list,
	                               taper_distance + displacement + load_width + well_location_list[-1])
	well_location_list = np.append(well_location_list, section_length - 1)
	well_location_list = np.insert(well_location_list, 0, 0)

	# Extend sediment thickness array accordingly
	sedth = np.append(sedth, 0)
	sedth = np.append(sedth, 0)
	sedth = np.insert(sedth, 0, 0)

	# Add additional control point for positive displacement
	if displacement > 0:
		# insert a point at mountain belt load front
		well_location_list = np.insert(well_location_list, 1, load_width)
		# sediment thickness recovery
		sedth = np.insert(sedth, 1, sedth[1] + extra_thickness)
		if load_shape == 'triangle':
			pass
		points_to_be_removed = [0, 1]  # remove added functional points
	else:
		points_to_be_removed = [0, -2, -1]  # remove added functional points

	# Interpolate sediment thickness along section
	func = interpolate.interp1d(well_location_list, sedth, kind=interp_method)
	x_new = np.arange(0, section_length, 1)
	y_new = func(x_new)

	# Zero out certain segments based on thrust belt inclusion flag
	if include_wells_involved_in_thrust_belt == 1:
		y_new[0:int(well_location_list[1])] = 0
	elif include_wells_involved_in_thrust_belt == 0:
		y_new[0:int(load_width)] = 0
		y_new[int(well_location_list[-4]) + 1:] = 0

	# Clean up well location list
	well_location_list_shifted = np.delete(well_location_list, points_to_be_removed)
	well_location_list_shifted1_2 = well_location_list_shifted[0:-4]
	return y_new, well_location_list_shifted1_2


def build_mountain_belt_geometry(load_height, load_height2, load_width, load_width2, load_geometry, load_geometry2, modeling_section_profile,
								 well_location_list_shifted, displacement2):
	"""
	Build mountain belt geometry based on specified shape parameters.

	Parameters:
	    load_height (float): Height of the first mountain load (km).
	    load_height2 (float): Height of the second mountain load (km).
	    load_width (float): Width of the first mountain load (km).
	    load_width2 (float): Width of the second mountain load (km).
	    load_geometry (str): Geometry type of first mountain load.
	    load_geometry2 (str): Geometry type of second mountain load.
	    modeling_section_profile (np.array): X-coordinates of the model section.
	    well_location_list_shifted (np.array): Updated well positions.
	    displacement2 (float): Displacement of secondary mountain belt.

	Returns:
	    tuple: Load distribution arrays and surface arrays.
	"""
	mountain_load_distribution = []
	mountain_load_distribution2 = []
	xSurface_triangle2 = []
	xSurface_triangle2_print = []

	# Build mountain load distribution
	if load_geometry == 'Triangle':
		surface_slope = load_height / load_width  # in tan
		xSurface_triangle = load_width - modeling_section_profile
		xSurface_triangle[xSurface_triangle < 0] = 0  # For places outboard mountain belt, load height is equal to zero.
		mountain_load_distribution = surface_slope * xSurface_triangle

	elif load_geometry2 == 'Triangle2':
		surface_slope2 = load_height2 / load_width2
		load_width_km2 = int(load_width2 / 1000)
		modeling_section_profile_km = modeling_section_profile / 1000
		xSurface_triangle2 = int(well_location_list_shifted[-1] + displacement2 + load_width_km2) - modeling_section_profile_km
		xSurface_triangle2[xSurface_triangle2 < 0] = 0  # For places outboard mountain belt, load height is equal to zero.
		xSurface_triangle2[xSurface_triangle2 > load_width_km2] = 0
		mountain_load_distribution2 = surface_slope2 * xSurface_triangle2 *1000
		xSurface_triangle2_print_ = xSurface_triangle2[int(well_location_list_shifted[-1]) + displacement2:int(well_location_list_shifted[-1]) + displacement2 + load_width_km2 + 1]
		xSurface_triangle2_print =np.flip(xSurface_triangle2_print_)
		mountain_load_distribution2[int(well_location_list_shifted[-1]) + displacement2:
									int(well_location_list_shifted[-1]) + displacement2 + load_width_km2 + 1] = surface_slope2 * xSurface_triangle2_print * 1000

	elif load_geometry == 'double_sided_triangle':
		surface_slope = load_height / load_width  # in tan
		xSurface_triangle = load_width - modeling_section_profile
		xSurface_triangle[xSurface_triangle < 0] = 0  # For places outboard mountain belt, load height is equal to zero.
		mountain_load_distribution = surface_slope * xSurface_triangle
		half_load_width = int(load_width / 2000)  # convert m to km because we care about index.
		mountain_load_distribution[0:half_load_width] = np.flip(mountain_load_distribution[half_load_width + 1:
                                                         half_load_width * 2 + 1])

	elif load_geometry == 'Rectangle':
		xSurface_rectangle = load_width - modeling_section_profile
		xSurface_rectangle[xSurface_rectangle >= 0] = 1
		xSurface_rectangle[xSurface_rectangle < 0] = 0
		mountain_load_distribution = load_height * xSurface_rectangle

	elif load_geometry2 == 'Rectangle2':
		load_width_km2 = int(load_width2 / 1000)
		modeling_section_profile_km = modeling_section_profile / 1000
		mountain_load_distribution2 = np.zeros(len(modeling_section_profile_km))
		mountain_load_distribution2[int(well_location_list_shifted[-1] + displacement2):
									int(well_location_list_shifted[-1] + displacement2 + load_width_km2) + 1] = load_height2 *1

	elif load_geometry == 'Rectangle_plus_triangle_at_left':  # rectangle plus triangle at its top
		surface_slope = load_height / load_width
		xSurface_triangle = load_width - modeling_section_profile
		xSurface_triangle[xSurface_triangle < 0] = 0  # For places outboard mountain belt, load height is equal to zero.
		mountain_load_distribution = surface_slope * xSurface_triangle * (1 / 3)
		load_width_km = int(load_width / 1000)
		mountain_load_distribution[0: load_width_km + 1] = mountain_load_distribution[0: load_width_km + 1] + int(2 / 3 * load_height)

	elif load_geometry == 'Rectangle_plus_triangle_at_right':  # rectangle plus triangle at its top
		surface_slope = (1 / 3) * load_height / load_width
		xSurface_triangle = load_width - modeling_section_profile
		xSurface_triangle[xSurface_triangle < 0] = 0  # For places outboard mountain belt, load height is equal to zero.
		mountain_load_distribution = surface_slope * xSurface_triangle
		load_width_km = int(load_width / 1000)
		mountain_load_distribution[0:load_width_km + 1] = np.flip(mountain_load_distribution[0:load_width_km + 1])
		mountain_load_distribution[0:load_width_km + 1] = mountain_load_distribution[0: load_width_km + 1] + int(
			2 / 3 * load_height)

	elif load_geometry2 == 'Rectangle_plus_triangle_at_right2':  # rectangle plus triangle at its top
		surface_slope2 = load_height2 / load_width2
		load_width_km2 = int(load_width2 / 1000)
		modeling_section_profile_km = modeling_section_profile / 1000
		xSurface_triangle2 = int(well_location_list_shifted[-1] + displacement2 + load_width_km2) - modeling_section_profile_km
		xSurface_triangle2[xSurface_triangle2 < 0] = 0  # For places outboard mountain belt, load height is equal to zero.
		xSurface_triangle2[xSurface_triangle2 > load_width_km2] = 0
		mountain_load_distribution2 = surface_slope2 * xSurface_triangle2 * 1000 * (1/3)
		xSurface_triangle2_print_ = xSurface_triangle2[int(well_location_list_shifted[-1]) + displacement2:int(
			well_location_list_shifted[-1]) + displacement2 + load_width_km2 + 1]
		xSurface_triangle2_print = np.flip(xSurface_triangle2_print_)
		mountain_load_distribution2[int(well_location_list_shifted[-1]) + displacement2:
									int(well_location_list_shifted[-1]) + displacement2 + load_width_km2 + 1] \
			= surface_slope2 * xSurface_triangle2_print * 1000 * (1/3) + int(load_height2 * (2/3))

	elif load_geometry == 'Isosceles_Trapezoid':
		surface_slope = load_height / (load_width * (1/3))
		xSurface_triangle = load_width - modeling_section_profile
		xSurface_triangle[xSurface_triangle < 0] = 0
		load_width_km = int(load_width / 1000)
		modeling_section_profile_km = modeling_section_profile / 1000
		mountain_load_distribution = surface_slope * xSurface_triangle
		one_third_load_width = int(load_width_km * (1/3))
		mountain_load_distribution[0: one_third_load_width] = np.flip(mountain_load_distribution[one_third_load_width * 2 + 1:
																								 one_third_load_width * 3 + 1])
		mountain_load_distribution[one_third_load_width: one_third_load_width * 2 + 1] = load_height

	elif load_geometry2 == 'Isosceles_Trapezoid2':
		surface_slope2 = load_height2 / (load_width2 * (1/3))
		load_width_km2 = int(load_width2 / 1000)
		modeling_section_profile_km = modeling_section_profile / 1000
		xSurface_triangle2 = int(well_location_list_shifted[-1] + displacement2 + load_width_km2) - modeling_section_profile_km
		xSurface_triangle2[xSurface_triangle2 < 0] = 0
		xSurface_triangle2[xSurface_triangle2 > load_width_km2] = 0
		mountain_load_distribution2 = surface_slope2 * xSurface_triangle2 * 1000
		one_third_load_width2 = int(load_width_km2 * (1 / 3))
		xSurface_triangle2_print_ = xSurface_triangle2[int(well_location_list_shifted[-1]) + displacement2 + one_third_load_width2 *2:
									int(well_location_list_shifted[-1]) + displacement2 + one_third_load_width2 *3 + 1]
		xSurface_triangle2_print =np.flip(xSurface_triangle2_print_)
		mountain_load_distribution2[int(well_location_list_shifted[-1]) + displacement2:
									int(well_location_list_shifted[-1]) + displacement2 + one_third_load_width2 + 1] \
			= surface_slope2 * xSurface_triangle2_print * 1000
		mountain_load_distribution2[int(well_location_list_shifted[-1]) + displacement2 + one_third_load_width2:
									int(well_location_list_shifted[-1]) + displacement2 + one_third_load_width2*2 + 1] = load_height2
		mountain_load_distribution2[int(well_location_list_shifted[-1]) + displacement2 + one_third_load_width2*2:
									int(well_location_list_shifted[-1]) + displacement2 + one_third_load_width2*3 + 1] \
			= surface_slope2 * xSurface_triangle2_print_ * 1000
	return mountain_load_distribution, mountain_load_distribution2, xSurface_triangle2_print


def cal_flexure(load_height, load_height2, load_width, load_width2, displacement2, well_location_list_shifted, alpha, alpha2, rho_c, rho_s_list, delta_rho, sedTh, modeling_section_profile,
                precomputed_const_paramter_list, precomputed_const_paramter_list2, flexure_modeling_mode='broken', load_geometry='triangle', load_geometry2='triangle',
                add_water_load_mode='off', rho_w=1000, water_depth_list=None):
	"""
	   Calculate flexural deflection induced by tectonic and sedimentary loads.

	   Parameters:
	       load_height (float): Height of primary mountain load (km).
	       load_height2 (float): Height of secondary mountain load (km).
	       load_width (float): Width of primary mountain load (km).
	       load_width2 (float): Width of secondary mountain load (km).
	       alpha (float): Flexural rigidity parameter (m).
	       rho_c (float): Crustal density (kg/m³).
	       rho_s_list (np.array): Sediment densities (kg/m³).
	       delta_rho (float): Density contrast between mantle and crust.
	       sedTh (np.array): Sediment thickness along profile.
	       modeling_section_profile (np.array): X-coordinates of model section.
	       precomputed_const_paramter_list (list): Precomputed constants for calculation.
	       flexure_modeling_mode (str): Method to compute deflection ('broken' or others).
	       add_water_load_mode (str): Whether to include water load ('on' or 'off').
	       rho_w (float): Water density (kg/m³).
	       water_depth_list (np.array): Water depth along profile.

	   Returns:
	       tuple: Various deflection and topography arrays.
	   """
	# Build mountain load distribution
	load_height = load_height * 1000  # convert unit km to meter
	load_width = load_width * 1000  # convert unit km to meter
	load_height2 = load_height2 * 1000  # convert unit km to meter
	load_width2 = load_width2 * 1000  # convert unit km to meter
	modeling_section_profile = modeling_section_profile * 1000

	g = 9.81  # gravity acceleration

	# Define flexure modeling profile
	section_length = len(modeling_section_profile)
	dx = modeling_section_profile[2] - modeling_section_profile[1]  # define point interval
	x_dim = modeling_section_profile[:]

	# Build_mountain_belt_geometry by calling build_mountain_belt_geometry function
	mountain_load_distribution, _, xSurface_triangle = build_mountain_belt_geometry(load_height=load_height, load_height2=None, load_width=load_width, load_width2=None, load_geometry=load_geometry, load_geometry2=None,
																modeling_section_profile=modeling_section_profile, well_location_list_shifted=None, displacement2=None)
	_, mountain_load_distribution2, xSurface_triangle2_print = build_mountain_belt_geometry(load_height=None, load_height2=load_height2, load_width = None, load_width2=load_width2, load_geometry=None, load_geometry2=load_geometry2,
																  modeling_section_profile=modeling_section_profile, well_location_list_shifted=well_location_list_shifted, displacement2=displacement2)

	mountain_load_distribution2_test = mountain_load_distribution2
	load_distribution_along_profile = mountain_load_distribution + mountain_load_distribution2 + sedTh

	# Define a numpy array to hold cumulative deflection by mountain load
	summed_deflection_by_mountain_load = np.zeros(section_length)
	summed_deflection_by_mountain_load2 = np.zeros(section_length)

	# Define a numpy array to hold cumulative deflection by sediment load
	summed_deflection_by_sediment = np.zeros(section_length)

	# Calculate deflection induced by mountain load and sediment separately
	# Slice load into slices at specified interval, both mountain belts and sediments
	for i in range(0, len(mountain_load_distribution), 1):
		A = precomputed_const_paramter_list[0]
		B = precomputed_const_paramter_list[1]
		precomputed_const_parameter_C = precomputed_const_paramter_list[2]
		if i == len(mountain_load_distribution) - 1:
			mountain_load_elevation = mountain_load_distribution[i]
			sediment_load_elevation = sedTh[i]
		else:
			mountain_load_elevation = (mountain_load_distribution[i] + mountain_load_distribution[i + 1]) * 0.5
			sediment_load_elevation = (sedTh[i] + sedTh[i + 1]) * 0.5
		load_location = modeling_section_profile[i]
		precomputed_item_A = abs(load_location - x_dim) / alpha

		# For broken model. We determine plate deflection using the solution given by Hetényi (1946)
		if flexure_modeling_mode == 'broken':
			precomputed_item_B = load_location / alpha
			C = np.exp(-precomputed_item_B) * (np.cos(precomputed_item_B) - np.sin(precomputed_item_B))
			D = np.exp(-precomputed_item_B) * np.cos(precomputed_item_B)
			E = np.exp(-precomputed_item_A) * (np.cos(precomputed_item_A) + np.sin(precomputed_item_A))
			w = (C + 2 * D) * A - 2 * (C + D) * B + E

		# Calculate deflection induced by mountain load (thrust load)
		mountain_slice_load = rho_c * g * mountain_load_elevation
		vo_mountain_load = mountain_slice_load * dx  # volume of mountain slice load

		# Calculate deflection along the profile
		w_max_by_mountain_load = vo_mountain_load / precomputed_const_parameter_C * (-1)
		w_dim_mountain_load = w_max_by_mountain_load * w
		summed_deflection_by_mountain_load = w_dim_mountain_load + summed_deflection_by_mountain_load
		sediment_slice_load = rho_s_list[i] * g * sediment_load_elevation

		if add_water_load_mode == 'on':
			sediment_slice_load = sediment_slice_load + water_depth_list[i] * g * rho_w

		vo_sediment_load = sediment_slice_load * dx  # volume of sediment slice load

		# Calculate deflection along the profile
		w_max_by_sediment_load = vo_sediment_load / precomputed_const_parameter_C * (-1)
		w_dim_sediment_load = w_max_by_sediment_load * w
		summed_deflection_by_sediment = summed_deflection_by_sediment + w_dim_sediment_load  #summed_deflection_by_sediment 单位为负数

	for j in range(0, len(mountain_load_distribution2), 1):
		A2 = precomputed_const_paramter_list2[0]
		B2 = precomputed_const_paramter_list2[1]
		precomputed_const_parameter_C2 = precomputed_const_paramter_list2[2]
		if j == len(mountain_load_distribution2) - 1 :
			mountain_load_elevation2 = mountain_load_distribution2[j]
		else:
			mountain_load_elevation2 = (mountain_load_distribution2[j] + mountain_load_distribution2[j + 1]) * 0.5

		load_location2 = modeling_section_profile[j]
		precomputed_item_A2 = abs(load_location2 - x_dim) / alpha2

		if flexure_modeling_mode == 'broken':
			# Numpy supports list operation element by element (element-wise operation).
			# These lines below are modified to improve computing efficiency
			precomputed_item_B2 = load_location2 / alpha2
			C2 = np.exp(-precomputed_item_B2) * (np.cos(precomputed_item_B2) - np.sin(precomputed_item_B2))
			D2 = np.exp(-precomputed_item_B2) * np.cos(precomputed_item_B2)
			E2 = np.exp(-precomputed_item_A2) * (np.cos(precomputed_item_A2) + np.sin(precomputed_item_A2))
			w2 = (C2 + 2 * D2) * A2 - 2 * (C2 + D2) * B2 + E2

		# Calculate deflection induced by mountain load (thrust load)
		mountain_slice_load2 = rho_c * g * mountain_load_elevation2
		vo_mountain_load2 = mountain_slice_load2 * dx  # volume of mountain slice load

		# Calculate deflection along the profile
		w_max_by_mountain_load2 = vo_mountain_load2 / precomputed_const_parameter_C2 * (-1)
		w_dim_mountain_load2 = w_max_by_mountain_load2 * w2
		summed_deflection_by_mountain_load2 = w_dim_mountain_load2 + summed_deflection_by_mountain_load2

	summed_deflection_by_mountain_load_only1_2 = summed_deflection_by_mountain_load + summed_deflection_by_mountain_load2
	summed_deflection = summed_deflection_by_sediment + summed_deflection_by_mountain_load_only1_2
	combined_topography = summed_deflection + load_distribution_along_profile
	final_mountain_topography = summed_deflection + mountain_load_distribution + mountain_load_distribution2
	forebulge_location = 0
	basin_width = 0
	RS_residual_subsidence = -(summed_deflection + sedTh)
	return final_mountain_topography, load_distribution_along_profile, combined_topography, summed_deflection_by_mountain_load_only1_2, summed_deflection_by_mountain_load, \
		summed_deflection_by_mountain_load2, summed_deflection, mountain_load_distribution, mountain_load_distribution2, forebulge_location, basin_width, mountain_load_distribution2_test, xSurface_triangle2_print, RS_residual_subsidence


def plot_flexure(basin_name, section_profile, summed_w, summed_deflection_by_mountain_load_only1_2, modified_sedth,
                 original_sedth, mountain_load_topography, load_wt_value, load_wt_value2, section_number,
                 displacement_value, displacement_value2, modeling_age, well_location_list_shifted,
                 flexure_ytick_range_min=-4000, flexure_ytick_interval=1000,
                 figure_format_list=None, ytick_mode='manually',
                 plot_font_size=15, well_marker_size=4,
                 well_location_shift_value=2, xtick_interval= 50, load_display_fraction=1, load_display_fraction2 =1,
                 number_of_wells_to_show=0, starting_well_number=0, font_type='English',
                 RS_shift_value=0, section_profile_base_unit=100, well_number_excluded_list=None, xtick_maximum_value=None):
	"""
	    Plot flexural subsidence and topography profile with well symbols.

	    Parameters:
	        basin_name (str): Name of the basin being modeled.
	        section_profile (np.array): X-coordinates of model section.
	        summed_w (np.array): Total deflection array.
	        modified_sedth (np.array): Adjusted sediment thickness.
	        well_location_list_shifted (np.array): Updated well positions.
	        figure_format_list (list): Output formats for saved figures (e.g., ['png', 'pdf']).
	        plot_font_size (int): Font size for plot labels.
	        well_marker_size (int): Size of well markers.
	        RS_shift_value (float): Residual subsidence adjustment.
	        ...

	    Returns:
	        None. Saves figure to disk.
	    """
	# Set default image format if not provided
	if figure_format_list is None:
		figure_format_list = ['png']

	# Initialize figure with specified size (inches)
	fig = plt.figure(1)
	width_cm = 30
	height_cm = 15
	fig = plt.figure(figsize=(width_cm/2.54, height_cm/2.54))
	axes1 = fig.add_subplot(2, 1, 1)

	# Convert well location list to numpy array for plotting
	well_location_x_coordinates = np.array(well_location_list_shifted)  # x coordinates of wells
	well_location_y_height = 1000
	well_location_y_coordinates = np.zeros(len(well_location_list_shifted)) + well_location_y_height

	# Calculate length of the well-defined profile
	well_section_profile_length = well_location_x_coordinates[-1] - well_location_x_coordinates[0]

	# Determine appropriate lower limit for Y-axis based on key deflection points
	displayed_load_length = int(load_wt_value - load_display_fraction * load_wt_value)
	displayed_load_length2 = int(
		well_location_list_shifted[-1] + displacement_value2 + load_wt_value2 - load_display_fraction2 * load_wt_value2)
	deflection_at_load_front = summed_w[load_wt_value]
	deflection_at_load_front2 = summed_w[int(well_location_list_shifted[-1]) + displacement_value2]
	deflection_at_figure_westernmost = summed_w[displayed_load_length]
	deflection_at_figure_esternmost = summed_w[displayed_load_length2]
	y_axis_bottom_limit = min(deflection_at_load_front, deflection_at_load_front2, deflection_at_figure_westernmost,
							  deflection_at_figure_esternmost)
	y_axis_bottom_limit = math.ceil(y_axis_bottom_limit / 1000) * 1000

	# Plot reference line at zero subsidence
	reference_line_value = 0
	x_reference_line_value = section_profile[
							 0:displayed_load_length2]
	plt.plot(x_reference_line_value, np.zeros(len(x_reference_line_value)) + reference_line_value, color='#808183', linestyle='--', alpha=0.5)
	plt.plot(section_profile[0:displayed_load_length2], summed_w[0:displayed_load_length2], 'black')
	plt.fill_between(section_profile[load_wt_value:int(well_location_x_coordinates[- 1]+ displacement_value2) + 1],
					 summed_deflection_by_mountain_load_only1_2[load_wt_value:int(well_location_x_coordinates[-1] + displacement_value2) + 1],
	                 summed_w[load_wt_value:int(well_location_x_coordinates[- 1]+ displacement_value2) + 1], facecolor='#EEC49E', alpha=1)

	plt.plot(section_profile[int(well_location_x_coordinates[starting_well_number]):int(well_location_x_coordinates[-1]) + 1],
	         -modified_sedth[int(well_location_x_coordinates[starting_well_number]):int(well_location_x_coordinates[-1])+1], color = '#33348E', linestyle='-')
	plt.plot(section_profile[0:displayed_load_length2],
			 summed_deflection_by_mountain_load_only1_2[0:displayed_load_length2], '#6B4E2F')

	# Plot well symbols along the section profile.
	if number_of_wells_to_show == 0:  # By default, we show all wells
		plt.plot(well_location_x_coordinates, well_location_y_coordinates, 'k|', markersize=well_marker_size)
	else:
		plt.plot(well_location_x_coordinates[starting_well_number:number_of_wells_to_show],
		         well_location_y_coordinates[starting_well_number:number_of_wells_to_show], 'k|',
		         markersize=well_marker_size)

	# Plot numbers above well symbols
	for well_index in range(len(well_location_x_coordinates)):
		if well_number_excluded_list is None or well_index + 1 not in well_number_excluded_list:  # then we plot well number
			plt.text(well_location_x_coordinates[well_index], well_location_y_coordinates[well_index] + 300,
			         well_index + 1, horizontalalignment='center', fontsize=plot_font_size - 3)

	# Fill area from thrust front to the last well or certain well
	plt.fill_between(section_profile[load_wt_value:int(well_location_x_coordinates[- 1]+ displacement_value2) + 1],
	                 summed_w[load_wt_value:int(well_location_x_coordinates[- 1]+ displacement_value2) + 1],
	                 -modified_sedth[load_wt_value:int(well_location_x_coordinates[- 1]+ displacement_value2) + 1],
	                 facecolor='grey', alpha=0.4)

	# Fill area between mountain load topography and total deflection
	plt.fill_between(section_profile[0:displayed_load_length2], summed_w[0:displayed_load_length2], mountain_load_topography[0:displayed_load_length2], facecolor='#7F5E40', alpha=0.8)

	# Set appropriate axis range for figure
	x_axis_right_limit = math.ceil((int(well_location_x_coordinates[-1] + displacement_value2 + load_wt_value2)+ 1)/
								   section_profile_base_unit) * section_profile_base_unit
	y_axis_top_limit = abs(max(mountain_load_topography[displayed_load_length:-1])) + 1000
	axes1.axis([displayed_load_length, displayed_load_length2, y_axis_bottom_limit, y_axis_top_limit])

	# Set X-axis ticks
	if xtick_maximum_value is not None:
		xtick_range = xtick_maximum_value
	else:
		xtick_range = math.ceil((x_axis_right_limit - well_location_x_coordinates[0]) / xtick_interval) * xtick_interval
	xticks_list = [x for x in range(0, xtick_range, xtick_interval)]
	plt.xticks(xticks_list + well_location_x_coordinates[0], xticks_list, fontsize=plot_font_size)

	# Set Y-axis ticks
	if ytick_mode == 'manually':
		axes1.set_yticks([y for y in range(flexure_ytick_range_min, int(y_axis_top_limit), flexure_ytick_interval)])
	else:
		axes1.set_yticks([y for y in range(y_axis_bottom_limit, int(y_axis_top_limit), 1000)])
	plt.yticks(fontsize=plot_font_size)

	# Compute residual subsidence
	residual_subsidence = modified_sedth + summed_w + RS_shift_value
	RS_computation_start_range = int(well_location_x_coordinates[starting_well_number])
	residual_subsidence_tailored = residual_subsidence[RS_computation_start_range:
	                                int(well_location_x_coordinates[- 1]) + 1]
	average_residual_subsidence = np.mean(residual_subsidence_tailored)

	# Set title and Y-axis label
	if font_type == 'English':
		axes1.set_title('ARS = %0.1f m Disp=%skm Disp2=%skm' %
		                (average_residual_subsidence, displacement_value, displacement_value2), fontsize=plot_font_size - 1)
		axes1.set_ylabel('Subsidence (m)', fontsize=plot_font_size)

	# Save figure into local disk. If the directory doesn't exist, then we create it
	figure_path = './data/flexural_modeling/%s/section_%s/figure/%s/' % (
		basin_name, section_number, modeling_age)
	if not os.path.exists(figure_path):
		os.makedirs(figure_path)

	# Save figure in specified formats
	for figure_format in figure_format_list:
		flexure_figure_path = figure_path
		if not os.path.exists(flexure_figure_path):
			os.makedirs(flexure_figure_path)

		flexure_figure_file = flexure_figure_path + '%sMa.%s' % \
		                      (modeling_age, figure_format)
		plt.savefig(flexure_figure_file, dpi=500)

	plt.close()
