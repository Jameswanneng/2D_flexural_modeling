#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# ======================================================
#    (c) China University of Geosciences (Beijing)
#            ALL RIGHTS RESERVED
#      Authors: Neng Wan, Xueyan Li
#		 Oct, 05, 2024


import math
import sys
import numpy as np
import configparser
import Function
from Function import get_well_locations
from Function import cal_flexure as cal_flexure
from Function import sediment_thickness_interpolation_along_profile as sediment_interp
from Function import read_age_dependant_data


# Configure parameters
if len(sys.argv) == 2:
	path_config_file = sys.argv[1]
	config_obj = configparser.ConfigParser()
	config_obj.read(path_config_file)

	# Section parameters.
	# Define and input general section parameters
	section_info_dict = config_obj['section_parameters']
	basin_name = section_info_dict['basin_name']  # datatype: str
	section_number = int(section_info_dict['section_number'])  # datatype: str
	plot_flexural_figure_mode = section_info_dict['plot_flexural_figure_mode']

	# Age sequence parameters
	modeling_age_sequence = eval(section_info_dict['modeling_age_sequence'])
	model_baseline_age = int(section_info_dict['model_baseline_age'])  # datatype: str
	basin_bottom_age_sequence = eval(section_info_dict['basin_bottom_age_sequence'])
	basin_top_age_sequence = eval(section_info_dict['basin_top_age_sequence'])
	modeling_section_length = int(section_info_dict['modeling_section_length'])
	sediment_taper_distance = int(section_info_dict['sediment_taper_distance'])

	# Geometry parameters: the first mountain load
	model_parameter_dict = config_obj['model_parameters']
	load_shape_list = eval(model_parameter_dict['load_shape'])

	# Elastic thickness (Te) parameters
	Te_value_lower_limit = eval((model_parameter_dict['Te_value_lower_limit']))
	Te_value_upper_limit = eval((model_parameter_dict['Te_value_upper_limit']))
	Te_value_increment = eval((model_parameter_dict['Te_value_increment']))
	Te_parameter_dict = dict(zip(basin_bottom_age_sequence, zip(Te_value_lower_limit, Te_value_upper_limit,
	                                                        Te_value_increment)))
	# Load width and height parameters
	load_width_lower_limit = eval(model_parameter_dict['load_width_value_lower_limit'])
	load_width_upper_limit = eval(model_parameter_dict['load_width_value_upper_limit'])
	load_width_increment = eval(model_parameter_dict['load_width_value_increment'])
	load_width_parameter_dict = dict(zip(basin_bottom_age_sequence, zip(load_width_lower_limit, load_width_upper_limit,
	                                                                    load_width_increment)))
	load_height_lower_limit = eval(model_parameter_dict['load_height_value_lower_limit'])
	load_height_upper_limit = eval(model_parameter_dict['load_height_value_upper_limit'])
	load_height_increment = eval(model_parameter_dict['load_height_value_increment'])
	load_height_parameter_dict = dict(zip(basin_bottom_age_sequence, zip(load_height_lower_limit, load_height_upper_limit,
	                                                                     load_height_increment)))

	# Displacement parameters
	displacement_lower_limit = eval(model_parameter_dict['displacement_value_lower_limit'])
	displacement_upper_limit = eval(model_parameter_dict['displacement_value_upper_limit'])
	displacement_increment = eval(model_parameter_dict['displacement_value_increment'])
	displacement_parameter_dict = dict(zip(basin_bottom_age_sequence, zip(displacement_lower_limit, displacement_upper_limit,
	                                                                      displacement_increment)))

	# Geometry parameters:: the second mountain load
	load_shape_list2 = eval(model_parameter_dict['load_shape2'])

	# Elastic thickness (Te) parameters
	Te_value_lower_limit2 = eval((model_parameter_dict['Te_value_lower_limit2']))
	Te_value_upper_limit2 = eval((model_parameter_dict['Te_value_upper_limit2']))
	Te_value_increment2 = eval((model_parameter_dict['Te_value_increment2']))
	Te_parameter_dict2 = dict(zip(basin_bottom_age_sequence, zip(Te_value_lower_limit2, Te_value_upper_limit2,
																Te_value_increment2)))

	# Load width and height parameters
	load_width_lower_limit2 = eval(model_parameter_dict['load_width_value_lower_limit2'])
	load_width_upper_limit2 = eval(model_parameter_dict['load_width_value_upper_limit2'])
	load_width_increment2 = eval(model_parameter_dict['load_width_value_increment2'])
	load_width_parameter_dict2 = dict(zip(basin_bottom_age_sequence, zip(load_width_lower_limit2, load_width_upper_limit2,
																		load_width_increment2)))
	load_height_lower_limit2 = eval(model_parameter_dict['load_height_value_lower_limit2'])
	load_height_upper_limit2 = eval(model_parameter_dict['load_height_value_upper_limit2'])
	load_height_increment2 = eval(model_parameter_dict['load_height_value_increment2'])
	load_height_parameter_dict2 = dict(
		zip(basin_bottom_age_sequence, zip(load_height_lower_limit2, load_height_upper_limit2,
										   load_height_increment2)))

	# Displacement parameters
	displacement_lower_limit2 = eval(model_parameter_dict['displacement_value_lower_limit2'])
	displacement_upper_limit2 = eval(model_parameter_dict['displacement_value_upper_limit2'])
	displacement_increment2 = eval(model_parameter_dict['displacement_value_increment2'])
	displacement_parameter_dict2 = dict(zip(basin_bottom_age_sequence, zip(displacement_lower_limit2, displacement_upper_limit2, displacement_increment2)))

	# Read sediment load parameters
	sediment_load_parameter_dict = config_obj['sediment_load_parameters']
	sediment_load_density_mode = sediment_load_parameter_dict['sediment_load_density_mode']
	sediment_thickness_interpolation_switch = sediment_load_parameter_dict['sediment_thickness_interpolation_switch']

	if sediment_load_density_mode == 'Given_variable_densities':
		sediment_load_average_density_list = eval(sediment_load_parameter_dict['average_density_of_sediment_load'])
		sediment_load_average_density_dict = dict(zip(basin_bottom_age_sequence, sediment_load_average_density_list))

	# Erosion recovery
	erosion_recovery_dict = {}
	for option in config_obj['erosion_recovery']:
		if option != 'extra_thickness_scaling_factor' and option != 'extra_thickness_scaling_factor2':
			erosion_recovery_dict[float(option)] = eval(config_obj['erosion_recovery'][option])
	extra_thickness_scaling_factor_list = eval(config_obj['erosion_recovery']['extra_thickness_scaling_factor'])
	extra_thickness_scaling_factor_dict = dict(zip(basin_bottom_age_sequence, extra_thickness_scaling_factor_list))
	extra_thickness_scaling_factor_list2 = eval(config_obj['erosion_recovery']['extra_thickness_scaling_factor2'])
	extra_thickness_scaling_factor_dict2 = dict(zip(basin_bottom_age_sequence, extra_thickness_scaling_factor_list2))

	# Figure parameters
	figure_parameter_dict = config_obj['figure_parameters']
	output_figure_format = eval(figure_parameter_dict['output_figure_format'])
	flexure_ytick_range_min_list = eval(figure_parameter_dict['flexure_ytick_range_min'])
	flexure_ytick_interval_list = eval(figure_parameter_dict['flexure_ytick_interval'])
	RS_shift_value = int(figure_parameter_dict['RS_shift_value'])
	flexure_xtick_interval = int(figure_parameter_dict['flexure_xtick_interval'])
	section_profile_base_unit = int(figure_parameter_dict['section_profile_base_unit'])
	well_location_shift_value = float(figure_parameter_dict['well_location_shift_value'])
	load_display_fraction_list = eval(figure_parameter_dict['load_display_fraction'])
	load_display_fraction_list2 = eval(figure_parameter_dict['load_display_fraction2'])
	number_of_wells_to_show_list = eval(figure_parameter_dict['number_of_wells_to_show'])
	starting_well_number_list = eval(figure_parameter_dict['starting_well_number'])
	well_number_excluded_list = eval(figure_parameter_dict['well_number_excluded_list'])
	y_axis_top_limit_value_list = eval(figure_parameter_dict['y_axis_top_limit_value'])
	xtick_maximum_value_list = eval(figure_parameter_dict['xtick_maximum_value'])
	font_type = figure_parameter_dict['font_type']  # By default, we use English. NSimSun:新宋体，SimHei:黑体, SimSun:宋体

else:
	sys.exit()

# Data input
# Specify path to files and read sediment thickness data at each time period
decompacted_total_subsidence_file = './data/%s/section_%s/%s.xls' \
                                    % (basin_name, section_number, basin_name)
basin_layer_thickness_dict = read_age_dependant_data(basin_bottom_age_sequence, model_baseline_age,
                                                     decompacted_total_subsidence_file, column_index=1)

# Read well location data
well_list_file = './data/%s/%s.txt' % (basin_name, basin_name)
well_distance_dict = get_well_locations(well_list_file)
well_location_list = np.array(list(well_distance_dict.values()))

# Define model constants
sea_level = 0  # paleo sea level above modern (meters)
rho_s_constant = 2400  # density of sediment infill (kg/m^3)
rho_c = 2850  # density of crustal load (kg/m^3)
rho_m = 3300  # density of displaced mantle (kg/m^3)
rho_w = 1000  # density of water (kg/m^3)
rho_infill = 0  # since we are going to model the deflection induced by thrust belt and sediments, we must set it to 0
delta_rho = rho_m - rho_infill  # density difference between infill material and mantle.
g = 9.8  # gravitational constant (m/s^2)
v = 0.25  # poisson's ratio
E = 70000000000  # Young's Modulus (70 GPa)

# Load shape dictionary (1-10 correspond to different geometries)
load_shape_dict = {1: 'Triangle', 2: 'Triangle2',  3: 'double_sided_triangle', 4: 'Rectangle', 5: 'Rectangle2',
				   6: 'Rectangle_plus_triangle_at_left', 7: 'Rectangle_plus_triangle_at_right', 8: 'Rectangle_plus_triangle_at_right2',
				   9: 'Isosceles_Trapezoid'	, 10: 'Isosceles_Trapezoid2'}

# Set output figure format
output_figure_format_dict = {1: 'png', 2: 'pdf', 3: 'svg'}
output_figure_format_list = [output_figure_format_dict[x] for x in output_figure_format]

# Run flexural modeling
model_index = 0
modeling_section_profile = np.arange(0, modeling_section_length, 1)

# Initialize model parameters and results
model_parameter_result_dict = {}
model_parameter_result_broken = {}
model_parameter_result_dict['broken'] = model_parameter_result_broken
model_parameter_result_dict2 = {}
model_parameter_result_broken2 = {}
model_parameter_result_dict2['broken'] = model_parameter_result_broken2
Te_value_dict = {}
displacement_value_dict = {}
load_height_value_dict = {}
load_width_value_dict = {}
Te_value_dict2 = {}
displacement_value_dict2 = {}
load_height_value_dict2 = {}
load_width_value_dict2 = {}

# Initialize parameter ranges for each modeling age
for modeling_age in modeling_age_sequence:
	# First load parameters
	Te_value_dict[modeling_age] = np.arange(Te_parameter_dict[modeling_age][0], Te_parameter_dict[modeling_age][1],
		                                        Te_parameter_dict[modeling_age][2])
	displacement_value_dict[modeling_age] = np.arange(displacement_parameter_dict[modeling_age][0],
		                                                  displacement_parameter_dict[modeling_age][1],
		                                                  displacement_parameter_dict[modeling_age][2])
	load_height_value_dict[modeling_age] = np.arange(load_height_parameter_dict[modeling_age][0],
		                                                 load_height_parameter_dict[modeling_age][1],
		                                                 load_height_parameter_dict[modeling_age][2])
	load_width_value_dict[modeling_age] = np.arange(load_width_parameter_dict[modeling_age][0],
		                                                load_width_parameter_dict[modeling_age][1],
		                                                load_width_parameter_dict[modeling_age][2])

	# Second load parameters
	Te_value_dict2[modeling_age] = np.arange(Te_parameter_dict2[modeling_age][0], Te_parameter_dict2[modeling_age][1],
		                                        Te_parameter_dict2[modeling_age][2])
	displacement_value_dict2[modeling_age] = np.arange(displacement_parameter_dict2[modeling_age][0],
		                                                  displacement_parameter_dict2[modeling_age][1],
		                                                  displacement_parameter_dict2[modeling_age][2])
	load_height_value_dict2[modeling_age] = np.arange(load_height_parameter_dict2[modeling_age][0],
		                                                 load_height_parameter_dict2[modeling_age][1],
		                                                 load_height_parameter_dict2[modeling_age][2])
	load_width_value_dict2[modeling_age] = np.arange(load_width_parameter_dict2[modeling_age][0],
		                                                load_width_parameter_dict2[modeling_age][1],
		                                                load_width_parameter_dict2[modeling_age][2])

# Main modeling loop
# Iterate through all modeling ages.
for modeling_age in modeling_age_sequence:
	age_index = basin_bottom_age_sequence.index(modeling_age)
	load_shape_value = load_shape_list[age_index]
	load_shape_value2 = load_shape_list2[age_index]
	load_shape = load_shape_dict[load_shape_value]
	load_shape2 = load_shape_dict[load_shape_value2]
	model_top_age = basin_top_age_sequence[age_index]
	flexure_ytick_range_min = flexure_ytick_range_min_list[age_index]
	flexure_ytick_interval = flexure_ytick_interval_list[age_index]
	load_display_fraction = load_display_fraction_list[age_index]
	load_display_fraction2 = load_display_fraction_list2[age_index]
	number_of_wells_to_show = number_of_wells_to_show_list[age_index]
	starting_well_number = starting_well_number_list[age_index]
	y_axis_top_limit_value = y_axis_top_limit_value_list[age_index]
	xtick_maximum_value = xtick_maximum_value_list[age_index]
	model_parameter_result_broken[modeling_age] = {}
	model_parameter_result_broken2[modeling_age] = {}

	# Loop through displacement values for each age
	for displacement_value in displacement_value_dict[modeling_age]:
		model_parameter_result_broken[modeling_age][displacement_value] = {}
		for displacement_value2 in displacement_value_dict2[modeling_age]:
			model_parameter_result_broken2[modeling_age][displacement_value2] = {}

			# Loop through Te values
			for ii in range(len(Te_value_dict[modeling_age])):
				Te_value = Te_value_dict[modeling_age][ii]
				Te_value2 = Te_value_dict2[modeling_age][ii]

				# Calculate D and alpha values
				D = (E * Te_value * Te_value * Te_value * math.pow(10, 9)) / (12 * (1 - v * v))
				D2 = (E * Te_value2 * Te_value2 * Te_value2 * math.pow(10, 9)) / (12 * (1 - v * v))
				alpha_value = np.power((4 * D / (delta_rho * g)), 0.25)
				alpha_value2 = np.power((4 * D2 / (delta_rho * g)), 0.25)

				# Precompute flexural response terms
				precomputed_item_A = modeling_section_profile * 1000 / alpha_value
				precomputed_item_A2 = modeling_section_profile * 1000 / alpha_value2
				precomputed_const_A = np.exp(-precomputed_item_A) * (
							np.cos(precomputed_item_A) + np.sin(precomputed_item_A))
				precomputed_const_A2 = np.exp(-precomputed_item_A2) * (
						np.cos(precomputed_item_A2) + np.sin(precomputed_item_A2))
				precomputed_const_B = np.exp(-precomputed_item_A) * np.sin(precomputed_item_A)
				precomputed_const_B2 = np.exp(-precomputed_item_A2) * np.sin(precomputed_item_A2)
				precomputed_const_parameter_C = (2 * delta_rho * g * alpha_value)
				precomputed_const_parameter_C2 = (2 * delta_rho * g * alpha_value2)
				precomputed_const_parameter_list = [precomputed_const_A, precomputed_const_B,
													precomputed_const_parameter_C]
				precomputed_const_parameter_list2 = [precomputed_const_A2, precomputed_const_B2,
													precomputed_const_parameter_C2]

				# Use dictionaries to hold model parameters and results
				model_parameter_result_broken[modeling_age][displacement_value][Te_value] = np.zeros(
					(9, len(load_height_value_dict[modeling_age]), len(load_width_value_dict[modeling_age])))
				model_parameter_result_broken2[modeling_age][displacement_value2][Te_value2] = np.zeros(
					(9, len(load_height_value_dict2[modeling_age]), len(load_width_value_dict2[modeling_age])))

				# Loop through load width and load height values
				for kk in range(len(load_width_value_dict[modeling_age])):
					load_wt_value = load_width_value_dict[modeling_age][kk]
					for tt in range(len(load_height_value_dict2[modeling_age])):
						load_wt_value2 = load_width_value_dict2[modeling_age][kk]
						for qq in range(len(load_height_value_dict[modeling_age])):
							load_height_value = load_height_value_dict[modeling_age][qq]
							for uu in range(len(load_height_value_dict2[modeling_age])):
								load_height_value2 = load_height_value_dict2[modeling_age][qq]

								# Process sediment thickness
								if sediment_thickness_interpolation_switch == 'on':
									original_sedth = basin_layer_thickness_dict[modeling_age]
									unmodified_sed_thickness, well_location_list_shifted1_2 = sediment_interp(well_location_list,
																									   original_sedth,
																									   displacement_value, displacement_value2,
																									   load_wt_value, load_wt_value2,
																									   section_length=modeling_section_length,
																									   interp_method="linear",
																									   extra_thickness=0,
																									   erosion_recovery_list=None,
																									   taper_distance=sediment_taper_distance)
									extra_thickness = displacement_value * extra_thickness_scaling_factor_dict[modeling_age]
									extra_thickness2 = displacement_value2 * extra_thickness_scaling_factor_dict2[modeling_age]
									modified_sed_thickness, _ = sediment_interp(well_location_list, original_sedth,
																			displacement_value, displacement_value2,
																			load_wt_value, load_wt_value2,
																			section_length=modeling_section_length,
																			interp_method="linear",
																			extra_thickness=extra_thickness, extra_thickness2=extra_thickness2,
																			erosion_recovery_list=erosion_recovery_dict[modeling_age],
																			taper_distance=sediment_taper_distance)
								# Set sediment density
								if sediment_load_density_mode == 'constant':
									rho_s = np.ones(modeling_section_length) * rho_s_constant
								elif sediment_load_density_mode == 'Given_variable_densities':
									rho_s = sediment_load_average_density_dict[modeling_age] * 1000 * np.ones(modeling_section_length)

								# Run flexural models. Use methods by Saylor et al., (2018) and Liu et al., (2004)
								for model_type in ['broken']:
									[mountain_load_topography, load_distribution, final_topography, summed_deflection_by_mountain_load_only1_2,
									deflection_by_mountain_only, deflection_by_mountain_only2, summed_w, mountain_load_configuration, mountain_load_configuration2,
									forebulge_location, basin_width, mountain_load_distribution2_test, xSurface_triangle2_print, RS_residual_subsidence] = cal_flexure(load_height_value, load_height_value2,
																				load_wt_value, load_wt_value2, displacement_value2, well_location_list_shifted1_2,
																				alpha_value, alpha_value2, rho_c, rho_s, delta_rho,
																				modified_sed_thickness,
																				modeling_section_profile,
																				precomputed_const_parameter_list, precomputed_const_parameter_list2,
																				flexure_modeling_mode=model_type,
																				load_geometry=load_shape, load_geometry2=load_shape2)
									# Generate figures
									if plot_flexural_figure_mode == 'on':
										Function.plot_flexure(basin_name, modeling_section_profile, summed_w,
																 summed_deflection_by_mountain_load_only1_2,
																 modified_sed_thickness, unmodified_sed_thickness,
																 mountain_load_topography,  load_wt_value, load_wt_value2,
																 section_number, displacement_value, displacement_value2,
																 modeling_age, well_location_list_shifted1_2,
																 flexure_ytick_range_min=flexure_ytick_range_min,
																 flexure_ytick_interval=flexure_ytick_interval,
																 figure_format_list=output_figure_format_list,
																 well_location_shift_value=well_location_shift_value,
																 xtick_interval=flexure_xtick_interval,
																 font_type=font_type,
																 load_display_fraction=load_display_fraction, load_display_fraction2=load_display_fraction2,
																 number_of_wells_to_show=number_of_wells_to_show,
																 starting_well_number=starting_well_number,
																 RS_shift_value=RS_shift_value,
																 section_profile_base_unit=section_profile_base_unit,
																 well_number_excluded_list=well_number_excluded_list,
																 xtick_maximum_value= xtick_maximum_value)

					model_index += 1
