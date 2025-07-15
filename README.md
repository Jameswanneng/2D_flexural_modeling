# 2D_flexural_modeling
A set of python scripts for running 2D flexural models.

Flexural Backstripping Modeling Toolkit

This repository contains a Python toolkit for 2D flexural backstripping modeling developed by researchers at China University of Geosciences (Beijing). 
The code was originally created by Neng Wan and later modified by Xueyan Li.

Key Features:
	1.Implements flexural isostatic modeling of sedimentary basins
	2.Handles multiple tectonic loads with various geometries (triangles, rectangles, trapezoids)
	3.Visualizes modeling results with publication-quality figures

Main Components:
	1.2D_flexural_modeling.py - Main driver script for flexural modeling
	2.Flexural_modeling_2D.py - Handles well location and stratigraphic data
	3.Core_basin_analysis.py - Core functions for flexure calculations and visualization
	4.__init__.py - Package initialization file

Usage:
	1.Configure parameters in the config file
	2.Run 2D_flexural_modeling.py with the config file as argument
	3.View output figures in the specified directory

Dependencies:
	1.Python 3.x
	2.NumPy
	3.SciPy
	4.Matplotlib
	5.xlrd

Citation:
When use this software, please cite:
	1. Wan, N., Liu, S. F., Li, X. Y., Zhang, B., Ren, R., & Chen, Z. X. (2023). Limited flexural control of fold-thrust belts on the Jurassic Sichuan Basin, South China. Frontiers in Earth Science, 11, 1276832. https://doi.org/10.3389/feart.2023.1276832 
  2. Liu, S. F., Zhang, B., Ma, P. F., Williams, S., Lin, C. F., Wan, N., Ran, C. L., & Gurnis, M. (2024). Craton deformation from flat-slab subduction and rollback. Nature Geoscience, 17, 1-8. https://doi.org/10.1038/s41561-024-01513-2
