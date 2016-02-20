# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 18:58:41 2012

@author: romain
"""
# LIBRARIES
import os as os
from copy import deepcopy

# CONSTANTS
Tr = True
Fl = False


# SCRIPT SETTINGS

load_data = Fl   #True if the data have to be loaded
load_library = Fl #True if the library has to be loaded
set_extent = Fl
resize = Fl
Offset = Fl
filter_sig = Fl
filtered = Tr
# Plot parameters
sensitivity = 0.05  	# when mode set to other
mode = "other"		# "other" or "max"
color_trace =["blue", "red"]
color_retrace = ["red", "blue"]

# SCRIPT BEGIN

if load_library :
	os.chdir("/home/romain")
	from setproc.data_plot.classes.map import Map

if load_data :
	os.chdir("/media/Iomega_HDD/Measures/Tb2C/sample3/7R/J9/Hysteresis_vs__Vg/")
	Map1t = Map("G_B_trace.json")
	Map1r = Map("G_B_retrace.json")
	os.chdir("/media/Iomega_HDD/Measures/Tb2C/sample3/7R/J9/Hysteresis_vs__Vg2/")
	Map2t = Map("G_B_trace.json")
	Map2r = Map("G_B_retrace.json")

if set_extent :
	# Set extent
	Map1t.SetExtent([-0.8,-0.4,-0.4,0.4])
	Map1r.SetExtent([-0.8,-0.4,-0.4,0.4])
	Map2t.SetExtent([-0.4,-0.2,-0.4,0.4])
	Map2r.SetExtent([-0.4,-0.2,-0.4,0.4])

if resize :
	# Resize
	Map1t.sweeps.Resize(-0.2,0.2)
	Map1r.sweeps.Resize(0.2,-0.2)
	Map2t.sweeps.Resize(-0.2,0.2)
	Map2r.sweeps.Resize(0.2,-0.2)
if Offset :
	Map1t.sweeps.SetOffset(0.0308)
	Map1r.sweeps.SetOffset(-0.0308)
	Map2t.sweeps.SetOffset(0.0308)
	Map2r.sweeps.SetOffset(-0.0308)

if filter_sig :
	# Copy
	Map1tf = deepcopy(Map1t)
	Map1rf = deepcopy(Map1r)
	Map2tf = deepcopy(Map2t)
	Map2rf = deepcopy(Map2r)
	# Filter 
	Map1tf.sweeps.GetFilter(0,4,1,1)
	Map1rf.sweeps.GetFilter(0,4,1,1)
	Map2tf.sweeps.GetFilter(0,4,1,1)
	Map2rf.sweeps.GetFilter(0,4,1,1)

if plot :
	fig= figure()
	ax1 = fig.add_subplot(311)
	ax2 = fig.add_subplot(312)
	ax3 = fig.add_subplot(313)
	if filtered :
		Map1tf.PlotFilteredExtrema(newplot=False, other_axes = ax2,
								   sensitivity = sensitivity, mode = mode,
								   color = color_trace)
		Map1rf.PlotFilteredExtrema(newplot=False, other_axes = ax3,
								   sensitivity = sensitivity, mode = mode,
								   color = color_retrace)
		Map2tf.PlotFilteredExtrema(newplot=False, other_axes = ax2,
								   sensitivity = sensitivity, mode = mode,
								   color = color_trace)
		Map2rf.PlotFilteredExtrema(newplot=False, other_axes = ax3,
								   sensitivity = sensitivity, mode = mode,
								   color = color_retrace)
	else :
		Map1tf.PlotExtrema(newplot=False, other_axes = ax2, color = color_trace)
		Map1rf.PlotExtrema(newplot=False, other_axes = ax3, color = color_retrace)
		Map2tf.PlotExtrema(newplot=False, other_axes = ax2, color = color_trace)
		Map2rf.PlotExtrema(newplot=False, other_axes = ax3, color = color_retrace)

	ax1.plot(P1t[0], P1t[1], color = "red")
	ax1.plot(P2t[0], P2t[1], color = "red")

	
	ax1.set_xlim(-0.7,-0.2)
	ax2.set_xlim(-0.7,-0.2)
	ax3.set_xlim(-0.7,-0.2)

	ax2.set_ylim(-0.06,0.06)
	ax3.set_ylim(-0.06,0.06)


