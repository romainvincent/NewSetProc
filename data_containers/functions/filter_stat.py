#external library
import numpy as np

def filter_stat_threshold(stat, mean, sensitivity) :
	"""
	Filter a statitic given a sensitivity

	Keyword arguments :
	stat        -- array containing the statistic to filter_stat
	mean        -- mean value of the sweep from which is extracted the statistic
	sensitivity -- sensitivity of the detection
	"""
	#init variable
	stat_len = len(stat) #should be equal to 4
	min_nbr = np.size(stat[0]) 
	max_nbr = np.size(stat[2])
	threshold = np.abs(mean * sensitivity)
	xmin = []
	vmin = []
	xmax = []
	vmax = []
	#start with min
	if min_nbr < 2  and stat[0] < - threshold :
		xmin.append(stat[1])
		vmin.append(stat[0])
	elif min_nbr > 1 :
		for i in range(min_nbr) :
			temp_vmin = stat[0][i]
			temp_xmin = stat[1][i]
			if temp_vmin < - threshold :
				xmin.append(temp_xmin)
				vmin.append(temp_vmin)


	#start with max
	if max_nbr < 2  and stat[2] > threshold :
		xmax.append(stat[3])
		vmax.append(stat[2])
	elif min_nbr > 1 :
		for i in range(max_nbr) :
			temp_vmax = stat[2][i]
			temp_xmax = stat[3][i]
			if temp_vmax > threshold :
				xmax.append(temp_xmax)
				vmax.append(temp_vmax)

	xmin = np.array(xmin)
	xmax = np.array(xmax)
	vmin = np.array(vmin)
	vmax = np.array(vmax)

	return [vmin, xmin, vmax, xmax]

def filter_stat_max(stat) :
	"""
	Filter a statitic given a sensitivity

	Keyword arguments :
	stat        -- array containing the statistic to filter_stat
	mean        -- mean value of the sweep from which is extracted the statistic
	sensitivity -- sensitivity of the detection
	"""
	#init variable
	xmin = []
	vmin = []
	xmax = []
	vmax = []
	if np.abs(stat[0]) > np.abs(stat[2]) :
		xmin.append(stat[1])
		vmin.append(stat[0])
	else :
		xmax.append(stat[3])
		vmax.append(stat[2])

	xmin = np.array(xmin)
	xmax = np.array(xmax)
	vmin = np.array(vmin)
	vmax = np.array(vmax)

	return [vmin, xmin, vmax, xmax]