import numpy as np


def get_arguments(sweep, begin, stop):
	si = np.size(sweep)
	if sweep[0] < sweep[si-1] :
		i_start, i_end = _increase_get_argument(sweep, begin, stop)
	else :
		i_start, i_end = _decrease_get_argument(sweep, begin, stop)

	return i_start, i_end

def _decrease_get_argument(sweep, begin, stop):
	si = np.size(sweep)
	detect_min = True
	detect_max = True
	i_start = 0
	i_end = si -1
	for i in range(si):
		if sweep[i] <= begin and detect_min :
			i_start = i
			detect_min = False

		if sweep[i] < stop and detect_max :
			i_end = i - 1
			detect_max = False
	return i_start, i_end

def _increase_get_argument(sweep, begin, stop):
	si = np.size(sweep)
	detect_min = True
	detect_max = True
	i_start = 0
	i_end = si -1
	for i in range(si):
		if sweep[i] >= begin and detect_min :
			i_start = i
			detect_min = False

		if sweep[i] > stop and detect_max :
			i_end = i - 1
			detect_max = False
	return i_start, i_end
