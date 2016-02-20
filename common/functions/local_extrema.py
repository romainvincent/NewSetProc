from numpy import size, array, shape, resize, inf, asarray, nonzero
from scipy.ndimage import maximum_filter1d,minimum_filter1d


def increment_array_size(to_extend):
    """ Increment the size of a numpy array by one"""
    shape_a = array(shape(to_extend))
    shape_a[0] = shape_a[0] + 1
    shape_a = tuple(shape_a)
    to_extend = resize(to_extend, shape_a)
    return to_extend

def local_extrema(sweep, span) :
    """    
    Find the extrema of a function

    Keyword arguments :
    sweep  -- array of data
    span   -- size of the windows for the local extrema

    Return value :
    four numpy array [down, arg_down, up, arg_up]

    Note :
    Problem with large span.. peaks at the end of the sweep are discarded.. to be corrected!!
    A Cython version may be available (or will be). Check in cfunctions !

    """

    #temporary elements
    temp_min = inf
    temp_min_arg = -1
    temp_max = -inf
    temp_max_arg = -1

    #This holds the result
    up_i = 0
    down_i = 0
    up = array([])
    arg_up = array([])
    down = array([])
    arg_down = array([])
    #init the writing bolean
    sweep_size = size(sweep)

    for i in range(1, sweep_size-1) :
        value = sweep[i]
        #check if it is max or min
        if value < temp_min :
            if value <= sweep[i-1] and value < sweep[i+1] :
                temp_min_arg = i
                temp_min = value

        if value > temp_max :
            if value >= sweep[i-1] and value > sweep[i+1] :
                temp_max_arg = i
                temp_max = value

        if (i - temp_min_arg) > span and temp_min_arg != -1 :
            #Reshape the array
            down = increment_array_size(down)
            arg_down = increment_array_size(arg_down)
            #Assign values
            down[down_i] = temp_min
            arg_down[down_i] = temp_min_arg
            down_i = down_i + 1
            temp_min = inf
            temp_min_arg = -1

        if (i - temp_max_arg) > span and temp_max_arg != -1 :
            #Reshape the array
            up = increment_array_size(up)
            arg_up = increment_array_size(arg_up)
            #Assign values
            up[up_i] = temp_max
            arg_up[up_i] = temp_max_arg
            up_i = up_i + 1
            temp_max = -inf
            temp_max_arg = -1
    return [down, arg_down, up, arg_up]
  
    
 
def local_extrema_2(arr, min_distance = 20):
    """Find all local maxima of the array, separated by at least min_distance."""
    cval = 0
    mode = 'constant'
    cval = arr.max()+1
    max_points = arr == maximum_filter1d(arr, min_distance)
    min_points = arr == minimum_filter1d(arr, min_distance)
    

    return [arr[nonzero(min_points==True)[0]],nonzero(min_points==True)[0], 
            arr[nonzero(max_points==True)[0]],nonzero(max_points==True)[0]]
