from numpy import size, array, shape, resize, inf, abs


def local_extrema_seuil(sweep, seuil1, seuil2, span) :
    """
    Return the local extrema

    Return the local extrema located between the threshold (seuil) given in argument
   
    Keyword arguments:
    sweep    -- data (numpy array)
    seuil1   -- minimum value for an extrema to be extracted
    seuil2   -- maximum value for an extrema to be extracted
    span     -- mimimum separation between two extrema
   
    Returned value:
    [min, argmin, max, argmax]   
    """

    #temporary elements
    temp_min = 0
    temp_min_arg = -1
    temp_max = 0
    temp_max_arg = -1

    #This holds the result
    up_i = 0
    down_i = 0
    up = array([])
    arg_up = array([])
    down = array([])
    arg_down = array([])
    #init the writing bolean
    min_write = True
    max_write = True
    sweep_size = size(sweep)

    for i in range(sweep_size) :
        value = sweep[i]
        #check if we are below the threshold, if yes, next point
        if abs(value) < seuil1 :
            max_write = True
            min_write = True
            if temp_max_arg != -1 :
                #Reshape the array
                s_up = array(shape(up))
                s_up[0] = s_up[0] + 1
                s_up = tuple(s_up)
                up = resize(up,s_up)
                arg_up = resize(arg_up,s_up)
                #Assign values
                up[up_i] = temp_max
                arg_up[up_i] = temp_max_arg
                up_i = up_i + 1
                temp_max = 0
                temp_max_arg = -1

            if temp_min_arg != -1 :
                #Reshape the array
                s_down = array(shape(down))
                s_down[0] = s_down[0] + 1
                s_down = tuple(s_down)
                down = resize(down,s_down)
                arg_down = resize(arg_down,s_down)
                #Assign values
                down[down_i] = temp_min
                arg_down[down_i] = temp_min_arg
                down_i = down_i + 1
                temp_min = 0
                temp_min_arg = -1

            continue


        #if we are in beetween the two threshold
        if abs(value) > seuil1 and abs(value) < seuil2 :
            if value < temp_min and min_write :
                temp_min = value
                temp_min_arg = i
            if value > temp_max and max_write:
                temp_max = value
                temp_max_arg = i

        #if we are above the threshold
        if abs(value) > seuil2 :
            #Make sure than min and max cannot be accessed before going back below seuil1
            if value < - seuil2 :
                min_write = False
                if(temp_min_arg + span > i) :
                    temp_min = 0
                    temp_min_arg = -1
            if value > seuil2 :
                max_write = False
                if(temp_max_arg + span > i) :
                    temp_max = 0
                    temp_max_arg = -1

    return [down, arg_down, up, arg_up]