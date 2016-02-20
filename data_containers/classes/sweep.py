#General library
import numpy as np

#Home made library
from setproc.common.cfunctions.moving_average import moving_average_c as mva
from setproc.active_sweep.functions.derivate import derivate
from setproc.sweep_set.functions.filter import filter
from setproc.common.functions.local_extrema import local_extrema
from setproc.common.functions.local_extrema_seuil import local_extrema_seuil
from setproc.data_containers.functions.filter_stat import filter_stat_threshold, filter_stat_max
from setproc.data_containers.functions.get_argument import get_arguments

class Sweep():

    """
    Return a Sweep object

    It contains a set of data on the form y = f(x).

    """

    def __init__(self):
        self.x = "None"
        self.y = "None"
        self.sweep_size = 0
        self.sweep_mean = "None"

    ####################
    #statistic functions

    def Filter(self, i_start, width, power, width2, normalized = True):
        """
        Apply a filter which convert steps into peaks.

        Keywords arguments:
        i_start    -- number of points truncated to the sweep before applying the filter
        width      -- width in points of the first mowing average windows
        power      -- power applied after filtering (should be soon removed)
        width2     -- width of the post filtering moving average
        normalized -- if set True, the signal is renormalized (default True)

        Returned value :
        Sweep object corresponding to the filtered signal.

        """
        bias = self.x.copy()
        bias = bias[i_start:]
        sweep = self.y.copy()
        sweep = sweep[i_start:]
        si = np.size(sweep)
        #apply the filter and resize the bias
        temp_x = bias[(width-1)+width2/2 : si - (width-1) -width2/2 +1 -width2%2]
        temp_y = filter(sweep, width, power, width2)
        if normalized :
                temp_y = si * temp_y / np.sum(temp_y)
        result = Sweep()
        result.LoadArrays(temp_x, temp_y)
        return result

    def GetExtrema(self):
        """Return the min and max of the sweep."""
        vmin = self.y.min()
        vmax = self.y.max()
        argmin = self.y.argmin()
        argmax = self.y.argmax()
        toA = np.array
        return [toA(vmin), toA(self.x[argmin]),toA(vmax), toA(self.x[argmax])]

    def GetFilteredExtrema(self, sensitivity=10, mode="max"):
        """
        Return a filtered min and max of a sweep

        Keyword argument :
        sensitivity -- number defining of many times an extrema has to exceed the mean value to be conserved
        """
        tofilter = self.GetExtrema()
        if mode == "max" :
            result = filter_stat_max(tofilter)
        else :
            result = filter_stat_threshold(tofilter, self.sweep_mean, sensitivity)
        return result


    def GetLocalExtrema(self,span):
        """
        Return the local extrema if the sweep.

        Keywrod arguments:
        span -- minimum distance between two extrema

        Returned value
        [min, argmin, max, argmax]
        """
        vmin, argmin, vmax, argmax = local_extrema(self.y,span)
        xmin = []
        xmax = []
        argmin = argmin.tolist()
        argmax = argmax.tolist()
        for i in argmin :
            xmin.append(self.x[i])
        for i in argmax :
            xmax.append(self.x[i])
        xmax = np.array(xmax)
        xmin = np.array(xmin)

        return [vmin,xmin,vmax,xmax]

    def GetFilteredLocalExtrema(self, span, sensitivity):
        """
        Return a filtered local min and max of a sweep

        Keyword argument :
        span        -- minimum distance between two extrema
        sensitivity -- number defining of many times an extrema has to exceed the mean value to be conserved

        """
        tofilter = self.GetLocalExtrema(self, span)
        result = filter_stat(tofilter, self.sweep_mean, sensitivity)
        return result


    def GetLocalExtremaSeuil(self, span, seuil1, seuil2):
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
        vmin, argmin, vmax, argmax = local_extrema_seuil(self.y, seuil1, seuil2, span)
        xmin = []
        xmax = []
        argmin = argmin.tolist()
        argmax = argmax.tolist()
        for i in argmin :
            xmin.append(self.x[i])
        for i in argmax :
            xmax.append(self.x[i])
        xmax = np.array(xmax)
        xmin = np.array(xmin)

        return [vmin,xmin,vmax,xmax]


    ##########################
    #modification of the sweep

    def SetOffset(self, offset):
        """Add the offset to the bias """
        self.x = self.x - offset
        return True

    def Resize(self,xmin,xmax):
        """Resize the sweep"""
        i_start, i_end = get_arguments(self.x, xmin, xmax)
        self.x = self.x[i_start:i_end]
        self.y = self.y[i_start:i_end]
        return True

    #######################
    #mathematical functions

    def Function(self, function):
        """
        Return the result of a function apllied to a sweep.

        Keywords arguments:
        function -- a function which can be applied to a sweep.

        Returned values:
        Depends on the function given as argument

        Notes:
        Given a Sweep Sw, the function as to take Sw.x as first argument
        and Sw.y as second one.

        """
        return function(self.x, self.y)

    def Derivative(self):
        """
        Return the derivate the Sweep

        Return values :
        A Sweep corresponding to the derivative of the original Sweep.

        """
        temp_x , temp_y = self.Function(derivate)
        result = Sweep()
        result.LoadArrays(temp_x, temp_y)
        return result

    def MovingAverage(self, width):
        """
        Return the moving average of the Sweep

        Keyword arguments :
        widht -- width of the moving average window.

        Returned value:
        A Sweep containing the moving average of the original Sweep.

        """
        temp_y = mva(self.y, width)
        temp_x = self.x[(width-1)/2 : self.sweep_size - (width-1)/2 - (width-1)%2]
        result =  Sweep()
        result.LoadArrays(temp_x, temp_y)
        return result

    ######################
    #Elementary operations

    def __sub__(self, b):
        """
        This function return a Sweep with the same bias but the difference between the two set of data
        """
        temp_y = self.y - b.y
        result = Sweep()
        result.LoadArrays(self.x, temp_y)
        return result

    def __add__(self, b):
        """
        This function return a Sweep with the same bias but the sum of the two set of data
        """
        temp_y = self.y - b.y
        result = Sweep()
        result.LoadArrays(self.x, temp_y)
        return result

    ########
    ##Saving

    #saving in npz format
    def SaveNpz(self, filename):
        """
        Save the Sweep in the npz format

        Keywords arguments :
        filename -- filename in string format without extension.

        """
        np.savez(filename, x=self.x, y=self.y)
        return True

    #saving in txt format
    def SaveTxt(self, filename, delimiter ='\t', newline = '\n'):
        """
        Save the Sweep in txt format

        Keywords arguments:
        delimiter -- delimiter between the values in the file (default '\t')
        newline   -- newline format (default '\n')

        Notes:
        The sweep is saved in a two columns format
        """
        np.savetxt(filename, np.array([self.x, self.y]).transpose(), delimiter = delimiter, newline = newline)

    ########
    #Loading

    #load from arrays
    def LoadArrays(self, x, y):
        """
        Load arrays into the Sweep

        Keywords arguments:
        x -- array (or list) corresponding to the x axis
        y -- array (or list) corresponding to the y axis

        Notes :
        Initiliaze the sweep_size parameter

        """
        self.x = np.array(x)
        self.y = np.array(y)
        self.sweep_size =  np.size(self.x)
        self.sweep_mean = np.mean(self.y)

    #loading from npz file
    def LoadNpz(self, filename):
        """
        Load data from a npz file

        Keywords arguments:
        filename -- name of the file containing the data

        Notes :
        Initialiaze the sweep_size parameter
        """
        temp = np.load(filename)
        self.x = temp["x"]
        self.y = temp["y"]
        self.sweep_size = np.size(self.x)
        self.sweep_mean = np.mean(self.y)
        return True

    #loadind from text file
    def LoadTxt(self, filename, delimiter ='\t') :
        """
        Load data from a txt file

        Keywords arguments:
        filename  -- name of the file containing the data
        delimiter -- delimter between the value in the text file (default '\t')

        Notes :
        Initialiaze the sweep_size parameter
        """
        temp = np.load(filename, delimiter =  delimiter)
        temp = temp.transpose()
        self.x = temp[0]
        self.y = temp[1]
        self.sweep_size = np.size(self.x)
        self.sweep_mean = np.mean(self.y)
