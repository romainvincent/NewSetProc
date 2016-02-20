#external library
import numpy as np
from copy import deepcopy
#internal library
from setproc.common.classes.json_to_sweeps import JsonToSweep
from setproc.data_containers.classes.sweep import Sweep

class SweepSet(dict):
    """
    Create an SweepSet object

    This object helps to handle a set of Sweep. It provides an interface to the Sweep method.

    Keyworg argument :
    filename -- name of the file containing the set of sweep to load
    filetype -- extension of the file containing the data (*.npz or *.json)

    Notes : for the moment, handle only json and npz files

    """
    def __init__(self, filename, filetype = "json") :
        dict.__init__(self)
        if filetype == "json":
            self._LoadJson_(filename)
        elif filetype == "npz":
            self._LoadNpz_(filename)
        elif filetype == "obj":
            self._LoadSweepSet_(filename)
        else:
            print("filetype not valid")

    ###############
    #stat functions

    def GetExtrema(self) :
        """ Extract the extrema of all the sweep"""
        result = []
        for i in range(self["sweep_number"]):
            temp_sweep = self["sweeps"][i]
            temp_result = temp_sweep.GetExtrema()
            result.append(temp_result)
        return result


    def GetFilteredExtrema(self, sensitivity=10, mode="max") :
        """
        Return a filtered min and max of all the sweeps

        Keyword argument :
        sensitivity -- number defining of many times an extrema has to exceed the mean value to be conserved
        """
        result = []
        for i in range(self["sweep_number"]):
            temp_sweep = self["sweeps"][i]
            temp_result = temp_sweep.GetFilteredExtrema(sensitivity,mode)
            result.append(temp_result)
        return result
    
    def GetLocalExtrema(self, span) :
        """ 
        Extract the extrema of all the sweep

        Keyword argument:
        span -- minimum distance between two extrema

        """
        result = []
        for i in range(self["sweep_number"]):
            temp_sweep = self["sweeps"][i]
            temp_result = temp_sweep.GetLocalExtrema(span)
            result.append(temp_result)
        return result


    def GetFilteredLocalExtrema(self, span , sensitivity) :
        """
        Return a filtered min and max of all the sweeps

        Keyword argument :
        span        -- minimum distance between two extrema
        sensitivity -- number defining of many times an extrema has to exceed the mean value to be conserved
        """
        result = []
        for i in range(self["sweep_number"]):
            temp_sweep = self["sweeps"][i]
            temp_result = temp_sweep.GetFilteredLocalExtrema(span, sensitivity)
            result.append(result)
        return result


    def GetFilter(self, i_start, width, width2, power=1, normalized=True):
        """
        Apply a filter to all the sweep which convert steps into peaks.

        Keywords arguments:
        i_start    -- number of points truncated to the sweep before applying the filter
        width      -- width in points of the first mowing average windows
        power      -- power applied after filtering (should be soon removed)
        width2     -- width of the post filtering moving average
        normalized -- if set True, the signal is renormalized (default True)

        """

        for i in range(self["sweep_number"]):
            tp_sw = self["sweeps"][i]
            self["sweeps"][i] = tp_sw.Filter(i_start, width, power, width2, normalized)
        return True


    #################
    #basic operations

    def __sub__(self, sweep_set):
        """
        Implement the substraction between SweepSet object

        Keyword argument:
        sweep_set -- SweepSet object

        Return value :
        SweepSet object

        Exemple :
        a.__sub__(b) --> a - b

        """
        temp_sweep_set = SweepSet(self, filetype="obj")
        si = temp_sweep_set.GetSize()
        for i in range(si):
            temp_sweep1 = temp_sweep_set["sweeps"][i]
            temp_sweep2 = sweep_set["sweeps"][i]
            temp_sweep_set["sweeps"][i] = temp_sweep1 - temp_sweep2
        return temp_sweep_set


    def GetSize(self):
        """Return the number of sweep """
        return self["sweep_number"]

    ###########################
    #modification of the sweeps

    def SetOffset(self, offset) :
        """
        Set a bias offset for all the sweeps

        Keyword argument :
        offset -- value to be subtracted to the sweep bias

        """
        si = self.GetSize()
        for i in range(si) :
            self["sweeps"][i].SetOffset(offset)
        return True

    def Resize(self, xmin, xmax):
        """Resize the sweeps """
        si = self.GetSize()
        for i in range(si) :
            self["sweeps"][i].Resize(xmin,xmax)

            
    ########
    #Loading

    def _LoadJson_(self, filename):
        """
        Load a Json file

        Keyword argument:
        filename -- name of the json file

        Note :
        Should not be used directly
        """
        temp = JsonToSweep(filename)
        for x in temp :
            self[x] = deepcopy(temp[x])
        del(temp)

        return True

    def _LoadNpz_(self, filename):
        """
        Load a npz file

        Keyword argument:
        filename -- name of the npz file

        Note :
        Should not be used directly !
        """
        self["sweeps"] = []
        temp = np.load(filename)
        temp_x = temp["x"]
        temp_y = temp["y"]
        si = len(temp_x)
        for i in range(si):
            temp_sweep = Sweep()
            temp_sweep.LoadArray(temp_x[i], temp_y[i])
            self["sweeps"].append(temp_sweep)

        return True

    def _LoadSweepSet_(self, sweep_set):
        """Copy an existing SweepSet object"""
        for x in sweep_set :
            self[x] = deepcopy(sweep_set[x])
        return True

    #######
    #saving

    def SaveNpz(self, filename):
        """
        Save in npz format

        The npz format is reallt efficient for loading and space saving. Should be always used!

        Keyword argument :
        filename -- name of the file

        Notes :
        Data can be reload using filetype = "npz" when calling SweepSet constructor.
        """
        temp_x = []
        temp_y = []
        for i in range(self["sweep_number"]) :
            temp_sweep = self["sweeps"][i]
            temp_x.append(temp_sweep.x)
            temp_y.append(temp_sweep.y)
        temp_x = np.array(temp_x)
        temp_y = np.array(temp_y)
        np.savez(filename, x=temp_x, y=temp_y)
        return True



    ############
    # DEPRECATED

    def GetStat(self, i_start, width1, width2, 
        normalized=True, power=1, stat_type ="global", 
        span=50, seuil=False, seuil1=0, seuil2=1):
        """
        Extract the statistics

        This function extract the abrupt changes in a sweep. Only the greatest is kept (one for positive change and one for negative one)

        Keyword arguments :
        i_start    -- number of points truncated to the sweep before applying the filter
        width1     -- width in points of the first mowing average windows
        width2     -- width of the post filtering moving average
        power      -- power applied after filtering (should be soon removed)
        normalized -- if set True, the signal is renormalized (default True)
        stat_type  -- if set "global" uses the global extrema, if set "local" use the local one (default "global")
        span       -- if stat_type set to "local", define the span between extrema

        Returned values :
        list of the change value and position in bias in the form [min,argmin,max,argmax]
        """

        print("deprecated!!! You should used ont of the GetExtrema method\n")
        
        result = []
        for i in range(self["sweep_number"]) :
            temp_sweep = self["sweeps"][i]
            temp_sweep = temp_sweep.Filter(i_start, width1, power, width2, normalized)
            if stat_type == "global" :
                temp_result = temp_sweep.GetExtrema() #return value : [min,argmin,max,argmax]
            elif stat_type == "local" and seuil == False :
                temp_result = temp_sweep.GetLocalExtrema(span) #return value : [min,argmin,max,argmax]
            elif stat_type == "local" and seuil == True :
                temp_result = temp_sweep.GetLocalExtremaSeuil(span, seuil1, seuil2) #return value : [min, argmin, max, argmax]
            else :
                print("unknown stat_type !(can be only \"global\" or \"local\" ")
            result.append(temp_result)
        return result
