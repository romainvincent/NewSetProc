#external library
from numpy import size, array
#internal library
from setproc.common.functions.json_handling import json_data
from setproc.common.functions.open_files import get_json
from setproc.data_containers.classes.sweep import Sweep

class JsonToSweep(dict):
    """
    Load Json files

    Parse a Json file and convert each measure to a Sweep object. It returns and object containing the information related to the measurements.

    Keyword argument :
    filename -- name of the json file (string)

    Returned value :
    An object subclass of dict which contains all the data and the informations concerning the measurement.

    Notes :
    Should not be used directly. Users should prefer SweepSet.

    """

    def __init__(self, filename) :
        dict.__init__(self)
        self.LoadJson(filename)

    def LoadJson(self, filename) :
        """
        Load Json file

        Keyworg argument
        filename -- name of the json file

        Returned value
        done -- True if the file has been properly loaded, False otherwise
        """

        done = True
        try :
            monjson = get_json(filename)
        except IOError :
            print "Problem loading the file"
            done = False

        if(done) :
            done = self._ExtractJsonData_(monjson)

        return done

    def _ExtractJsonData_(self, monjson) :
        """
        Extract data from the json object

        NOT TO BE USED DIRECTLY!!!
        """
        done = True
        try :
            K = monjson.keys()
            K.pop(K.index("vim_modeline")) #eliminate vim mode parameter

            #the function json_data(json_object,i,j) return the jth column of the ith measurement of a json_object
            self["sweep_number"] = size(monjson["measures"])
            self["bias"] = json_data(monjson, 0, 1)
            self["date"] = []
            self["sweeps"] = []
            sweep_number = self["sweep_number"]
            for i in range(sweep_number):
                self["date"].append(monjson["measures"][i]["start"])
                temp_sweep = Sweep()
                temp_x = array(self["bias"])
                temp_y = array(json_data(monjson, i, 2))
                temp_sweep.LoadArrays(temp_x, temp_y)
                self["sweeps"].append(temp_sweep)

            K.pop(K.index("measures"))
            if (size(K) > 0) :
                self["metadata"] = dict([])
                for x in K :
                    self["metadata"][x]  = monjson[x]
        except KeyError :
            print "Problem while loading the file. If an object is howerver loaded, it can be incomplete"
            done = False

        return done
