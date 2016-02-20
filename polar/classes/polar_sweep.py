#External library
import numpy as np
from matplotlib.pyplot import figure, subplot
#Internal library
from setproc.common.classes.measure import Measure
from setproc.sweep_set.classes.stat_point import StatPoint
from setproc.sweep_set.functions.filter import filter


#short-cut
size = np.size
linspace = np.linspace
array = np.array
matrix = np.matrix


class PolarSweep(Measure) :
    """
    This class is used to handle data stored in a json file. It takes a jsonfile as argument ang generate an object that makes data easier to manipulate
    """
    def __init__(self, filename, mode):
        Measure.__init__(self, filename, mode)
        self.f = "None"
        self.ax = "None"
        self.p = "None"
        self.cb = "None"

        try :
            self["theta"] = [self["metadata"]["theta_min"], self["metadata"]["theta_max"]]
        except :
            print "This file does not contain any information about the angle range\n Please, enter these informations"
            self["theta"] = []
            self["theta"].append(input("theta_min : "))
            self["theta"].append(input("theta_max : "))
        self.__CheckDataSize__()


    def __CheckDataSize__(self):
        """
        This method makes sure that all the sweep have the same number of points than the bias
        """
        sib = size(self["bias"])
        dsi = []
        for i in range(self["sweep_number"]) :
            sit = size(self["data"][i])
            dsi.append(abs(sit-sib))
            if (sit < sib ) :
                for j in range(sib - sit) :
                    #recopie la derniere valeur pour completer
                    self["data"][i].append(self["data"][i][sit-j-1])
            if (sit > sib ) :
                #recopie la derniere valeur pour completer
                self["data"][i] = self["data"][i][0:sib]
        print "Maximum points modified ----->  " , max(dsi)


    def MapPhase(self, delta, offset):
        #Si = size(self["bias"])
        th = linspace(self["theta"][0]+delta, self["theta"][1]+delta, self["sweep_number"])
        r = array(self["bias"])
        self.f = figure()
        self.ax = subplot(111, projection = 'polar')
        self.p = self.ax.imshow(th, r, array(matrix(array(self["data"])-offset).transpose()))
        self.p.set_cmap("jet")
        self.ax.set_rmax(abs(max(self["bias"])))
        self.cb = self.f.colorbar(self.p)

    def GetJump(self, nbr, i_start, w=4, power=1, sw=1, si = "None") :
        up = StatPoint()
        down = StatPoint()
        if si == "None" :
            si = np.size(self["bias"])
        bsweep = self["bias"][i_start+(w-1)+sw/2:si-(w-1)-sw/2+1-sw%2]
        temp_array = np.array(self["data"][nbr][i_start:], dtype = np.float)
        jump = filter(temp_array, w, power, sw)
        max_jump = jump.max()
        min_jump = jump.min()
        #construct up
        up.field = bsweep[jump.argmax()]
        up.value = max_jump
        up.up = True
        up.sweep_nbr = nbr
        #construct down
        down.field = bsweep[jump.argmin()]
        down.value = min_jump
        down.up = False
        down.sweep_nbr = nbr
        return down, up
