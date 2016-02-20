#External library
from numpy import size, array, sign, linspace, matrix
from matplotlib.pyplot import figure, subplot, polar
#Internal library
from setproc.polar.classes.polar_sweep import PolarSweep

class FullPolar(dict) :
    """
    This function handles the Polar plot
    """
    def __init__(self, trace, retrace, mode) :
        dict.__init__(self)
        self.trace = PolarSweep(trace, mode)
        self.retrace = PolarSweep(retrace, mode)
        self.p = "None"
        self.f = "None"
        self.cb = "None"
        self.ax = "None"
        self.__DoDifference__()
        self["detection"] = []

    def __DoDifference__(self) :
        sit = size(self.trace["bias"])
        sir = size(self.retrace["bias"])
        if(sit != sir) :
            print "Size problem !! Are you sure that the data are compatible"
        else :
            self["data"] = []
            for i in range(self.trace["sweep_number"]) :
                YT = array(self.trace["data"][i])
                YR = self.retrace["data"][i]
                YR.reverse() #trace starts at min and retrace et max so it has to be reversed
                self["data"].append(YT -array(YR))
                YR.reverse() #put the sweeps back to normal
                for j in range(sit) :
                    self["data"][i][j] = self["data"][i][j]*sign(self.trace["bias"][j])
        return True


    def MapPhase(self, delta, offset):
        th = linspace(self.trace["theta"][0]+delta, self.trace["theta"][1]+delta, self.trace["sweep_number"])
        r = array(self.trace["bias"])
        self.f = figure()
        self.ax = subplot(111, projection = 'polar')
        self.p = self.ax.pcolormesh(th, r, array(matrix(array(self["data"])-offset).transpose()))
        self.p.set_cmap("jet")
        self.ax.set_rmax(abs(max(self.trace["bias"])))
        #self.ax.yaxis.set_major_locator(MaxNLocator(3))
        self.cb = self.f.colorbar(self.p)

    def GetStat(self, i_start, w, power=1, sw=1):
        self["detection"] = []
        self["th"] = []
        th = linspace(self.trace["theta"][0], self.trace["theta"][1], self.trace["sweep_number"])
        sweep_number = max(self.trace['sweep_number'], self.retrace['sweep_number'])
        for i in range(sweep_number):
            #trace
            Down, Up = self.trace.GetJump(i, i_start, w, power, sw)
            self["detection"].append(Up.field)
            self["th"].append(th[i])
            self["detection"].append(Down.field)
            self["th"].append(th[i])
            #retrace
            Down, Up = self.retrace.GetJump(i, i_start, w, power, sw)
            self["detection"].append(Up.field)
            self["th"].append(th[i])
            self["detection"].append(Down.field)
            self["th"].append(th[i])

    def PlotStat(self, marker_type, marker_size, color = "black", delta = 0):
        polar(array(self["th"]) + delta, self["detection"], marker_type, markersize = marker_size, color = color)
