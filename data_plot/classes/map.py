#external library
import numpy as np
from matplotlib.pyplot import imshow, figure, draw
#internal library
from setproc.map.functions.map_functions import plot_profile_h, plot_profile
from setproc.data_containers.classes.sweep_set import SweepSet



class Map(dict):
    """
    Create a Map object

    This object is used whenever a SweepSet can be represented by a colormap

    Keyword arguments
    data -- either a filename or a SweepSet object
    type -- either the filename extension or "obj" in case of SweepSet object

    Note : for the moment only *.npz and *.json files are supported

    """

    def __init__(self, data, data_type="json"):
        dict.__init__(self)

        self.fig = "None" #contains the figure related to the plot
        self.im = "None"  #contains the image related to the plot
        self.ax = "None"  #contains the axes related to the plot

        if data_type == "obj" :
            self.sweeps = data
        else :
            self.sweeps = SweepSet(data, data_type)

    def SetExtent(self, extent):
        """
        Define the extent of the map

        Keyword argument :
        extent -- array of the form [xmin, xmax, ymin, ymax]

        Notes :
        ymin and ymax can be different from the one of the sweeps.. be carreful !!
        """
        self["extent"] = extent
        try :
            self.im.set_extent(extent)
            self.ax.set_aspect("auto")
            draw()
        except AttributeError :
            print("You habe to create a plot first")

    def SetXLabel(self, label):
        """ Set the xlabel of the plot """
        self["xlabel"] =  label
        try :
            self.ax.set_xlabel(label)
            draw()
        except AttributeError:
            print("You habe to create a plot first")

    def SetYLabel(self, label):
        """ Set the ylabel of the plot """
        self["ylabel"] =  label
        try :
            self.ax.set_ylabel(label)
            draw()
        except AttributeError:
            print("You habe to create a plot first")

    def PlotMap(self):
        """Plot the colormap """
        #etxract data
        data = []
        si = self.sweeps.GetSize()
        for i in range(si) :
            temp_sweep = self.sweeps["sweeps"][i]
            data.append(temp_sweep.y)
        data = np.matrix(data)
        data = data.transpose()
        data = data.tolist()
        #plot data
        self.fig = figure()
        self.ax = self.fig.add_subplot(111)
        self.im = imshow(data, interpolation="nearest", origin="lower")


    def PlotExtrema(self, newplot=True,other_axes="None", color = ["blue", "red"]):
        """Plot the extream position"""
        toplot = self.sweeps.GetExtrema()
        si = self.sweeps.GetSize()
        proceed = False
        if newplot :
            fig = figure()
            ax = fig.add_subplot(111)
        elif other_axes != "None":
            ax = other_axes
        else :
            ax = self.ax

        try :
            xmin = self["extent"][0]
            xmax = self["extent"][1]
            proceed = True
        except KeyError:
            print("you have to defined extent first")

        if(proceed):
            bias = np.linspace(xmin, xmax, si)
            x =[]
            ymin = []
            ymax = []
            for i in range(si):
                x.append(bias[i])
                ymin.append(toplot[i][1])
                ymax.append(toplot[i][3])
            ax.plot(x, ymin, ".", color = color[0], markersize=1)
            ax.plot(x, ymax, ".", color = color[1], markersize=1)

    def PlotFilteredExtrema(self, sensitivity=10, newplot=True, 
                            other_axes="None",mode ='max',
                            color =["blue","red"]):
        """
        Plot the filtered extrema
        """
        toplot = self.sweeps.GetFilteredExtrema(sensitivity, mode)
        si = self.sweeps.GetSize()
        proceed = False
        if newplot :
            fig = figure()
            ax = fig.add_subplot(111)
        elif other_axes != "None":
            ax = other_axes
        else :
            ax = self.ax

        try :
            xmin = self["extent"][0]
            xmax = self["extent"][1]
            proceed = True
        except KeyError:
            print("you have to defined extent first")

        if(proceed):
            bias = np.linspace(xmin, xmax, si)
            bmin = []
            ymin = []
            bmax = []
            ymax = []
            for i in range(si):
                #min first
                if np.size(toplot[i][1]) != 0:
                    bmin.append(bias[i])
                    ymin.append(toplot[i][1])
                #then max
                if np.size(toplot[i][3]) != 0:
                    bmax.append(bias[i])
                    ymax.append(toplot[i][3])
            ax.plot(bmin, ymin, ".", color = color[0], markersize=1)
            ax.plot(bmax, ymax, ".", color = color[1], markersize=1)



    def PlotSteps(self,i_start, width1, width2, stat_type="global",
                 span=50, normalized=True, power=1, newplot=True, 
                 other_axes="None",seuil=False,seuil1=0,seuil2=1):
        """
        Plot the step positions in the set of sweep

        Keyword arguments :
        i_start    -- number of points truncated to the sweep before applying the filter
        width1     -- width in points of the first mowing average windows
        width2     -- width of the post filtering moving average
        normalized -- if set True, the signal is renormalized (default True)
        power      -- power applied after filtering (default 1 -- should be soon removed)
        newplot    -- if True plotted in a different plot (default True)
        other_axes -- axes, if newplot set to false, the curve will be plot using these axes
        stat_type  -- if set "global", use the global extrema, if set to local, use the local ones
        span       -- if stat_type set to "local", define the span between two extrema
        seuil      -- if True, the local minimum use the seuil parameters (default False)
        seuil1     -- minimum value for a jump to be taken (default 0)
        seuil2     -- maximum value for a jump to be taken (default 1)
        """
        toplot = self.sweeps.GetStat(i_start, width1, width2,
                                    normalized, power,
                                    stat_type, span,
                                    seuil, seuil1, seuil2)
        si = self.sweeps.GetSize()
        proceed = False

        if  newplot :
            fig = figure()
            ax = fig.add_subplot(111)
        elif other_axes != "None" :
            ax = other_axes
        else :
            ax = self.ax

        try :
            xmin = self["extent"][0]
            xmax = self["extent"][1]
            proceed = True
        except KeyError:
            print("You have to define the exetent first")

        if(proceed and stat_type == "global"):
            bias = np.linspace(xmin, xmax, si)
            x = []
            ymin = []
            ymax = []
            for i in range(si) :
                x.append(bias[i])
                ymin.append(toplot[i][1])
                ymax.append(toplot[i][3])
            ax.plot(x, ymin, ".b", markersize=1)
            ax.plot(x, ymax, ".r", markersize=1)

        elif proceed and stat_type == "local" :
            bias = np.linspace(xmin, xmax, si)
            xmin = []
            xmax = []
            ymin = []
            ymax = []
            for i in range(si):
                si_min = np.size(toplot[i][1])
                si_max = np.size(toplot[i][3])
                for j in range(si_min):
                    ymin.append(toplot[i][1][j])
                    xmin.append(bias[i])
                for j in range(si_max):
                    ymax.append(toplot[i][3][j])
                    xmax.append(bias[i])
            ax.plot(xmin, ymin, ".b", markersize=1)
            ax.plot(xmax, ymax, ".r", markersize=1)



    def PlotVProfile(self):
        """
        Plot a vertical profile of the Map

        Select a point on the colormap for which you want a profile

        """
        print('select a point on the colormap')
        return plot_profile(self.im)

    def PlotHProfile(self):
        """
        Plot a horizontal profile of the Map

        Select a point on the colormap for which you want a profile

        """
        print('select a point on the colormap')
        return plot_profile_h(self.im)        