#General library
import numpy as np

#Home made library
from setproc.common.cfunctions.moving_average import moving_average_c as mva
from setproc.active_sweep.functions.derivate import derivate

class Sweep():
    """
    This object as the aim to load data and perform operations on them. It should be able to load data from
    objects or files according to the need. It should be able to release the data after doing is task.
    """

    def __init__(self):
        self.x = "None"
        self.y = "None"
        self.sweep_size = 0

    def Load(self, x, y):
        """
        Load the Sweep if the nummpy arrays given in arguments.
        """
        self.x = np.array(x)
        self.y = np.array(y)
        self.sweep_size =  np.size(self.x)

    def Function(self, function):
        """
        This method allows to easely apply a function to the sweep. The function as to take no argument
        and has to be defined as follows f(x,y).
        """
        return function(self.x, self.y)

    def Release(self):
        """
        This function release the element stored in the sweep. Should be used each time the content of
        Sweep object is no more needed.
        """
        del(self.x)
        del(self.y)
        self.x = "None"
        self.y = "None"

    def Derivative(self):
        """
        This function returns the derivative of the sweep in the from of [bias,data].
        """
        return self.function(derivate)

    def Derivate(self):
        """
        This function replace a Sweep by its derivative.
        """
        temp = self.derivative()
        self.x = temp[0]
        self.y = temp[1]

    def MovingAverage(self, width):
        """
        Compute the moving average of a sweep
        """
        temp_y = mva(self.y, width)
        temp_x = self.x[(width-1)/2 : self.sweep_size - (width-1)/2 - (width-1)%2]
        return [temp_x, temp_y]

    def MovingAveraging(self, width):
        """
        Replace a sweep by its corresponding moving average
        """
        temp = self.MovingAverage(width)
        self.x = temp[0]
        self.y = temp[1]
