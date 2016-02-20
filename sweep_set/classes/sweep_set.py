from copy import deepcopy
from setproc.common.functions.local_extrema import local_extrema
from setproc.common.classes.open_json import OpenJson
from setproc.common.classes.open_bin import OpenBin
from setproc.sweep_set.classes.stat_point import StatPoint
import setproc.common.cfunctions.extract_stat as exst
from setproc.sweep_set.functions.merge_gb import merge_gb
from setproc.sweep_set.functions.filter import filter
from matplotlib.pyplot import figure, draw, get_fignums,title, hist, ginput, plot, scatter, xlabel, ylabel,xlim, ylim,text,tick_params, show

import numpy as np
import cPickle



class SweepSetOpen(dict) :
    """
    It generated an sweep_set_open object. Due to the huge data it can contains, the usual saving systeme can be replaced by npz format faster and more compact at a price of less versatility in what is saved (only "bias", "sweep_nbr" and "data")
    """
    def __init__(self, filename, mode) :
        dict.__init__(self)
        if mode == "npz" :
            temp = np.load(filename)
            for x in temp.files :
                self[x] = deepcopy(temp[x])
            temp.close()
            del(temp)
        else :
            if mode == "Json" :
                temp = OpenJson(filename)
            elif mode == "Bin" :
                temp = OpenBin(filename)
            for x in temp :
                self[x] = deepcopy(temp[x])
            del(temp)

    def GetLocalMin(self, nbr, i_start, span, w = 4, power = 1, sw = 1):
        si = np.size(self["bias"])
        bsweep = self["bias"][i_start+(w-1)+sw/2 : si - (w-1) -sw/2 +1 -sw%2]
        temp_array = np.array(self["data"][nbr][i_start:],dtype = np.float)
        temp_array = temp_array / np.mean(temp_array)
        jump = filter(temp_array,w,power,sw)
        jump = si * jump / np.sum(jump)
        [up, arg_up, down, arg_down] = local_extrema(jump, span)
        result_up = [[],[]]
        result_down = [[],[]]
        up = up.tolist()
        i = 0
        for x in up :
            result_up[0].append(x)
            result_up[1].append(bsweep[np.int(arg_up[i])])
            i = i+1
        down = down.tolist()
        i = 0
        for x in down :
            result_down[0].append(x)
            result_down[1].append(bsweep[np.int(arg_down[i])])
            i = i+1

        return result_up, result_down


    def GetJump(self, nbr, i_start, w=4, power=1, sw=1, si = "None", 
                mode ="classic", seuil1 = 0, seuil2 =0, span = 50) :
        up = StatPoint()
        down = StatPoint()
        if si == "None" :
            si = np.size(self["bias"])
        bsweep = self["bias"][i_start+(w-1)+sw/2:si-(w-1)-sw/2+1-sw%2]
        temp_array = np.array(self["data"][nbr][i_start:], dtype = np.float)
        #temp_array = temp_array / np.mean(temp_array)
        jump = filter(temp_array, w, power, sw)
        #jump = si * jump / np.sum(jump)
        if mode == "classic" :
            if( nbr == 0) :
                print("In classic mode")
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
        else :
            down.value, min_arg , up.value, max_arg = exst.extract_stat(jump, seuil1, seuil2, span)
            up.up = True
            up.sweep_nbr = nbr
            if max_arg < 0 :
                up.field = 0
            else :
                up.field = bsweep[max_arg]
                down.up = False
                up.sweep_nbr = nbr
            if min_arg < 0 :
                down.field = 0
            else :
                down.field = bsweep[min_arg]
        return down, up

    

    def SanityCheck(self) :
        size_bias = np.size(self["bias"])
        nbr = None
        temp = 0
        for i in range(self["sweep_number"]) :
            size_data = np.size(self["data"][i])
            if (size_data < size_bias ) :
                nbr = i
                for j in range(size_bias - size_data) :
                    #recopie la derniere valeur pour completer
                    self["data"][i].append(self["data"][i][size_data-j-1])

            if (size_data > size_bias ) :
                #recopie la derniere valeur pour completer
                self["data"][i] = self["data"][i][0:size_bias]
                nbr = i
            if (abs(size_data - size_bias) > temp ) :
                temp = abs(size_data - size_bias)

        if(nbr == None) :
            print "\tPerfect matching"
        else :
            print "\tThe maximum difference is ", temp, " points with the sweep ", nbr

        return True


    def __add__(self, other) :
        return merge_gb([self, other])

    ##To plot and acess
    def PlotCurve(self, nbr, i_start, w=4, pw = 1, sw = 1) :
        try :
            self.fig.clear()
            self.ax1.clear()
            self.ax2.clear()
        except :
            self.fig = figure()
        X = self["bias"][i_start:]
        Y = self["data"][nbr][i_start:]
        si = np.size(self["bias"][i_start:])
        self.ax1 = self.fig.add_subplot(211)
        self.ax1.plot(X, Y)
        self.ax1.set_xlim(X[0], X[si-1])
        self.ax2 = self.fig.add_subplot(212)
        self.ax2.plot(X[(w-1)+sw/2:si-(w-1)-sw/2 +1 -sw%2], filter(np.array(Y), w, pw, sw))
        self.ax2.set_xlim(X[0], X[si-1])
        for i in get_fignums() :
            figure(i)
            draw()
        show()
        return True             
    
    def GetSweep(self, nbr) :
        return [self["bias"], self["data"][nbr]]


    def Save(self, savename) :
        """
        Save all the keys field of a dictionnary in savename file. The "*.bin" extension should be used
        """
        done = True
        try :
            stream = open(savename, "w")
        except IOError :
            print "Problem while saving the file"
            done = False
        l = []
        temp = self.keys()
        l.append(temp)
        for x in self:
            l.append(self[x])
        cPickle.dump(l, stream, 1)
        return done

    def Savez(self, filename) :
        np.savez(filename, data = self["data"], bias= self["bias"], sweep_number = self["sweep_number"])

    def GetPeaks(self,thres,offset=0,mode='max',x_rge=False,i_start=0,w=4,pw=1,sw=17,span=20):
        """
        if mode = max : GetPeaks returns the max up or down jump 
                        if they are above a certain threshold
        """
        if i_start >= 0:
            i_stop = -1       
        else:
            i_stop = i_start
            i_start = 0;
            
        if x_rge == False:
            x_rge=[0,self["sweep_number"]-1]
        
        if (self["bias"][-1]-self["bias"][0])>0:
            trace = True
        else:
            trace = False
        
        si = np.size(self["bias"])
        swp = np.array(self["bias"][i_start+(w-1)+sw/2:si-(w-1)-sw/2+1-sw%2+i_stop])
        swp_x  = np.array(np.linspace(x_rge[0],x_rge[1],self['sweep_number']))
       
        jump_up = []
        jump_down = []

        for nbr in range(self['sweep_number']):
            temp_array = np.array(self["data"][nbr][i_start:i_stop], dtype = float)
            
            if mode == 'all':
                extr = np.array(local_extrema(filter(temp_array, w, pw, sw),span)) 
                j_up   = extr[3][:,abs(extr[2])>thres] #indices of up jumps
                j_down = extr[1][:,abs(extr[0])>thres] #indices of down jumps                          
                
                for jj in range(np.size(j_up)):  
                    jump_up.append([ swp_x[nbr],swp[j_up[jj]]+offset])
                
                for kk in range(np.size(j_down)):    
                    jump_down.append([ swp_x[nbr],swp[j_down[kk]]+offset])
                 
            elif mode =='max':
                jump = filter(temp_array, w, pw, sw)
                max_jump = max(jump)
                min_jump = min(jump)
                
                if abs(max_jump) > abs(min_jump):
                    if abs(max_jump)>thres:
                        #print(nbr,max_jump)
                        jump_up.append([ swp_x[nbr],swp[jump.argmax()]+offset,max_jump])
                    
                else:
                    if abs(min_jump) > thres:
                        #print(nbr,min_jump)
                        jump_down.append([ swp_x[nbr],swp[jump.argmin()]+offset,min_jump])
            
            elif mode =='number':
                extr = np.array(local_extrema(filter(temp_array, w, pw, sw),span)) 
                j_up   = extr[3][:,abs(extr[2])>thres] #indices of up jumps
                j_down = extr[1][:,abs(extr[0])>thres] #indices of down jumps                          
                
                jump_up.append(np.size(j_up)+np.size(j_down))
                     
                
        jump_down = np.array(jump_down)
        jump_up = np.array(jump_up)
        
        return(jump_up,jump_down)


