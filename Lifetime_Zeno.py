# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 19:07:03 2012

@author: stefan
"""
import sys
import matplotlib.pyplot as plt
import matplotlib.mpl
import numpy as np
from numpy import size, exp, pi, nonzero, zeros, ones, sin, cos
import pylab as pl
sys.path.append("/home/stefan/1PhD/")
from setproc.common.classes.to_save_object import ToSaveObject
from setproc.cycle_process.classes.cycle_process import CycleProcess
from setproc.sweep_set.functions.filter import filter
from setproc.common.functions.local_extrema import local_extrema
from scipy.optimize.minpack import curve_fit
from scipy.interpolate import interp1d
import os
import json
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import qutip as qt
import subprocess
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
            
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx        
        
class Lifetime(CycleProcess):
    def __init__(self,delta_t=2.5):    
        
        self.trace =[]
        self.retrace =[]
        self.delta_t = delta_t
        self.raw = []
        self.raw2 = []
        self.data = []
        self.miss = []
        self.dmi1 = []
        self.dmi2 = []
        self.dmi3 = []
        self.dmi4 = []
        self.cnt = 0;
        self.trm = []
        
    def GetOffset(self,istart_tr,istart_rt,w,pw,sw,span,thres,off_tr,off_rt,strech=1): 
        sw_nbr = self.trace['sweep_number']
        up_tr,down_tr = self.trace.GetPeaks(istart_tr,w,pw,sw,span,thres,
                                            off_tr,[0,sw_nbr],'max')
        up_rt,down_rt = self.retrace.GetPeaks(istart_rt,w,pw,sw,span,thres,
                                              off_rt,[0,sw_nbr],'max')
        up_tr[:,1]   *= strech
        down_tr[:,1] *= strech
        up_rt[:,1]   *= strech
        down_rt[:,1] *= strech
        f = plt.figure()
        ax1 = f.add_subplot(211);   ax2 = f.add_subplot(212)
        ax1.set_title('trace');     ax2.set_title('retrace')
        ax1.set_xlabel('Bp (T)');   ax2.set_xlabel('Bp (T)');
        ax1.hist(up_tr[:,1],200);   ax2.hist(up_rt[:,1],200); 
        
            
    def GetJumps(self,istart_tr,istart_rt,w,pw,sw,span,thres,off_tr,off_rt,strech=1): 
        sw_nbr = self.trace['sweep_number']
        dt = self.delta_t
        up_tr,down_tr = self.trace.GetPeaks(istart_tr,w,pw,sw,span,thres,
                                            off_tr,[0,2*sw_nbr*dt-2*dt],'max')
        up_rt,down_rt = self.retrace.GetPeaks(istart_rt,w,pw,sw,span,thres,
                                              off_rt,[dt,2*sw_nbr*dt-dt],'max')
        up_tr[:,1]   *= strech
        down_tr[:,1] *= strech
        up_rt[:,1]   *= strech
        down_rt[:,1] *= strech
        
        jmp = np.row_stack((up_tr,down_tr,up_rt,down_rt))
        ind = np.lexsort((jmp[:,1],jmp[:,0])) 
        jmp = jmp[ind]
        temp = []
        temp2 = []
        for ii in range(size(jmp[:,0])):
            if ((jmp[ii,1] > -0.047) and (jmp[ii,1]<-0.033)):
                temp.append([jmp[ii,0],  1.5])
                temp2.append([jmp[ii,0], jmp[ii,1]])
            elif ((jmp[ii,1] > -0.020) and (jmp[ii,1]<-0.006)):
                temp.append([jmp[ii,0],  0.5])
                temp2.append([jmp[ii,0], jmp[ii,1]])
            elif ((jmp[ii,1] > 0.006) and (jmp[ii,1]<0.020)):
                temp.append([jmp[ii,0], -0.5])
                temp2.append([jmp[ii,0], jmp[ii,1]])
            elif ((jmp[ii,1] > 0.033) and (jmp[ii,1]<0.047)):
                temp.append([jmp[ii,0], -1.5])
                temp2.append([jmp[ii,0], jmp[ii,1]])
            elif ((jmp[ii,1] > -0.074) and (jmp[ii,1]<-0.064)):
                temp.append([jmp[ii,0], -1.5])
                temp2.append([jmp[ii,0], jmp[ii,1]])
            
        
        self.raw = jmp
        self.raw2 = np.array(temp2)
        self.data = np.array(temp)
        return[up_tr,down_tr,up_rt,down_rt]
    
    def GetJumpsSat(self,istart_tr,istart_rt,w,pw,sw,span,thres,off_tr,off_rt,strech=1): 
        sw_nbr = self.trace['sweep_number']
        dt = self.delta_t
        up_tr,down_tr = self.trace.GetPeaks(istart_tr,w,pw,sw,span,thres,
                                            off_tr,[0,2*sw_nbr*dt-2*dt],'max')
        up_rt,down_rt = self.retrace.GetPeaks(istart_rt,w,pw,sw,span,thres,
                                              off_rt,[dt,2*sw_nbr*dt-dt],'max')
        up_tr[:,1]   *= strech
        down_tr[:,1] *= strech
        up_rt[:,1]   *= strech
        down_rt[:,1] *= strech
        
        jmp = np.row_stack((up_tr,down_tr,up_rt,down_rt))
        ind = np.lexsort((jmp[:,1],jmp[:,0])) 
        jmp = jmp[ind]
        temp = []
        for ii in range(size(jmp[:,0])):
            if ((jmp[ii,1] > -0.047) and (jmp[ii,1]<-0.033)):
                temp.append([jmp[ii,0],  1.5])
            elif ((jmp[ii,1] > -0.020) and (jmp[ii,1]<-0.006)):
                temp.append([jmp[ii,0],  0.5])
            elif ((jmp[ii,1] > 0.006) and (jmp[ii,1]<0.020)):
                temp.append([jmp[ii,0], -0.5])
            elif ((jmp[ii,1] > 0.033) and (jmp[ii,1]<0.047)):
                temp.append([jmp[ii,0], -1.5])
            elif ((jmp[ii,1] > 0.070) and (jmp[ii,1]<0.084)):
                temp.append([jmp[ii,0],  -1.5])
            elif ((jmp[ii,1] > -0.081) and (jmp[ii,1]<-0.065)):
                temp.append([jmp[ii,0],  1.5])
        
        self.raw = jmp
        self.data = np.array(temp)
    
    def JumpFreq(self):
        #only one peak per cycle        
        dt1 = []
        dt2 = []
        for ii in range(np.size(self.data[:,0])-1):
            if self.data[ii,1] == self.data[ii+1,1]:
                if self.data[ii,1] == 0.5:
                    dt1.append(self.data[ii+1,0]-self.data[ii,0])
                elif self.data[ii,1] == -0.5:
                    dt2.append(self.data[ii+1,0]-self.data[ii,0])
                
        dt1 = np.array(dt1)
        dt2 = np.array(dt2)
        
        fig1 = pl.figure()
        ax1 = fig1.add_subplot(211)
        ax1.set_xlim((0.2,2))
        ax1.set_xlabel(r'$\Delta t \ (\mathsf{s}))$')
        ax1.set_ylabel('occurence')
        ax1.set_title('mi = 0.5')
        ax1.hist(np.log10(dt1),200,log=False)      
        
        ax2 = fig1.add_subplot(212)
        ax2.set_xlim((0.2,2))
        ax2.set_xlabel(r'$\Delta t \ (\mathsf{s}))$')
        ax2.set_ylabel('occurence')
        ax2.set_title(r"m_{\mathsf{I}} = -1/2$")
        ax2.hist(np.log10(dt2),200,log=False)  
        fig1.tight_layout()        
            
    def Lifet(self,dt_max,t_max,save = False):
        delta_t = self.delta_t
        temp1 = []
        temp2 = []        
        tau = 0
        for ii in range(size(self.data[:,0])-1):           
            if (self.data[ii,1] == self.data[ii+1,1]):
                if(self.data[ii+1,0]-self.data[ii,0]<dt_max):
                    tau += self.data[ii+1,0]-self.data[ii,0]
                else:
                    if (self.data[ii,1] == +0.5): temp1.append(tau)
                    if (self.data[ii,1] == -0.5): temp2.append(tau)
                    tau = 0
            else:
                if (self.data[ii,1] == +0.5): temp1.append(tau)
                if (self.data[ii,1] == -0.5): temp2.append(tau)
                tau = 0
        
        temp = [temp1,temp2]
                            
        exp_fit = lambda t,tau,a : a*exp(-t/tau)    
        time = np.linspace(0,120,100)
        
        f1 = plt.figure()
        
        #for i in range(2):
        ax = f1.add_subplot(111)
        ax.set_yscale('log')
        ax.set_ylim(0.005,1)
        ax.set_yticks((0.01,0.1,1))
        ax.set_xlabel(r'$t \ (\mathsf{s})$')
        ax.set_ylabel(r"$\langle m_{\mathsf{I}} \ = \ -1/2 \rangle $");
        ax.set_xlim((0,t_max));
        ax.set_xticks((0,20,40,60,80,100,120))
        #create Histogramm of lifetime distribution
        H1 = np.histogram(np.array(temp[0]),
                          bins=np.linspace(delta_t,120,int(120/delta_t)))
    
        #extract all nonzero elements  #
        lft = (np.reshape(np.concatenate((H1[1][nonzero(H1[0]!=0)],
                                          H1[0][nonzero(H1[0]!=0)])),
                         (size(H1[1][nonzero(H1[0]!=0)]),2),order='F'))
                       
        params, cov = curve_fit(exp_fit,lft[:,0],lft[:,1],[10,100])
        
        fit = exp(-time/params[0])
    
        #normalize                     
        lft[:,1] = lft[:,1]/params[1]
        ax.scatter(lft[:,0],lft[:,1],c='k')  
        ax.plot(time,fit,'r--',linewidth=3)         
        ax.text(0.6,0.8,r'$\tau$'+' = '+str(round(params[0]*100)/100)+'s',
                                         fontsize = 'x-large',
                                         transform = ax.transAxes)      
                                             
        f1.tight_layout()
        
        if save == True:
            f1.savefig('lifetime.png',dpi=300,format='png')
            f1.savefig('lifetime.pdf',dpi=300,format='pdf')
        

    def MissingLZ(self,his,lim=[0,100],legloc=3, save = False):
        delta_t = self.delta_t
        #create missing events        
        miss = np.zeros(4)
        miss_d = []
        for ii in range(size(self.data[:,0])-1):  
            t_diff = self.data[ii+1,0]-self.data[ii,0]
            if( t_diff > delta_t):
                n = int(t_diff/delta_t)-1
                idx1 = int(1.5-self.data[ii,1])
                idx2 = int(1.5-self.data[ii+1,1])                 
                miss[idx1] += n/2.0
                miss[idx2] += n/2.0
                for jj in range(n):
                    miss_d.append([self.data[ii,0]+(jj+1)*self.delta_t,self.data[ii,1]])
                    miss_d.append([self.data[ii,0]+(jj+1)*self.delta_t,self.data[ii+1,1]])
        
        self.miss = np.array(miss_d)        
        f = plt.figure(figsize=(8,4))
        ax = f.add_axes([0.05,0.2,0.42,0.78])
        ax2 = f.add_axes([0.60,0.2,0.37,0.78])
        ax.set_yticklabels((""))
        ax2.set_ylim((-2,2))
        ax2.set_yticks((1.5,0.5,-0.5,-1.5))
        ax2.set_xlabel("occurence "+r"$\times 10^3$")
        ax.set_xlabel(r"$t \ (\mathsf{s}) \ \times 10^3$")
        ax.set_yticks((1.5,0.5,-0.5,-1.5))
        ax.set_ylim((-2,2))
        ax2.set_yticklabels([r"$|+\frac{3}{2}\rangle$",
                            r"$|+\frac{1}{2}\rangle$",
                            r"$|-\frac{1}{2}\rangle$",
                            r"$|-\frac{3}{2}\rangle$"])
        ax.set_xlim(lim)
        ax.scatter(self.data[:,0]/1000.,self.data[:,1],s=36,
                   marker='o',color='r',label="QTM transition")
        ax.scatter(self.miss[:,0]/1000.,self.miss[:,1],s=36,
                   marker='o',color='white',
                   edgecolor="black",label="no transition")
        ax.plot(self.data[:,0]/1000.,self.data[:,1],c="k",linewidth=2)
        ax.legend(loc=legloc,prop={'size':14})
        dy = 0.15
        ax2.barh([1.5+dy,0.5+dy,-0.5+dy,-1.5+dy],his/1000.,color=("blue","green","red","black"),
                 height=0.25,align="center",alpha=0.5)
        ax2.barh([1.5-dy,0.5-dy,-0.5-dy,-1.5-dy],miss/1000.,color="white",
                 edgecolor="black",height=0.25,align="center")
        
        x = max(his)*1.1/1000.         
        ax2.set_xlim((0,x*1.4))
        
        for i in range(4):
            ax2.text(x,1.5-i,str(round(100.*his[i]/
            (his[i]+miss[i])))+"%",fontsize="x-large")            
            
        if save == True:
            f.savefig('missing_LZ.png',dpi=300)
            f.savefig('missing_LZ.pdf')
            f.savefig('missing_LZ.eps')
        return miss


if __name__ == '__main__':
    os.chdir("/home/stefan/1PhD/EMig/Diluette/14R/Lifetime_Bt/Lifetime_dB80mT_dt2.5s_Bt=800mT/")    
    lt = Lifetime(2.5)
    temp = CycleProcess('Stat.bin') 
    temp.LoadSweeps()
    lt.trace = temp.trace
    lt.retrace = temp.retrace    
    lt.GetJumps(0,0,4,1,17,10,10**-30.85,-0.0065,0.021)
    #lt.GetHist_Gauss()
    #lt.Fidelity()
    #lt.DeltaMI()
    #lt.RelEx()    
    #lt.Lifet(True)
    #lt.Trans_Matrix(True)
    lt.Movie()