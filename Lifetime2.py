# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 19:07:03 2012

@author: stefan
"""
import sys
import matplotlib.pyplot as plt
import matplotlib.mpl
import numpy as np
from numpy import size, exp, pi, nonzero, zeros, ones
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
    def __init__(self,trace=None,retrace=None,interval=None,mode=None, dt=2.5):    
        #temp = CycleProcess("0G_B_trace.json","0G_B_retrace.json",range(10),"Json")
        #temp.GetStat(0,4,1,17)
        #temp.SaveAll("trace","retrace","Stat.bin")
        
        self.trace =[]
        self.retrace =[]
        self.dt = dt
        self.raw = []
        self.data = []
        self.dmi1 = []
        self.dmi2 = []
        self.dmi3 = []
        self.dmi4 = []
        self.cnt = 0;
        self.trm = []
    
    def GetJumps(self,istart_tr,istart_rt,w,pw,sw,span,thres,off_tr,off_rt,dt): 
        sw_nbr = self.trace['sweep_number']
        speed = (max(self.trace['bias'])-min(self.trace['bias']))/self.dt
        #up_tr[Bp,time,ind]
        up_tr,down_tr = self.trace.GetPeaks(istart_tr,w,pw,sw,span,thres,
                                            off_tr,[0,2*sw_nbr*dt-2*dt],'max')
        
        up_tr[:,0]=up_tr[:,0]+abs(up_tr[:,1]-self.trace['bias'][0])/speed
        down_tr[:,0]=down_tr[:,0]+abs(down_tr[:,1]-self.trace['bias'][0])/speed
        
        up_rt,down_rt = self.retrace.GetPeaks(istart_rt,w,pw,sw,span,thres,
                                              off_rt,[dt,2*sw_nbr*dt-dt],'max')
        
        up_rt[:,0]=up_rt[:,0]+abs(up_rt[:,1]-self.retrace['bias'][0])/speed
        down_rt[:,0]=down_rt[:,0]+abs(down_rt[:,1]-self.retrace['bias'][0])/speed
                                      
        jmp = np.row_stack((up_tr,down_tr,up_rt,down_rt))
        ind = np.lexsort((jmp[:,1],jmp[:,0])) 
        jmp = jmp[ind]
        temp = []
        for ii in range(size(jmp[:,0])):
            if ((jmp[ii,1] > -0.047) and (jmp[ii,1]<-0.033)):
                temp.append([jmp[ii,0],-1.5])
            elif ((jmp[ii,1] > -0.020) and (jmp[ii,1]<-0.006)):
                temp.append([jmp[ii,0],-0.5])
            elif ((jmp[ii,1] > 0.006) and (jmp[ii,1]<0.020)):
                temp.append([jmp[ii,0], 0.5])
            elif ((jmp[ii,1] > 0.033) and (jmp[ii,1]<0.047)):
                temp.append([jmp[ii,0], 1.5])
        
        self.raw = jmp
        self.data = np.array(temp)
        
    def GetHist_Gauss(self):
        
        fit_gauss = lambda x, a, x0, sigma: a*exp(-np.square(x-x0)/
                                                  (2*np.square(sigma)))    
        
        f1 = plt.figure(figsize=(8,12))
        ax = f1.add_subplot(211); ax3 = f1.add_subplot(212)
        ax2 = ax.twinx(); ax4 = ax3.twinx();
        ax2.set_yticklabels(("")); ax4.set_yticklabels((""))
        #ax.set_ylim((0,500))
        ax.set_title('trace'); ax3.set_title('retrace')
        ax.set_xlabel(r'$\mu_0H_{||} \ (\mathsf{T})$')
        ax3.set_xlabel(r'$\mu_0H_{||} \ (\mathsf{T})$')
        ax.set_ylabel('occurence'); ax3.set_ylabel('occurence')
        ax.set_xlim((-0.06,0.06)); ax3.set_xlim((-0.06,0.06))
        #trace
        hist1 = np.histogram(self.raw[nonzero(self.raw[:,0]%(2*self.dt)<self.dt)][:,1],
                                      200,range=[-0.06,0.06]) 
        #retrace
        hist2 = np.histogram(self.raw[nonzero(self.raw[:,0]%(2*self.dt)>=self.dt)][:,1],
                                      200,range=[-0.06,0.06]) 
        
        cl = ['grey','blue','green','red']
        c = [-0.04,-0.013,0.013,0.04]
        
        for i in range(4):    
            params, cov = curve_fit(fit_gauss,hist1[1][i*50:(i+1)*50],
                                    hist1[0][i*50:(i+1)*50],[1000,c[i],0.005])
            fit = fit_gauss(hist1[1][i*50:(i+1)*50],params[0],params[1],params[2])
            ax2.bar(hist1[1][i*50:(i+1)*50]-3e-4,hist1[0][i*50:(i+1)*50],
                    width=6e-4,color=cl[i],edgecolor='None')
            ax2.plot(hist1[1][i*50:(i+1)*50],fit,'k--')
            ax.bar(params[1],np.sqrt(2*pi)*params[0]*params[2],align='center',
                alpha=0.3,width=6*params[2],color='grey',edgecolor='black',
                linewidth=2)            
        
        for i in range(4):    
            params, cov = curve_fit(fit_gauss,hist2[1][i*50:(i+1)*50],
                                    hist2[0][i*50:(i+1)*50],[1000,c[i],0.005])
            fit = fit_gauss(hist2[1][i*50:(i+1)*50],params[0],params[1],params[2])
            ax4.bar(hist2[1][i*50:(i+1)*50]-3e-4,hist2[0][i*50:(i+1)*50],
                    width=6e-4,color=cl[i],edgecolor='None')
            ax4.plot(hist2[1][i*50:(i+1)*50],fit,'k--')
            ax3.bar(params[1],np.sqrt(2*pi)*params[0]*params[2],align='center',
                alpha=0.3,width=6*params[2],color='grey',edgecolor='black',
                linewidth=2)  
            
        f1.tight_layout() 
        f1.savefig('Fidelity.png',dpi=300,format='png')
        f1.savefig('Fidelity.pdf',dpi=300,format='pdf')
    
    def Movie(self):
        
        def add_circle(theta,res):
            phi = np.linspace(0,2*pi,res)
            x = sin(theta)*cos(phi)
            y = sin(theta)*sin(phi)
            z = cos(theta)*np.ones(res)
            return[x,y,z]

        os.chdir('/home/stefan/1PhD/EMig/Diluette/14R/Lifetime_dB80mT_4/movie3/')
        fig = plt.figure(figsize=(16,9))
        rec1 = [0.1,0.15,0.4,0.7]
        rec2 = [0.60,0.125,0.4,0.725]   
        theta = np.array([39.3,75,104.9,140.6])*pi/180   
        colors = np.array(['k','b','g','r'])             
        res = [50,75,75,50]                     
        for jj in range(20000):
            
            ax1 = fig.add_axes(rec1)
            ax1.set_xlabel(r'$t \ (\mathsf{s})$')
            ax1.set_ylabel(r'$B_{||} \ (\mathsf{T})$')
            ax3 = ax1.twinx()
            ax3.set_ylabel(r'$m_{\mathsf{I}}$')
            ax3.set_yticks((-1.5,-0.5,0.5,1.5))
            ax3.set_yticklabels((r'$|+\frac{3}{2}>$',r'$|+\frac{1}{2}>$',
                                         r'$|-\frac{1}{2}>$',r'$|-\frac{3}{2}>$'))

            #ax1.plot([0.5,0.5],[0,1],'k.-.',transform=ax1.transAxes)
            ax1.set_ylim((0.06,-0.06))
            ax3.set_ylim((2.25,-2.25))    
            ax1.scatter(self.raw[:,0],self.raw[:,1],s=30,c='grey')
            ax3.plot(self.data[:,0],self.data[:,1],color='r',linewidth=4)
            ax1.bar(0.5,1,0.05,0.0,align='center',color='grey',edgecolor='black',alpha=0.5,transform=ax1.transAxes)
                      
            rge = 100            
            ii = jj;
            ax1.set_xlim((ii,ii+2*rge))        
            #ax1.bar(ii+rge,0.12,5,-0.06,align='center',color='grey',edgecolor='black',alpha=0.5)
            idx = int(round(self.data[find_nearest(self.data[:,0],ii+rge),1]+1.5))      
            
            
            ax2 = p3.Axes3D(fig,rec2)
            B = qt.Bloch(fig,ax2)
            B.font_size = 30
            B.zlabel = [r'$z$','']
            B.xlpos = [1.7,-1.7]
            B.ylpos = [1.55,-1.55]
            B.zlpos = [1.35,-1.35]
            B.point_color =  colors[idx]
            #B.vector_color = colors[idx] 
            #B.add_vectors([0,sin(theta[idx]),cos(theta[idx])])
            B.add_points(add_circle(theta[idx],res[idx]))
            B.show()
            
            #draw a vector
            vz = Arrow3D([0,0],[0,0],[0,1.3], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
            vx = Arrow3D([0,0],[0,-1.5],[0,0], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
            vy = Arrow3D([0,1.5],[0,0],[0,0], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
            vec = Arrow3D([0,1.05*sin(theta[idx])],[0,0],[0,1.05*cos(theta[idx])], mutation_scale=20, lw=6, arrowstyle="-|>", color="k") 
            ax2.add_artist(vz)
            ax2.add_artist(vy)
            ax2.add_artist(vx)
            ax2.add_artist(vec)
            ax2.axis("off")
            ax2.view_init(10,-45)
            ax1.ticklabel_format(style='plain')
            #fig.show()
            fig.savefig(str('%05d' % jj)+'.png',dpi=100,format='png')
            print('create figure '+str(jj))
            fig.clear()
            
        command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=1600:h=900:fps=25',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4:vbitrate=3000',
           '-oac',
           'copy',
           '-o',
           'output.avi')

        #os.spawnvp(os.P_WAIT, 'mencoder', command)
    
        print "\n\nabout to execute:\n%s\n\n" % ' '.join(command)
        subprocess.check_call(command)
    
        print "\n\n The movie was written to 'output.avi'"
    
        print "\n\n You may want to delete *.png now.\n\n"
            
        
    def DeltaMI(self):
        temp = []           
        for ii in range(size(self.data[:,0])-1):
            temp.append([self.data[ii+1,1]-self.data[ii,1]])
        f1 = plt.figure()
        ax1 = f1.add_subplot(111)
        self.dmi = np.array(temp)
        ax1.hist(self.dmi,
                 bins=(-3,-2,-1,0,1.0,2,3),
                 align='left',rwidth=0.5,normed=True,color='r')
        ax1.set_xlabel(r'$\Delta m_{\mathsf{I}}$')
        ax1.set_ylabel('probability')
        ax1.set_xlim((-3,3))
        f1.savefig('delta_m.png',dpi=300,format='png')
        f1.savefig('delta_m.pdf',dpi=300,format='pdf')
        
    def RelEx(self):
        temp1 = []
        temp2 = []
        temp3 = []
        temp4 = []
        
        for ii in range(size(self.data[:,0])-1):
            if ((self.data[ii,0] % 5) == 0): #trace
                if((self.data[ii+1,0]-self.data[ii,0]) == 2.5):
                    dm = self.data[ii+1,1]-self.data[ii,1]                    
                    if (self.data[ii,1] == 1.5):
                        if (dm == 0):
                            temp1.append([0,dm])#nothing happend
                        else:
                            temp1.append([1,dm])#excitation between tr and rt
   
                    if (self.data[ii,1] == 0.5):
                        if dm == 0:
                            temp2.append([0,dm])#nothing happend
                        elif dm < 0:
                            temp2.append([1,dm])#excitation between tr and rt
                        else:
                            temp2.append([-1,dm])#relaxation between tr and rt

                    if self.data[ii,1] == -0.5:
                        if dm == 0:
                            temp3.append([0,dm])#nothing happend
                        elif dm < 0:
                            temp3.append([1,dm])#excitation between tr and rt
                        else:
                            temp3.append([-1,dm])#relaxation between tr and rt

                    if self.data[ii,1] == -1.5:
                        if dm == 0:
                            temp4.append([0,dm])#nothing happend
                        else:
                            temp4.append([-1,dm])#relaxation between tr and rt
            
            else: #retrace
                if((self.data[ii+1,0]-self.data[ii,0]) == 2.5) :
                    dm = self.data[ii+1,1]-self.data[ii,1]                    
                    if self.data[ii,1] == 1.5:
                        if dm == 0:
                            temp1.append([0,dm])#nothing happend
                        else:
                            temp1.append([-1,dm])#relaxation between tr and rt
   
                    if self.data[ii,1] == 0.5:
                        if dm == 0:
                            temp2.append([0,dm])#nothing happend
                        elif dm < 0:
                            temp2.append([-1,dm])#relaxation between tr and rt
                        else:
                            temp2.append([+1,dm])#excitation between tr and rt

                    if self.data[ii,1] == -0.5:
                        if dm == 0:
                            temp3.append([0,dm])#nothing happend
                        elif dm < 0:
                            temp3.append([-1,dm])#relaxation between tr and rt
                        else:
                            temp3.append([+1,dm])#excitation between tr and rt

                    if self.data[ii,1] == -1.5:
                        if dm == 0:
                            temp4.append([0,dm])#nothing happend
                        else:
                            temp4.append([1,dm])#excitation between tr and rt
        
        f1 = plt.figure(figsize=(12,8))
        ax1 = f1.add_subplot(221)
        ax2 = f1.add_subplot(222)
        ax3 = f1.add_subplot(223)
        ax4 = f1.add_subplot(224)
        ax1.set_ylim((0,1));ax2.set_ylim((0,1));ax3.set_ylim((0,1));ax4.set_ylim((0,1));
        self.dmi1 = np.array(temp1)
        self.dmi2 = np.array(temp2)
        self.dmi3 = np.array(temp3)
        self.dmi4 = np.array(temp4)
        
        rel_m1 = size(nonzero(self.dmi1[nonzero(self.dmi1[:,0]==-1),1][0]== -1))
        exc_m1 = size(nonzero(self.dmi1[nonzero(self.dmi1[:,0]== 1),1][0]== -1))
        stay   = size(self.dmi1[nonzero(self.dmi1[:,0]== 0),1][0])
        norm = size(self.dmi1[:,0])*1.0        
        rel_m1 = rel_m1*1.0/norm; exc_m1 = exc_m1*1.0/norm; stay = stay*1.0/norm;
        ax1.bar(-1.3,rel_m1,width=0.3,color='red',label='relaxation')
        ax1.bar(-1.0,exc_m1,width=0.3,color='blue',label='excitation')
        ax1.bar(-0.25,stay,width=0.5,color='black')
        ax1.set_xlabel(r'$\Delta m_{\mathsf{I}}$')
        ax1.set_ylabel('probability')
        ax1.set_title(r'$m_{\mathsf{I}} \  = \  \frac{3}{2}$')
        ax1.set_xlim((-2,2))
        ax1.set_xticks([-2,-1,0,1,2])
        ax1.plot([-2,2],[0.85,0.85],'k.-.',label='theory')
        ax1.legend(loc='upper right',prop={'size':14})
        
        rel_m1 = size(nonzero(self.dmi2[nonzero(self.dmi2[:,0]==-1),1][0]== -1))
        rel_p1 = size(nonzero(self.dmi2[nonzero(self.dmi2[:,0]==-1),1][0]== +1))
        exc_m1 = size(nonzero(self.dmi2[nonzero(self.dmi2[:,0]==+1),1][0]== -1))
        exc_p1 = size(nonzero(self.dmi2[nonzero(self.dmi2[:,0]==+1),1][0]== +1))
        stay   = size(self.dmi2[nonzero(self.dmi2[:,0]== 0),1][0])
        norm = size(self.dmi2[:,0])*1.0
        rel_m1 = rel_m1/norm; rel_p1 = rel_p1/norm
        exc_m1 = exc_m1/norm; exc_p1 = exc_p1/norm
        stay = stay/norm
        ax2.bar(-1.3,rel_m1,width=0.3,color='red',label='relaxation')
        ax2.bar(-1.0,exc_m1,width=0.3,color='blue',label='excitation')
        ax2.bar( 0.7,rel_p1,width=0.3,color='red',label='relaxation')
        ax2.bar( 1.0,exc_p1,width=0.3,color='blue',label='excitation')
        ax2.bar(-0.25,stay,width=0.5,color='black')
        ax2.set_xlabel(r'$\Delta m_{\mathsf{I}}$')
        ax2.set_ylabel('probability')
        ax2.set_title(r'$m_{\mathsf{I}} \ = \ \frac{1}{2}$')
        ax2.set_xlim((-2,2))
        ax2.set_xticks([-2,-1,0,1,2])
        ax2.plot([-2,2],[0.67,0.67],'k.-.',label='theory')
        
        rel_m1 = size(nonzero(self.dmi3[nonzero(self.dmi3[:,0]==-1),1][0]== -1))
        rel_p1 = size(nonzero(self.dmi3[nonzero(self.dmi3[:,0]==-1),1][0]== +1))
        exc_m1 = size(nonzero(self.dmi3[nonzero(self.dmi3[:,0]==+1),1][0]== -1))
        exc_p1 = size(nonzero(self.dmi3[nonzero(self.dmi3[:,0]==+1),1][0]== +1))
        stay   = size(self.dmi3[nonzero(self.dmi3[:,0]== 0),1][0])
        norm = size(self.dmi3[:,0])*1.0
        rel_m1 = rel_m1/norm; rel_p1 = rel_p1/norm
        exc_m1 = exc_m1/norm; exc_p1 = exc_p1/norm
        stay = stay/norm        
        ax3.bar(-1.3,rel_m1,width=0.3,color='red',label='relaxation')
        ax3.bar(-1.0,exc_m1,width=0.3,color='blue',label='excitation')
        ax3.bar( 0.7,rel_p1,width=0.3,color='red',label='relaxation')
        ax3.bar( 1.0,exc_p1,width=0.3,color='blue',label='excitation')
        ax3.bar(-0.25,stay,width=0.5,color='black')
        ax3.set_xlabel(r'$\Delta m_{\mathsf{I}}$')
        ax3.set_ylabel('probability')
        ax3.set_title(r'$m_{\mathsf{I}} \ = \ -\frac{1}{2}$')
        ax3.set_xlim((-2,2))
        ax3.set_xticks([-2,-1,0,1,2])
        ax3.plot([-2,2],[0.701,0.701],'k.-.',label='theory')
        
        rel_p1 = size(nonzero(self.dmi4[nonzero(self.dmi4[:,0]==-1),1][0]== 1))
        exc_p1 = size(nonzero(self.dmi4[nonzero(self.dmi4[:,0]== 1),1][0]== 1))
        stay   = size(self.dmi4[nonzero(self.dmi4[:,0]== 0),1][0])
        norm = size(self.dmi4[:,0])*1.0    
        rel_p1 = rel_p1/norm; exc_p1 = exc_p1/norm; stay = stay/norm;
        ax4.bar( 0.7,rel_p1,width=0.3,color='red',label='relaxation')
        ax4.bar( 1.0,exc_p1,width=0.3,color='blue',label='excitation')
        ax4.bar(-0.25,stay,width=0.5,color='black')         
        ax4.set_xlabel(r'$\Delta m_{\mathsf{I}}$')
        ax4.set_ylabel('probability')
        ax4.set_title(r'$m_{\mathsf{I}} \ = \  -\frac{3}{2}$')
        ax4.set_xlim((-2,2))
        ax4.set_xticks([-2,-1,0,1,2])
        ax4.plot([-2,2],[0.833,0.833],'k.-.',label='theory')
        
        f1.tight_layout()
        f1.savefig('delta_m.png',dpi=300,format='png')
        f1.savefig('delta_m.pdf',dpi=300,format='pdf')
        
    def Trans_Matrix(self,dt=2.5):
        rel_tr = np.zeros(3)
        exc_tr = np.zeros(3)
        rel_rt = np.zeros(3)
        exc_rt = np.zeros(3)
        
        for ii in range(size(self.data[:,0])-1):
            
            delta_t = self.data[ii+1,0]-self.data[ii,0]           
            if delta_t == dt:            
                dm = self.data[ii+1,1]-self.data[ii,1] 
                if(abs(dm) == 1):
                    if ((self.data[ii,0] % (2*dt)) == 0): #trace
                        #dm=-1 relaxation dm=+1 excitation                    
                        if dm == -1:
                            rel_tr[int(self.data[ii,1]+0.5)]+=1
                        elif dm == 1:
                            exc_tr[int(self.data[ii,1]+1.5)]+=1
                    
                    else: #retrace
                        #dm=+1 relaxation dm=-1 excitation                    
                        if dm == 1:
                            rel_rt[int(0.5-self.data[ii,1])]+=1
                        elif dm == -1:
                            exc_rt[int(1.5-self.data[ii,1])]+=1
                        
        self.trm = [rel_tr,exc_tr,rel_rt,exc_rt]

        #Relaxation        
        f = plt.figure(figsize=(12,8))
        label = ['relaxation','excitation']
        xtks =[[r"$ - \frac{3}{2} \leftarrow - \frac{1}{2} $",
                r"$ - \frac{1}{2} \leftarrow   \frac{1}{2} $",
                r"$   \frac{1}{2} \leftarrow   \frac{3}{2} $"],
               [r"$ - \frac{3}{2} \rightarrow - \frac{1}{2} $",
                r"$ - \frac{1}{2} \rightarrow   \frac{1}{2} $",
                r"$   \frac{1}{2} \rightarrow   \frac{3}{2} $"],
               [r"$  \frac{3}{2} \leftarrow  \frac{1}{2} $",
                r"$  \frac{1}{2} \leftarrow  -\frac{1}{2} $",
                r"$  -\frac{1}{2} \leftarrow  -\frac{3}{2} $"],
               [r"$  \frac{3}{2} \rightarrow  \frac{1}{2} $",
                r"$  \frac{1}{2} \rightarrow  -\frac{1}{2} $",
                r"$  -\frac{1}{2} \rightarrow  -\frac{3}{2} $"]]
        cl = ['red','blue']
        title = ['trace','retrace']

        for i in range(4):
            ax = f.add_subplot(2,2,i+1)
            ax.set_title(title[int(i/2)])
            ax.set_xlabel(label[i%2])
            ax.set_ylabel("occurence")
            ax.set_xlim((0.5,3.5))
            ax.set_xticks((1,2,3))        
            ax.set_xticklabels(xtks[i]) 
            ax.bar([1,2,3],self.trm[i],width=0.3,
                align='center',alpha=1,color=cl[i%2])
        
        f.subplots_adjust(wspace=0.4, hspace=0.6)
        f.savefig('trans_matrix.png',dpi=300)
        f.savefig('trans_matrix.pdf')

    def Lifet(self):
       
        temp1 = []
        temp2 = []
        temp3 = []
        temp4 = []
        
        tau = 0
        for ii in range(size(self.data[:,0])-1):           
            if (self.data[ii,1] == self.data[ii+1,1]):
                tau += self.data[ii+1,0]-self.data[ii,0]
            else:
                if (self.data[ii,1] == -1.5): temp1.append(tau)
                if (self.data[ii,1] == -0.5): temp2.append(tau)
                if (self.data[ii,1] ==  0.5): temp3.append(tau)
                if (self.data[ii,1] ==  1.5): temp4.append(tau)                          
                tau = 0
        
        temp = [temp1,temp2,temp3,temp4]
                            
        exp_fit = lambda t,tau,a : a*exp(-t/tau)    
        time = np.linspace(0,120,100)
        
        f1 = plt.figure(figsize=(12,8))
        
        for i in range(4):
            ax = f1.add_subplot(2,2,i+1)
            ax.set_xlabel(r'$t \ (s)$')
            ax.set_ylabel(r'$<m_{\mathsf{I}} \ = \ $'+str(+1.5-i)+r"$>$");
            ax.set_xlim((0,120));
            ax.set_xticks((0,30,60,90,120))
            ax.set_yticks((0,0.5,1))
            ax.set_ylim((-0.1,1.1)); 
            #create Histogramm of lifetime distribution
            H1 = np.histogram(np.array(temp[i]),bins=np.linspace(1,48,120))
        
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
            ax.text(0.3,0.8,r'$\tau$'+' = '+str(round(params[0]*100)/100)+'s',
                                             fontsize = 'x-large',
                                             transform = ax.transAxes)      
                                             
        f1.tight_layout()
        f1.savefig('lifetime.png',dpi=300,format='png')
        f1.savefig('lifetime.pdf',dpi=300,format='pdf')

if __name__ == '__main__':
    os.chdir("/home/stefan/1PhD/EMig/Diluette/14R/Lifetime_Bt/Lifetime_dB80mT_dt2.5s_Bt=800mT")    
    lt = Lifetime()
    temp = CycleProcess('Stat.bin') 
    temp.LoadSweeps()
    lt.trace = temp.trace
    lt.retrace = temp.retrace    
    lt.GetJumps(0,0,4,1,17,10,10**-31.5,-0.0065,0.021,2.5)
    lt.GetHist_Gauss()
    #lt.DeltaMI()
    #lt.RelEx()    
    lt.Lifet()
    #lt.Trans_Matrix(2.5)
    #lt.Movie()