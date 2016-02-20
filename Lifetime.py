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
        up_tr,down_tr = self.trace.GetPeaks(thres,off_tr,'max')
        up_rt,down_rt = self.retrace.GetPeaks(thres,off_rt,'max')
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
        up_tr,down_tr = self.trace.GetPeaks(thres,off_tr,'max',[0,2*sw_nbr*dt-2*dt])
        up_rt,down_rt = self.retrace.GetPeaks(thres,off_rt,'max',[dt,2*sw_nbr*dt-dt])
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
        up_tr,down_tr = self.trace.GetPeaks(thres,off_tr,'max',[0,2*sw_nbr*dt-2*dt])
        up_rt,down_rt = self.retrace.GetPeaks(thres,off_rt,'max',[dt,2*sw_nbr*dt-dt])
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
    
 
    def Fidelity(self,bins = 360,save = False):
        delta_t = self.delta_t
        fit_gauss = lambda x, a, x0, sigma: a*exp(-np.square(x-x0)/
                                                  (2*np.square(sigma)))    


        f1 = plt.figure()
        ax = f1.add_subplot(111); 
        ax2 = ax.twinx(); 
        ax2.set_ylabel("occurance "+r"$\times 10^3$")
        #ax.set_ylim((0,500))
        ax.set_xlabel(r'$\mu_0H_{||} \ (\mathsf{T})$')
        ax.set_ylabel("probability "); 
        ax.set_xlim((-0.0535,0.0535)); 
        #ax.set_ylim((0,0.5))
        ax3 = ax.twiny();
        ax3.set_xlim((-2,2))
        ax3.set_xticks((-1.5,-0.5,0.5,1.5))
        ax3.invert_xaxis()
        ax3.set_xticklabels([r"$|-\frac{3}{2} \rangle $",
                            r"$|-\frac{1}{2} \rangle$",
                            r"$|+\frac{1}{2} \rangle$",
                            r"$|+\frac{3}{2} \rangle $"])
        hist1 = np.histogram(self.raw[:,1],bins,range=[-0.0535,0.0535]) 
        
        cl = ['blue','green','red',"grey"]
        c = [-0.04,-0.013,0.013,0.04]
        occ = np.zeros(4)
        z = int(bins/4.0)
        
        for i in range(4):
            occ[i] += size(self.data[nonzero(self.data[:,1]==(1.5-i))][:,1])
        
        y_max = max(occ)/sum(occ)
        ax.set_ylim(0,1.2*y_max)
        ax.set_yticks((0,0.1,0.2,0.3))
        y2_max = max(hist1[0])/1000.0
        ax2.set_ylim(0,1.3*y2_max)
        for i in range(4):    
            params, cov = curve_fit(fit_gauss,hist1[1][i*z:(i+1)*z],
                                    hist1[0][i*z:(i+1)*z],[1000,c[i],0.005])
            fit = fit_gauss(hist1[1][i*z:(i+1)*z],
                                params[0],params[1],params[2])
            
            ax2.bar(hist1[1][i*z:(i+1)*z],hist1[0][i*z:(i+1)*z]/1000.0,
                    width=4*0.0535/bins,color=(0.5,0.5,0.5),edgecolor=(0.5,0.5,0.5))
            ax2.plot(hist1[1][i*z:(i+1)*z],fit/1000.0,color="black",linestyle='--')  
            ax.bar(params[1],occ[i]/sum(occ),align='center',alpha=0.5,
                   width=7*params[2],color=cl[i],edgecolor='black',linewidth=2)
             
            idx1 = find_nearest(hist1[1],params[1]-3.5*params[2])
            idx2 = find_nearest(hist1[1],params[1]+3.5*params[2])
            noise = sum(hist1[0][i*z:idx1])+sum(hist1[0][idx2:(i+1)*z])    
            fid1 = occ[i]/(noise+occ[i])
            ax.text(params[1],y_max*1.1,str(round(fid1*1000)/10)+'%',
                 horizontalalignment='center',fontsize='x-large',
                 backgroundcolor='white')
            print "tolerance: +-",3.5*params[2]
        
        f1.subplots_adjust(left=0.15, top=0.86,right=0.85)
        if save == True:        
            f1.savefig('Fidelity.png',dpi=300,)
            f1.savefig('Fidelity.pdf',)
        
        return occ
        
        
    def Movie(self):
        
        def add_circle(theta,res):
            phi = np.linspace(0,2*pi,res)
            x = sin(theta)*cos(phi)
            y = sin(theta)*sin(phi)
            z = cos(theta)*np.ones(res)
            return[x,y,z]

        fig = plt.figure(figsize=(10,5.625))
        rec1 = [0.13,0.2,0.37,0.7]
        rec2 = [0.60,0.125,0.4,0.725]   
        theta = np.array([39.3,75,104.9,140.6])*pi/180   
        colors = np.array(['b','g','r','k'])             
        res = [50,75,75,50]     
        shift = 0                
        for jj in range(1000):
            
            ax1 = fig.add_axes(rec1)
            ax1.set_xlabel(r'$t \ (\mathsf{s})$')
            ax1.set_ylabel(r'$B_{||} \ (\mathsf{T})$')
            ax1b = ax1.twinx()
            ax1b.set_ylabel(r'$m_{\mathsf{I}}$')
            ax1.set_ylim((0.06,-0.06))
            ax1b.set_ylim((0.06,-0.06))  
            ax1b.set_yticks((-0.04,-0.013,0.013,0.04))
            ax1b.set_yticklabels((r'$|+\frac{3}{2}>$',r'$|+\frac{1}{2}>$',
                                         r'$|-\frac{1}{2}>$',r'$|-\frac{3}{2}>$'))
                                         
            ax1.plot(self.data[:,0]-shift,self.data[:,1]/(-1.5)*0.04,color='r',linewidth=3)
            ax1b.scatter(self.raw2[:,0]-shift,self.raw2[:,1],s=120,c="white",edgecolor="k")
            ax1b.scatter(self.raw2[:,0]-shift,self.raw2[:,1],s=25,c=(0.3,0.3,0.3),edgecolor=(0.3,0.3,0.3))
            ax1b.bar(0.5,1,0.05,0.0,align='center',color='grey',edgecolor='black',alpha=0.5,transform=ax1b.transAxes)
                      
            rge = 100            
            ii = jj;
            ax1.set_xlim((ii,ii+2*rge))        
            #ax1.bar(ii+rge,0.12,5,-0.06,align='center',color='grey',edgecolor='black',alpha=0.5)
            idx = int(round(1.5-self.data[find_nearest(self.data[:,0]-shift,ii+rge),1]))      
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
            fig.savefig(str('%05d' % jj)+'.png',dpi=150)
            print('create figure '+str(jj))
            fig.clear()
            
        command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=576:h=325:fps=25',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4:vbitrate=4000',
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
        ax1.hist(self.dmi,bins=(-3,-2,-1,0,1.0,2,3),align='left',rwidth=0.5,normed=True,color='r')
        ax1.set_xlabel(r'$\Delta m_{\mathsf{I}}$')
        ax1.set_ylabel('probability')
        ax1.set_xlim((-3,3))
        f1.savefig('delta_m.png',dpi=300,format='png')
        f1.savefig('delta_m.pdf',dpi=300,format='pdf')
        
       
    def Trans_Matrix(self,save = False):
        delta_t = self.delta_t
        rel_tr = np.zeros(3)
        exc_tr = np.zeros(3)
        rel_rt = np.zeros(3)
        exc_rt = np.zeros(3)
        stay_tr = np.zeros(4)
        stay_rt = np.zeros(4)
        
        
        for ii in range(size(self.data[:,0])-1):
            
            dt = self.data[ii+1,0]-self.data[ii,0]           
            if dt == delta_t:            
                dm = self.data[ii+1,1]-self.data[ii,1] 
                if(abs(dm) == 1):
                    if ((self.data[ii,0] % (2*delta_t)) == 0): #trace
                        #dm= +1 relaxation dm=-1 excitation                    
                        if dm == 1:
                            rel_tr[int(0.5-self.data[ii,1])]+=1
                        else:
                            exc_tr[int(1.5-self.data[ii,1])]+=1
                    
                    else: #retrace
                        #dm=-1 relaxation dm=+1 excitation                    
                        if dm == -1:
                            rel_rt[int(0.5+self.data[ii,1])]+=1
                        else:
                            exc_rt[int(1.5+self.data[ii,1])]+=1
                elif (dm == 0): #stay in the same state
                    mj = self.data[ii,1]
                    if ((self.data[ii,0] % (2*delta_t)) == 0): #trace                        
                        stay_tr[int(1.5-mj)] += 1
                    else: #retrace
                        stay_rt[int(1.5+mj)] += 1                         
        
        
        rel_tr /= (stay_tr[1:] +rel_tr+ np.array([exc_tr[1],exc_tr[2],0]))
        exc_tr /= (stay_tr[:-1]+exc_tr+ np.array([0,rel_tr[0],rel_tr[1]]))
        rel_rt /= (stay_rt[1:] +rel_tr+ np.array([exc_rt[1],exc_rt[2],0]))
        exc_rt /= (stay_rt[:-1]+exc_rt+ np.array([0,rel_rt[0],rel_rt[1]]))
                       
        self.trm = [rel_tr,exc_tr,rel_rt,exc_rt]

        #Relaxation        
        f = plt.figure(figsize=(12,8))
        label = ['relaxation','excitation']
        xtks =[[r"$  \frac{3}{2} \leftarrow   \frac{1}{2} $",
                r"$  \frac{1}{2} \leftarrow  -\frac{1}{2} $",
                r"$ -\frac{1}{2} \leftarrow  -\frac{3}{2} $"],
               [r"$  \frac{3}{2} \rightarrow  \frac{1}{2} $",
                r"$  \frac{1}{2} \rightarrow -\frac{1}{2} $",
                r"$ -\frac{1}{2} \rightarrow -\frac{3}{2} $"],
               [r"$ -\frac{3}{2} \leftarrow  -\frac{1}{2} $",
                r"$ -\frac{1}{2} \leftarrow   \frac{1}{2} $",
                r"$  \frac{1}{2} \leftarrow   \frac{3}{2} $"],
               [r"$ -\frac{3}{2} \rightarrow -\frac{1}{2} $",
                r"$ -\frac{1}{2} \rightarrow  \frac{1}{2} $",
                r"$  \frac{1}{2} \rightarrow  \frac{3}{2} $"]]
        cl = ['red','blue']
        title = ['trace','retrace']

        for i in range(4):
            ax = f.add_subplot(2,2,i+1)
            ax.set_title(title[int(i/2)])
            ax.set_xlabel(label[i%2])
            ax.set_ylabel("probability")
            ax.set_xlim((0.5,3.5))
            ax.set_xticks((1,2,3))        
            ax.set_xticklabels(xtks[i]) 
            ax.bar([1,2,3],self.trm[i],width=0.3,
                align='center',alpha=1,color=cl[i%2])
        
        f.subplots_adjust(wspace=0.4, hspace=0.6)
        
        if save == True:
            f.savefig('trans_matrix.png',dpi=300)
            f.savefig('trans_matrix.pdf')
    
    def Trans_Matrix2(self,save = False):
        delta_t = self.delta_t
        rel = np.zeros(3)
        exc = np.zeros(3)
        sty = np.zeros(4)
    
        for ii in range(size(self.data[:,0])-1):
            
            dt = self.data[ii+1,0]-self.data[ii,0]           
            if dt == delta_t:            
                dm = self.data[ii+1,1]-self.data[ii,1] 
                if(abs(dm) == 1):
                    if ((self.data[ii,0] % (2*delta_t)) == 0): #trace
                        #dm= +1 relaxation dm=-1 excitation                    
                        if dm == 1:
                            rel[int(0.5-self.data[ii,1])]+=1
                        else:
                            exc[int(1.5-self.data[ii,1])]+=1
                    
                    else: #retrace
                        #dm=-1 relaxation dm=+1 excitation                    
                        if dm == -1:
                            rel[int(0.5+self.data[ii,1])]+=1
                        else:
                            exc[int(1.5+self.data[ii,1])]+=1
                elif (dm == 0): #stay in the same state
                    mj = self.data[ii,1]
                    if ((self.data[ii,0] % (2*delta_t)) == 0): #trace                        
                        sty[int(1.5-mj)] += 1
                    else: #retrace
                        sty[int(1.5+mj)] += 1                         
        
        
        rel /= (sty[1:] +rel+ np.array([exc[1],exc[2],0]))
        exc /= (sty[:-1]+exc+ np.array([0,rel[0],rel[1]]))
        self.trm = [rel,exc]

        #Relaxation        
        f = plt.figure(figsize=(7,11))
        xtks =[[r"$ 0 \leftarrow   1 $",
                r"$ 1 \leftarrow   2 $",
                r"$ 2 \leftarrow   3 $"],
               [r"$ 0 \rightarrow  1 $",
                r"$ 1 \rightarrow  2 $",
                r"$ 2 \rightarrow  3 $"]]
                
        cl = ['red','blue']
        title = ['relaxation','excitation']

        for i in range(2):
            ax = f.add_subplot(2,1,i+1)
            ax.set_title(title[i],fontsize='xx-large')
            ax.set_xlabel("transitions")
            ax.set_ylabel("probability")
            ax.set_xlim((0.5,3.5))
            ax.set_xticks((1,2,3))        
            ax.set_xticklabels(xtks[i]) 
            ax.bar([1,2,3],self.trm[i],width=0.3,
                align='center',alpha=1,color=cl[i])
        
        f.subplots_adjust(hspace=0.6)
        
        if save == True:
            f.savefig('trans_matrix.png',dpi=300)
            f.savefig('trans_matrix.pdf')        
            
    def Lifet(self,xmax,save = False):
        delta_t = self.delta_t
        temp1 = []
        temp2 = []
        temp3 = []
        temp4 = []
        
        tau = 0
        for ii in range(size(self.data[:,0])-1):           
            if (self.data[ii,1] == self.data[ii+1,1]):
                tau += self.data[ii+1,0]-self.data[ii,0]
            else:
                if (self.data[ii,1] == +1.5): temp1.append(tau)
                if (self.data[ii,1] == +0.5): temp2.append(tau)
                if (self.data[ii,1] == -0.5): temp3.append(tau)
                if (self.data[ii,1] == -1.5): temp4.append(tau)                          
                tau = 0
        
        temp = [temp1,temp2,temp3,temp4]
                            
        exp_fit = lambda t,tau,a : a*exp(-t/tau)    
        time = np.linspace(0,120,100)
        
        f1 = plt.figure(figsize=(12,8))
        
        for i in range(4):
            ax = f1.add_subplot(2,2,i+1)
            ax.set_yscale('log')
            ax.set_ylim(0.005,1)
            ax.set_yticks((0.01,0.1,1))
            ax.set_xlabel(r'$t \ (\mathsf{s})$')
            ax.set_ylabel(r'$\langle m_{\mathsf{I}} \ = \ $'+str(1.5-i)+r"$\rangle $");
            ax.set_xlim((0,xmax));
            ax.set_xticks((0,20,40,60,80,100,120))
            #create Histogramm of lifetime distribution
            H1 = np.histogram(np.array(temp[i]),
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
        ax.plot(self.data[:,0]/1000.,self.data[:,1],c="k",linewidth=2)
        ax.scatter(self.data[:,0]/1000.,self.data[:,1],s=36,
                   marker='o',color='r',label="QTM transition")
        ax.scatter(self.miss[:,0]/1000.,self.miss[:,1],s=36,
                   marker='o',color='white',
                   edgecolor="black",label="no transition")
        
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