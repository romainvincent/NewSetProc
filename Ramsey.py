# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 10:27:02 2013

@author: stefan
"""

import matplotlib.pylab as pl
import numpy as np
from setproc.cycle_process.classes.cycle_process import CycleProcess
from scipy.optimize.minpack import curve_fit
import pickle
from scipy.interpolate import interp1d

pi = np.pi       
thres = 10**-31.5 
off_tr = 0.123
off_rt = 0.148
bins=[-0.06,-0.050,-0.030,-0.023,-0.003,0.003,0.023,0.030,0.050,0.06]


def run(excl=0, plot_hist=False, pre = "0_"):
    His1 = []
    His2 = []
    His3 = []
    His4 = []
    g_min = []
    g_max = []
    for f in width:
        try:
            temp = CycleProcess(pre+str(round(f,2))+"us Stat.bin")
            temp.LoadSweeps()
        except:
            temp = CycleProcess("tau:"+str(round(f,2))+"us_trace.json",
                                "tau:"+str(round(f,2))+"us_retrace.json",
                                [pre],"Json")
            temp.GetStat(0,4,1,17)
            temp.SaveAll(pre+str(round(f,2))+"us trace",
                         pre+str(round(f,2))+"us retrace",
                         pre+str(round(f,2))+"us Stat.bin")
                         
        temp.SetPlotCalibration([thres,10**-25,[-0.06,0.06],0])   
        temp.GetAR()   
        
        if plot_hist == True:
            temp.GetValueStat()
            up_tr,down_tr = temp.trace.GetPeaks(thres,off_tr,'max',
                                                [-0.06,0.06],0)
            up_rt,down_rt = temp.retrace.GetPeaks(thres,off_rt,'max',
                                                  [-0.06,0.06],0)
                                               
            tr = np.row_stack((up_tr,down_tr))
            tr = tr[np.lexsort((tr[:,1],tr[:,0]))]
            rt = np.row_stack((up_rt,down_rt))
            rt = rt[np.lexsort((rt[:,1],rt[:,0]))]
            
            f3 = pl.figure()
            ax3 = f3.add_subplot(111)
            ax3b = ax3.twinx()
            ax3.hist(tr[:,1],[-0.048,-0.024,0,0.024,0.048],rwidth=0.95,
                     range=(-0.06,0.06),alpha=0.5)
            ax3b.hist(tr[:,1],150,range=(-0.08,0.08),color="black")
            f4 = pl.figure()
            ax4 = f4.add_subplot(111)
            ax4b = ax4.twinx()
            ax4.hist(rt[:,1],[-0.048,-0.024,0,0.024,0.048],rwidth=0.95,
                     range=(-0.06,0.06),alpha=0.5)
            ax4b.hist(rt[:,1],150,range=(-0.08,0.08),color="black")    
            pl.draw()
            
            f2 = pl.figure()
            ax1 = f2.add_subplot(111)
            H, xedges, yedges = np.histogram2d(np.array(temp["AvsR"][0])+off_tr,
                                            np.array(temp["AvsR"][1])+off_rt,
                                            bins=[bins,bins])
            extentim = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            im1 = ax1.imshow(H.T,extent = extentim, origin = 'lower',vmin=0,vmax=35)
            ax1.set_title(f)    
            ax1.set_xlabel("trace")
            ax1.set_ylabel("retrace")
            f2.colorbar(im1) 
            pl.draw()

        g_m = []        
        sz = int(np.size(temp.trace["bias"])/2)
        for i in range(temp.trace["sweep_number"]):
            g_m.append(temp.trace["data"][i][sz])
        g_m = np.array(g_m)
        g_min.append(g_m.min())
        g_max.append(g_m.max())
        
        nbr_itv = 1
        itv = np.round(1.0*np.size(temp["AvsR"][0])/nbr_itv)
        off1 = []
        off2 = []
        for kk in range(nbr_itv):
            if kk == nbr_itv-1:
                H1 , xedges, yedges = np.histogram2d(np.array(temp["AvsR"][0])[kk*itv:]+off_tr,
                           np.array(temp["AvsR"][1])[kk*itv:]+off_rt,
                           bins=[bins,bins])
            else:
                H1 , xedges, yedges = np.histogram2d(np.array(temp["AvsR"][0])[kk*itv:(kk+1)*itv]+off_tr,
                           np.array(temp["AvsR"][1])[kk*itv:(kk+1)*itv]+off_rt,
                           bins=[bins,bins])
                
            off1.append(H1[1][3]+H1[3][1])
            off2.append(H1[5][7]+H1[7][5])
        
        #delete the smallest and largest values        
        off1.sort()   
        off2.sort() 
        for ii in range(excl):
            off1.pop(0)
            off1.pop(-1)
            off2.pop(0)
            off2.pop(-1)
        
        His1.append(sum(off1))
        His2.append(sum(off2))
            
        
        temp.GetRA()
        itv = np.round(1.0*np.size(temp["RvsA"][0])/nbr_itv)
        off3 = []
        off4 = []
        for kk in range(nbr_itv):
            if kk == nbr_itv-1:
                H2 , xedges, yedges = np.histogram2d(np.array(temp["RvsA"][0])[kk*itv:]+off_tr,
                           np.array(temp["RvsA"][1])[kk*itv:]+off_rt,
                           bins=[bins,bins])
            else:
                H2 , xedges, yedges = np.histogram2d(np.array(temp["RvsA"][0])[kk*itv:(kk+1)*itv]+off_tr,
                           np.array(temp["RvsA"][1])[kk*itv:(kk+1)*itv]+off_rt,
                           bins=[bins,bins])
                
            off3.append(H2[1][3]+H2[3][1])#-H2[1][1]+H2[3][3])
            off4.append(H2[5][7]+H2[7][5])#-H2[5][5]-H2[7][7])
        
        #delete the two smallest and largest values        
        off3.sort()
        off4.sort()             
        for ii in range(excl):        
            off3.pop(0)
            off3.pop(-1)
            off4.pop(0)
            off4.pop(-1)
            
        His3.append(sum(off3))
        His4.append(sum(off4))
        
    norm = 1.0*temp.trace["sweep_number"]/8.0
    print norm    
    His1 = np.array(His1)/norm
    His2 = np.array(His2)/norm
    His3 = np.array(His3)/norm
    His4 = np.array(His4)/norm
    
    #rab_AR = His1[0::2]+His1[1::2]+His4[0::2]+His4[1::2]
    rab_AR = His1#+His4

    #rab_RA = His2[0::2]+His2[1::2]+His3[0::2]+His3[1::2]
    rab_RA = His2#+His3    
    
    return [His1,His2,His3,His4,np.array(g_min)*-1e9,np.array(g_max)*-1e9]
 
if __name__ == '__main__':
    width = np.arange(0.01,6.5,0.5) 
    pre = "0_"
    His1,His2,His3,His4, g_min, g_max = run(0,False,pre)
    
    His = np.array([His1,His4,His1+His4])    
    omega_s = 5
    
    try:
        tosave = pickle.load(open("analysis","rb"))
    except:
        tosave = []
        
    print tosave
    
    
    f1 = pl.figure(figsize=(5,9))
    f = lambda x,a,b,d: a*np.cos((x+0)/b*2*np.pi)+d
    
    for i in range(3):    
        ax = f1.add_subplot(3,1,i+1)    
        ax.set_xlabel(r"$width \ \mathsf{(\mu s)}$")
        ax.set_ylabel(r"$\mathsf{counts}$")
        ax.set_xlim((width[0],width[-1]))
        #ax.set_ylim((0,1.5))
        #ax1.set_yticklabels(())        
        params,dev = curve_fit(f,width,His[i],[-0.1,omega_s,0.2]) 
        time = np.linspace(width[0],width[-1],100)
        fit = f(time,params[0],params[1],params[2])
        trans = np.fft.fft(His[i])
        ax.text(0.05,0.8,r"$\Omega_{R} = $"+str(round(100/params[1])/100)+" MHz",
                 fontsize=24)
        print("1: ",round(100*params[0])/100,round(100/params[1])/100,
              round(100*np.real(trans[1]))/100,round(100*np.real(trans[2]))/100)
        tosave.append(round(100*params[0])/100);
        tosave.append(round(100/params[1])/100);
        tosave.append(round(100*np.real(trans[1]))/100);    
        ax.scatter(width,His[i],s=100,c="r",marker="o")
        ax.plot(time,fit,'r--',linewidth=3)
        
    
    """
    ax1c = ax1.twiny()
    ax1c.set_xticks((np.arange(0,2,0.5)*np.pi))
    ax1c.set_xticklabels((0,r"$\frac{\pi}{2}$",r"$\pi$",r"$\frac{3 \pi}{2}$",
                          r"$2 \pi$",r"$\frac{5\pi}{2}$",r"$2\pi$"))
    ax1c.set_xlim((0,params[1]))
    ax1.plot([params[1]/4,params[1]/4],[0,rab_AR.max()*1.1)],"k--")    
    ax1.text(0.15,0,r"$\tau \left(\frac{\pi}{2}\right) = $"+str(round(params[1]/4*10000)/10)+"ns",size="x-large")
    """

    
    ax.set_title(r"$g = ($"+str(round((g_max.max()+g_min.min()))/2)+
                  r"$ \pm$" +str(round((g_max.max()-g_min.min()))/2)+") nS")
    
    f1.subplots_adjust(left=0.15,top=0.95,bottom=0.05,hspace=0)
    
    f1.savefig(pre+"Rabi.png",format="png",dpi=300)
    f1.savefig(pre+"Rabi.pdf",format="pdf",)
    
        
    f4  = pl.figure()
    ax4 = f4.add_subplot(111)
    ax4.set_xlabel(r"$width \ \mathsf{(\mu s)}$")
    ax4.set_ylabel(r"$g \ (nS)$")
    ax4.scatter(width,g_min,s=50,c="blue",label="g_max")
    ax4.scatter(width,g_max,s=50,c="red",label="g_min")
    ax4.legend(loc=7)
    
    Ram = (His4[0::2]+His4[1::2])/2
    tau = (width[0::2]+width[1::2])/2    
    f5  = pl.figure()
    ax5 = f5.add_subplot(111)
    ax5.set_xlim((tau[0],tau[-1]))
    ax5.set_xlabel(r"$width \ \mathsf{(\mu s)}$")
    ax5.set_ylabel(r"$\mathsf{counts\ (a.u.)}$")
    ax5.scatter(tau,Ram,s=100,c="red")
    ax5.plot(tau,f(tau,0.4,5.15,0.55),"k--",linewidth=3)
    
    """
    f2 = plt.figure()
    ax2 = f2.add_subplot(111)
    ax2.set_xlabel(r"$width \ \mathsf{(\mu s)}$")
    ax2.set_ylabel(r"$\mathsf{counts}$")
    ax2.set_title("0dBm")
    ax2.set_xlim((0,width[-1]))
    ax2.set_ylim((0,1.0))
    ax2.scatter(width,rab_RA,s=100,c="r",marker="o")
    params,dev = curve_fit(f,width,rab_RA,[-0.1,0.3,0.3])
    print("2: ",params[0],params[1])
    time = np.linspace(0,width[-1],100)
    fit = f(time,params[0],params[1],params[2])
    ax2.plot(time,fit,'r--',linewidth=3)
    ax2.set_title(r"$g = ($"+str(round((g_max.max()+g_min.min()))/2)+
                  r"$ \pm$" +str(round((g_max.max()-g_min.min()))/2)+") nS")
    ax2.text(0.3,0.05,r"$\Omega_R = $"+str(round(100/params[1])/100)+" MHz",
             fontsize=24)
    f2.subplots_adjust(left=0.15,top=0.87)
    """
    
    