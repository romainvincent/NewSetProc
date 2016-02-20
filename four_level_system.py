# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 10:27:02 2013

@author: stefan
"""

import matplotlib.pylab as pl
import numpy as np
from setproc.cycle_process.classes.cycle_process import CycleProcess
from scipy.optimize.minpack import curve_fit

pre = "2.207V_"
pi = np.pi
width = np.arange(0.05,0.75,0.05)        
thres = 10**-31.5 
off_tr = 0.123
off_rt = 0.148
bins=[-0.06,-0.050,-0.030,-0.023,-0.003,0.003,0.023,0.030,0.050,0.06]


def run(excl=0, plot_hist=False):
    His1 = []
    His2 = []
    His3 = []
    His4 = []
    for f in width:
        try:
            temp = CycleProcess(pre+str(round(f,2))+"us Stat.bin")
            temp.LoadSweeps()
        except:
            temp = CycleProcess("width:"+str(round(f,2))+"us_trace.json",
                                "width:"+str(round(f,2))+"us_retrace.json",
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
            
            f3 = plt.figure()
            ax3 = f3.add_subplot(111)
            ax3b = ax3.twinx()
            ax3.hist(tr[:,1],[-0.048,-0.024,0,0.024,0.048],rwidth=0.95,
                     range=(-0.06,0.06),alpha=0.5)
            ax3b.hist(tr[:,1],150,range=(-0.08,0.08),color="black")
            f4 = plt.figure()
            ax4 = f4.add_subplot(111)
            ax4b = ax4.twinx()
            ax4.hist(rt[:,1],[-0.048,-0.024,0,0.024,0.048],rwidth=0.95,
                     range=(-0.06,0.06),alpha=0.5)
            ax4b.hist(rt[:,1],150,range=(-0.08,0.08),color="black")    
            draw()
            
            f2 = plt.figure()
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
            draw()
        
        nbr_itv = 10
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
                
            off1.append(H1[1][3]+H1[3][1])#-H1[1][1]-H1[3][3])
            off2.append(H1[5][7]+H1[7][5])#-H1[5][5]-H1[7][7])
        
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
    His1 = np.array(His1)/norm
    His2 = np.array(His2)/norm
    His3 = np.array(His3)/norm
    His4 = np.array(His4)/norm
    
    #rab_AR = His1[0::2]+His1[1::2]+His4[0::2]+His4[1::2]
    rab_AR = His1+His4

    #rab_RA = His2[0::2]+His2[1::2]+His3[0::2]+His3[1::2]
    rab_RA = His2+His3
    
    return [rab_AR,rab_RA]
 
if __name__ == '__main__':
    
    rab_AR, rab_RA = run(0,False)

    f1 = pl.figure()
    ax1 = f1.add_subplot(111)    
    ax1.set_xlabel(r"$width \ \mathsf{(\mu s)}$")
    ax1.set_ylabel(r"$\mathsf{counts}$")
    ax1.set_xlim((0,width[-1]))
    ax1.set_ylim((0,rab_AR.max()*1.1))
    #ax1.set_yticklabels(())
    f = lambda x,a,b,c: a*np.cos(x/b*2*np.pi)+c
    params,dev = curve_fit(f,width,rab_AR,[-0.1,0.5,0.25]) 
    time = np.linspace(0,width[-1],100)
    fit = f(time,-0.1,0.55,0.26)
    """    
    ax1b = ax1.twiny()
    ax1b.set_xticks((np.arange(0,2,0.5)*np.pi))
    ax1b.set_xticklabels((0,r"$\frac{\pi}{2}$",r"$\pi$",r"$\frac{3 \pi}{2}$",
                          r"$2 \pi$",r"$\frac{5\pi}{2}$",r"$2\pi$"))
    ax1b.set_xlim((0,params[1]))
    ax1.plot([params[1]/4,params[1]/4],[0,rab_AR.max()*1.1)],"k--")    
    ax1.text(0.15,0,r"$\tau \left(\frac{\pi}{2}\right) = $"+str(round(params[1]/4*10000)/10)+"ns",size="x-large")
    """
    ax1.scatter(width,rab_AR,s=100,c="r",marker="o",label="f =2.4512GHz")
    ax1.plot(time,fit,'r--',linewidth=3)
    ax1.set_title(r"$g = (243\pm50)nS$")
    f1.subplots_adjust(left=0.15,top=0.87)
    
    f1.savefig(pre+"Rabi.png",format="png",dpi=300)
    f1.savefig(pre+"Rabi.pdf",format="pdf",)
    """
    f2 = plt.figure()
    ax2 = f2.add_subplot(111)
    ax2.set_xlabel(r"$width \ \mathsf{(\mu s)}$")
    ax2.set_ylabel(r"$\mathsf{counts}$")
    ax2.set_title("0dBm")
    ax2.set_xlim((0,width[-1]))
    ax2.scatter(width,rab_RA,s=100,c="b",marker="o",label="off01+10")
    params,dev = curve_fit(f,time,rab_RA,[-15,0.25,15])
    fit = f(width,params[0],params[1],params[2])
    ax2.plot(width,fit,'b--',linewidth=3)
    draw()
    """