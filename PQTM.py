# -*- coding: utf-8 -*-
'''
Created on Thu Feb  6 10:20:25 2014

@author: stefan
'''

import pylab 
from setproc.cycle_process.classes.cycle_process import CycleProcess
import os
from scipy.optimize import curve_fit
from setproc.sweep_set.functions.filter import filter
from setproc.common.functions.local_extrema import local_extrema
import pickle
from copy import deepcopy

class QTM:
    
    def __init__(self):
        
        self.off = 166
        
        self.vg = 2.1
        self.temp = []
        self.fig =[]
        self.ax1 = [];self.ax2 = [];self.ax3 = [];self.ax4 = [];
        self.ax5 = [];self.ax6 = [];self.ax7 = [];self.ax8 = [];
        
        self.jmp_thres = []
        self.th_tr = []
        self.th_rt = []
        self.maxl = []
        self.maxr =[]
        self.swp_nbr = 1
        self.ymax_ax5 = 0
        self.ymax_ax8 = 0
        
    def MakePlot(self):
        self.fig = figure(figsize=(16,8))
        
        self.ax1 = self.fig.add_subplot(241)
        self.ax2 = self.fig.add_subplot(242)
        self.ax3 = self.fig.add_subplot(243)
        self.ax4 = self.fig.add_subplot(244)
        self.ax5 = self.fig.add_subplot(245)
        self.ax6 = self.fig.add_subplot(246)
        self.ax7 = self.fig.add_subplot(247)
        self.ax8 = self.fig.add_subplot(248)
        self.fig.subplots_adjust(left=0.05,wspace=0.3,hspace=0.5)

    
    def ClearPlot(self):
        self.ax1.cla()
        self.ax2.cla() 
        self.ax3.cla()
        self.ax4.cla()
        self.ax5.cla()
        self.ax6.cla() 
        self.ax7.cla()
        self.ax8.cla()
        
        
    def OpenFile(self,filename,v):
        
        os.chdir(filename)
        try:
            self.temp = CycleProcess(str(round(v,3))+'V_Stat.bin')
            self.temp.LoadSweeps()
        except:
            self.temp = CycleProcess('V_trace.json','V_retrace.json',[v],'Json') 
            self.temp.GetStat(0,4,1,17)
            self.temp.SaveAll(str(round(v,3))+'V_trace',str(round(v,3))+'V_retrace',
                     str(round(v,3))+'V_Stat.bin')
        self.swp_nbr = self.temp.trace['sweep_number']
        
        
    def PlotStat1(self):
        
        self.ax1.set_title('jump stat trace')
        self.ax1.set_xlabel('log10(jump hight)')
        self.ax1.set_ylabel('occurence')
        self.ax4.set_title('jump stat retrace')
        self.ax4.set_xlabel('lo10(jump hight)')
        self.ax4.set_ylabel('occurence')
        
        'jump statistic'
        
        value_stat = [[], []]
        
        #filter the top and bottom 3%        
        off = int(self.swp_nbr*0.03)
        
        for i in range(size(self.temp['detection'])) :
            topush = abs(self.temp['detection'][i].value)
            towhere = self.temp['detection'][i].trace
            if(towhere) :
                if topush > 0 :
                    value_stat[0].append(topush)
            else :
                if topush > 0 :
                    value_stat[1].append(topush)
        
        value_stat[0] = sort(array(value_stat[0]))
        value_stat[1] = sort(array(value_stat[1]))
        
        self.ax1.hist(log10(value_stat[0][off:-off]), 150)
        self.ax4.hist(log10(value_stat[1][off:-off]), 150)
        draw()
        H1 = histogram(log10(value_stat[0][off:-off]), 150)
        H2 = histogram(log10(value_stat[0][off:-off]), 150)
        
        print 'click once in -jump stat trace & retrace- to define thresholds'
        # get threshold
        g = ginput(1)
        thres1 = g[0][0]
        self.ax1.plot([g[0][0],g[0][0]],[0,max(H1[0])],'r',lw=3)
        draw()
        g = ginput(1)
        thres2 = g[0][0]
        self.ax4.plot([g[0][0],g[0][0]],[0,max(H2[0])],'r',lw=3)
        draw()
        
        self.jmp_thres = 10**(0.5*(thres1+thres2))
        
        
    def PlotStat2(self):
        
        self.ax5.set_title('stat diff trace')
        self.ax5.set_xlabel('jump hight')
        self.ax5.set_ylabel('occurence')
        self.ax8.set_title('stat diff retrace')
        self.ax8.set_xlabel('jump hight')
        self.ax8.set_ylabel('occurence')
    
        diff_tr = []
        diff_rt = []
        for i in range(1,self.swp_nbr):
            data = self.temp.trace['data'][i][self.off:]
            diff_tr.append(data[0]-data[-1])
            data = self.temp.retrace['data'][i][self.off:]
            diff_rt.append(data[0]-data[-1])
            
        h_tr = histogram(diff_tr,150)
        h_rt = histogram(diff_rt,150) 
        self.ymax_ax5 = max(h_tr[0])
        self.ymax_ax8 = max(h_rt[0])
        
        #plot histogram jump hight
        self.ax5.hist(diff_tr,150)    
        self.ax8.hist(diff_rt,150)
        draw()
        #get threshold
        print 'click twice in -stat diff trace & retrace- to define thresholds'
        g = ginput(1)
        thres_tr1 = g[0][0]
        self.ax5.plot([thres_tr1,thres_tr1],[0,self.ymax_ax5],'r',lw=3)
        draw()
        g = ginput(1)
        thres_tr2 = g[0][0]
        self.ax5.plot([thres_tr2,thres_tr2],[0,self.ymax_ax5],'r',lw=3)
        draw()
        g = ginput(1)
        thres_rt1 = g[0][0]
        self.ax8.plot([thres_rt1,thres_rt1],[0,self.ymax_ax8],'r',lw=3) 
        draw()
        g = ginput(1)
        thres_rt2 = g[0][0]
        self.ax8.plot([thres_rt2,thres_rt2],[0,self.ymax_ax8],'r',lw=3) 
        draw()
        
        self.th_tr=sort(array([thres_tr1,thres_tr2]))
        self.th_rt=sort(array([thres_rt1,thres_rt2]))
        
    def PlotStat3(self):
        
        self.ax6.set_title('stat data left')
        self.ax6.set_xlabel('g (S)')
        self.ax6.set_ylabel('occurance')
        self.ax7.set_title('stat data right')
        self.ax7.set_xlabel('g (S)')
        self.ax7.set_ylabel('occurance')
        
        l=[]; r=[]
        
        for i in range(self.swp_nbr):
            r.append(self.temp.trace['data'][i][-1])
            r.append(self.temp.retrace['data'][i][self.off])
            l.append(self.temp.trace['data'][i][self.off]) 
            l.append(self.temp.retrace['data'][i][-1]) 
            
        self.ax6.hist(l,150)
        self.ax7.hist(r,150)
        draw()   
        
        #find maxima
        print 'click on the two maxima in -stat data left-'
        g = ginput(1)
        l1 = g[0][0]
        ymax = max(histogram(l,150)[0])
        self.ax6.plot([l1,l1],[0,ymax],'r',lw=3)
        draw()
        g = ginput(1)
        l2 = g[0][0]
        self.ax6.plot([l2,l2],[0,ymax],'r',lw=3)
        draw()
        
        print 'click on the two maxima in -stat data right-'
        g = ginput(1)
        r1 = g[0][0]
        ymax = max(histogram(r,150)[0])
        self.ax7.plot([r1,r1],[0,ymax],'r',lw=3)
        draw()
        g = ginput(1)
        r2 = g[0][0]
        self.ax7.plot([r2,r2],[0,ymax],'r',lw=3)
        draw()
        
        self.maxl = sort(array([l1,l2]))
        self.maxr = sort(array([r1,r2]))
        
        
    def AnalyzeDataVis(self):
        
        self.ax2.set_title('data trace')
        self.ax2.set_xlabel('B (T)')
        self.ax2.set_ylabel('g (S)')
        self.ax3.set_title('data retrace')
        self.ax3.set_xlabel('B (T)')
        self.ax3.set_ylabel('g (S)')
        
        self.ax2.plot([self.temp.trace['bias'][self.off],self.temp.trace['bias'][-1]],
                      [self.maxl[0],self.maxr[0]],'r',lw=5)
        self.ax2.plot([self.temp.trace['bias'][self.off],self.temp.trace['bias'][-1]],
                      [self.maxl[1],self.maxr[1]],'r',lw=5)
        self.ax3.plot([self.temp.trace['bias'][self.off],self.temp.trace['bias'][-1]],
                      [self.maxl[0],self.maxr[0]],'r',lw=5)
        self.ax3.plot([self.temp.trace['bias'][self.off],self.temp.trace['bias'][-1]],
                      [self.maxl[1],self.maxr[1]],'r',lw=5)
        draw()
        
        #init state
        diff_gs = abs(self.maxl[1]-self.temp.trace['data'][0][self.off])
        diff_es = abs(self.maxl[0]-self.temp.trace['data'][0][self.off])
        if diff_gs < diff_es:
            state = 'gs'
        else:
            state = 'es'
        
        #init last data point
        last_dp = self.temp.trace['data'][0][self.off]
        
        #cycle through data
        for i in range(1,self.swp_nbr):
            
            bias = self.temp.trace['bias'][self.off:]
            data = self.temp.trace['data'][i][self.off:]
            diff = data[0]-data[-1]
                 
            jmps = np.array(local_extrema(abs(filter(data, 4, 1, 17)),20))[2] 
            nbr = size(jmps[jmps>self.jmp_thres])
            self.ax2.plot(bias,data)
            self.ax5.plot([diff,diff],[0,self.ymax_ax5],'g',lw=3)
            draw()
            
            if (last_dp-data[0]) < (-0.5*(self.maxl[1]-self.maxl[0])):
                state = 'gs'
                print 'inelastic rel','es->gs'
            elif (last_dp-data[0]) > (0.5*(self.maxl[1]-self.maxl[0])):
                state = 'es'
                print 'inelastic exc','gs->es'
                
            if diff < self.th_tr[0]:
                state = 'es'
                print 'jump up',nbr,'es->es'
                
            elif diff > self.th_tr[1]:
                state = 'gs'
                print 'jump down',nbr,'gs->gs'
                
            else:
                if state == 'gs':
                    state = 'es'
                    if nbr != 0:
                        print 'double jump',nbr,'gs->es'
                    else:
                        print 'no jump',nbr,'gs->es'
                        
                elif state == 'es':
                    state = 'gs'
                    if nbr != 0:
                        print 'double jump',nbr,'es->gs'
                    else:
                        print 'no jump',nbr,'es->gs'
                        
            last_dp = data[-1]
            
            
            data = self.temp.retrace['data'][i][self.off:]
            diff = data[0]-data[-1]
        
            jmps = np.array(local_extrema(abs(filter(data, 4, 1, 17)),20))[2] 
            nbr = size(jmps[jmps>self.jmp_thres])
            self.ax3.plot(bias[::-1],data)
            self.ax8.plot([diff,diff],[0,self.ymax_ax8],'g',lw=3)
            draw()
            
            if (last_dp-data[0]) > (0.5*(self.maxr[1]-self.maxr[0])):
                state = 'gs'
                print 'inelastic rel','es->gs'
            elif (last_dp-data[0]) < (-0.5*(self.maxr[1]-self.maxr[0])):
                state = 'es'
                print 'inelastic exc','gs->es'
             
            if diff < self.th_rt[0]:
                state = 'gs'
                print 'jump up',nbr,'gs->gs'
                
            elif diff > self.th_rt[1]:
                state = 'es'
                print 'jump down',nbr,'es->es'
    
            else:
                if state == 'gs':
                    state = 'es'
                    if nbr != 0:
                        print 'double jump',nbr,'gs->es'
                    else:
                        print 'no jump',nbr,'gs->es'
                    
                elif state == 'es':
                    state = 'gs'
                    if nbr != 0:
                        print 'double jump',nbr,'es->gs'
                    else:
                        print 'no jump',nbr,'es->gs'
            
            last_dp = data[-1]
            
            ginput(1)
            del(self.ax2.lines[-1])
            del(self.ax3.lines[-1])
            del(self.ax5.lines[-1])
            del(self.ax8.lines[-1])
            
    def AnalyzeData(self):
        
        print 'start analysis'
        gs = 0
        es = 0
        gs_es = 0
        gs_gs = 0
        es_gs = 0
        es_es = 0        
        djmp = 0
        rel = 0
        exc = 0
        
        #init state
        diff_gs = abs(self.maxl[1]-self.temp.trace['data'][0][self.off])
        diff_es = abs(self.maxl[0]-self.temp.trace['data'][0][self.off])
        if diff_gs < diff_es:
            state = 'gs'
        else:
            state = 'es'
        
        #init last data point
        last_dp = self.temp.trace['data'][0][self.off]
        
        #cycle through data
        for i in range(1,self.swp_nbr):
            
            if (i-1)%int(self.swp_nbr/20.)==0:
                print 'completed:',int(100.*i/self.swp_nbr),'%'
            
            bias = self.temp.trace['bias'][self.off:]
            data = self.temp.trace['data'][i][self.off:]
            diff = data[0]-data[-1]
                 
            jmps = np.array(local_extrema(abs(filter(data, 4, 1, 17)),20))[2] 
            nbr = size(jmps[jmps>self.jmp_thres])
            
            if (last_dp-data[0]) < (-0.5*(self.maxl[1]-self.maxl[0])):
                rel += 1                
                state = 'gs'
            elif (last_dp-data[0]) > (0.5*(self.maxl[1]-self.maxl[0])):
                exc += 1                
                state = 'es'
                
            if diff < self.th_tr[0]:
                state = 'es'
                es += 1
                es_es += 1
                #print('jump up',nbr,'es->es')
                
            elif diff > self.th_tr[1]:
                state = 'gs'
                gs += 1
                gs_gs += 1
                #print('jump down',nbr,'gs->gs')
                
            else:
                if state == 'gs':
                    state = 'es'
                    if nbr != 0:
                        djmp += 1                        
                        #print('double jump',nbr,'gs->es')
                    else:
                        gs += 1
                        gs_es += 1
                        #print('no jump',nbr,'gs->es')  
                        
                elif state == 'es':
                    state = 'gs'
                    if nbr != 0:
                        djmp += 1
                        #print('double jump',nbr,'es->gs')
                    else:
                        es += 1
                        es_gs += 1
                        #print('no jump',nbr,'es->gs')
                        
            last_dp = data[-1]
            
            
            data = self.temp.retrace['data'][i][self.off:]
            diff = data[0]-data[-1]
        
            jmps = np.array(local_extrema(abs(filter(data, 4, 1, 17)),20))[2] 
            nbr = size(jmps[jmps>self.jmp_thres])
    
            if (last_dp-data[0]) > (0.5*(self.maxr[1]-self.maxr[0])):
                rel += 1                
                state = 'gs'
                #print('inelastic rel','es->gs')
            elif (last_dp-data[0]) < (-0.5*(self.maxr[1]-self.maxr[0])):
                exc += 1                
                state = 'es'
                #print('inelastic exc','gs->es')
             
            if diff < self.th_rt[0]:
                gs += 1
                gs_gs += 1
                state = 'gs'
                #print('jump up',nbr,'gs->gs')
                
            elif diff > self.th_rt[1]:
                es += 1
                es_es +=1
                state = 'es'
                #print('jump down',nbr,'es->es')
    
            else:
                if state == 'gs':
                    state = 'es'
                    if nbr != 0:
                        djmp += 1
                        #print('double jump',nbr,'gs->es')
                    else:
                        gs += 1
                        gs_es += 1
                        #print('no jump',nbr,'gs->es')
                    
                elif state == 'es':
                    state = 'gs'
                    if nbr != 0:
                        djmp += 1
                        #print('double jump',nbr,'es->gs')
                    else:
                        es += 1
                        es_gs += 1
                        #print('no jump',nbr,'es->gs')
            
            last_dp = data[-1]
            
        return[gs,gs_es,gs_gs,es,es_gs,es_es,djmp,rel,exc]
        
        
if __name__ == '__main__':
    
    res = []
       
    vg = arange(2.14,2.205,0.005)
    
    
    directory = '/home/stefan/Sionludi Remote/QTM/Stat_Vg/B_speed 0.05'
    try:
        res = pickle.load( open( directory+'analysis2.p', "rb" ) )
        
    except:
        p = QTM()
        p.MakePlot()
        for v in vg:
            p.ClearPlot()        
            p.OpenFile(directory,v)
            #p.swp_nbr = 10000                
            p.PlotStat1() 
            p.PlotStat2()
            p.PlotStat3()
            res.append(p.AnalyzeData())
            
        pickle.dump( res, open( "analysis.p", "wb" ) )
            
    res = reshape(array(res),(int(size(res)/9),9))
    """
    temp = deepcopy(res[1]) 
    res[1][0] = temp[3]; res[1][1] = temp[4]; res[1][2] = temp[5]; 
    res[1][3] = temp[0]; res[1][4] = temp[1]; res[1][5] = temp[2];
    res[1][7] = temp[8]; res[1][8] = temp[7];
    temp = deepcopy(res[2]) 
    res[2][0] = temp[3]; res[2][1] = temp[4]; res[2][2] = temp[5]; 
    res[2][3] = temp[0]; res[2][4] = temp[1]; res[2][5] = temp[2];
    res[2][7] = temp[8]; res[2][8] = temp[7];
    """
    fig = figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlabel("Vg")
    ax1.set_ylabel("probability")
    ax1.plot(vg,1.*res.T[2]/res.T[0],label=(r"$P_{LZ}$"))
    ax1.plot(vg,1.*res.T[5]/res.T[3],label=(r"$P_{LZ}^{\sf{inv}}$"))
    ax1.plot(vg,1.*res.T[6]/p.swp_nbr,label=(r"$\sf{2\ jumps}$")) 
    ax1.plot(vg,1.*res.T[7]/p.swp_nbr,label=(r"$\sf{relaxation}$")) 
    ax1.plot(vg,1.*res.T[8]/p.swp_nbr,label=(r"$\sf{excitation}$")) 
    legend(loc=4)
    os.chdir(directory)
    fig.savefig("PLZ.pdf")