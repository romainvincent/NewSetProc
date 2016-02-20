# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:18:58 2014

@author: stefan
"""

import json
from os import chdir
from setproc.cycle_process.classes.cycle_process import CycleProcess

verz = '/home/stefan/1PhD/Python/test/' 
filename ='0_'
chdir(verz)
sweeps = 100

# cut trace file
with open(filename+"trace.json",'r') as json_data:
    data = json.load(json_data)

meas = data.pop("measures")

swp_nbr = size(meas)
loop_max = int(1.*swp_nbr/sweeps)

for i in range(loop_max):
    data['measures']=meas[i*sweeps:(i+1)*sweeps]
    if i<loop_max:
        json.dump(data,open(str(i)+filename+"trace.json",'w'))

if (loop_max*sweeps)%swp_nbr != 0:
    data['measures']=meas[loop_max*sweeps:]
    json.dump(data,open(str(loop_max)+filename+"trace.json",'w'))

del(data)        

# cut retrace file
with open(filename+"retrace.json",'r') as json_data:
    data = json.load(json_data)

meas = data.pop("measures")

swp_nbr = size(meas)
loop_max = int(1.*swp_nbr/sweeps)

for i in range(loop_max):
    data['measures']=meas[i*sweeps:(i+1)*sweeps]
    if i<loop_max:
        json.dump(data,open(str(i)+filename+"retrace.json",'w'))

if (loop_max*sweeps)%swp_nbr != 0:
    data['measures']=meas[loop_max*sweeps:]
    json.dump(data,open(str(loop_max)+filename+"retrace.json",'w'))
    loop_max += 1
    
del(data)



try:
    temp = CycleProcess(filename+"trace.json",filename+"retrace.json",range(loop_max),'Json') 
    temp.GetStat(0,4,1,17)
    temp.SaveAll('trace','retrace','Stat.bin')
except:
    print 'process failes'
