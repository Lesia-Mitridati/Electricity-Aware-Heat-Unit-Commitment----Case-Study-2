# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 08:57:30 2016

@author: lemitri
"""

#%%
            
import os
import pandas as pd
import scipy.stats as sp
#import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
sb.set_style('ticks')

import gurobipy as gb
import itertools as it

import numpy as np

os.chdir("/Users/lmitridati3/Documents/Postdoc projects/Electricity-aware unit commitment/Python/Electricity Aware Unit Commitment copy")

 
#%% initialize/indexes
    
D=365
T=24
day = ['d{0:03d}'.format(t+1) for t in range(D)]
time = ['h{0:02d}'.format(t) for t in range(T)]
time_year =  ['t{0:04d}'.format(t) for t in range(T*D)]

N=14+9+6
KBH = ['KBH_{0}'.format(n+1) for n in range(14)]
AAR = ['AAR_{0}'.format(n+1) for n in range(9)]
TVIS = ['TVIS_{0}'.format(n+1) for n in range(6)]
node = TVIS + AAR + KBH
zone_DK = ['DK1','DK2']
#zone = ['DK1','DK2','DE','NO2','SE3','SE4']
zone = ['DK1','DK2']

prod_type=['gen','CHP','incinerator','boiler','HP','heat storage']
gen=['G{0}'.format(x+1) for x in range(117)]
CHP = ['CHP{0}'.format(x+1) for x in range(11)]
HP =['HP{0}'.format(x+1) for x in range(6)]
incinerator = ['IS{0}'.format(x+1) for x in range(6)]
boiler = ['HO{0}'.format(x+1) for x in range(20)]
heat_storage = ['HS{0}'.format(x+1) for x in range(3)]

Nb = 2
    
model_type = ['HUC','EAHUC','CHPUC']

#%% import time series for whole year (wind, solar, loads)

wind_prod_csv = pd.read_csv('wind_prod.csv',index_col='hiy')
solar_prod_csv = pd.read_csv('solar_prod.csv',index_col='hiy')
heat_load_csv = pd.read_csv('heat_load.csv',index_col='hiy')
elec_load_csv = pd.read_csv('elec_load.csv',index_col='hiy')

hiy = {}
for t in range(len(time_year)):
    hiy['d{0:03d}'.format(np.int(heat_load_csv['diy'][t])),'h{0:02d}'.format(np.int(heat_load_csv['hid'][t]))]=time_year[t]

dhiy = {}
for t in range(len(time_year)):
    dhiy[time_year[t]]=['d{0:03d}'.format(np.int(heat_load_csv['diy'][t])),'h{0:02d}'.format(np.int(heat_load_csv['hid'][t]))]

    
heat_load = {}
for t in range(len(time_year)):
    for n in node:
        heat_load[n,time_year[t]]=heat_load_csv[n][t]
                
elec_load = {}
for t in range(len(time_year)):
    for n in zone_DK:
        elec_load[n,time_year[t]]=elec_load_csv[n][t]

       
wind_prod = {}
for t in range(len(time_year)):
    for n in zone_DK:
        wind_prod[n,time_year[t]]=wind_prod_csv[n][t]
    


solar_prod = {}
for t in range(len(time_year)):
    for n in zone_DK:
        solar_prod[n,time_year[t]]=solar_prod_csv[n][t]



#%% load generators data


generators_csv = pd.read_csv('generators.csv',index_col = 'UNIT ID')
  
prod_node={n:[] for n in node}
incinerator_node={n:[] for n in node}
HP_node={n:[] for n in node}
heat_storage_node={n:[] for n in node}
CHP_node={n:[] for n in node}

for n in node:
    for g in HP:
        if generators_csv['Node'][g] == n:
            prod_node[n].append(g)
            HP_node[n].append(g)
            
    for g in incinerator:
        if generators_csv['Node'][g] == n:
            prod_node[n].append(g)
            incinerator_node[n].append(g)
            
    for g in boiler:
        if generators_csv['Node'][g] == n:
            prod_node[n].append(g)

    
    for g in heat_storage:
        if generators_csv['Node'][g] == n:
            heat_storage_node[n].append(g)
        
    for g in CHP:
        if generators_csv['Node'][g] == n:
            CHP_node[n].append(g)
            
CHP_zone={n:[] for n in zone_DK}
gen_zone={n:[] for n in zone_DK}
HP_zone={n:[] for n in zone_DK}
incinerator_zone={n:[] for n in zone_DK}
prod_zone={n:[] for n in zone_DK}

for n in zone_DK:
    for g in CHP+gen+incinerator:
        if generators_csv['Zone'][g] == n:
            prod_zone[n].append(g)

    for g in HP:
        if generators_csv['Zone'][g] == n:
            HP_zone[n].append(g)

    for g in gen:
        if generators_csv['Zone'][g] == n:
            gen_zone[n].append(g)


    for g in CHP:
        if generators_csv['Zone'][g] == n:
            CHP_zone[n].append(g)

    for g in incinerator:
        if generators_csv['Zone'][g] == n:
            incinerator_zone[n].append(g)
            
            
zone_prod = {}
for g in HP+CHP+gen+incinerator:
    zone_prod[g] = generators_csv['Zone'][g]
        
node_prod = {}
for g in HP+CHP+incinerator+heat_storage+boiler:
    node_prod[g] = generators_csv['Node'][g]


#%%
            
pmax = {}
pmin = {}
qmax = {}
qmin={}
fmax={}
fmin={}
c_elec={}
c_heat={}
c_0={}
c_start={}
emax={}
rho_elec={}
rho_heat={}
rmin = {}
rho_charge = {}
rho_discharge = {}
heat_loss = {}
time_start = {}
time_on_min = {}
time_off_min = {}
q_bid = {g:[] for g in CHP}

for g in gen:
    pmax[g] = generators_csv['Pmax'][g]           
    pmin[g] = generators_csv['Pmin'][g]
    c_elec[g] = generators_csv['Ce'][g]
    c_0[g] = generators_csv['C0'][g]
    time_start[g] = np.int(generators_csv['Ns'][g])
    time_on_min[g] = max(1,time_start[g])
    time_off_min[g] = max(0,time_start[g]-1)
    
for g in CHP:
    qmax[g] = generators_csv['Qmax'][g]
    qmin[g] = generators_csv['Qmin'][g]
    fmax[g] = generators_csv['Fmax'][g]
    fmin[g] = generators_csv['Fmin'][g]
    rho_elec[g] = generators_csv['Relec'][g]
    rho_heat[g] = generators_csv['Rheat'][g]
    rmin[g] = generators_csv['Rmin'][g]
    c_0[g] = generators_csv['C0'][g]
    c_elec[g] = generators_csv['Ce'][g]
    c_heat[g] = generators_csv['Ch'][g]
    time_start[g] = 3
    time_on_min[g] = 3
    time_off_min[g] = 2
    q_bid[g].append(fmin[g]/(rmin[g]*rho_elec[g]+rho_heat[g]))
    q_bid[g].append(qmax[g]-q_bid[g][0])
    
for g in incinerator:
    pmax[g] = generators_csv['Pmax'][g]           
    pmin[g] = generators_csv['Pmin'][g]
    qmax[g] = generators_csv['Qmax'][g]           
    qmin[g] = generators_csv['Qmin'][g]    
    rmin[g] = generators_csv['Rmin'][g]
    c_0[g] = generators_csv['C0'][g]
    c_elec[g] = generators_csv['Ce'][g]
    c_heat[g] = generators_csv['Ch'][g]
    time_start[g] = 3
    time_on_min[g] = 3
    time_off_min[g] = 2

for g in boiler:
    qmax[g] = generators_csv['Qmax'][g]           
    qmin[g] = generators_csv['Qmin'][g]
    c_0[g] = generators_csv['C0'][g]
    c_heat[g] = generators_csv['Ch'][g]
    time_start[g] = 1
    time_on_min[g] = 0
    time_off_min[g] = 0   
    
for g in HP:
    qmax[g] = generators_csv['Qmax'][g]           
    qmin[g] = generators_csv['Qmin'][g]    
    rmin[g] = generators_csv['Rmin'][g]
    c_0[g] = generators_csv['C0'][g] 
    c_heat[g] = generators_csv['Ch'][g]
    time_start[g] = 1
    time_on_min[g] = 1
    time_off_min[g] = 1   

for g in heat_storage:
    qmax[g] = generators_csv['Qmax'][g]
    emax[g] = generators_csv['Emax'][g]
    rho_charge[g] = 0.95
    rho_discharge[g] = 1.05
    heat_loss[g] =  emax[g]/100


time_start_range = {g:[x+1 for x in range(time_start[g])] for g in gen+CHP+boiler+HP+incinerator}

c_start = {}
for g in CHP+HP+boiler+incinerator+gen:
    for x in time_start_range[g]:
        c_start[g,x] = generators_csv['Cs'][g] + generators_csv['Cs'][g]/2*x
        
        
        
#%% build heat and elec networks topology / and link 2 networks

heat_network_topology_csv = pd.read_csv('heat_network_topology.csv',index_col='Node')
elec_network_topology_csv = pd.read_csv('elec_network_topology.csv',index_col='Zone')
        
pipeline = []
pipeline_connexion = {}
node_node = {n:[] for n in node}
for n in node:
    for m in node:
        if heat_network_topology_csv[n][m] ==1:
            node_node[n].append(m)
            pipeline.append('P{0:02d}'.format(len(pipeline)+1))
            pipeline_connexion[pipeline[-1]]=[n,m]

pipeline_end = {z:[] for z in node}
pipeline_start = {z:[] for z in node}

for z in node:
    for l in pipeline:
        if pipeline_connexion[l][0]==z:
           pipeline_start[z].append(l) 
        if pipeline_connexion[l][1]==z:
           pipeline_end[z].append(l) 
           
line = []
line_connexion = {}
zone_zone ={z:[] for z in zone}
for n in range(len(zone)):
    for m in zone[:n]:
        if elec_network_topology_csv[zone[n]][m] ==1:
            zone_zone[zone[n]].append(m) 
    for m in zone[n:]:
        if elec_network_topology_csv[zone[n]][m] ==1:
            line.append('L{0:02d}'.format(len(line)+1))
            line_connexion[line[-1]]=[zone[n],m]
            zone_zone[zone[n]].append(m)

line_end = {z:[] for z in zone}
line_start = {z:[] for z in zone}

for z in zone:
    for l in line:
        if line_connexion[l][0]==z:
           line_start[z].append(l) 
        if line_connexion[l][1]==z:
           line_end[z].append(l) 

    
zone_DK_node={'DK1':AAR+TVIS,'DK2':KBH}

node_zone_DK={}
for n in node:
    for z in zone_DK:
        if n in zone_DK_node[z]:
            node_zone_DK[n]=z
            
ATC_upper_csv = pd.read_csv('ATC_upper.csv',index_col='hiy')
ATC_lower_csv = pd.read_csv('ATC_lower.csv',index_col='hiy')

ATC_upper = {}
for t in range(len(time_year)):
    for l in line:
        ATC_upper[l,time_year[t]]=ATC_upper_csv[l][t]

ATC_lower = {}
for t in range(len(time_year)):
    for l in line:
        ATC_lower[l,time_year[t]]=ATC_lower_csv[l][t]

time_delay_max = 0
        
f_upper = {}
for n1 in node:
    for n2 in node_node[n1]:
        f_upper[n1,n2] = min(1000,max(sum(qmax[g] for g in prod_node[n1]),sum(qmax[g] for g in prod_node[n2])))


#%% initial states and bids
T_init = max(time_delay_max,3)
   
time_init=['t00-{0}'.format(T_init-k) for k in range(T_init)]   
time_final=['t+{0}'.format(k+1) for k in range(time_delay_max)]   
time_year_extended= time_init+time_year
time_extended= time_init+time

u_init = {(m1,g,D_day,t):0 for g in CHP+HP+boiler+incinerator+gen for D_day in day for t in time_init for m1 in model_type}
v_on_init = {(m1,g,D_day,t):0 for g in CHP+HP+boiler+incinerator+gen for D_day in day for t in time_init for m1 in model_type}
v_off_init = {(m1,g,D_day,t):0 for g in CHP+HP+boiler+incinerator+gen  for D_day in day for t in time_init for m1 in model_type}
energy_stored_init = {(m1,g,D_day,t):emax[g]/2 for g in heat_storage  for D_day in day for t in time_init for m1 in model_type}
HF_in_init = {(m1,n1,n2,D_day,t):0 for n1 in node for n2 in node_node[n1] for D_day in day for t in time_init for m1 in model_type}
HF_direction_init = {(m1,n1,n2,D_day,t):1 for n1 in node for n2 in node_node[n1] for D_day in day for t in time_init for m1 in model_type}
time_delay_init = {}
for m1 in model_type:
    for n1 in node:
        for n2 in node_node[n1]:
            for D_day in day:
                for t in time_init:
                    time_delay_init[m1,n1,n2,D_day,t,0] = 1
                    for x in range(1,time_delay_max+1):
                        time_delay_init[m1,n1,n2,D_day,t,x] = 0 
                    
u_set = {}
v_on_set = {}
v_off_set = {}
Q_set = {}
elec_price_estimate = {}
c_heat_bid = {}
pmin_bid = {}
pmax_bid = {}
HF_in_set = {}
HF_direction_set = {}
time_delay_set = {}
energy_stored_set = {}
elec_price_set = {}
P_set = {}   
wind_prod_set = {}
solar_prod_set = {}

heat_cost = {}
elec_cost = {}
total_cost = {}

