#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 16:26:12 2020

@author: lmitridati3
"""
def solve_markets(m1,D):


    if m1 == 'HUC' or m1 == 'EAHUC':        
        compute_bids(m1,D)
        
    solve_heat(m1,D)
    solve_elec(m1,D)

#%%
for m1 in model_type:
    for D in range(len(day)):
        solve_markets(m1,D)

#%%
model_type = ['HUC','EAHUC']
DHN_zone_DK = {'DK1':['TVIS','AAR'],'DK2':['KBH']}

bid_selected_set = {(m,g,d,t):[u_set[m,g,d,t],0] for m in model_type for g in CHP for d in day for t in time}

for m in model_type:
    for g in CHP: 
        for d in day: 
            for t in time:
                if Q_set[m,g,d,t]>q_bid[g][0]:
                    bid_selected_set[m,g,d,t][1] = 1
    

DHN_dic = {'AAR':AAR,'KBH':KBH,'TVIS':TVIS}  
DHN_name = ['AAR','KBH','TVIS']
           
heat_price_estimate = {}

for m in model_type:
    for DHN in DHN_name:
        for d in day: 
            for t in time:
                lll = [c_heat_bid[m,g,d,t]*u_set[m,g,d,t] for g in prod_node[n] for n in DHN_dic[DHN] ]+[c_heat_bid[m,g,d,t][x]*bid_selected_set[m,g,d,t][x] for g in CHP_node[n] for n in DHN_dic[DHN] for x in range(2)]  
                heat_price_estimate[m,DHN,d,t] = max(lll)

#%%
RE_utilization={}
for m in ['HUC','EAHUC']:
    RE_utilization[m]= 100-100*(sum(wind_prod[z,t] + solar_prod[z,t] for z in zone_DK for t in time_year) -  sum(solar_prod_set[m,z,d,t] + wind_prod_set[m,z,d,t] for z in zone_DK for d in day for t in time))/sum(wind_prod[z,t] + solar_prod[z,t] for z in zone_DK for t in time_year)
    
#%%
                
nb_hours={g:0 for g in CHP+HP+incinerator}

loss_hours={}
                
for d in day:
    for t in time:
        for z in zone_DK:
            for DHN in DHN_zone_DK[z]:
                for n in DHN_dic[DHN]:

                    for g in HP_node[n]+incinerator_node[n]:
                        
                        if rmin[g]*elec_price_set['HUC',z,d,t]*u_set['HUC',g,d,t] > c_heat_bid['HUC',g,d,t]*u_set['HUC',g,d,t]:
                           nb_hours[g] = nb_hours[g] + 1
                        
                        loss_hours[g,d,t]= min(0, heat_price_estimate['HUC',DHN,d,t] - rmin[g]*elec_price_set['HUC',z,d,t] )*Q_set['HUC',g,d,t]
                           
                    for g in CHP_node[n]:
                        
                        if rho_heat[g]/rho_elec[g]*elec_price_set['HUC',z,d,t]*bid_selected_set['HUC',g,d,t][0] > c_heat_bid['HUC',g,d,t][0]*bid_selected_set['HUC',g,d,t][0] or (c_elec[g]*rmin[g] + c_heat[g] - elec_price_set['HUC',z,d,t]*rmin[g])*bid_selected_set['HUC',g,d,t][1] > c_heat_bid['HUC',g,d,t][1]*bid_selected_set['HUC',g,d,t][1] :
                           nb_hours[g] = nb_hours[g] + 1
                        
                        loss_hours[g,d,t]= min(0, heat_price_estimate['HUC',DHN,d,t] - max( rho_heat[g]/rho_elec[g]*elec_price_set['HUC',z,d,t]*bid_selected_set['HUC',g,d,t][0] , (c_elec[g]*rmin[g] + c_heat[g] - elec_price_set['HUC',z,d,t]*rmin[g])*bid_selected_set['HUC',g,d,t][1] ) )*Q_set['HUC',g,d,t]
 
            
loss_total = {g:sum(loss_hours[g,d,t] for d in day for t in time)/1000 for g in CHP+incinerator+HP}
sum(loss_total[g] for g in CHP+HP+incinerator)
                    
