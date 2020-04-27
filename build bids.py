# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 23:56:37 2020

@author: lmitridati3
"""
      

#%% solve UC for elec system (CHP can offer all elec production)


class expando(object):
    '''   
        A small class which can have attributes set
    '''
    pass


class elec_UC:
    def __init__(self,m1,D_day):
        #self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        #self._load_data()
        self._build_model(m1,D_day)

    
    def optimize(self):
        self.model.optimize()

    #def computeIIS(self):
        #self.model.computeIIS()
    
    #def _load_data(self):
   
    def _build_model(self,m1,D_day):
        
        self.model = gb.Model()
        self._build_variables(m1,D_day)
        self._build_objective()
        self._build_constraints(m1,D_day)
    
    def _build_variables(self,m1,D_day):
        m = self.model

        self.variables.L = {} #elec consumption of heat pump
        for t in time:
            for h in HP:                
                self.variables.L[h,t] = m.addVar(lb=0,ub=rmin[h]*qmax[h],name='L({0},{1})'.format(h,t))

        self.variables.P = {} #elec production
        for t in time:
            for h in CHP+gen+incinerator:                
                self.variables.P[h,t] = m.addVar(lb=0,name='P({0},{1})'.format(h,t))


        self.variables.solar_prod = {} #elec production
        for t in time:
            for h in zone_DK:                
                self.variables.solar_prod[h,t] = m.addVar(lb=0,ub=max(solar_prod[h,hiy[D_day,t]],0),name='solar prod({0},{1})'.format(h,t))


                
        self.variables.wind_prod = {} #elec production
        for t in time:
            for h in zone_DK:                
                self.variables.wind_prod[h,t] = m.addVar(lb=0,ub=max(0,wind_prod[h,hiy[D_day,t]]),name='wind prod({0},{1})'.format(h,t))

        
        self.variables.S = {} #load shedding
        for t in time:
            for h in zone_DK:                
                self.variables.S[h,t] = m.addVar(lb=0,ub=elec_load[h,hiy[D_day,t]],name='S({0},{1})'.format(h,t))
                

        self.variables.F = {} #load shedding
        for t in time:
            for h in line:                
                self.variables.F[h,t] = m.addVar(lb=-ATC_lower[h,hiy[D_day,t]],ub=ATC_upper[h,hiy[D_day,t]],name='F({0},{1})'.format(h,t))

        self.variables.u = {} # status on/off
        
        for h in CHP+gen+incinerator:
           for t in time:                            
                self.variables.u[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='u({0},{1})'.format(h,t))
           for t in time_init:                       
                self.variables.u[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=u_init[m1,h,D_day,t],ub=u_init[m1,h,D_day,t],name='u({0},{1})'.format(h,t))

    

        self.variables.v_on = {} # status turn on
        
        for h in CHP+gen+incinerator:
           for t in time:                            
                self.variables.v_on[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='v on({0},{1})'.format(h,t))
           for t in time_init:                           
                self.variables.v_on[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=v_on_init[m1,g,D_day,t],ub=v_on_init[m1,g,D_day,t],name='v on({0},{1})'.format(h,t))


        self.variables.v_off = {} # status turn off
        
        for h in CHP+gen+incinerator:
           for t in time:                           
                self.variables.v_off[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='v off({0},{1})'.format(h,t))
           for t in time_init:                            
                self.variables.v_off[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=v_off_init[m1,g,D_day,t],ub=v_off_init[m1,g,D_day,t],name='v off({0},{1})'.format(h,t))


        self.variables.r = {} # start up cost
        for t in time:
            for h in CHP+gen+incinerator:
                self.variables.r[h,t] = m.addVar(lb=0,name='start up cost({0},{1})'.format(h,t)) # dispatch of electricity generators
                
        m.update()
    
    def _build_objective(self): # building the objective function for the heat maret clearing
  
        m = self.model    
              
        m.setObjective(gb.quicksum(c_elec[g]*self.variables.P[g,t] +  c_0[g]*self.variables.u[g,t] + self.variables.r[g,t] for t in time for g in CHP+gen+incinerator) + gb.quicksum(500*self.variables.S[n,t] for t in time for n in zone_DK) - gb.quicksum(max(c_heat[k] for k in CHP+boiler+incinerator)*rmin[g]*self.variables.L[g,t] for t in time for g in HP),   
            gb.GRB.MINIMIZE)
            
        
    def _build_constraints(self,m1,D_day):
    
        m = self.model
        

                
        self.constraints.start_cost = {}

        for g in CHP+gen+incinerator:        
            for x in range(len(time)):
                for h in time_start_range[g]:
                    self.constraints.start_cost[g,time[x],h] = m.addConstr(
                            self.variables.r[g,time[x]],
                            gb.GRB.GREATER_EQUAL,
                            c_start[g,h]*(self.variables.u[g,time[x]] - gb.quicksum(self.variables.u[g,time_extended[len(time_init)+x-k-1]] for k in range(h))),name='start up cost compute({0},{1},{2})'.format(g,time[x],h))


        self.constraints.time_on_min = {}
        self.constraints.time_off_min = {}

        for g in CHP+gen+incinerator:
            for x in range(len(time)):
                
                self.constraints.time_on_min[g,time[x]] = m.addConstr(
                        gb.quicksum(self.variables.v_on[g,time_extended[len(time_init)+x-k]] for k in range(time_on_min[g])),
                        gb.GRB.LESS_EQUAL,
                        self.variables.u[g,time_extended[len(time_init)+x]],name='min time on({0},{1})'.format(g,time[x]))

                self.constraints.time_off_min[g,time[x]] = m.addConstr(
                        gb.quicksum(self.variables.v_on[g,time_extended[len(time_init)+x-k]] for k in range(time_off_min[g])),
                        gb.GRB.LESS_EQUAL,
                        1-self.variables.u[g,time_extended[len(time_init)+x-time_off_min[g]]],name='min time off({0},{1})'.format(g,time[x]))
                
                
        self.constraints.uc = {}
        for g in CHP+gen+incinerator:
            for x in range(len(time)):
                self.constraints.uc[g,time[x]] = m.addConstr(
                        self.variables.u[g,time_extended[len(time_init)+x]]-self.variables.u[g,time_extended[len(time_init)+x-1]],
                        gb.GRB.EQUAL,      
                        self.variables.v_on[g,time_extended[len(time_init)+x]]-self.variables.v_off[g,time_extended[len(time_init)+x]],name='uc binary variables ({0},{1})'.format(g,time[x]))



        #4) ELEC GEN MAX PROD and total production variable
        
                    
        self.constraints.elec_maxprod = {}
        self.constraints.elec_minprod = {}
        
        for t in time:

            for g in gen+incinerator:
                
                self.constraints.elec_maxprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    pmax[g]*self.variables.u[g,t],name='elec maxprod({0},{1})'.format(g,t))
                
                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.GREATER_EQUAL,
                    pmin[g]*self.variables.u[g,t],name='elec minprod({0},{1})'.format(g,t))

            for g in CHP:
                
                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.GREATER_EQUAL,
                    max(0,(rmin[g]*fmin[g])/(rmin[g]*rho_elec[g]+rho_heat[g]))*self.variables.u[g,t],name='elec minprod({0},{1})'.format(g,t))                   

                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    fmax[g]/rho_elec[g]*self.variables.u[g,t],name='elec minprod({0},{1})'.format(g,t))                   

       
        # elec balance

        self.constraints.elec_balance = {}
        
        for t in time:
                for z in zone_DK: 
                    
                    self.constraints.elec_balance[z,t] = m.addConstr(
                        gb.quicksum(self.variables.P[g,t] for g in prod_zone[z])+self.variables.wind_prod[z,t]+self.variables.solar_prod[z,t]+self.variables.S[z,t],
                        gb.GRB.GREATER_EQUAL,
                        elec_load[z,hiy[D_day,t]] + gb.quicksum(self.variables.L[h,t] for h in HP_zone[z]) + gb.quicksum(self.variables.F[l,t] for l in line_start[z]) - gb.quicksum(self.variables.F[l,t] for l in line_end[z]),name='elec balance({0},{1})'.format(z,t))





class elec_prices_estimator:
    def __init__(self,m1,D_day):
        #self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        #self._load_data()
        self._build_model(m1,D_day)

    
    def optimize(self):
        self.model.optimize()

    #def computeIIS(self):
        #self.model.computeIIS()
    
    #def _load_data(self):
   
    def _build_model(self,m1,D_day):
        
        self.model = gb.Model()
        self._build_variables(D_day)
        self._build_objective()
        self._build_constraints(m1,D_day)
    
    def _build_variables(self,D_day):
        m = self.model

        self.variables.L = {} #elec consumption of heat pump
        for t in time:
            for h in HP:                
                self.variables.L[h,t] = m.addVar(lb=0,ub=rmin[h]*qmax[h],name='L({0},{1})'.format(h,t))

        self.variables.P = {} #elec production
        for t in time:
            for h in CHP+gen+incinerator:                
                self.variables.P[h,t] = m.addVar(lb=0,name='P({0},{1})'.format(h,t))


        self.variables.solar_prod = {} #elec production
        for t in time:
            for h in zone_DK:                
                self.variables.solar_prod[h,t] = m.addVar(lb=0,ub=max(solar_prod[h,hiy[D_day,t]],0),name='solar prod({0},{1})'.format(h,t))


                
        self.variables.wind_prod = {} #elec production
        for t in time:
            for h in zone_DK:                
                self.variables.wind_prod[h,t] = m.addVar(lb=0,ub=max(0,wind_prod[h,hiy[D_day,t]]),name='wind prod({0},{1})'.format(h,t))

        
        self.variables.S = {} #load shedding
        for t in time:
            for h in zone_DK:                
                self.variables.S[h,t] = m.addVar(lb=0,ub=elec_load[h,hiy[D_day,t]],name='S({0},{1})'.format(h,t))
                

        self.variables.F = {} #load shedding
        for t in time:
            for h in line:                
                self.variables.F[h,t] = m.addVar(lb=-ATC_lower[h,hiy[D_day,t]],ub=ATC_upper[h,hiy[D_day,t]],name='F({0},{1})'.format(h,t))
          
        m.update()
    
    def _build_objective(self): # building the objective function for the heat maret clearing
  
        m = self.model    
              
        m.setObjective(gb.quicksum(c_elec[g]*self.variables.P[g,t] for t in time for g in CHP+gen+incinerator) + gb.quicksum(500*self.variables.S[n,t] for t in time for n in zone_DK) - gb.quicksum(max(c_heat[k] for k in CHP+boiler+incinerator)*rmin[g]*self.variables.L[g,t] for t in time for g in HP),   
            gb.GRB.MINIMIZE)
            
        
    def _build_constraints(self,m1,D_day):
    
        m = self.model
        
        #4) ELEC GEN MAX PROD and total production variable
        
                    
        self.constraints.elec_maxprod = {}
        self.constraints.elec_minprod = {}
        
        for t in time:

            for g in gen+incinerator:
                
                self.constraints.elec_maxprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    pmax[g]*u_set[m1,g,D_day,t],name='elec maxprod({0},{1})'.format(g,t))
                
                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.GREATER_EQUAL,
                    pmin[g]*u_set[m1,g,D_day,t],name='elec minprod({0},{1})'.format(g,t))

            for g in CHP:
                
                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.GREATER_EQUAL,
                    max(0,(rmin[g]*fmin[g])/(rmin[g]*rho_elec[g]+rho_heat[g]))*u_set[m1,g,D_day,t],name='elec minprod({0},{1})'.format(g,t))                   

                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    fmax[g]/rho_elec[g]*u_set[m1,g,D_day,t],name='elec minprod({0},{1})'.format(g,t))                   

       
        # elec balance

        self.constraints.elec_balance = {}
        
        for t in time:
                for z in zone_DK: 
                    
                    self.constraints.elec_balance[z,t] = m.addConstr(
                        gb.quicksum(self.variables.P[g,t] for g in prod_zone[z])+self.variables.wind_prod[z,t]+self.variables.solar_prod[z,t]+self.variables.S[z,t],
                        gb.GRB.GREATER_EQUAL,
                        elec_load[z,hiy[D_day,t]] + gb.quicksum(self.variables.L[h,t] for h in HP_zone[z]) + gb.quicksum(self.variables.F[l,t] for l in line_start[z]) - gb.quicksum(self.variables.F[l,t] for l in line_end[z]),name='elec balance({0},{1})'.format(z,t))



#%%

   
#for D_day in [day[0],day[1]]:

def compute_bids(m1,D):
    
    global CHP
    global HP
    global incinerator
    global gen
    global zone_DK
    global time
    global time_init
    global day  
    
    global pmin
    global pmax
    global u_init
    global v_on_init
    global v_off_init
    global u_set
    global v_on_set
    global v_off_set
    global pmin_bid
    global pmax_bid
    global c_heat_bid
    global elec_price_estimate
    global zone_prod
    global rho_elec
    global rho_heat
    global rmin
    global c_elec
    global c_heat
    global boiler
    
    D_day =  day[D]
    
    UC = elec_UC(m1,day[D])
    UC.optimize()                       

    for g in gen+CHP+incinerator:
        
        for t in time:
            
            u_set[m1,g,day[D],t] = UC.variables.u[g,t].x
            v_on_set[m1,g,day[D],t] = UC.variables.v_on[g,t].x
            v_off_set[m1,g,day[D],t] = UC.variables.v_off[g,t].x
            

        if D<len(day)-1:
            for x in range(len(time_init)):
                
                u_init[m1,g,day[D+1],time_init[len(time_init)-x-1]] = u_set[m1,g,day[D],time[len(time)-x-1]]
                v_on_init[m1,g,day[D+1],time_init[len(time_init)-x-1]] = v_on_set[m1,g,day[D],time[len(time)-x-1]]
                v_off_init[m1,g,day[D+1],time_init[len(time_init)-x-1]] = v_off_set[m1,g,day[D],time[len(time)-x-1]]

    for g in gen:
        for t in time: 
            
            pmin_bid[m1,g,day[D],t] = u_set[m1,g,day[D],t]*pmin[g]
            pmax_bid[m1,g,day[D],t] = u_set[m1,g,day[D],t]*pmax[g]


    elec = elec_prices_estimator(m1,day[D])
    elec.optimize()
    
    for t in time:
        for z in zone_DK:
            elec_price_estimate[m1,z,day[D],t] = elec.constraints.elec_balance[z,t].Pi
            
        for g in CHP:
            c_heat_bid[m1,g,day[D],t] = []
            c_heat_bid[m1,g,day[D],t].append(elec_price_estimate[m1,zone_prod[g],day[D],t]*rho_heat[g]/rho_elec[g])
    
            if elec_price_estimate[m1,zone_prod[g],day[D],t] < c_elec[g]:
                c_heat_bid[m1,g,day[D],t].append(c_elec[g]*rmin[g] + c_heat[g] - elec_price_estimate[m1,zone_prod[g],day[D],t]*rmin[g])
            else:
                c_heat_bid[m1,g,day[D],t].append(elec_price_estimate[m1,zone_prod[g],day[D],t]*rho_heat[g]/rho_elec[g])
            
            
            
            
        for g in HP:
            c_heat_bid[m1,g,day[D],t] = max(0,elec_price_estimate[m1,zone_prod[g],day[D],t]*rmin[g])
        for g in incinerator:    
            c_heat_bid[m1,g,day[D],t] = max(0,(c_elec[g]*rmin[g] + c_heat[g]) - elec_price_estimate[m1,zone_prod[g],day[D],t]*rmin[g]) 
        for g in boiler:    
            c_heat_bid[m1,g,day[D],t] = c_heat[g]  
