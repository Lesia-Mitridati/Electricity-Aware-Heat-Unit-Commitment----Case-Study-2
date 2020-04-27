#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 19:07:41 2020

@author: lmitridati3
"""


class EM:
    def __init__(self,m1,D_day):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._build_model(m1,D_day)

    
    def optimize(self):
        self.model.optimize()
        
    def _build_model(self,m1,D_day):
        
        self.model = gb.Model()
        self._build_variables(m1,D_day)
        self._build_objective()
        self._build_constraints(m1,D_day)
    
    def _build_variables(self,m1,D_day):
        m = self.model

        # EM variables

        self.variables.P = {} #elec production
        for t in time:
            for h in CHP+gen:                
                self.variables.P[h,t] = m.addVar(lb=0,name='P({0},{1})'.format(h,t))

        self.variables.solar_prod = {} #elec production
        for t in time:
            for h in zone_DK:                
                self.variables.solar_prod[h,t] = m.addVar(lb=0,ub=max(solar_prod[h,hiy[D_day,t]],0),name='solar prod({0},{1})'.format(h,t))

                
        self.variables.wind_prod = {} #elec production
        for t in time:
            for h in zone_DK:                
                self.variables.wind_prod[h,t] = m.addVar(lb=0,ub=max(0,wind_prod[h,hiy[D_day,t]]),name='wind prod({0},{1})'.format(h,t))

        
        self.variables.ES = {} # elec load shedding
        for t in time:
            for h in zone_DK:
                self.variables.ES[h,t] = m.addVar(lb=0,ub=elec_load[h,hiy[D_day,t]],name='S({0},{1})'.format(h,t))

        self.variables.EF = {} #load shedding
        for t in time:
            for h in line:                
                self.variables.EF[h,t] = m.addVar(lb=-ATC_lower[h,hiy[D_day,t]],ub=ATC_upper[h,hiy[D_day,t]],name='elec Flow({0},{1})'.format(h,t))
        m.update()
    
    def _build_objective(self): # building the objective function for the EAHUC
  
        m = self.model    
              

        m.setObjective( gb.quicksum(c_elec[g]*self.variables.P[g,t] for t in time for g in CHP+gen) + gb.quicksum(500*self.variables.ES[n,t] for t in time for n in zone_DK) ,   
            gb.GRB.MINIMIZE)


        
    def _build_constraints(self,m1,D_day):
     
        m = self.model
                
        ## LL: electricity market constraints

                
        #4) ELEC GEN MAX PROD and total production variable

        
        self.constraints.elec_maxprod = {} # c18b
        self.constraints.elec_minprod = {} # c18a
        
        for t in time:

            for g in gen:

                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.GREATER_EQUAL,
                    pmin_bid[m1,g,D_day,t],name='elec minprod({0},{1})'.format(g,t))
                
                self.constraints.elec_maxprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    pmax_bid[m1,g,D_day,t],name='elec maxprod({0},{1})'.format(g,t))

            for g in CHP:

                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.GREATER_EQUAL,
                    pmin_bid[m1,g,D_day,t],name='elec minprod({0},{1})'.format(g,t))
                    
                self.constraints.elec_maxprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    pmax_bid[m1,g,D_day,t],name='elec maxprod({0},{1})'.format(g,t))                               

        # c19a and c19b = bounds on elec load shedding
        
        #c20a and c20b : bounds on solar and wind prod
        
        #c21a and c21b : bounds on solar and wind prod
        
        # c22a and c22b = bounds on elec flows      
        
        # elec balance
        
        self.constraints.elec_balance = {} # c23
        
        for t in time:
            for z in zone_DK:
                
                self.constraints.elec_balance[z,t] = m.addConstr(
                    gb.quicksum(self.variables.P[g,t] for g in gen_zone[z]+CHP_zone[z]) + gb.quicksum(rmin[h]*Q_set[m1,h,D_day,t] for h in incinerator_zone[z]) + self.variables.wind_prod[z,t]+self.variables.solar_prod[z,t]+self.variables.ES[z,t],
                    gb.GRB.EQUAL,
                    elec_load[z,hiy[D_day,t]] + gb.quicksum(rmin[h]*Q_set[m1,h,D_day,t] for h in HP_zone[z]) + gb.quicksum(self.variables.EF[l,t] for l in line_start[z]) - gb.quicksum(self.variables.EF[l,t] for l in line_end[z]),name='elec balance({0},{1})'.format(z,t))

#%%
                
def solve_elec(m1,D):
    
    global CHP
    global HP
    global incinerator
    global gen
    global zone_DK
    global time
    global time_init
    global day 
    global node
    global node_node
    
    global pmin
    global pmax
    global pmin_bid
    global pmax_bid
    global u_init
    global v_on_init
    global v_off_init
    global u_set
    global v_on_set
    global v_off_set
    global Q_set
    global c_heat_bid
    global elec_price_estimate
    global zone_prod
    global rho_elec
    global rho_heat
    global rmin
    global c_elec
    global c_heat
    global boiler
    global wind_prod_set
    global solar_prod_set
    
    global solve_it
    
    global elec_cost
    
    solve_it = EM(m1,day[D])
        
    solve_it.optimize()
    elec_cost[m1,day[D]] = solve_it.model.getObjective().getValue()

    for g in CHP+gen:           
        for t in time:

            P_set[m1,g,day[D],t] = solve_it.variables.P[g,t].x
            
    for z in zone_DK:
            
        for t in time:
            
            elec_price_set[m1,z,day[D],t] = solve_it.constraints.elec_balance[z,t].Pi
            wind_prod_set[m1,z,day[D],t] = solve_it.variables.wind_prod[z,t].x
            solar_prod_set[m1,z,day[D],t] = solve_it.variables.solar_prod[z,t].x

