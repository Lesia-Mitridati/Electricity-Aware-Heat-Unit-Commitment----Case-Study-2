#%%
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:43:46 2016

@author: lemitri
"""

#TODO: 1) bid validity constraint 2) objective value 3) dual constraints :( 4) EM with values of heat commitment fixed 5) LOOP over all days

#%% HUC


class HUC:
    def __init__(self,D_day):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._build_model(D_day)

    
    def optimize(self):
        self.model.optimize()
        
    def _build_model(self,D_day):
        
        self.model = gb.Model()
        self.model.Params.TimeLimit = 60
        self._build_variables(D_day)
        self._build_objective(D_day)
        self._build_constraints(D_day)
    
    def _build_variables(self,D_day):
        m = self.model


        ## UL: binary variables and start up cost
        
        self.variables.u = {} # status on/off
        
        for h in CHP+boiler+incinerator+HP:
           for t in time:                            
                self.variables.u[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='u({0},{1})'.format(h,t))
           for t in time_init:                       
                self.variables.u[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=u_init['HUC',h,D_day,t],ub=u_init['HUC',h,D_day,t],name='u({0},{1})'.format(h,t))


        self.variables.bid_accepted = {} # status on/off
        
        for t in time: 
            for h in CHP:      
                for x in range(2):                     
                    self.variables.bid_accepted[h,t,x] = m.addVar(vtype=gb.GRB.BINARY,name='bid accepted({0},{1},{2})'.format(h,t,x))
    

        self.variables.v_on = {} # status turn on
        
        for h in CHP+boiler+incinerator+HP:
           for t in time:                            
                self.variables.v_on[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='v on({0},{1})'.format(h,t))
           for t in time_init:                           
                self.variables.v_on[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=v_on_init['HUC',h,D_day,t],ub=v_on_init['HUC',h,D_day,t],name='v on({0},{1})'.format(h,t))


        self.variables.v_off = {} # status turn off
        
        for h in CHP+boiler+incinerator+HP:
           for t in time:                           
                self.variables.v_off[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='v off({0},{1})'.format(h,t))
           for t in time_init:                            
                self.variables.v_off[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=v_off_init['HUC',h,D_day,t],ub=v_off_init['HUC',h,D_day,t],name='v off({0},{1})'.format(h,t))


        self.variables.r = {} # start up cost
        for t in time:
            for h in CHP+boiler+incinerator+HP:
                self.variables.r[h,t] = m.addVar(lb=0,name='start up cost({0},{1})'.format(h,t)) # dispatch of electricity generators


        self.variables.time_delay = {} # time delay auxiliary variables self.variables.time_delay[n1,n2,t,x] = 1 iif time delay of heat flow entering pipeline (node n1) at time t, exits pipeline (node n2) at time t+x (time_delay[n1,n2,t]=sum_{x=0}^{time_delay_max} x*self.variables.time_delay[n1,n2,t,x] )
        
        for n1 in node:
            for n2 in node_node[n1]:
                    for t in time:
                        for x in range(time_delay_max+1):
                            self.variables.time_delay[n1,n2,t,x] = m.addVar(vtype = gb.GRB.BINARY,name='time delay auxiliary variables({0},{1},{2},{3}])'.format(n1,n2,t,x)) # dispatch of electricity generators

                    for t in time_init:
                        for x in range(time_delay_max+1):
                            self.variables.time_delay[n1,n2,t,x] = m.addVar(vtype = gb.GRB.BINARY,lb=time_delay_init['HUC',n1,n2,D_day,t,x],ub=time_delay_init['HUC',n1,n2,D_day,t,x],name='time delay auxiliary variables({0},{1},{2},{3}])'.format(n1,n2,t,x)) # dispatch of electricity generators

    
        self.variables.HF_direction = {} # direction of flow
        for n1 in node:
            for n2 in node_node[n1]:
                for t in time:
                    self.variables.HF_direction[n1,n2,t] = m.addVar(vtype = gb.GRB.BINARY,name='HF direction({0},{1},{2})'.format(n1,n2,t)) # dispatch of electricity generators

                for t in time_init:
                    self.variables.HF_direction[n1,n2,t] = m.addVar(vtype = gb.GRB.BINARY,lb=HF_direction_init['HUC',n1,n2,D_day,t],ub=HF_direction_init['HUC',n1,n2,D_day,t],name='HF direction({0},{1},{2})'.format(n1,n2,t)) # dispatch of electricity generators

        # ML: HM variables
        
        self.variables.Q = {} #elec production
        
        for t in time:
            
            for h in HP+incinerator+boiler:                
                self.variables.Q[h,t] = m.addVar(lb=0,ub=qmax[h],name='Q({0},{1})'.format(h,t))
            
            for h in CHP:
                for x in range(2):                
                    self.variables.Q[h,t,x] = m.addVar(lb=0,ub=q_bid[h][x],name='Q({0},{1},{2})'.format(h,t,x))

        self.variables.Q_charge = {} #discharging heat storage
        self.variables.Q_discharge = {} #charging heat storage
        
        for t in time:
            for h in heat_storage:  
                
                self.variables.Q_charge[h,t] = m.addVar(lb=0,ub=qmax[h],name='Q charge({0},{1})'.format(h,t))
                self.variables.Q_discharge[h,t] = m.addVar(lb=0,ub=qmax[h],name='Q discharge({0},{1})'.format(h,t))

        self.variables.energy_stored = {} #charging heat storage
        
        for h in heat_storage:
            for t in time:
                              
                self.variables.energy_stored[h,t] = m.addVar(lb=0,ub=emax[h],name='energy stored({0},{1})'.format(h,t))
                

            for t in time_init:
                              
                self.variables.energy_stored[h,t] = m.addVar(lb=energy_stored_init['HUC',h,D_day,t],ub=energy_stored_init['HUC',h,D_day,t],name='energy stored({0},{1})'.format(h,t))
                

        self.variables.HS = {} #heat load shedding

        for n1 in node:              
            for t in time:
                self.variables.HS[n1,t] = m.addVar(lb=0,ub=heat_load[n1,hiy[D_day,t]],name='heat load shedding({0},{1})'.format(n1,t))

                
        self.variables.HF_in = {} # heat flow entering pipeline (node n1) at time t

        for n1 in node:
            for n2 in node_node[n1]:               
                for t in time:
                    self.variables.HF_in[n1,n2,t] = m.addVar(lb=0,ub=f_upper[n1,n2],name='heat flow in({0},{1},{2})'.format(n1,n2,t))

                for t in time_init:
                    self.variables.HF_in[n1,n2,t] = m.addVar(lb=HF_in_init['HUC',n1,n2,D_day,t],ub=HF_in_init['HUC',n1,n2,D_day,t],name='heat flow in({0},{1},{2})'.format(n1,n2,t))


        self.variables.HF_in_delay = {} # heat entering pipeline (node n1) at time t-x AND exiting (node n2) at time t

        for n1 in node:
            for n2 in node_node[n1]:               
                for t in time:
                    for x in range(time_delay_max+1):
                        self.variables.HF_in_delay[n1,n2,t,x] = m.addVar(lb=0,ub=f_upper[n1,n2],name='heat flow out({0},{1},{2},{3})'.format(n1,n2,t,x))
                


        self.variables.HF_out = {} # heat flow exiting pipeline (node n2) at time t

        for n1 in node:
            for n2 in node_node[n1]:               
                for t in time:
                    self.variables.HF_out[n1,n2,t] = m.addVar(lb=0,ub=f_upper[n1,n2],name='heat flow out({0},{1},{2})'.format(n1,n2,t))
                
                
        self.variables.pmin_bid = {} # ADJUSTED BIDS of CHPs and incinerators

        self.variables.pmax_bid = {} # ADJUSTED BIDS of CHPs and incinerators
        
        for g in CHP:
            for t in time:
                self.variables.pmin_bid[g,t] = m.addVar(lb=0,name='elec min bid ({0},{1})'.format(g,t))
                self.variables.pmax_bid[g,t] = m.addVar(lb=0,name='elec max bid({0},{1})'.format(g,t))

        m.update()
    
    def _build_objective(self,D_day): # building the objective function for the HUC
  
        m = self.model    
              

        m.setObjective(    - (sum(c_heat_bid['HUC',g,D_day,t] for g in HP+incinerator+boiler for t in time) + sum(c_heat_bid['HUC',g,D_day,t][0] for g in CHP for t in time) )/(len(time)*len(HP+CHP+CHP+incinerator+boiler))*gb.quicksum( self.variables.HF_in[n1,n2,t] - self.variables.HF_out[n1,n2,t] for n1 in node for n2 in node_node[n1] for t in time)    - (sum(c_heat_bid['HUC',g,D_day,t] for g in HP+incinerator+boiler for t in time) + sum(c_heat_bid['HUC',g,D_day,t][0] for g in CHP for t in time) )/(len(time)*len(HP+CHP+CHP+incinerator+boiler))*gb.quicksum(self.variables.energy_stored[g,time[-1]] for g in heat_storage)   + gb.quicksum(500*self.variables.HS[n,t] for n in node for t in time) + gb.quicksum(c_heat_bid['HUC',g,D_day,t][x]*self.variables.Q[g,t,x] for t in time for g in CHP for x in range(2)) + gb.quicksum(c_heat_bid['HUC',g,D_day,t]*self.variables.Q[g,t] for t in time for g in incinerator+HP+boiler) + gb.quicksum(self.variables.r[h,t] + c_0[h]*self.variables.u[h,t] for t in time for h in CHP+incinerator+HP+boiler) ,   
            gb.GRB.MINIMIZE)


        
    def _build_constraints(self,D_day):
     
        m = self.model
        
        # UL : UC, accepted (valid bids), time delays and flow direction in pipelines
                
        self.constraints.start_cost = {}

        for g in CHP+boiler+incinerator+HP:        
            for x in range(len(time)):
                for h in time_start_range[g]:
                    self.constraints.start_cost[g,time[x],h] = m.addConstr(
                            self.variables.r[g,time[x]],
                            gb.GRB.GREATER_EQUAL,
                            c_start[g,h]*(self.variables.u[g,time[x]] - gb.quicksum(self.variables.u[g,time_extended[len(time_init)+x-k-1]] for k in range(h))),name='start up cost compute({0},{1},{2})'.format(g,time[x],h))


        self.constraints.time_on_min = {}
        self.constraints.time_off_min = {}

        for g in CHP+boiler+incinerator+HP:
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
        
        for g in CHP+boiler+incinerator+HP:
            for x in range(len(time)):
                self.constraints.uc[g,time[x]] = m.addConstr(
                        self.variables.u[g,time_extended[len(time_init)+x]]-self.variables.u[g,time_extended[len(time_init)+x-1]],
                        gb.GRB.EQUAL,      
                        self.variables.v_on[g,time_extended[len(time_init)+x]]-self.variables.v_off[g,time_extended[len(time_init)+x]],name='uc binary variables ({0},{1})'.format(g,time[x]))


        self.constraints.bid_accepted = {}
        
        for g in CHP:
            for t in time:
                self.constraints.bid_accepted[g,t,0] = m.addConstr(
                        self.variables.bid_accepted[g,t,0],
                        gb.GRB.LESS_EQUAL,      
                        self.variables.u[g,t],name='bid accepted'.format(g,t,0))

                self.constraints.bid_accepted[g,t,1] = m.addConstr(
                        self.variables.bid_accepted[g,t,1],
                        gb.GRB.LESS_EQUAL,      
                        self.variables.bid_accepted[g,t,0],name='bid accepted'.format(g,t,1))



        # flow direction

        self.constraints.flow_direction = {}  

        for t in time:

            for n1 in node:
                for n2 in node_node[n1]:

                    self.constraints.flow_direction[n1,n2,t] = m.addConstr(
                        self.variables.HF_direction[n2,n1,t]+self.variables.HF_direction[n1,n2,t],
                        gb.GRB.EQUAL,
                        1)

        # time delays
        
        self.constraints.time_delay_sum = {}
        
        self.constraints.time_delay_direction = {}
        
        for t in time:
            for n1 in node:
                for n2 in node_node[n1]:
                        
                    self.constraints.time_delay_sum[n1,n2,t] = m.addConstr(
                        gb.quicksum( self.variables.time_delay[n1,n2,t,x] for x in range(time_delay_max+1) ),
                        gb.GRB.EQUAL,
                        1)
                        
                    self.constraints.time_delay_direction[n1,n2,t] = m.addConstr(
                        self.variables.time_delay[n1,n2,t,0],
                        gb.GRB.GREATER_EQUAL,
                        1-self.variables.HF_direction[n1,n2,t])

        # time delays bounds: first in first out
        
        self.constraints.first_in_first_out = {}
        
        for t in range(len(time)):
            for x in range(1,time_delay_max+1):
                for n1 in node:
                    for n2 in node_node[n1]:                            
                        self.constraints.first_in_first_out[n1,n2,time[t],time_extended[len(time_init)+t-x]] = m.addConstr(
                            t + gb.quicksum(k*self.variables.time_delay[n1,n2,time[t],k] for k in range(time_delay_max+1)),
                            gb.GRB.GREATER_EQUAL,
                            t - x + gb.quicksum(k*self.variables.time_delay[n1,n2,time_extended[len(time_init)+t-x],k] for k in range(time_delay_max+1)))

        
        ## ML: heat market constraints
        
        # min/max production depending on accepted bids and on/off status
        
        self.constraints.heat_minprod = {} #c1a                
        self.constraints.heat_maxprod = {} # c1b 
        
        for t in time:

            for g in boiler+incinerator+HP:
                               
                self.constraints.heat_minprod[g,t] = m.addConstr(
                    self.variables.Q[g,t],
                    gb.GRB.GREATER_EQUAL,
                    qmin[g]*self.variables.u[g,t],name='elec minprod({0},{1})'.format(g,t))

                self.constraints.heat_maxprod[g,t] = m.addConstr(
                    self.variables.Q[g,t],
                    gb.GRB.LESS_EQUAL,
                    qmax[g]*self.variables.u[g,t],name='elec maxprod({0},{1})'.format(g,t))
 
    
            for g in CHP:

                self.constraints.heat_minprod[g,t,0] = m.addConstr(
                    self.variables.Q[g,t,0],
                    gb.GRB.GREATER_EQUAL,
                    q_bid[g][0]*self.variables.bid_accepted[g,t,1],name='elec minprod({0},{1},{2})'.format(g,t,0))
                
                for x in range(2):

                    self.constraints.heat_maxprod[g,t,x] = m.addConstr( 
                        self.variables.Q[g,t,x],
                        gb.GRB.LESS_EQUAL,
                        q_bid[g][x]*self.variables.bid_accepted[g,t,x],name='elec maxprod({0},{1},{2})'.format(g,t,x))

        # c2a and c2b: bounds on load shedding
        
        # flow bounds based on flow direction
                  
        self.constraints.flow_direction_in = {} # c3 (implicitly lower bounded in c5 and c6)
        self.constraints.flow_direction_out = {} # c4 (implicitly lower bounded in c7)
        
        for t in time:
            for n1 in node:
                for n2 in node_node[n1]:
                    
                    self.constraints.flow_direction_in[n1,n2,t] = m.addConstr(
                        self.variables.HF_in[n1,n2,t],
                        gb.GRB.LESS_EQUAL,
                        f_upper[n1,n2]*self.variables.HF_direction[n1,n2,t])
                    
                    # I THINK THIS MAY BE REDUNDANT?        
                    self.constraints.flow_direction_out[n1,n2,t] = m.addConstr(
                        self.variables.HF_out[n1,n2,t],
                        gb.GRB.LESS_EQUAL,
                        f_upper[n1,n2]*self.variables.HF_direction[n1,n2,t]) 
        
        # in/out flows based on time delays


        self.constraints.HF_delay_upper_1 = {} #c5b
        self.constraints.HF_delay_lower_1 = {} #c5a

        self.constraints.HF_delay_upper_2 = {} #c6b
        self.constraints.HF_delay_lower_2 = {} #c6a
        # these constraints also provide upper and lower bounds for HF_in and HF_in_delay
        
        for t in range(len(time)):
            for x in range(time_delay_max+1):
                for n1 in node:
                    for n2 in node_node[n1]:
                            
                        self.constraints.HF_delay_upper_1[n1,n2,time[t],x] = m.addConstr(
                            - self.variables.HF_in_delay[n1,n2,time[t],x] + self.variables.HF_in[n1,n2,time_extended[len(time_init)+t-x]],
                            gb.GRB.LESS_EQUAL,
                            f_upper[n1,n2]*(1-self.variables.time_delay[n1,n2,time_extended[len(time_init)+t-x],x]))
                        
                        self.constraints.HF_delay_lower_1[n1,n2,time[t],x] = m.addConstr(
                            - self.variables.HF_in_delay[n1,n2,time[t],x] + self.variables.HF_in[n1,n2,time_extended[len(time_init)+t-x]],
                            gb.GRB.GREATER_EQUAL,
                            0)          
    
                            
                        self.constraints.HF_delay_upper_2[n1,n2,time[t],x] = m.addConstr(
                            self.variables.HF_in_delay[n1,n2,time[t],x],
                            gb.GRB.LESS_EQUAL,
                            f_upper[n1,n2]*self.variables.time_delay[n1,n2,time_extended[len(time_init)+t-x],x])
                        
                        self.constraints.HF_delay_lower_2[n1,n2,time[t],x] = m.addConstr(
                            self.variables.HF_in_delay[n1,n2,time[t],x],
                            gb.GRB.GREATER_EQUAL,
                            0) 
                        
        self.constraints.HF_out = {} #c7
        
        for t in time:
            for n1 in node:
                for n2 in node_node[n1]:
                    
                    self.constraints.HF_out[n1,n2,t] = m.addConstr(
                        self.variables.HF_out[n1,n2,t],
                        gb.GRB.EQUAL,
                        gb.quicksum(self.variables.HF_in_delay[n1,n2,t,x] for x in range(time_delay_max+1)) )
                        
        # heat balance at each node

        self.constraints.heat_balance = {} # c8
        
        for t in time:
            for n1 in node:
                        
                self.constraints.heat_balance[n1,t] = m.addConstr(
                    gb.quicksum(self.variables.Q_discharge[g,t] - self.variables.Q_charge[g,t] for g in heat_storage_node[n1]) + gb.quicksum(self.variables.Q[g,t] for g in prod_node[n1]) + gb.quicksum(self.variables.Q[g,t,x] for g in CHP_node[n1] for x in range(2)) + gb.quicksum( self.variables.HF_out[n2,n1,t] - self.variables.HF_in[n1,n2,t] for n2 in node_node[n1]),
                    gb.GRB.EQUAL,
                    heat_load[n1,hiy[D_day,t]]-self.variables.HS[n1,t])        
                        

        # heat storage constraints
        
        self.constraints.heat_storage_update = {} # c9
         
        for t in range(len(time)):
            for g in heat_storage:
                        
                self.constraints.heat_storage_update[g,t] = m.addConstr(
                    self.variables.energy_stored[g,time_extended[len(time_init)+t]],
                    gb.GRB.EQUAL,
                    self.variables.energy_stored[g,time_extended[len(time_init)+t-1]] - rho_discharge[g]*self.variables.Q_discharge[g,time_extended[len(time_init)+t]] + rho_charge[g]*self.variables.Q_charge[g,time_extended[len(time_init)+t]])           

        # c10 = energy bounds , c11 = charge bounds, c12 = discharge bounds
        
        # adjsuted bids of CHPs for electricity market

#        self.constraints.pmin_bid_lower_1 = {} #c13a        
#        self.constraints.pmin_bid_upper_1 = {} #c13b
#        self.constraints.pmin_bid_upper_2 = {} #c14
#        self.constraints.pmin_bid_upper_3 = {} #c15
#
#
#        self.constraints.pmax_bid_lower_1 = {} #c16a        
#        self.constraints.pmax_bid_upper_1 = {} #c16b
#        self.constraints.pmax_bid_lower_2 = {}#c17a 
#        self.constraints.pmax_bid_upper_2 = {} #c17b
#        
#        for t in time:
#            for g in CHP:
#
#
#                self.constraints.pmin_bid_lower_1[g,t] = m.addConstr(
#                    self.variables.pmin_bid[g,t],
#                    gb.GRB.GREATER_EQUAL,
#                    0,name='pmin bid lower 1({0},{1})'.format(g,t))
#
#                self.constraints.pmin_bid_upper_1[g,t] = m.addConstr(
#                    self.variables.pmin_bid[g,t],
#                    gb.GRB.LESS_EQUAL,
#                    fmax[g]/rho_elec[g]*self.variables.u[g,t],name='pmin bid upper 1({0},{1})'.format(g,t))
#
#
#
#
#                                
#                self.constraints.pmin_bid_upper_2[g,t] = m.addConstr(
#                    - self.variables.pmin_bid[g,t] + (fmin[g]-rho_heat[g]*gb.quicksum(self.variables.Q[g,t,x] for x in range(2)))/rho_elec[g],
#                    gb.GRB.LESS_EQUAL,
#                    fmax[g]/rho_elec[g]*(1-self.variables.u[g,t]),name='pmin bid upper 2({0},{1})'.format(g,t))
#
#
#
#
#
#                self.constraints.pmin_bid_upper_3[g,t] = m.addConstr(
#                    - self.variables.pmin_bid[g,t] + rmin[g]*gb.quicksum(self.variables.Q[g,t,x] for x in range(2)),
#                    gb.GRB.LESS_EQUAL,
#                    fmax[g]/rho_elec[g]*(1-self.variables.u[g,t]),name='pmin bid upper 3({0},{1})'.format(g,t))                
#
#
#
#                
#                self.constraints.pmax_bid_lower_1[g,t] = m.addConstr(
#                    - self.variables.pmax_bid[g,t] + (fmax[g]-rho_heat[g]*gb.quicksum(self.variables.Q[g,t,x] for x in range(2)))/rho_elec[g],
#                    gb.GRB.GREATER_EQUAL,
#                    0,name='pmax bid lower 1({0},{1})'.format(g,t))
#                
#                self.constraints.pmax_bid_upper_1[g,t] = m.addConstr(
#                    - self.variables.pmax_bid[g,t] + (fmax[g]-rho_heat[g]*gb.quicksum(self.variables.Q[g,t,x] for x in range(2)))/rho_elec[g],
#                    gb.GRB.LESS_EQUAL,
#                    fmax[g]/rho_elec[g]*(1-self.variables.u[g,t]),name='pmax bid upper 1({0},{1})'.format(g,t))
#
#
#
#
#                self.constraints.pmax_bid_lower_2[g,t] = m.addConstr(
#                    self.variables.pmax_bid[g,t],
#                    gb.GRB.GREATER_EQUAL,
#                    0,name='pmax bid lower 2({0},{1})'.format(g,t))
#                
#                self.constraints.pmax_bid_upper_2[g,t] = m.addConstr(
#                    self.variables.pmax_bid[g,t],
#                    gb.GRB.LESS_EQUAL,
#                    fmax[g]/rho_elec[g]*self.variables.u[g,t],name='pmax bid upper 2({0},{1})'.format(g,t))


#%% combined
        
class CHPUC:
    
    def __init__(self,D_day):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._build_model(D_day)

    
    def optimize(self):
        self.model.optimize()
        
    def _build_model(self,D_day):
        
        self.model = gb.Model()
        self.model.Params.TimeLimit = 60
        self._build_variables(D_day)
        self._build_objective(D_day)
        self._build_constraints(D_day)
    
    def _build_variables(self,D_day):
        m = self.model


        ## UL: binary variables and start up cost
        
        self.variables.u = {} # status on/off
        
        for h in CHP+boiler+incinerator+HP+gen:
           for t in time:                            
                self.variables.u[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='u({0},{1})'.format(h,t))
           for t in time_init:                       
                self.variables.u[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=u_init['CHPUC',h,D_day,t],ub=u_init['CHPUC',h,D_day,t],name='u({0},{1})'.format(h,t))

        self.variables.v_on = {} # status turn on
        
        for h in CHP+boiler+incinerator+HP+gen:
           for t in time:                            
                self.variables.v_on[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='v on({0},{1})'.format(h,t))
           for t in time_init:                           
                self.variables.v_on[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=v_on_init['CHPUC',h,D_day,t],ub=v_on_init['CHPUC',h,D_day,t],name='v on({0},{1})'.format(h,t))


        self.variables.v_off = {} # status turn off
        
        for h in CHP+boiler+incinerator+HP+gen:
           for t in time:                           
                self.variables.v_off[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='v off({0},{1})'.format(h,t))
           for t in time_init:                            
                self.variables.v_off[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=v_off_init['CHPUC',h,D_day,t],ub=v_off_init['CHPUC',h,D_day,t],name='v off({0},{1})'.format(h,t))


        self.variables.r = {} # start up cost
        for t in time:
            for h in CHP+boiler+incinerator+HP+gen:
                self.variables.r[h,t] = m.addVar(lb=0,name='start up cost({0},{1})'.format(h,t)) # dispatch of electricity generators


        self.variables.time_delay = {} # time delay auxiliary variables self.variables.time_delay[n1,n2,t,x] = 1 iif time delay of heat flow entering pipeline (node n1) at time t, exits pipeline (node n2) at time t+x (time_delay[n1,n2,t]=sum_{x=0}^{time_delay_max} x*self.variables.time_delay[n1,n2,t,x] )
        
        for n1 in node:
            for n2 in node_node[n1]:
                    for t in time:
                        for x in range(time_delay_max+1):
                            self.variables.time_delay[n1,n2,t,x] = m.addVar(vtype = gb.GRB.BINARY,name='time delay auxiliary variables({0},{1},{2},{3}])'.format(n1,n2,t,x)) # dispatch of electricity generators

                    for t in time_init:
                        for x in range(time_delay_max+1):
                            self.variables.time_delay[n1,n2,t,x] = m.addVar(vtype = gb.GRB.BINARY,lb=time_delay_init['CHPUC',n1,n2,D_day,t,x],ub=time_delay_init['CHPUC',n1,n2,D_day,t,x],name='time delay auxiliary variables({0},{1},{2},{3}])'.format(n1,n2,t,x)) # dispatch of electricity generators

    
        self.variables.HF_direction = {} # direction of flow
        for n1 in node:
            for n2 in node_node[n1]:
                for t in time:
                    self.variables.HF_direction[n1,n2,t] = m.addVar(vtype = gb.GRB.BINARY,name='HF direction({0},{1},{2})'.format(n1,n2,t)) # dispatch of electricity generators

                for t in time_init:
                    self.variables.HF_direction[n1,n2,t] = m.addVar(vtype = gb.GRB.BINARY,lb=HF_direction_init['CHPUC',n1,n2,D_day,t],ub=HF_direction_init['CHPUC',n1,n2,D_day,t],name='HF direction({0},{1},{2})'.format(n1,n2,t)) # dispatch of electricity generators

        # ML: HM variables
        
        self.variables.Q = {} #elec production
        
        for t in time:
            
            for h in HP+incinerator+boiler+CHP:                
                self.variables.Q[h,t] = m.addVar(lb=0,ub=qmax[h],name='Q({0},{1})'.format(h,t))


        self.variables.Q_charge = {} #discharging heat storage
        self.variables.Q_discharge = {} #charging heat storage
        
        for t in time:
            for h in heat_storage:  
                
                self.variables.Q_charge[h,t] = m.addVar(lb=0,ub=qmax[h],name='Q charge({0},{1})'.format(h,t))
                self.variables.Q_discharge[h,t] = m.addVar(lb=0,ub=qmax[h],name='Q discharge({0},{1})'.format(h,t))

        self.variables.energy_stored = {} #charging heat storage
        
        for h in heat_storage:
            for t in time:
                              
                self.variables.energy_stored[h,t] = m.addVar(lb=0,ub=emax[h],name='energy stored({0},{1})'.format(h,t))
                

            for t in time_init:
                              
                self.variables.energy_stored[h,t] = m.addVar(lb=energy_stored_init['CHPUC',h,D_day,t],ub=energy_stored_init['CHPUC',h,D_day,t],name='energy stored({0},{1})'.format(h,t))
                

        self.variables.HS = {} #heat load shedding

        for n1 in node:              
            for t in time:
                self.variables.HS[n1,t] = m.addVar(lb=0,ub=heat_load[n1,hiy[D_day,t]],name='heat load shedding({0},{1})'.format(n1,t))

                
        self.variables.HF_in = {} # heat flow entering pipeline (node n1) at time t

        for n1 in node:
            for n2 in node_node[n1]:               
                for t in time:
                    self.variables.HF_in[n1,n2,t] = m.addVar(lb=0,ub=f_upper[n1,n2],name='heat flow in({0},{1},{2})'.format(n1,n2,t))

                for t in time_init:
                    self.variables.HF_in[n1,n2,t] = m.addVar(lb=HF_in_init['CHPUC',n1,n2,D_day,t],ub=HF_in_init['CHPUC',n1,n2,D_day,t],name='heat flow in({0},{1},{2})'.format(n1,n2,t))


        self.variables.HF_in_delay = {} # heat entering pipeline (node n1) at time t-x AND exiting (node n2) at time t

        for n1 in node:
            for n2 in node_node[n1]:               
                for t in time:
                    for x in range(time_delay_max+1):
                        self.variables.HF_in_delay[n1,n2,t,x] = m.addVar(lb=0,ub=f_upper[n1,n2],name='heat flow out({0},{1},{2},{3})'.format(n1,n2,t,x))
                


        self.variables.HF_out = {} # heat flow exiting pipeline (node n2) at time t

        for n1 in node:
            for n2 in node_node[n1]:               
                for t in time:
                    self.variables.HF_out[n1,n2,t] = m.addVar(lb=0,ub=f_upper[n1,n2],name='heat flow out({0},{1},{2})'.format(n1,n2,t))

        # LL: EM variables
        
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
    
    def _build_objective(self,D_day): # building the objective function for the EAHUC
  
        m = self.model    
              

        m.setObjective( - sum(c_heat[g] for g in HP+incinerator+boiler+CHP)/len(HP+CHP+incinerator+boiler)*gb.quicksum( self.variables.HF_in[n1,n2,t] - self.variables.HF_out[n1,n2,t] for n1 in node for n2 in node_node[n1] for t in time)    - - sum(c_heat[g] for g in HP+incinerator+boiler+CHP)/len(HP+CHP+incinerator+boiler)*gb.quicksum(self.variables.energy_stored[g,time[-1]] for g in heat_storage) + gb.quicksum(500*self.variables.HS[n,t] for n in node for t in time) + gb.quicksum(c_heat[g]*self.variables.Q[g,t] for t in time for g in CHP+HP+incinerator+boiler) + gb.quicksum(self.variables.r[h,t] + c_0[h]*self.variables.u[h,t] for t in time for h in CHP+incinerator+HP+boiler+gen) + gb.quicksum(c_elec[g]*self.variables.P[g,t] for t in time for g in CHP+gen) + gb.quicksum(500*self.variables.ES[n,t] for t in time for n in zone_DK) ,   
            gb.GRB.MINIMIZE)


        
    def _build_constraints(self,D_day):
     
        m = self.model
        
        # UL : UC, accepted (valid bids), time delays and flow direction in pipelines
                
        self.constraints.start_cost = {}

        for g in CHP+boiler+incinerator+HP+gen:        
            for x in range(len(time)):
                for h in time_start_range[g]:
                    self.constraints.start_cost[g,time[x],h] = m.addConstr(
                            self.variables.r[g,time[x]],
                            gb.GRB.GREATER_EQUAL,
                            c_start[g,h]*(self.variables.u[g,time[x]] - gb.quicksum(self.variables.u[g,time_extended[len(time_init)+x-k-1]] for k in range(h))),name='start up cost compute({0},{1},{2})'.format(g,time[x],h))


        self.constraints.time_on_min = {}
        self.constraints.time_off_min = {}

        for g in CHP+boiler+incinerator+HP+gen:
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
        
        for g in CHP+boiler+incinerator+HP+gen:
            for x in range(len(time)):
                self.constraints.uc[g,time[x]] = m.addConstr(
                        self.variables.u[g,time_extended[len(time_init)+x]]-self.variables.u[g,time_extended[len(time_init)+x-1]],
                        gb.GRB.EQUAL,      
                        self.variables.v_on[g,time_extended[len(time_init)+x]]-self.variables.v_off[g,time_extended[len(time_init)+x]],name='uc binary variables ({0},{1})'.format(g,time[x]))

        # flow direction

        self.constraints.flow_direction = {}

        for t in time:

            for n1 in node:
                for n2 in node_node[n1]:

                    self.constraints.flow_direction[n1,n2,t] = m.addConstr(
                        self.variables.HF_direction[n2,n1,t]+self.variables.HF_direction[n1,n2,t],
                        gb.GRB.EQUAL,
                        1)

        # time delays
        
        self.constraints.time_delay_sum = {}
        
        self.constraints.time_delay_direction = {}
        
        for t in time:
            for n1 in node:
                for n2 in node_node[n1]:
                        
                    self.constraints.time_delay_sum[n1,n2,t] = m.addConstr(
                        gb.quicksum( self.variables.time_delay[n1,n2,t,x] for x in range(time_delay_max+1) ),
                        gb.GRB.EQUAL,
                        1)
                        
                    self.constraints.time_delay_direction[n1,n2,t] = m.addConstr(
                        self.variables.time_delay[n1,n2,t,0],
                        gb.GRB.GREATER_EQUAL,
                        1-self.variables.HF_direction[n1,n2,t])

        # time delays bounds: first in first out
        
        self.constraints.first_in_first_out = {}
        
        for t in range(len(time)):
            for x in range(1,time_delay_max+1):
                for n1 in node:
                    for n2 in node_node[n1]:                            
                        self.constraints.first_in_first_out[n1,n2,time[t],time_extended[len(time_init)+t-x]] = m.addConstr(
                            t + gb.quicksum(k*self.variables.time_delay[n1,n2,time[t],k] for k in range(time_delay_max+1)),
                            gb.GRB.GREATER_EQUAL,
                            t - x + gb.quicksum(k*self.variables.time_delay[n1,n2,time_extended[len(time_init)+t-x],k] for k in range(time_delay_max+1)))
     
        
        ## ML: heat market constraints
        
        # min/max production depending on accepted bids and on/off status
        
        self.constraints.heat_minprod = {} #c1a                
        self.constraints.heat_maxprod = {} # c1b 
        
        for t in time:

            for g in boiler+incinerator+HP+CHP:
                               
                self.constraints.heat_minprod[g,t] = m.addConstr(
                    self.variables.Q[g,t],
                    gb.GRB.GREATER_EQUAL,
                    qmin[g]*self.variables.u[g,t],name='elec minprod({0},{1})'.format(g,t))

                self.constraints.heat_maxprod[g,t] = m.addConstr(
                    self.variables.Q[g,t],
                    gb.GRB.LESS_EQUAL,
                    qmax[g]*self.variables.u[g,t],name='elec maxprod({0},{1})'.format(g,t))

        # c2a and c2b: bounds on load shedding
        
        # flow bounds based on flow direction
                  
        self.constraints.flow_direction_in = {} # c3 (implicitly lower bounded in c5 and c6)
        self.constraints.flow_direction_out = {} # c4 (implicitly lower bounded in c7)
        
        for t in time:
            for n1 in node:
                for n2 in node_node[n1]:
                    
                    self.constraints.flow_direction_in[n1,n2,t] = m.addConstr(
                        self.variables.HF_in[n1,n2,t],
                        gb.GRB.LESS_EQUAL,
                        f_upper[n1,n2]*self.variables.HF_direction[n1,n2,t])
                    
                    # I THINK THIS MAY BE REDUNDANT?        
                    self.constraints.flow_direction_out[n1,n2,t] = m.addConstr(
                        self.variables.HF_out[n1,n2,t],
                        gb.GRB.LESS_EQUAL,
                        f_upper[n1,n2]*self.variables.HF_direction[n1,n2,t]) 
        
        # in/out flows based on time delays


        self.constraints.HF_delay_upper_1 = {} #c5b
        self.constraints.HF_delay_lower_1 = {} #c5a

        self.constraints.HF_delay_upper_2 = {} #c6b
        self.constraints.HF_delay_lower_2 = {} #c6a
        # these constraints also provide upper and lower bounds for HF_in and HF_in_delay
        
        for t in range(len(time)):
            for x in range(time_delay_max+1):
                for n1 in node:
                    for n2 in node_node[n1]:
                            
                        self.constraints.HF_delay_upper_1[n1,n2,time[t],x] = m.addConstr(
                            - self.variables.HF_in_delay[n1,n2,time[t],x] + self.variables.HF_in[n1,n2,time_extended[len(time_init)+t-x]],
                            gb.GRB.LESS_EQUAL,
                            f_upper[n1,n2]*(1-self.variables.time_delay[n1,n2,time_extended[len(time_init)+t-x],x]))
                        
                        self.constraints.HF_delay_lower_1[n1,n2,time[t],x] = m.addConstr(
                            - self.variables.HF_in_delay[n1,n2,time[t],x] + self.variables.HF_in[n1,n2,time_extended[len(time_init)+t-x]],
                            gb.GRB.GREATER_EQUAL,
                            0)          
    
                            
                        self.constraints.HF_delay_upper_2[n1,n2,time[t],x] = m.addConstr(
                            self.variables.HF_in_delay[n1,n2,time[t],x],
                            gb.GRB.LESS_EQUAL,
                            f_upper[n1,n2]*self.variables.time_delay[n1,n2,time_extended[len(time_init)+t-x],x])
                        
                        self.constraints.HF_delay_lower_2[n1,n2,time[t],x] = m.addConstr(
                            self.variables.HF_in_delay[n1,n2,time[t],x],
                            gb.GRB.GREATER_EQUAL,
                            0) 
                        
        self.constraints.HF_out = {} #c7
        
        for t in time:
            for n1 in node:
                for n2 in node_node[n1]:
                    
                    self.constraints.HF_out[n1,n2,t] = m.addConstr(
                        self.variables.HF_out[n1,n2,t],
                        gb.GRB.EQUAL,
                        gb.quicksum(self.variables.HF_in_delay[n1,n2,t,x] for x in range(time_delay_max+1)) )
                        
        # heat balance at each node

        self.constraints.heat_balance = {} # c8
        
        for t in time:
            for n1 in node:
                        
                self.constraints.heat_balance[n1,t] = m.addConstr(
                    gb.quicksum(self.variables.Q_discharge[g,t] - self.variables.Q_charge[g,t] for g in heat_storage_node[n1]) + gb.quicksum(self.variables.Q[g,t] for g in prod_node[n1]+CHP_node[n1]) + gb.quicksum(self.variables.HF_out[n2,n1,t] - self.variables.HF_in[n1,n2,t] for n2 in node_node[n1]),
                    gb.GRB.EQUAL,
                    heat_load[n1,hiy[D_day,t]]-self.variables.HS[n1,t])        
                        

        # heat storage constraints
        
        self.constraints.heat_storage_update = {} # c9
         
        for t in range(len(time)):
            for g in heat_storage:
                        
                self.constraints.heat_storage_update[g,t] = m.addConstr(
                    self.variables.energy_stored[g,time_extended[len(time_init)+t]],
                    gb.GRB.EQUAL,
                    self.variables.energy_stored[g,time_extended[len(time_init)+t-1]] - rho_discharge[g]*self.variables.Q_discharge[g,time_extended[len(time_init)+t]] + rho_charge[g]*self.variables.Q_charge[g,time_extended[len(time_init)+t]] - heat_loss[g]/10)           

        # c10 = energy bounds , c11 = charge bounds, c12 = discharge bounds
        
        # adjsuted bids of CHPs for electricity market


        self.constraints.f_upper = {} #c14
        self.constraints.f_lower = {} #c15
        self.constraints.rmin = {} #c16

        
        for t in time:
            for g in CHP:

                                
                self.constraints.f_lower[g,t] = m.addConstr(
                    rho_elec[g]*self.variables.P[g,t] + rho_heat[g]*self.variables.Q[g,t],
                    gb.GRB.GREATER_EQUAL,
                    fmin[g]*self.variables.u[g,t],name='f lower ({0},{1})'.format(g,t))
                
                self.constraints.f_upper[g,t] = m.addConstr(
                    rho_elec[g]*self.variables.P[g,t] + rho_heat[g]*self.variables.Q[g,t],
                    gb.GRB.LESS_EQUAL,
                    fmax[g]*self.variables.u[g,t],name='f upper ({0},{1})'.format(g,t))
                
                self.constraints.rmin[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.GREATER_EQUAL,
                    rmin[g]*self.variables.Q[g,t],name='rmin ({0},{1})'.format(g,t))
                
        ## LL: electricity market constraints

                
        #4) ELEC GEN MAX PROD and total production variable

        self.constraints.elec_minprod = {} # c18a        
        self.constraints.elec_maxprod = {} # c18b

        
        for t in time:

            for g in gen:

                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.GREATER_EQUAL,
                    pmin[g]*self.variables.u[g,t],name='elec minprod({0},{1})'.format(g,t))
                
                self.constraints.elec_maxprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    pmax[g]*self.variables.u[g,t],name='elec maxprod({0},{1})'.format(g,t))


        # c19a and c19b = bounds on elec load shedding
        
        #c20a and c20b : bounds on solar prod
        
        #c21a and c21b : bounds on wind prod
        
        # c22a and c22b = bounds on elec flows      
        
        # elec balance
        
        self.constraints.elec_balance = {} # c23
        
        for t in time:
            for z in zone_DK:
                
                self.constraints.elec_balance[z,t] = m.addConstr(
                    gb.quicksum(self.variables.P[g,t] for g in gen_zone[z]+CHP_zone[z]) + gb.quicksum(rmin[h]*self.variables.Q[h,t] for h in incinerator_zone[z]) + self.variables.wind_prod[z,t]+self.variables.solar_prod[z,t]+self.variables.ES[z,t],
                    gb.GRB.EQUAL,
                    elec_load[z,hiy[D_day,t]] + gb.quicksum(rmin[h]*self.variables.Q[h,t] for h in HP_zone[z]) + gb.quicksum(self.variables.EF[l,t] for l in line_start[z]) - gb.quicksum(self.variables.EF[l,t] for l in line_end[z]),name='elec balance({0},{1})'.format(z,t))       

        
#%% EAHUC

alpha = 0.99
M = 5000

class EAHUC:
    def __init__(self,D_day):
        self.data = expando()
        self.variables = expando()
        self.constraints = expando()
        self._build_model(D_day)

    
    def optimize(self):
        self.model.optimize()
        
    def _build_model(self,D_day):
        
        self.model = gb.Model()
        self.model.Params.TimeLimit = 60
        self._build_variables(D_day)
        self._build_objective(D_day)
        self._build_constraints(D_day)
    
    def _build_variables(self,D_day):
        m = self.model


        ## UL: binary variables and start up cost
        
        self.variables.u = {} # status on/off
        
        for h in CHP+boiler+incinerator+HP:
           for t in time:                            
                self.variables.u[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='u({0},{1})'.format(h,t))
           for t in time_init:                       
                self.variables.u[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=u_init['EAHUC',h,D_day,t],ub=u_init['EAHUC',h,D_day,t],name='u({0},{1})'.format(h,t))


        self.variables.bid_accepted = {} # status on/off
        
        for t in time: 
            for h in CHP:      
                for x in range(2):                     
                    self.variables.bid_accepted[h,t,x] = m.addVar(vtype=gb.GRB.BINARY,name='bid accepted({0},{1},{2})'.format(h,t,x))
    

        self.variables.v_on = {} # status turn on
        
        for h in CHP+boiler+incinerator+HP:
           for t in time:                            
                self.variables.v_on[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='v on({0},{1})'.format(h,t))
           for t in time_init:                           
                self.variables.v_on[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=v_on_init['EAHUC',h,D_day,t],ub=v_on_init['EAHUC',h,D_day,t],name='v on({0},{1})'.format(h,t))


        self.variables.v_off = {} # status turn off
        
        for h in CHP+boiler+incinerator+HP:
           for t in time:                           
                self.variables.v_off[h,t] = m.addVar(vtype=gb.GRB.BINARY,name='v off({0},{1})'.format(h,t))
           for t in time_init:                            
                self.variables.v_off[h,t] = m.addVar(vtype=gb.GRB.BINARY,lb=v_off_init['EAHUC',h,D_day,t],ub=v_off_init['EAHUC',h,D_day,t],name='v off({0},{1})'.format(h,t))


        self.variables.r = {} # start up cost
        for t in time:
            for h in CHP+boiler+incinerator+HP:
                self.variables.r[h,t] = m.addVar(lb=0,name='start up cost({0},{1})'.format(h,t)) # dispatch of electricity generators


        self.variables.time_delay = {} # time delay auxiliary variables self.variables.time_delay[n1,n2,t,x] = 1 iif time delay of heat flow entering pipeline (node n1) at time t, exits pipeline (node n2) at time t+x (time_delay[n1,n2,t]=sum_{x=0}^{time_delay_max} x*self.variables.time_delay[n1,n2,t,x] )
        
        for n1 in node:
            for n2 in node_node[n1]:
                    for t in time:
                        for x in range(time_delay_max+1):
                            self.variables.time_delay[n1,n2,t,x] = m.addVar(vtype = gb.GRB.BINARY,name='time delay auxiliary variables({0},{1},{2},{3}])'.format(n1,n2,t,x)) # dispatch of electricity generators

                    for t in time_init:
                        for x in range(time_delay_max+1):
                            self.variables.time_delay[n1,n2,t,x] = m.addVar(vtype = gb.GRB.BINARY,lb=time_delay_init['EAHUC',n1,n2,D_day,t,x],ub=time_delay_init['EAHUC',n1,n2,D_day,t,x],name='time delay auxiliary variables({0},{1},{2},{3}])'.format(n1,n2,t,x)) # dispatch of electricity generators

    
        self.variables.HF_direction = {} # direction of flow
        for n1 in node:
            for n2 in node_node[n1]:
                for t in time:
                    self.variables.HF_direction[n1,n2,t] = m.addVar(vtype = gb.GRB.BINARY,name='HF direction({0},{1},{2})'.format(n1,n2,t)) # dispatch of electricity generators

                for t in time_init:
                    self.variables.HF_direction[n1,n2,t] = m.addVar(vtype = gb.GRB.BINARY,lb=HF_direction_init['EAHUC',n1,n2,D_day,t],ub=HF_direction_init['EAHUC',n1,n2,D_day,t],name='HF direction({0},{1},{2})'.format(n1,n2,t)) # dispatch of electricity generators

        # ML: HM variables
        
        self.variables.Q = {} #elec production
        
        for t in time:
            
            for h in HP+incinerator+boiler:                
                self.variables.Q[h,t] = m.addVar(lb=0,ub=qmax[h],name='Q({0},{1})'.format(h,t))
            
            for h in CHP:
                for x in range(2):                
                    self.variables.Q[h,t,x] = m.addVar(lb=0,ub=q_bid[h][x],name='Q({0},{1},{2})'.format(h,t,x))

        self.variables.Q_charge = {} #discharging heat storage
        self.variables.Q_discharge = {} #charging heat storage
        
        for t in time:
            for h in heat_storage:  
                
                self.variables.Q_charge[h,t] = m.addVar(lb=0,ub=qmax[h],name='Q charge({0},{1})'.format(h,t))
                self.variables.Q_discharge[h,t] = m.addVar(lb=0,ub=qmax[h],name='Q discharge({0},{1})'.format(h,t))

        self.variables.energy_stored = {} #charging heat storage
        
        for h in heat_storage:
            for t in time:
                              
                self.variables.energy_stored[h,t] = m.addVar(lb=0,ub=emax[h],name='energy stored({0},{1})'.format(h,t))
                

            for t in time_init:
                              
                self.variables.energy_stored[h,t] = m.addVar(lb=energy_stored_init['EAHUC',h,D_day,t],ub=energy_stored_init['EAHUC',h,D_day,t],name='energy stored({0},{1})'.format(h,t))
                

        self.variables.HS = {} #heat load shedding

        for n1 in node:              
            for t in time:
                self.variables.HS[n1,t] = m.addVar(lb=0,ub=heat_load[n1,hiy[D_day,t]],name='heat load shedding({0},{1})'.format(n1,t))

                
        self.variables.HF_in = {} # heat flow entering pipeline (node n1) at time t

        for n1 in node:
            for n2 in node_node[n1]:               
                for t in time:
                    self.variables.HF_in[n1,n2,t] = m.addVar(lb=0,ub=f_upper[n1,n2],name='heat flow in({0},{1},{2})'.format(n1,n2,t))

                for t in time_init:
                    self.variables.HF_in[n1,n2,t] = m.addVar(lb=HF_in_init['EAHUC',n1,n2,D_day,t],ub=HF_in_init['EAHUC',n1,n2,D_day,t],name='heat flow in({0},{1},{2})'.format(n1,n2,t))


        self.variables.HF_in_delay = {} # heat entering pipeline (node n1) at time t-x AND exiting (node n2) at time t

        for n1 in node:
            for n2 in node_node[n1]:               
                for t in time:
                    for x in range(time_delay_max+1):
                        self.variables.HF_in_delay[n1,n2,t,x] = m.addVar(lb=0,ub=f_upper[n1,n2],name='heat flow out({0},{1},{2},{3})'.format(n1,n2,t,x))
                


        self.variables.HF_out = {} # heat flow exiting pipeline (node n2) at time t

        for n1 in node:
            for n2 in node_node[n1]:               
                for t in time:
                    self.variables.HF_out[n1,n2,t] = m.addVar(lb=0,ub=f_upper[n1,n2],name='heat flow out({0},{1},{2})'.format(n1,n2,t))
                
                
        self.variables.pmin_bid = {} # ADJUSTED BIDS of CHPs and incinerators

        self.variables.pmax_bid = {} # ADJUSTED BIDS of CHPs and incinerators
        
        for g in CHP:
            for t in time:
                self.variables.pmin_bid[g,t] = m.addVar(lb=0,name='elec min bid ({0},{1})'.format(g,t))
                self.variables.pmax_bid[g,t] = m.addVar(lb=0,name='elec max bid({0},{1})'.format(g,t))


        # LL: EM variables
        
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


        # dual variables 

        self.variables.mu_c1a = {} #heat min
        self.variables.mu_c1b = {} #heat max

        for t in time:
            for g in incinerator+HP+boiler:
                self.variables.mu_c1a[g,t] = m.addVar(lb=0)
                self.variables.mu_c1b[g,t] = m.addVar(lb=0) 
            for g in CHP:
                self.variables.mu_c1a[g,t] = m.addVar(lb=0)
                for x in range(2):                    
                    self.variables.mu_c1b[g,t,x] = m.addVar(lb=0)

        self.variables.mu_c2a = {} # load shedding min  
        self.variables.mu_c2b = {} # load shedding max

        for t in time:
            for n in node:
                self.variables.mu_c2a[n,t] = m.addVar(lb=0)
                self.variables.mu_c2b[n,t] = m.addVar(lb=0) 


                
        self.variables.mu_c3 = {} # HF in min / direction
        self.variables.mu_c4 = {} # HF out min / direction

        self.variables.mu_c5a = {} # HF in delay lower 1
        self.variables.mu_c5b = {} # HF in delay upper 1
        self.variables.mu_c6a = {} # HF in delay lower 2
        self.variables.mu_c6b = {} # HF in delay upper 2
        
        self.variables.lambda_c7 = {} # HF out equality
        
        for n1 in node:
            for n2 in node_node[n1]:
                for t in time:
                    
                    self.variables.mu_c3[n1,n2,t] = m.addVar(lb=0)
                    self.variables.mu_c4[n1,n2,t] = m.addVar(lb=0)

                    for x in range(time_delay_max+1):
                        self.variables.mu_c5a[n1,n2,t,x] = m.addVar(lb=0)
                        self.variables.mu_c5b[n1,n2,t,x] = m.addVar(lb=0)
                        self.variables.mu_c6a[n1,n2,t,x] = m.addVar(lb=0)
                        self.variables.mu_c6b[n1,n2,t,x] = m.addVar(lb=0)

                    self.variables.lambda_c7[n1,n2,t] = m.addVar(lb=0)

        self.variables.heat_price = {} # heat prices (balance equation) <=> lambda c8
        for t in time:
            for n in node:
                self.variables.heat_price[n,t] = m.addVar(lb=-300,ub=500,name='elec price({0},{1})'.format(n,t))


        self.variables.lambda_c9 = {} # energy stored updates
        
        self.variables.mu_c10a = {} # energy stored min
        self.variables.mu_c10b = {} # energy stored max
        self.variables.mu_c11a = {} # energy charged min
        self.variables.mu_c11b = {} # energy charged max
        self.variables.mu_c12a = {} # energy discharged min
        self.variables.mu_c12b = {} # energy discharged max

        for t in time:
            for g in heat_storage:
                
                self.variables.lambda_c9 = m.addVar(lb=-gb.GRB.INFINITY)
                
                self.variables.mu_c10a[g,t] = m.addVar(lb=0)
                self.variables.mu_c10b[g,t] = m.addVar(lb=0) 
                self.variables.mu_c11a[g,t] = m.addVar(lb=0)
                self.variables.mu_c11b[g,t] = m.addVar(lb=0) 
                self.variables.mu_c12a[g,t] = m.addVar(lb=0)
                self.variables.mu_c12b[g,t] = m.addVar(lb=0) 



        self.variables.mu_c14 = {} # pmin bid upper 1
        self.variables.mu_c15 = {} # pmin bid upper 2
        self.variables.mu_c16a = {} # pmax bid lower 1
        self.variables.mu_c16b = {} # pmax bid upper 1
        
        for t in time:
            for g in CHP:
                
                self.variables.mu_c14[g,t] = m.addVar(lb=0)
                self.variables.mu_c15[g,t] = m.addVar(lb=0)
                self.variables.mu_c16a[g,t] = m.addVar(lb=0)
                self.variables.mu_c16b[g,t] = m.addVar(lb=0)                
      
        
        self.variables.mu_c18a = {} # pmax
        self.variables.mu_c18b = {} # pmin 
        
        for t in time:
            for n in gen+CHP:
                
                self.variables.mu_c18a[n,t] = m.addVar(lb=0)
                self.variables.mu_c18b[n,t] = m.addVar(lb=0)

        
        self.variables.mu_c19a = {} # min load shedding
        self.variables.mu_c19b = {} # max load shedding
        for t in time:
            for n in zone_DK:
                self.variables.mu_c19a[n,t] = m.addVar(lb=0)
                self.variables.mu_c19b[n,t] = m.addVar(lb=0)


        self.variables.mu_c20a = {} # min solar
        self.variables.mu_c20b = {} # max solar
        for t in time:
            for n in zone_DK:
                self.variables.mu_c20a[n,t] = m.addVar(lb=0)
                self.variables.mu_c20b[n,t] = m.addVar(lb=0)


        self.variables.mu_c21a = {} # min wind
        self.variables.mu_c21b = {} # max wind
        
        for t in time:
            for n in zone_DK:
                self.variables.mu_c21a[n,t] = m.addVar(lb=0)
                self.variables.mu_c21b[n,t] = m.addVar(lb=0)

        self.variables.mu_c22a = {} # min elec flow
        self.variables.mu_c22b = {} # max elec flow
        for t in time:
            for l in line:
                self.variables.mu_c22a[n,t] = m.addVar(lb=0)
                self.variables.mu_c22b[n,t] = m.addVar(lb=0)
                
        self.variables.elec_price = {} # elec prices <=> lambda c23
        for t in time:
            for n in zone_DK:
                self.variables.elec_price[n,t] = m.addVar(lb=-300,ub=500,name='elec price({0},{1})'.format(n,t))
                            
        m.update()
    
    def _build_objective(self,D_day): # building the objective function for the EAHUC
  
        m = self.model    
              

        m.setObjective( alpha*( - (sum(c_heat_bid['EAHUC',g,D_day,t] for g in HP+incinerator+boiler for t in time) + sum(c_heat_bid['EAHUC',g,D_day,t][0] for g in CHP for t in time) )/(len(time)*len(HP+CHP+CHP+incinerator+boiler))*gb.quicksum( self.variables.HF_in[n1,n2,t] - self.variables.HF_out[n1,n2,t] for n1 in node for n2 in node_node[n1] for t in time)    - (sum(c_heat_bid['EAHUC',g,D_day,t] for g in HP+incinerator+boiler for t in time) + sum(c_heat_bid['EAHUC',g,D_day,t][0] for g in CHP for t in time) )/(len(time)*len(HP+CHP+CHP+incinerator+boiler))*gb.quicksum(self.variables.energy_stored[g,time[-1]] for g in heat_storage) + gb.quicksum(500*self.variables.HS[n,t] for n in node for t in time) + gb.quicksum(c_heat_bid['EAHUC',g,D_day,t][x]*self.variables.Q[g,t,x] for t in time for g in CHP for x in range(2)) + gb.quicksum(c_heat_bid['EAHUC',g,D_day,t]*self.variables.Q[g,t] for t in time for g in incinerator+HP+boiler) + gb.quicksum(self.variables.r[h,t] + c_0[h]*self.variables.u[h,t] for t in time for h in CHP+incinerator+HP+boiler) ) +  (1-alpha)*( gb.quicksum(c_elec[g]*self.variables.P[g,t] for t in time for g in CHP+gen) + gb.quicksum(500*self.variables.ES[n,t] for t in time for n in zone_DK) ) ,   
            gb.GRB.MINIMIZE)


        
    def _build_constraints(self,D_day):
     
        m = self.model
        
        # UL : UC, accepted (valid bids), time delays and flow direction in pipelines
                
        self.constraints.start_cost = {}

        for g in CHP+boiler+incinerator+HP:        
            for x in range(len(time)):
                for h in time_start_range[g]:
                    self.constraints.start_cost[g,time[x],h] = m.addConstr(
                            self.variables.r[g,time[x]],
                            gb.GRB.GREATER_EQUAL,
                            c_start[g,h]*(self.variables.u[g,time[x]] - gb.quicksum(self.variables.u[g,time_extended[len(time_init)+x-k-1]] for k in range(h))),name='start up cost compute({0},{1},{2})'.format(g,time[x],h))


        self.constraints.time_on_min = {}
        self.constraints.time_off_min = {}

        for g in CHP+boiler+incinerator+HP:
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
        
        for g in CHP+boiler+incinerator+HP:
            for x in range(len(time)):
                self.constraints.uc[g,time[x]] = m.addConstr(
                        self.variables.u[g,time_extended[len(time_init)+x]]-self.variables.u[g,time_extended[len(time_init)+x-1]],
                        gb.GRB.EQUAL,      
                        self.variables.v_on[g,time_extended[len(time_init)+x]]-self.variables.v_off[g,time_extended[len(time_init)+x]],name='uc binary variables ({0},{1})'.format(g,time[x]))


        self.constraints.bid_accepted = {}
        
        for g in CHP:
            for t in time:
                self.constraints.bid_accepted[g,t,0] = m.addConstr(
                        self.variables.bid_accepted[g,t,0],
                        gb.GRB.LESS_EQUAL,      
                        self.variables.u[g,t],name='bid accepted'.format(g,t,0))

                self.constraints.bid_accepted[g,t,1] = m.addConstr(
                        self.variables.bid_accepted[g,t,1],
                        gb.GRB.LESS_EQUAL,      
                        self.variables.bid_accepted[g,t,0],name='bid accepted'.format(g,t,1))



        # flow direction

        self.constraints.flow_direction = {}  

        for t in time:

            for n1 in node:
                for n2 in node_node[n1]:

                    self.constraints.flow_direction[n1,n2,t] = m.addConstr(
                        self.variables.HF_direction[n2,n1,t]+self.variables.HF_direction[n1,n2,t],
                        gb.GRB.EQUAL,
                        1)

        # time delays
        
        self.constraints.time_delay_sum = {}
        
        self.constraints.time_delay_direction = {}
        
        for t in time:
            for n1 in node:
                for n2 in node_node[n1]:
                        
                    self.constraints.time_delay_sum[n1,n2,t] = m.addConstr(
                        gb.quicksum( self.variables.time_delay[n1,n2,t,x] for x in range(time_delay_max+1) ),
                        gb.GRB.EQUAL,
                        1)
                        
                    self.constraints.time_delay_direction[n1,n2,t] = m.addConstr(
                        self.variables.time_delay[n1,n2,t,0],
                        gb.GRB.GREATER_EQUAL,
                        1-self.variables.HF_direction[n1,n2,t])

        # time delays bounds: first in first out
        
        self.constraints.first_in_first_out = {}
        
        for t in range(len(time)):
            for x in range(1,time_delay_max+1):
                for n1 in node:
                    for n2 in node_node[n1]:                            
                        self.constraints.first_in_first_out[n1,n2,time[t],time_extended[len(time_init)+t-x]] = m.addConstr(
                            t + gb.quicksum(k*self.variables.time_delay[n1,n2,time[t],k] for k in range(time_delay_max+1)),
                            gb.GRB.GREATER_EQUAL,
                            t - x + gb.quicksum(k*self.variables.time_delay[n1,n2,time_extended[len(time_init)+t-x],k] for k in range(time_delay_max+1)))
     

        # bid validity
        
        
        self.constraints.bid_validity_1 = {}
        self.constraints.bid_validity_2 = {}
        
        for t in time:
            for g in CHP:
                self.constraints.bid_validity_1[g,t,0] = m.addConstr(
                    c_heat_bid['EAHUC',g,D_day,t][0] + 501*(1-self.variables.bid_accepted[g,t,0]),
                    gb.GRB.GREATER_EQUAL,
                    rho_heat[g]/rho_elec[g]*self.variables.elec_price[zone_prod[g],t])
          

                self.constraints.bid_validity_1[g,t,1] = m.addConstr(
                    c_heat_bid['EAHUC',g,D_day,t][1] + 501*(1-self.variables.bid_accepted[g,t,1]),
                    gb.GRB.GREATER_EQUAL,
                    rho_heat[g]/rho_elec[g]*self.variables.elec_price[zone_prod[g],t])
                
                self.constraints.bid_validity_2[g,t] = m.addConstr(
                    c_heat_bid['EAHUC',g,D_day,t][1] + 501*(1-self.variables.bid_accepted[g,t,1]),
                    gb.GRB.GREATER_EQUAL,
                    (c_elec[g] - self.variables.elec_price[zone_prod[g],t])*rmin[g] + c_heat[g])

            for g in HP:
                self.constraints.bid_validity_2[g,t] = m.addConstr(
                    c_heat_bid['EAHUC',g,D_day,t] + 501*(1-self.variables.u[g,t]),
                    gb.GRB.GREATER_EQUAL,
                    self.variables.elec_price[zone_prod[g],t]*rmin[g])

            for g in incinerator:
                self.constraints.bid_validity_2[g,t] = m.addConstr(
                    c_heat_bid['EAHUC',g,D_day,t] + 501*(1-self.variables.u[g,t]),
                    gb.GRB.GREATER_EQUAL,
                    (c_elec[g] - self.variables.elec_price[zone_prod[g],t])*rmin[g] + c_heat[g])

        
        ## ML: heat market constraints
        
        # min/max production depending on accepted bids and on/off status
        
        self.constraints.heat_minprod = {} #c1a                
        self.constraints.heat_maxprod = {} # c1b 
        
        for t in time:

            for g in boiler+incinerator+HP:
                               
                self.constraints.heat_minprod[g,t] = m.addConstr(
                    self.variables.Q[g,t],
                    gb.GRB.GREATER_EQUAL,
                    qmin[g]*self.variables.u[g,t],name='elec minprod({0},{1})'.format(g,t))

                self.constraints.heat_maxprod[g,t] = m.addConstr(
                    self.variables.Q[g,t],
                    gb.GRB.LESS_EQUAL,
                    qmax[g]*self.variables.u[g,t],name='elec maxprod({0},{1})'.format(g,t))
 
    
            for g in CHP:

                self.constraints.heat_minprod[g,t] = m.addConstr(
                    self.variables.Q[g,t,0],
                    gb.GRB.GREATER_EQUAL,
                    q_bid[g][0]*self.variables.bid_accepted[g,t,1],name='elec minprod({0},{1},{2})'.format(g,t,0))
                
                for x in range(2):

                    self.constraints.heat_maxprod[g,t,x] = m.addConstr( 
                        self.variables.Q[g,t,x],
                        gb.GRB.LESS_EQUAL,
                        q_bid[g][x]*self.variables.bid_accepted[g,t,x],name='elec maxprod({0},{1},{2})'.format(g,t,x))

        # c2a and c2b: bounds on load shedding
        
        # flow bounds based on flow direction
                  
        self.constraints.flow_direction_in = {} # c3 (implicitly lower bounded in c5 and c6)
        self.constraints.flow_direction_out = {} # c4 (implicitly lower bounded in c7)
        
        for t in time:
            for n1 in node:
                for n2 in node_node[n1]:
                    
                    self.constraints.flow_direction_in[n1,n2,t] = m.addConstr(
                        self.variables.HF_in[n1,n2,t],
                        gb.GRB.LESS_EQUAL,
                        f_upper[n1,n2]*self.variables.HF_direction[n1,n2,t])
                    
                    # I THINK THIS MAY BE REDUNDANT?        
                    self.constraints.flow_direction_out[n1,n2,t] = m.addConstr(
                        self.variables.HF_out[n1,n2,t],
                        gb.GRB.LESS_EQUAL,
                        f_upper[n1,n2]*self.variables.HF_direction[n1,n2,t]) 
        
        # in/out flows based on time delays


        self.constraints.HF_delay_upper_1 = {} #c5b
        self.constraints.HF_delay_lower_1 = {} #c5a

        self.constraints.HF_delay_upper_2 = {} #c6b
        self.constraints.HF_delay_lower_2 = {} #c6a
        # these constraints also provide upper and lower bounds for HF_in and HF_in_delay
        
        for t in range(len(time)):
            for x in range(time_delay_max+1):
                for n1 in node:
                    for n2 in node_node[n1]:
                            
                        self.constraints.HF_delay_upper_1[n1,n2,time[t],x] = m.addConstr(
                            - self.variables.HF_in_delay[n1,n2,time[t],x] + self.variables.HF_in[n1,n2,time_extended[len(time_init)+t-x]],
                            gb.GRB.LESS_EQUAL,
                            f_upper[n1,n2]*(1-self.variables.time_delay[n1,n2,time_extended[len(time_init)+t-x],x]))
                        
                        self.constraints.HF_delay_lower_1[n1,n2,time[t],x] = m.addConstr(
                            - self.variables.HF_in_delay[n1,n2,time[t],x] + self.variables.HF_in[n1,n2,time_extended[len(time_init)+t-x]],
                            gb.GRB.GREATER_EQUAL,
                            0)          
    
                            
                        self.constraints.HF_delay_upper_2[n1,n2,time[t],x] = m.addConstr(
                            self.variables.HF_in_delay[n1,n2,time[t],x],
                            gb.GRB.LESS_EQUAL,
                            f_upper[n1,n2]*self.variables.time_delay[n1,n2,time_extended[len(time_init)+t-x],x])
                        
                        self.constraints.HF_delay_lower_2[n1,n2,time[t],x] = m.addConstr(
                            self.variables.HF_in_delay[n1,n2,time[t],x],
                            gb.GRB.GREATER_EQUAL,
                            0) 
                        
        self.constraints.HF_out = {} #c7
        
        for t in time:
            for n1 in node:
                for n2 in node_node[n1]:
                    
                    self.constraints.HF_out[n1,n2,t] = m.addConstr(
                        self.variables.HF_out[n1,n2,t],
                        gb.GRB.EQUAL,
                        gb.quicksum(self.variables.HF_in_delay[n1,n2,t,x] for x in range(time_delay_max+1)) )
                        
        # heat balance at each node

        self.constraints.heat_balance = {} # c8
        
        for t in time:
            for n1 in node:
                        
                self.constraints.heat_balance[n1,t] = m.addConstr(
                    gb.quicksum(self.variables.Q_discharge[g,t] - self.variables.Q_charge[g,t] for g in heat_storage_node[n1]) + gb.quicksum(self.variables.Q[g,t] for g in prod_node[n1]) + gb.quicksum(self.variables.Q[g,t,x] for g in CHP_node[n1] for x in range(2)) + gb.quicksum(self.variables.HF_out[n2,n1,t] - self.variables.HF_in[n1,n2,t] for n2 in node_node[n1]),
                    heat_load[n1,hiy[D_day,t]]-self.variables.HS[n1,t])        
                        

        # heat storage constraints
        
        self.constraints.heat_storage_update = {} # c9
         
        for t in range(len(time)):
            for g in heat_storage:
                        
                self.constraints.heat_storage_update[g,t] = m.addConstr(
                    self.variables.energy_stored[g,time_extended[len(time_init)+t]],
                    gb.GRB.EQUAL,
                    self.variables.energy_stored[g,time_extended[len(time_init)+t-1]] - rho_discharge[g]*self.variables.Q_discharge[g,time_extended[len(time_init)+t]] + rho_charge[g]*self.variables.Q_charge[g,time_extended[len(time_init)+t]] - heat_loss[g]/10)           

        # c10 = energy bounds , c11 = charge bounds, c12 = discharge bounds
        
        # adjsuted bids of CHPs for electricity market


        self.constraints.pmin_bid_upper_1 = {} #c14
        self.constraints.pmin_bid_upper_2 = {} #c15


        self.constraints.pmax_bid_lower = {} #c16a        
        self.constraints.pmax_bid_upper = {} #c16b
        
        for t in time:
            for g in CHP:

                                
                self.constraints.pmin_bid_upper_1[g,t] = m.addConstr(
                    - self.variables.pmin_bid[g,t] + (fmin[g]-rho_heat[g]*gb.quicksum(self.variables.Q[g,t,x] for x in range(2)))/rho_elec[g],
                    gb.GRB.LESS_EQUAL,
                    fmin[g]/rho_elec[g]*(1-self.variables.u[g,t]),name='pmin bid upper 2({0},{1})'.format(g,t))



                self.constraints.pmin_bid_upper_2[g,t] = m.addConstr(
                    - self.variables.pmin_bid[g,t] + rmin[g]*gb.quicksum(self.variables.Q[g,t,x] for x in range(2)),
                    gb.GRB.LESS_EQUAL,
                    0,name='pmin bid upper 3({0},{1})'.format(g,t))                



                
                self.constraints.pmax_bid_lower_1[g,t] = m.addConstr(
                    - self.variables.pmax_bid[g,t] + (fmax[g]-rho_heat[g]*gb.quicksum(self.variables.Q[g,t,x] for x in range(2)))/rho_elec[g],
                    gb.GRB.GREATER_EQUAL,
                    fmax[g]/rho_elec[g]*(1-self.variables.u[g,t]),name='pmax bid lower 1({0},{1})'.format(g,t))
                
                self.constraints.pmax_bid_upper_1[g,t] = m.addConstr(
                    - self.variables.pmax_bid[g,t] + (fmax[g]-rho_heat[g]*gb.quicksum(self.variables.Q[g,t,x] for x in range(2)))/rho_elec[g],
                    gb.GRB.LESS_EQUAL,
                    fmax[g]/rho_elec[g]*(1-self.variables.u[g,t]),name='pmax bid upper 1({0},{1})'.format(g,t))

                
        ## LL: electricity market constraints

                
        #4) ELEC GEN MAX PROD and total production variable

        self.constraints.elec_minprod = {} # c18a        
        self.constraints.elec_maxprod = {} # c18b

        
        for t in time:

            for g in gen:

                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.GREATER_EQUAL,
                    pmin[g]*u_set['EAHUC',g,D_day,t],name='elec minprod({0},{1})'.format(g,t))
                
                self.constraints.elec_maxprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    pmax[g]*u_set['EAHUC',g,D_day,t],name='elec maxprod({0},{1})'.format(g,t))

            for g in CHP:

                self.constraints.elec_minprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.GREATER_EQUAL,
                    self.variables.pmin_bid[g,t],name='elec minprod({0},{1})'.format(g,t))
                    
                self.constraints.elec_maxprod[g,t] = m.addConstr(
                    self.variables.P[g,t],
                    gb.GRB.LESS_EQUAL,
                    self.variables.pmax_bid[g,t],name='elec maxprod({0},{1})'.format(g,t))                               

        # c19a and c19b = bounds on elec load shedding
        
        #c20a and c20b : bounds on solar prod
        
        #c21a and c21b : bounds on wind prod
        
        # c22a and c22b = bounds on elec flows      
        
        # elec balance
        
        self.constraints.elec_balance = {} # c23
        
        for t in time:
            for z in zone_DK:
                
                self.constraints.elec_balance[z,t] = m.addConstr(
                    gb.quicksum(self.variables.P[g,t] for g in gen_zone[z]+CHP_zone[z]) + gb.quicksum(rmin[h]*self.variables.Q[h,t] for h in incinerator_zone[z]) + self.variables.wind_prod[z,t]+self.variables.solar_prod[z,t]+self.variables.ES[z,t],
                    gb.GRB.EQUAL,
                    elec_load[z,hiy[D_day,t]] + gb.quicksum(rmin[h]*self.variables.Q[h,t] for h in HP_zone[z]) + gb.quicksum(self.variables.EF[l,t] for l in line_start[z]) - gb.quicksum(self.variables.EF[l,t] for l in line_end[z]),name='elec balance({0},{1})'.format(z,t))


        #  DUAL CONSTRAINTS
        
        # stationarity constraints DL/dx = 0 
        
        self.constraints.dL_dQ = {} 
        
        for t in time:
            for g in boiler:
                
                self.constraints.dL_dQ[g,t] = m.addConstr(
                    alpha*c_heat_bid[g,D_day,t]-self.variables.mu_c1a[g,t]+self.variables.mu_c1b[g,t]-self.variables.heat_price[node_prod[g],t],
                    gb.GRB.EQUAL,
                    0)

            for g in HP:
                
                self.constraints.dL_dQ[g,t] = m.addConstr(
                    alpha*c_heat_bid[g,D_day,t]-self.variables.mu_c1a[g,t]+self.variables.mu_c1b[g,t]-self.variables.heat_price[node_prod[g],t]+self.variables.elec_price[zone_prod[g],t],
                    gb.GRB.EQUAL,
                    0)

            for g in incinerator:
                
                self.constraints.dL_dQ[g,t] = m.addConstr(
                    alpha*c_heat_bid[g,D_day,t]-self.variables.mu_c1a[g,t]+self.variables.mu_c1b[g,t]-self.variables.heat_price[node_prod[g],t]-self.variables.elec_price[zone_prod[g],t],
                    gb.GRB.EQUAL,
                    0)
                
            for g in CHP:
                
                self.constraints.dL_dQ[g,t,0] = m.addConstr(
                    alpha*c_heat_bid[g,D_day,t][0]-self.variables.mu_c1a[g,t]+self.variables.mu_c1b[g,t,0]-self.variables.heat_price[node_prod[g],t] - rho_heat[g]/rho_elec[g]*self.variables.mu_c14[g,t] + rmin[g]*self.variables.mu_c15[g,t] + rho_heat[g]/rho_elec[g]*self.variables.mu_c16a[g,t] - rho_heat[g]/rho_elec[g]*self.variables.mu_c16b[g,t],
                    gb.GRB.EQUAL,
                    0)        

                self.constraints.dL_dQ[g,t,1] = m.addConstr(
                    alpha*c_heat_bid[g,D_day,t][1]+self.variables.mu_c1b[g,t,1]-self.variables.heat_price[node_prod[g],t] - rho_heat[g]/rho_elec[g]*self.variables.mu_c14[g,t] + rmin[g]*self.variables.mu_c15[g,t] + rho_heat[g]/rho_elec[g]*self.variables.mu_c16a[g,t] - rho_heat[g]/rho_elec[g]*self.variables.mu_c16b[g,t],
                    gb.GRB.EQUAL,
                    0) 
                    
        self.constraints.dL_dHS = {} 
        
        for t in time:
            for n in node:
                
                self.constraints.dL_dHS[n,t] = m.addConstr(
                    alpha*500-self.variables.mu_c2a[n,t]+self.variables.mu_c2b[n,t] - self.variables.heat_price[n,t],
                    gb.GRB.EQUAL,
                    0)

        self.constraints.dL_dHF_in = {} 
        self.constraints.dL_dHF_in_delay = {} 
        self.constraints.dL_dHF_out = {} 
        
        for x in range(len(time)):
            for n1 in node:
                for n2 in node_node[n1]:
                
                    self.constraints.dL_dHF_in[n1,n2,time[x]] = m.addConstr(
                        alpha*(sum(c_heat_bid['EAHUC',g,D_day,t] for g in HP+incinerator+boiler for t in time) + sum(c_heat_bid['EAHUC',g,D_day,t][0] for g in CHP for t in time) )/(len(time)*len(HP+CHP+CHP+incinerator+boiler))  +  self.variables.mu_c3[n1,n2,time[x]] + gb.quicksum( - self.variables.mu_c5a[n1,n2,time[x+k],k]  + self.variables.mu_c5b[n1,n2,time[x+k],k] for k in range(min(time_delay_max+1,len(time)-x))) + self.variables.heat_price[n1,time[x]],
                        gb.GRB.EQUAL,
                        0)
                    
                    for x in range(time_delay_max+1):
                        
                        self.constraints.dL_dHF_in_delay[n1,n2,t,x] = m.addConstr(
                            self.variables.mu_c5a[n1,n2,t,x]-self.variables.mu_c5b[n1,n2,t,x]-self.variables.mu_c6a[n1,n2,t,x]+self.variables.mu_c6b[n1,n2,t,x]-self.variables.lambda_c7[n1,n2,t],
                            gb.GRB.EQUAL,
                            0)

                    self.constraints.dL_dHF_out[n1,n2,t] = m.addConstr(
                        - alpha*(sum(c_heat_bid['EAHUC',g,D_day,t] for g in HP+incinerator+boiler for t in time) + sum(c_heat_bid['EAHUC',g,D_day,t][0] for g in CHP for t in time) )/(len(time)*len(HP+CHP+CHP+incinerator+boiler)) + self.variables.mu_c4[n1,n2,t] + self.variables.lambda_c7[n1,n2,t] - self.variables.heat_price[n2,t],
                        gb.GRB.EQUAL,
                        0)
                    
        self.constraints.dL_denergy_stored = {}                     
        self.constraints.dL_dQ_charge = {} 
        self.constraints.dL_dQ_discharge = {} 

        for g in heat_storage:

            for (t,t2) in zip(time[:-1],time[1:]):

                self.constraints.dL_denergy_stored[g,t] = m.addConstr(
                    - self.variablers.mu_c10a[g,t] + self.variablers.mu_c10b[g,t] + self.variables.lambda_c9[g,t] - self.variables.lambda_c9[g,t2],
                    gb.GRB.EQUAL,
                    0)  

            self.constraints.dL_denergy_stored[g,time[-1]] = m.addConstr(
                - alpha*(sum(c_heat_bid['EAHUC',g,D_day,x] for g in HP+incinerator+boiler for x in time) + sum(c_heat_bid['EAHUC',g,D_day,x][0] for g in CHP for x in time) )/(len(time)*len(HP+CHP+CHP+incinerator+boiler)) - self.variablers.mu_c10a[g,time[-1]] + self.variablers.mu_c10b[g,time[-1]] + self.variables.lambda_c9[g,time[-1]],
                gb.GRB.EQUAL,
                0)
                
            for t in time:    
                
                self.constraints.dL_dQ_charge[g,t] = m.addConstr(
                    self.variables.heat_price[node_prod[g],t] - self.variables.mu_c11a[g,t] + self.variables.mu_c11b[g,t] - rho_charge[g]*self.variables.lambda_c9[g,t],
                    gb.GRB.EQUAL,
                    0)                    

                self.constraints.dL_dQ_discharge[g,t] = m.addConstr(
                    -self.variables.heat_price[node_prod[g],t] - self.variables.mu_c12a[g,t] + self.variables.mu_c12b[g,t] + rho_discharge[g]*self.variables.lambda_c9[g,t],
                    gb.GRB.EQUAL,
                    0)

        
        self.constraints.dL_dpmin_bid = {} 
        self.constraints.dL_dpmax_bid = {} 
        
        for t in time:
            for g in CHP:
                
                self.constraints.dL_dpmin_bid[g,t] = m.addConstr(
                    self.variables.mu_c18a[g,t] - self.variables.mu_c14[g,t] - self.variables.mu_c15[g,t],
                    gb.GRB.EQUAL,
                    0)
                
                self.constraints.dL_dpmax_bid[g,t] = m.addConstr(
                    -self.variables.mu_c18b[g,t] + self.variables.mu_c16a[g,t]- self.variables.mu_c16b[g,t],
                    gb.GRB.EQUAL,
                    0)

        self.constraints.dL_dP = {} 
        
        for t in time:
            for g in CHP+gen:
                
                self.constraints.dL_dP[g,t] = m.addConstr(
                    (1-alpha)*c_elec[g]-self.variables.elec_price[zone_prod[g],t]-self.variables.mu_c18a[g,t]+self.variables.mu_c18b[g,t],
                    gb.GRB.EQUAL,
                    0)
                
        self.constraints.dL_dES = {} 
        self.constraints.dL_dwind_prod = {} 
        self.constraints.dL_dsolar_prod = {} 
        
        for t in time:
            for z in zone_DK:
                
                self.constraints.dL_dES[z,t] = m.addConstr(
                    (1-alpha)*500 - self.variables.elec_price[z,t] - self.variables.mu_c19a[z,t] + self.variables.mu_c19b[z,t],
                    gb.GRB.EQUAL,
                    0)

                self.constraints.dL_dsolar_prod[z,t] = m.addConstr(
                    -self.variables.mu_c20a[z,t]+self.variables.mu_c20b[z,t] - self.variables.elec_price[z,t],
                    gb.GRB.EQUAL,
                    0)
                
                self.constraints.dL_dwind_prod[z,t] = m.addConstr(
                    -self.variables.mu_c21a[z,t]+self.variables.mu_c21b[z,t] - self.variables.elec_price[z,t],
                    gb.GRB.EQUAL,
                    0)

        self.constraints.dL_dEF = {} 
        
        for t in time:
            for l in line:
                
                self.constraints.dL_dEF[l,t] = m.addConstr(
                    - self.variables.mu_c22a[l,t] + self.variables.mu_c22b[l,t] + self.variables.elec_price[line_connexion[l][0],t] - self.variables.elec_price[line_connexion[l][1],t],
                    gb.GRB.EQUAL,
                    0)

        # strong duality or complementarity conditions???
        
        self.constraint.strong_duality = m.addConstr( gb.quicksum(0 for t in time),
                gb.GRB.EQUAL,
                alpha*( - (sum(c_heat_bid['EAHUC',g,D_day,t] for g in HP+incinerator+boiler for t in time) + sum(c_heat_bid['EAHUC',g,D_day,t][0] for g in CHP for t in time) )/(len(time)*len(HP+CHP+CHP+incinerator+boiler))*gb.quicksum( self.variables.HF_in[n1,n2,t] - self.variables.HF_out[n1,n2,t] for n1 in node for n2 in node_node[n1] for t in time)    - (sum(c_heat_bid['EAHUC',g,D_day,t] for g in HP+incinerator+boiler for t in time) + sum(c_heat_bid['EAHUC',g,D_day,t][0] for g in CHP for t in time) )/(len(time)*len(HP+CHP+CHP+incinerator+boiler))*gb.quicksum(self.variables.energy_stored[g,time[-1]] for g in heat_storage) + gb.quicksum(500*self.variables.HS[n,t] for n in node for t in time) + gb.quicksum(c_heat_bid['EAHUC',g,D_day,t][x]*self.variables.Q[g,t,x] for t in time for g in CHP for x in range(2)) + gb.quicksum(c_heat_bid['EAHUC',g,D_day,t]*self.variables.Q[g,t] for t in time for g in incinerator+HP+boiler) ) +  (1-alpha)*( gb.quicksum(c_elec[g]*self.variables.P[g,t] for t in time for g in CHP+gen) + gb.quicksum(500*self.variables.ES[n,t] for t in time for n in zone_DK) ) )

        # McCormick envelopes (EXACT) of bilinear term u*\mu in strong duality!                     

        self.constraints.MC1_upper_1 = {} 
        self.constraints.MC1_lower_1 = {} 
        self.constraints.MC1_upper_2 = {} 
        self.constraints.MC1_lower_2 = {} 
        
        for t in time:
            for g in incinerator+HP+boiler:
                
                self.constraints.MC1_lower_1[g,t] = m.addConstr(
                    self.variables.MC1[g,t],
                    gb.GRB.GREATER_EQUAL,
                    0)

                self.constraints.MC1_upper_1[g,t] = m.addConstr(
                    self.variables.MC1[g,t],
                    gb.GRB.LESS_EQUAL,
                    )

                self.constraints.MC1_lower_2[g,t] = m.addConstr(
                    self.variables.MC1[g,t],
                    gb.GRB.GREATER_EQUAL,
                    0)

                self.constraints.MC1_upper_2[g,t] = m.addConstr(
                    self.variables.MC1[g,t],
                    gb.GRB.LESS_EQUAL,
                    )
                    
#%% SOLVE                       

def solve_heat(m1,D):
    
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
    
    global solve_it
    
    global heat_cost
    global total_cost
    
    if m1 == 'EAHUC':
        solve_it = EAHUC(day[D])
        
    if m1 == 'HUC':
        solve_it = HUC(day[D])
        
        
    if m1 == 'CHPUC':
        solve_it = CHPUC(day[D])
        
        
    solve_it.optimize()
        
    for g in HP+incinerator+boiler+CHP:   
        
        for t in time:
            
            u_set[m1,g,day[D],t] = solve_it.variables.u[g,t].x
            v_on_set[m1,g,day[D],t] = solve_it.variables.v_on[g,t].x
            v_off_set[m1,g,day[D],t] = solve_it.variables.v_off[g,t].x

        if D<len(day)-1:
            
            for x in range(len(time_init)):
                
                u_init[m1,g,day[D+1],time_init[len(time_init)-x-1]] = u_set[m1,g,day[D],time[len(time)-x-1]]
                v_on_init[m1,g,day[D+1],time_init[len(time_init)-x-1]] = v_on_set[m1,g,day[D],time[len(time)-x-1]]
                v_off_init[m1,g,day[D+1],time_init[len(time_init)-x-1]] = v_off_set[m1,g,day[D],time[len(time)-x-1]]
    


    for g in HP+incinerator+boiler: 
        for t in time:

            Q_set[m1,g,day[D],t] = solve_it.variables.Q[g,t].x
            
    

    for g in heat_storage:
        
        for t in time:
            
            energy_stored_set[m1,g,day[D],t] = solve_it.variables.energy_stored[g,t].x 

        if D<len(day)-1:
            for x in range(len(time_init)):
                            
                energy_stored_init[m1,g,day[D+1],time_init[len(time_init)-x-1]] = energy_stored_set[m1,g,day[D],time[len(time)-x-1]]


    for n1 in node:
        for n2 in node_node[n1]:
            
            for t in time:
    
                HF_in_set[m1,n1,n2,day[D],t] = solve_it.variables.HF_in[n1,n2,t].x
                HF_direction_set[m1,n1,n2,day[D],t] = solve_it.variables.HF_direction[n1,n2,t].x            
            
                for k in range(time_delay_max+1):
                    time_delay_set[m1,n1,n2,day[D],t,k] = solve_it.variables.time_delay[n1,n2,t,k].x


            if D<len(day)-1:
                for x in range(len(time_init)):
                
                              
                    HF_in_init[m1,n1,n2,day[D+1],time_init[len(time_init)-x-1]] = HF_in_set[m1,n1,n2,day[D],time[len(time)-x-1]]
                    HF_direction_init[m1,n1,n2,day[D+1],time_init[len(time_init)-x-1]] = HF_direction_set[m1,n1,n2,day[D],time[len(time)-x-1]]
            
                    for k in range(time_delay_max+1):
                
                        time_delay_init[m1,n1,n2,day[D+1],time_init[len(time_init)-x-1],k] = time_delay_set[m1,n1,n2,day[D],time[len(time)-x-1],k]
            
    
        
    if m1 == 'CHPUC':
        
        total_cost[m1,day[D]] = solve_it.model.getObjective().getValue()        
        
        for g in gen: 
            
            for t in time:
            
                u_set[m1,g,day[D],t] = solve_it.variables.u[g,t].x
                v_on_set[m1,g,day[D],t] = solve_it.variables.v_on[g,t].x
                v_off_set[m1,g,day[D],t] = solve_it.variables.v_off[g,t].x

            if D<len(day)-1:
                
                for x in range(len(time_init)):
                    
                    u_init[m1,g,day[D+1],time_init[len(time_init)-x-1]] = u_set[m1,g,day[D],time[len(time)-x-1]]
                    v_on_init[m1,g,day[D+1],time_init[len(time_init)-x-1]] = v_on_set[m1,g,day[D],time[len(time)-x-1]]
                    v_off_init[m1,g,day[D+1],time_init[len(time_init)-x-1]] = v_off_set[m1,g,day[D],time[len(time)-x-1]]
        

        for g in CHP:    
            for t in time:            
                Q_set[m1,g,day[D],t] = solve_it.variables.Q[g,t].x

                pmin_bid[m1,g,day[D],t] = max(rmin[g]*Q_set[m1,g,day[D],t],(fmin[g]-rho_heat[g]*Q_set[m1,g,day[D],t])/rho_elec[g])
                pmax_bid[m1,g,day[D],t] = (fmax[g]-rho_heat[g]*Q_set[m1,g,day[D],t])/rho_elec[g]

        for g in gen:
            for t in time: 
                
                pmin_bid[m1,g,day[D],t] = u_set[m1,g,day[D],t]*pmin[g]
                pmax_bid[m1,g,day[D],t] = u_set[m1,g,day[D],t]*pmax[g]



    if m1 == 'HUC':
        
        heat_cost[m1,day[D]] = solve_it.model.getObjective().getValue()
        
        for g in CHP:
            for t in time:
                
                Q_set[m1,g,day[D],t] = solve_it.variables.Q[g,t,0].x + solve_it.variables.Q[g,t,1].x

                pmin_bid[m1,g,day[D],t] = max(rmin[g]*Q_set[m1,g,day[D],t],(fmin[g]-rho_heat[g]*Q_set[m1,g,day[D],t])/rho_elec[g])
                pmax_bid[m1,g,day[D],t] = (fmax[g]-rho_heat[g]*Q_set[m1,g,day[D],t])/rho_elec[g]

    if m1 == 'EAHUC':
        
        heat_cost[m1,day[D]] = solve_it.model.getObjective().getValue()
        
        for g in CHP:
            for t in time:
                
                Q_set[m1,g,day[D],t] = solve_it.variables.Q[g,t,0].x + solve_it.variables.Q[g,t,1].x

                pmin_bid[m1,g,day[D],t] = solve_it.variables.pmin_bid[g,t].x
                pmax_bid[m1,g,day[D],t] = solve_it.variables.pmax_bid[g,t].x

#%%

#if solve_it.model.status >= 2 :
##    
#solve_it.computeIIS()
#constraints = solve_it.model.getConstrs()
#variables = solve_it.model.getVars()
#
##Print the names of all of the constraints in the IIS set.
#print('constraint:')
#for c in constraints:
#    if c.IISConstr>0:
#        print(c.ConstrName)               
#
##Print the names of all of the variables in the IIS set.
#print('lower bound:')
#for v in variables:
#    if v.IISLB>0:
#        print(v.VarName)
#
#print('upper bound:')        
#for v in variables:
#    if v.IISUB>0:
#        print(v.VarName)