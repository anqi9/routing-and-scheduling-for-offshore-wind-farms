# -*- coding: utf-8 -*-

import numpy as np
from gurobipy import Model, GRB, quicksum
import itertools

class GB():
    def __init__(self, data):
        self.period = data['period']  # days: 3-7
        # available working hours based on the weather and vessels type (/hour)
        self.time_window = data['time_window']
        # max. num. of turbines can be visited at one period
        self.eta = data['eta']
        # indicate whether the base serves the wind farm
        self.lambda_bf = data['lambda']
        # technician cost (/day) corresponding to technician type
        self.tech_cost = data['technician']
        # num of wind farms
        self.B_coordinate = data['base']['coordinate']
        # num. of available technician at bases
        self.B_tech = data['base']['technician']
        # indicate available vessels at bases
        self.B_vessel = data['base']['vessel']
        # coordinate of turbines in wind farms
        self.WF_coordinate = data['wind_farm']['coordinate']
        # maintenance time (/hour) of corresponding turbines
        self.WF_mainTime = data['wind_farm']['maintenance_time']
        # num. of needed technician of each type
        self.WF_tech = data['wind_farm']['technician']
        # weight of needed spare parts by the turbine
        self.WF_partsWeight = data['wind_farm']['parts_weight'] 
        # turbines need vessels to be present
        self.WF_present = data['wind_farm']['present']
        # recommended period for maintenance
        self.WF_deadline = data['wind_farm']['deadline']
        # penalty cost if turbines not repaired before deadline
        self.WF_penCost = data['wind_farm']['penalty_cost']
        # indicate whether spare parts for the turbine are available at the period
        self.WF_availParts = data['wind_farm']['avail_parts']
        # spare parts capacity (/tonne) of vessel
        self.V_parts = data['vessel']['capacity']
        # indicate whether vessel is able to transfer spare parts needed by the turbine
        self.V_toTurbine = data['vessel']['to_turbine']
        # max num. of technician on board
        self.V_tech = data['vessel']['technician']
        # fuel cost (/hour) of vessels
        self.V_cost = data['vessel']['cost']
        # speed (km/hour) of vessels
        self.V_speed = data['vessel']['speed']
        # transfer time for technicians and equipment from vessel to turbine
        self.V_transTime = data['vessel']['trans_time']
        # indicate whether vessel available at the period
        self.V_available = data['vessel']['availability']
        # num of wind farms
        self.num_WF = 1#len(self.WF_tech)
        # num of bases
        self.num_B = 1#len(self.B_tech)
        # num of vessels
        self.num_V = 2#len(self.V_cost)

    # to calculate distance including a base
    def createDistanceMatrix(self, coordinate, base):
        Num = len(coordinate) + 1
        d = np.zeros((2*Num, 2*Num))
        for i in range(Num-1):
            p1 = np.array([coordinate[i][0], coordinate[i][1]])
            p3 = np.array([base[0], base[1]])
            d[0, i+1] = d[i+1, 0] = np.linalg.norm(p1 - p3)
            for j in range(i + 1, Num-1):
                p2 = np.array([coordinate[j][0], coordinate[j][1]])
                d[i+1, j+1] = d[j+1, i+1] = np.linalg.norm(p1 - p2)
        d[:Num, 2*Num-1] = d[:Num, 0]
        d[:Num,Num:-1] = d[:Num, 1:Num]
        d[Num:-1,:] = d[1:Num,:]
        d[2*Num-1,:] = d[0,:]
        # d = np.round(d, 2)
        # print (d)
        return d  # type: np.array

    def Solver(self):
        model = Model()
        n = len(self.WF_mainTime[0])
        dropoff = list(range(1, n+1))
        pickup = list(range(n+1, 2*n+1))
        nodes = list(range(2*n+2))
        B = list(range(self.num_B))
        WF = list(range(self.num_WF))
        D = list(range(self.period))
        V = list(range(self.num_V))
        P = list(range(len(self.B_tech[0])))
        waiting = [[i for i in dropoff if self.WF_present[wf][i-1] == 1] for wf in WF]
        
        # calculate distance
        dist = [[] for b in B]
        for b in B:
            for wf in WF:
                if self.lambda_bf[b][wf] == 1:
                    dist[b].append(self.createDistanceMatrix(self.WF_coordinate[wf], \
                    self.B_coordinate[b]))
                else: dist[b].append([])
        # print (dist)

        # calculate travel time
        travelTime = [[] for v in V]
        for v in V:
            for b in B:
                if self.B_vessel[b][v] == 1:
                    for wf in WF:
                        if self.lambda_bf[b][wf] == 1:
                            travelTime[v].append(dist[b][wf]/self.V_speed[v])
                        else: travelTime[v].append([])        
        # set decision variables
        AA = [(v,t,wf,i,j) for i in nodes for j in nodes for wf in WF \
            for t in D for v in V]
        BB = [(v,t,wf,i) for i in nodes for wf in WF for t in D for v in V]
        CC = [(v,t,wf,p,i) for i in nodes for p in P for wf in WF for t in D for v in V]
        DD = [(wf,i) for i in dropoff for wf in WF]
        y = model.addVars(AA, vtype = GRB.BINARY)
        T = model.addVars(BB, lb = 0, vtype = GRB.CONTINUOUS)
        Q = model.addVars(CC, lb = 0, vtype = GRB.INTEGER)
        x = model.addVars(DD, lb = 0, vtype = GRB.INTEGER)

################# Constraints ####################
        model.addConstrs(y[v,t,wf,i,i] == 0 for i in nodes for wf in WF \
            for t in D for v in V)
        model.addConstrs(quicksum(y[v,t,wf,0,j] for j in dropoff) <= 1  for wf in WF \
            for t in D for v in V)
        model.addConstrs(quicksum(y[v,t,wf,j,2*n+1] for j in pickup) <= 1 for wf in WF \
            for t in D for v in V)
        model.addConstrs(quicksum(y[v,t,wf,i,j] for i in nodes for t in D \
            for v in V) == 1 for j in pickup for wf in WF)
        model.addConstrs(quicksum(y[v,t,wf,i,j] for j in nodes) == \
            quicksum(y[v,t,wf,j,i] for j in nodes) for i in nodes for wf in WF \
            for t in D for v in V)
        model.addConstrs(quicksum(y[v,t,wf,i,j] for j in nodes) ==\
            quicksum(y[v,t,wf,n+i,j] for j in nodes) for i in dropoff for wf in WF \
            for t in D for v in V)
        model.addConstrs(quicksum(y[v,t,wf,i,n+i] for t in D for v in V) == 1 \
            for wf in WF for i in waiting[wf])
        
        model.addConstrs(T[v,t,wf,0] == 0 for wf in WF for t in D for v in V)
        model.addConstrs(T[v,t,wf,n+i] - T[v,t,wf,i] >= quicksum(y[v,t,wf,0,j] \
            for j in dropoff) * (self.WF_mainTime[wf][i-1] + self.V_transTime[v]) \
            for i in dropoff for v in V for t in D for wf in WF)
        model.addConstrs((y[v,t,wf,0,j]==1) >> (T[v,t,wf,j] >= travelTime[v][wf][0,j]) \
            for j in dropoff for wf in WF for t in D for v in V)
       
        model.addConstrs(quicksum(Q[v,t,wf,p,0] for p in P) <= self.V_tech[v] \
             for wf in WF for t in D for v in V)
        model.addConstrs(quicksum(Q[v,t,wf,p,0] * self.B_vessel[b][v] for v in V \
            for wf in WF) <= self.B_tech[b][p][t] for p in P for t in D for b in B)
        model.addConstrs(Q[v,t,wf,p,i] <= Q[v,t,wf,p,0] for p in P for i in nodes \
            for t in D for v in V for wf in WF)
        for i in nodes[1:-1]:
            for j in dropoff:
                model.addConstrs((y[v,t,wf,i,j] == 1) >> \
                    (Q[v,t,wf,p,i] - Q[v,t,wf,p,j] == self.WF_tech[wf][j-1][p]) \
                    for p in P for wf in WF for t in D for v in V)
            for j in pickup:
                model.addConstrs((y[v,t,wf,i,j] == 1) >> \
                    (Q[v,t,wf,p,j] - Q[v,t,wf,p,i] == self.WF_tech[wf][j-n-1][p]) \
                    for p in P for wf in WF for t in D for v in V)
            for j in nodes[1:]:
                if j == i+n:
                    continue    
                model.addConstrs((y[v,t,wf,i,j]==1) >> (T[v,t,wf,j] - T[v,t,wf,i] >=\
                    self.V_transTime[v] + travelTime[v][wf][i,j]) for wf in WF for t in D \
                    for v in V)

        model.addConstrs(quicksum(t * y[v,t,wf,i,j] - x[wf,j] for i in nodes for v in V \
            for t in D) <= self.WF_deadline[wf][j-1] for wf in WF for j in dropoff)


        model.addConstrs(y[v,t,wf,j,2*n+1] == 0 for j in dropoff for wf in WF \
            for t in D for v in V)
        model.addConstrs(y[v,t,wf,n+i,i] == 0 for i in dropoff for wf in WF \
            for t in D for v in V)
        model.addConstrs(y[v,t,wf,0,j] == 0 for j in pickup for wf in WF \
            for t in D for v in V)
        model.addConstrs(y[v,t,wf,i,j] <= self.lambda_bf[b][wf] for i in nodes \
            for j in nodes for t in D for v in V for wf in WF for b in B)
        model.addConstrs(y[v,t,wf,i,j] <= self.B_vessel[b][v] for i in nodes \
            for j in nodes for t in D for wf in WF for b in B for v in V)

        c_qr = quicksum(Q[v,t,wf,p,0] * self.tech_cost[p] for p in P \
            for wf in WF for t in D for v in V)
        c_tr = quicksum(self.V_cost[v] * travelTime[v][wf][i,j] * y[v,t,wf,i,j] \
            for i in nodes for j in nodes for wf in WF for t in D for v in V)
        c_lr = quicksum(x[wf,j] * self.WF_penCost[wf][j-1] for j in dropoff for wf in WF)
        
        model.modelSense = GRB.MINIMIZE
        model.setObjective(c_qr+c_tr+c_lr)
        model.optimize()
        if model.status == GRB.OPTIMAL:
            c_total = model.objVal
            print (c_total)
        else:
            print('infeasible')
        return
