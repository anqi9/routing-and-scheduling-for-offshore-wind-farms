# -*- coding: utf-8 -*-
import numpy as np
from sys import exit
from gurobipy import Model, GRB, quicksum
import copy
from itertools import combinations

''''
this class is to solve routing&scheduling using Dantzig-Wolf decomposition
'''
class MasterProblem():
    def __init__(self, data):
        self.period = data['period']  # days: 3-7
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
        # self.WF_availParts = data['wind_farm']['avail_parts']
        # spare parts capacity (/tonne) of vessel
        self.V_parts = data['vessel']['capacity']
        # indicate whether vessel is able to transfer spare parts needed by the turbine
        # self.V_toTurbine = data['vessel']['to_turbine']
        # max num. of technician on board
        self.V_tech = data['vessel']['technician']
        # fuel cost (/hour) of vessels
        self.V_cost = data['vessel']['cost']
        # speed (km/hour) of vessels
        self.V_speed = data['vessel']['speed']
        # transfer time for technicians and equipment from vessel to turbine
        self.V_transTime = data['vessel']['trans_time']
        # available working hours based on the weather and vessels type (/hour)
        self.time_window = data['time_window']
        # indicate whether vessel available at the period
        self.V_available = data['vessel']['availability']
        # num of wind farms
        self.num_WF = len(self.WF_tech)
        # num of bases
        self.num_B = len(self.B_tech)
        # num of vessels
        self.num_V = len(self.time_window)
        # types of available technicians
        self.types = list(range(len(self.tech_cost)))
        # possible subtours with max length eta
        self.subtours = self.findsubsets()
   
    # to calculate distance including a base
    def createDistanceMatrix(self, coordinate, base):
        Num = len(coordinate) + 1
        d = np.zeros((Num, Num))
        for i in range(Num-1):
            p1 = np.array([coordinate[i][0], coordinate[i][1]])
            p3 = np.array([base[0], base[1]])
            d[0, i+1] = d[i+1, 0] = np.linalg.norm(p1 - p3)
            for j in range(i + 1, Num-1):
                p2 = np.array([coordinate[j][0], coordinate[j][1]])
                d[i+1, j+1] = d[j+1, i+1] = np.linalg.norm(p1 - p2)
        # d = np.round(d, 2)
        return d  # type: np.array

    # to generate all possible subsets with length m
    def findsubsets(self):
        subset = [[] for i in range(self.num_WF)]
        for wf in range(self.num_WF):
            for i in range(2, self.eta + 1):
                S = list(range(len(self.WF_coordinate[wf])))
                sub = list(combinations(S, i))
                sub = [list(k) for k in sub]
                subset[wf].extend(sorted(sub))
        return subset

    def Mater_Problem(self):
        [R_turbine, routes, NumTechs] = self.Generating_Feasible_Routes()
        model = Model()
        
        # set parameters
        A = []
        for v in range(self.num_V):
            for t in range(self.period):
                for b in range(self.num_B):
                    for wf in range(self.num_WF):
                        if R_turbine[v][t][b][wf] == []:
                            continue
                        for r in range(np.size(R_turbine[v][t][b][wf][0],1)):
                            A.append((v,t,b,wf,r))
        x = model.addVars(A, vtype=GRB.BINARY)
        
        # add constraints
        for v in range(self.num_V):
            for t in range(self.period):
                tmp1 = 0
                for b in range(self.num_B):
                    for wf in range(self.num_WF):
                        if R_turbine[v][t][b][wf] == []:
                            continue
                        ntmp = np.size(R_turbine[v][t][b][wf][0], 1)
                        tmp1 += sum(x[v,t,b,wf,r] for r in range(ntmp))
                model.addConstr(tmp1 <= 1)

        # explore every turbines in the farm
        for b in range(self.num_B):
            for wf in range(self.num_WF):
                for j in range(len(self.WF_tech[wf])):
                    tmp = 0
                    for v in range(self.num_V):
                        for t in range(self.period):
                            if R_turbine[v][t][b][wf] == []:
                                continue
                            ntmp = np.size(R_turbine[v][t][b][wf][0], 1)
                            tmp += quicksum(x[v,t,b,wf,r] * R_turbine[v][t][b][wf][0][j,r] \
                                    for r in range(ntmp))
                    model.addConstr(tmp == 1)

        for b in range(self.num_B):
            for t in range(self.period):
                for p in self.types:
                    tmp = 0
                    for v in range(self.num_V):
                        for wf in range(self.num_WF):
                            if R_turbine[v][t][b][wf] == []:
                                continue
                            tmp += quicksum(x[v,t,b,wf,r] * NumTechs[v][t][b][wf][r][p] \
                                for r in range(np.size(R_turbine[v][t][b][wf][0],1)))
                    model.addConstr(tmp <= self.B_tech[b][p][t])
       
        ctmp = 0
        for v in range(self.num_V):
            for t in range(self.period):
                for b in range(self.num_B):
                    for wf in range(self.num_WF):
                        if R_turbine[v][t][b][wf] == []:
                                continue
                        ctmp += quicksum(R_turbine[v][t][b][wf][0][-1,r] * x[v,t,b,wf,r] \
                            for r in range(np.size(R_turbine[v][t][b][wf][0],1)))
        
        model.modelSense = GRB.MINIMIZE
        model.setParam('OutputFlag', 0)
        model.setObjective(ctmp)
        model.optimize()
        
        _routes = [[] for i in range(self.num_V)]
        if model.status == GRB.OPTIMAL:
            cost_opt = model.objVal

            for b in range(self.num_B):
                for wf in range(self.num_WF):
                    if R_turbine[v][t][b][wf] == []:
                            continue
                    turbines = []
                    for v in range(self.num_V):
                        for t in range(self.period):
                            # if R_turbine[v][t][b][wf] == []:
                            #     continue
                            if sum(x[v,t,b,wf,r].x for r in range(np.size(R_turbine[v][t][b][wf][0],1))) < 1:
                                _routes[v].append([])
                                continue
                            for r in range(np.size(R_turbine[v][t][b][wf][0],1)):
                                if x[v,t,b,wf,r].x == 1:
#                                     print ('(b,wf,v,t)[%d,%d,%d,%d]' % (b,wf,v,t),\
#                                         'routes: ', routes[v][t][b][wf][r])
                                    _routes[v].append(routes[v][t][b][wf][r])
            return [cost_opt, _routes]
        else: 
            print('Master Problem: Infeasible')
            return [None, None, None]

    def Generating_Feasible_Routes(self):
        np.set_printoptions(suppress=True)
        # store all feasible sets of turbines including cost at corresponding period
        Routes_turbine = [[[[[] for i in range(self.num_WF)] for j in range(self.num_B)] for k in range(self.period)] for m in range(self.num_V)]
        # store all feasible route at corresponding period
        Routes = [[[[[] for i in range(self.num_WF)] for j in range(self.num_B)] for k in range(self.period)] for m in range(self.num_V)]
        # store corresponding num of technicians required by the route
        NumTechs = [[[[[] for i in range(self.num_WF)] for j in range(self.num_B)] for k in range(self.period)] for m in range(self.num_V)]
        visited = [[[[[] for i in range(self.num_WF)] for j in range(self.num_B)] for k in range(self.period)] for m in range(self.num_V)]

        for vessel in range(self.num_V):
            for base in range(self.num_B):
                if self.B_vessel[base][vessel] == 0:
                    continue
                for farm in range(self.num_WF):
                    # if self.lambda_bf[base][farm] == 0:
                    #     continue
                    # calculate the distance
                    dist = self.createDistanceMatrix(self.WF_coordinate[farm], self.B_coordinate[base])
                    # travel time between turbines
                    travelTime = dist/self.V_speed[vessel]
                    numTurbine = len(self.WF_coordinate[farm])
                    for t in range(self.period):
                        # feasibility condition
                        if self.V_available[vessel][t] == 0: 
                            continue
                        ''' 
                        rows: turbines + cost, columns: feasible routes,
                        R_turbine[i,p] = 1 if route p visits i
                        last 4 rows: travel, personnel, penalty and total cost
                        to store feasible routes for R_turbine at one period
                        '''
                        R_turbine = np.zeros((numTurbine + 4, 1))
                        for turbine in range(numTurbine):
                            # feasibility conditions
                            if self.V_toTurbine[vessel][farm][turbine] == 0 or self.WF_availParts[farm][turbine][t] == 0: 
                            #     continue
                            # total travel time
                            _travle_time = 2*travelTime[0,turbine]+self.WF_mainTime[farm][turbine] + self.V_transTime[vessel]
                            if _travle_time >= self.time_window[vessel][t][farm]:
                                continue
                            # base -> turbine -> base
                            newRoute = np.zeros((numTurbine + 4, 1))
                            r_start = np.zeros((self.num_B, 1))
                            # travel cost
                            newRoute[-4] = 2*(self.V_cost[vessel]*travelTime[0,turbine])
                            # personnel cost
                            newRoute[-3] = round(sum(self.WF_tech[farm][turbine][p] * self.tech_cost[p] for p in self.types),2)
                            # penalty cost
                            newRoute[-2] = max(0, t - self.WF_deadline[farm][turbine] * self.WF_penCost[farm][turbine])
                            # total cost
                            newRoute[-1] = sum(newRoute[-4:-1,0])
                            newRoute[:-4,0] = np.round(newRoute[:-4,0])
                            newRoute[-4:,0] = np.round(newRoute[-4:,0],2)
                            newRoute[turbine] = 1
                            r_start[base] = 1
                            R_turbine = np.append(R_turbine, newRoute, axis = 1)
                            Routes[vessel][t][base][farm].append([turbine,turbine])
                            numTechs = [self.WF_tech[farm][turbine][p] for p in self.types]
                            NumTechs[vessel][t][base][farm].append(numTechs)
                            visited[vessel][t][base][farm].append([turbine])

                        for L in self.subtours[farm]:
                            if sum(self.WF_partsWeight[farm][tt] for tt in L) > \
                                    self.V_parts[vessel]:
                                    continue
                            if sum(self.WF_mainTime[farm][tt] for tt in L) >=\
                                    self.time_window[vessel][t][farm]:
                                    continue
                            newRoute = np.zeros((numTurbine + 4, 1))
                            # Use MILP solve subproblem to get the min cost
                            [newRoute[-4:,0], route_tmp, numTechs] = \
                            self.Solving_MILP(L, dist, travelTime, vessel, t, \
                                base, farm, Routes_turbine, Routes, NumTechs, visited)
                            if route_tmp == None:
                                continue
                            visited[vessel][t][base][farm].append(L)
                            newRoute[L] = 1
                            newRoute = np.round(newRoute)
                            R_turbine = np.append(R_turbine, newRoute, axis = 1)
                            r_start = np.zeros((self.num_B, 1))
                            r_start[base] = 1
                            Routes[vessel][t][base][farm].append(route_tmp)
                            NumTechs[vessel][t][base][farm].append(numTechs)

                        Routes_turbine[vessel][t][base][farm].append(R_turbine[:,1:])
        return [Routes_turbine, Routes, NumTechs]

    def Solving_MILP(self, turbines, dist, travelTime, vessel, t, base, farm, R, routes, Techs, visited):

        for t_1 in range(1, t):
            for b_1 in range(self.num_B):
                if turbines not in visited[vessel][t_1][b_1][farm]:
                    continue
                index = visited[vessel][t_1][b_1][farm].index(turbines)
                if self.time_window[vessel][t][farm] == self.time_window[vessel][t_1][farm]:
                    c_t = R[vessel][t_1][b_1][farm][0][-4,index]
                    c_p = R[vessel][t_1][b_1][farm][0][-3,index]
                    c_lr = sum(max(0, t - self.WF_deadline[farm][j]) \
                        * self.WF_penCost[farm][j] for j in turbines)
                    c_total = c_t + c_p + c_lr
                    path = routes[vessel][t_1][b_1][farm][index]
                    NumTechs = Techs[vessel][t_1][b_1][farm][index]
                    return [[c_t, c_p, c_lr, c_total], path, NumTechs]

        # total num of turbines in the wind farm
        n = len(self.WF_deadline[farm])
        # set of drop off nodes [1,...,n]
        dropoff = list(map(lambda qq: qq+1, turbines))
        d1 = list(range(n+1, 2*n+1))
        # set of pick up nodes [n+1, n+2,...,2n]
        pickup = [d1[i-1] for i in dropoff]
        # set of nodes need vessels to be present
        waiting = [i+1 for i in turbines if self.WF_present[farm][i] == 1]
        nodes = [0] + dropoff + pickup + [2*n+1]
        
        model = Model()
        T = model.addVars(nodes, lb=0.0, vtype=GRB.CONTINUOUS)
        A = [(i,j) for i in nodes for j in nodes if i!=j]
        B = [(p,i) for p in self.types for i in nodes]
        y = model.addVars(A, vtype=GRB.BINARY)
        Q = model.addVars(B, lb=0, vtype=GRB.INTEGER)

        model.addConstr(y[2*n+1,0] == 1)
        model.addConstrs(quicksum(y[i,j] for j in nodes if i!=j) == 1 for i in nodes) #9
        model.addConstr(quicksum(y[0,j] for j in dropoff) == 1) #10
        model.addConstr(quicksum(y[j,2*n+1] for j in pickup) == 1) #11
        model.addConstrs(quicksum(y[i,j] for j in nodes if i!=j) == \
            quicksum(y[j,i] for j in nodes if i!=j) for i in nodes) #12
        model.addConstrs(quicksum(y[i,j] for j in nodes if i!=j) ==\
            quicksum(y[n+i,j] for j in nodes if n+i!=j) for i in dropoff) #13
        model.addConstrs(y[i,n+i] == 1 for i in waiting) #14
        model.addConstrs(T[n+i] - T[i] >= self.WF_mainTime[farm][i-1] \
                + self.V_transTime[vessel] for i in dropoff) #15
        model.addConstr(T[2*n+1] <= self.time_window[vessel][t][farm]) #16
        model.addConstr(T[0] == 0) #17
        model.addConstrs(y[0,j]==0 for j in pickup)
        model.addConstrs(y[j,2*n+1]==0 for j in dropoff)
        model.addConstr(quicksum(Q[p,0] for p in self.types) <= self.V_tech[vessel]) #18
        model.addConstrs(Q[p,0] <= self.B_tech[base][p][t] for p in self.types)
        model.addConstrs(Q[p,0] >= Q[p,i] for p in self.types for i in nodes)

        for i in nodes[1:-1]:
            for j in dropoff:
                if i == j:
                    continue
                model.addConstrs((y[i,j]==1) >> \
                    (Q[p,i] - Q[p,j] == self.WF_tech[farm][j-1][p]) for p in self.types)
            for j in pickup:
                if i == j:
                    continue
                model.addConstrs((y[i,j]==1) >> \
                    (Q[p,j] - Q[p,i] == self.WF_tech[farm][j-n-1][p]) for p in self.types)
      
        model.addConstrs((y[0,j]==1) >> \
                    (T[0] + travelTime[0,j] <= T[j]) for j in dropoff)
        model.addConstrs((y[j,2*n+1]==1) >> (T[j] + self.V_transTime[vessel] \
            + travelTime[j-n,0] <= T[2*n+1]) for j in pickup)
        
        for i in nodes[1:-1]:
            for j in nodes[1:-1]:
                if j == i+n or i==j:
                    continue
                i1 = i if i<=n else i-n
                j1 = j if j<=n else j-n
                model.addConstr((y[i,j]==1) >> \
                    (T[i]+ self.V_transTime[vessel] + travelTime[i1,j1] <= T[j])) #19
            
        # cost for technicians
        c_qr = sum(Q[p,0] * self.tech_cost[p] for p in self.types)
        # travel cost (fuel cost)
        c_tr = 0
        for i in nodes:
            for j in nodes:
                if i != j:
                    if i == 2*n + 1:
                        i1 = 0
                    elif i <= n:
                        i1 = i
                    else:
                        i1 = i - n
                        i1 = i if i <= n else i - n
                    if j == 2*n + 1:
                        j1 = 0
                    elif j <= n:
                        j1 = j
                    else:
                        j1 = j - n
                    c_tr += self.V_cost[vessel] * travelTime[i1,j1] * y[i,j]


        # penalty cost
        c_lr = sum(max(0, t+1 - self.WF_deadline[farm][j]) \
            * self.WF_penCost[farm][j] for j in turbines)
        # ojective function
        model.modelSense = GRB.MINIMIZE
        model.setObjective(c_qr+c_tr+c_lr)
        model.setParam('OutputFlag', 0)
        model.optimize()

        if model.status == GRB.OPTIMAL:
            print('optimal')
            c_total = model.objVal
            # cost for technicians
            c_q = sum(Q[p,0].x * self.tech_cost[p] for p in self.types)
            # travel cost (fuel cost)
            c_t = 0
            for i in nodes:
                for j in nodes:
                    if i != j:
                        if i == 2*n + 1:
                            i1 = 0
                        elif i <= n:
                            i1 = i
                        else:
                            i1 = i - n
                            i1 = i if i <= n else i - n
                        if j == 2*n + 1:
                            j1 = 0
                        elif j <= n:
                            j1 = j
                        else:
                            j1 = j - n
                        c_t += self.V_cost[vessel] * travelTime[i1,j1] * y[i,j].x

            path = [] # store the optimal route without base
            TT = [0]
            rr = 0
            ii = 0
            while rr < len(nodes) - 2:
                for k in nodes:
                    if ii!= k and round(y[ii,k].x) == 1:
                        ii = k
                        path.append(k-1 if k <= n else k-n-1)
                        TT.append(round(T[k].x, 2))
                        break
                rr += 1
            TT.append(T[nodes[-1]].x)
            [c_t, c_q, c_lr, c_total] = [round(i,2) for i in [c_t, c_q, c_lr, c_total]]
            num_techs = [Q[p,0].x for p in self.types]
#             print('path: ', ['vessel[%d]' % vessel] + ['period[%d]' % t] + \
#                  ['wf[%d]' % farm] + ['b[%d]' % base] + path)
            return [[c_t, c_q, c_lr, c_total], path, num_techs]
        else:
            print('infeasible')
            return [[None], None, None]
