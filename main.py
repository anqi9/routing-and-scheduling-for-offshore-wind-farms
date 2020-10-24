# -*- coding: utf-8 -*-
import time
from GB import GB
from MasterProblem import MasterProblem

def main():
    data = create_data_model()
    '''using DW-decomposition to solve the routing & scheduling problem'''
    opt = MasterProblem(data)
    [cost, routes] = opt.Mater_Problem()
   
    ''' to solve the routing&scheduling problem directly '''
    # to create an optimization model
    opt_GB = GB(data)
    # to solve the model using Gurobi solver
    opt_GB.Solver()

if __name__ == "__main__":
    main()
