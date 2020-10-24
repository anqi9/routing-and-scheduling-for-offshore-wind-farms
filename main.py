# -*- coding: utf-8 -*-
import time
from GB import GB

def main():
    data = create_data_model()
    # to create an optimization model
    opt_GB = GB(data)
    # to solve the model using Gurobi solver
    opt_GB.Solver()

if __name__ == "__main__":
    main()
