# -*- coding: utf-8 -*-
import time
from GB2 import GB2

def create_data_model():
    """Stores the data for the problem."""
    data = {}
    data['period'] = 3
    # data['time_window'] = [
    #     [[6, 6, 0], [6, 6, 0], [12, 12, 0]],
    #     [[12, 12, 0], [12, 12, 0], [12, 12, 0]],
    #     [[0, 7, 7], [0, 7, 7], [0, 12, 12]],
    #     [[0, 12, 12], [0, 12, 12], [0, 12, 12]]
    # ]
    data['time_window'] = [
        [[6, 6, 0], [6, 6, 0], [12, 12, 0]],
        [[12, 12, 0], [0, 12, 12], [0, 12, 12]],
        [[12, 12, 0], [0, 12, 12], [0, 12, 12]],
        [[0, 9, 9], [0, 9, 9], [0, 12, 12]],
        [[0, 12, 12], [0, 12, 12], [0, 12, 12]]
    ]
    data['eta'] = 5
    data['lambda'] = [[1, 1, 0], [0, 1, 1]]
    # data['lambda'] = [[1]]
    data['technician'] = [300, 325, 350]
    data['base'] = {}
    # data['base']['technician'] = [
    #     [[6, 6, 6], [6, 6, 6], [6, 6, 6]],
    #     [[6, 6, 6], [6, 6, 6], [6, 6, 6]]
    # ]
    data['base']['technician'] = [
        [[12,12,12], [12,12,12], [12,12,12]],
        [[12,12,12], [12,12,12], [12,12,12]]
    ]
    data['base']['vessel'] = [[1, 1, 1, 0, 0], [0, 0, 0, 1, 1]]
    data['wind_farm'] = {}
    data['wind_farm']['maintenance_time'] = [
        [4, 3, 5, 2, 4, 5, 2, 2],
        [3, 2, 4, 2, 3, 4, 4, 3],
        [3, 4, 3, 2, 3, 2, 2, 5]
    ]
    data['wind_farm']['technician'] = [
        [[2,0,1], [0,1,1], [3,0,0], [1,0,2], [1,2,2], [3,0,0], [2,2,1], [3,0,1]],
        [[0,1,1], [0,1,1], [1,3,0], [0,2,1], [3,2,0], [1,2,0], [3,0,1], [3,1,0]],
        [[3,1,0], [3,0,1], [0,2,2], [0,3,0], [1,0,2], [0,1,1], [0,0,2], [3,2,0]]
    ]
    data['wind_farm']['parts_weight'] = [
        [700, 700, 300, 900, 600, 900, 900, 500],
        [400, 600, 800, 700, 600, 800, 400, 800],
        [700, 800, 800, 700, 600, 800, 700, 300]
    ]
    data['wind_farm']['present'] = [
        [0, 0, 1, 0, 0, 1, 0, 0],
        [0, 0, 1, 0, 0, 0, 1, 0],
        [0, 1, 0, 0, 0, 0, 0, 1]
    ]
    data['wind_farm']['deadline'] = [
        [3, 2, 4, 1, 1, 1, 4, 1],
        [1, 3, 2, 1, 2, 3, 3, 2],
        [1, 1, 4, 2, 1, 4, 3, 2]
    ]
    data['wind_farm']['penalty_cost'] = [
        [1900, 1500, 1600, 1900, 1200, 1600, 1800, 1100],
        [1300, 1500, 1400, 1900, 1900, 1600, 1500, 1800],
        [1300, 1500, 1900, 1200, 1600, 1900, 1500, 1900]
    ]
    data['wind_farm']['avail_parts'] = [
        [[1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1]],
        [[1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1]],
        [[1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1], [1,1,1]]
    ]
    data['vessel'] = {}
    # tonnes
    # data['vessel']['capacity'] = [1.5, 26, 2, 15]
    # data['vessel']['capacity'] = [12, 26, 26, 20, 15]
    data['vessel']['capacity'] = [20, 26, 20, 26, 26]
    data['vessel']['capacity'] = [i * 1000 for i in data['vessel']['capacity']]
    data['vessel']['to_turbine'] = [
        [[1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]],
        [[1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]],
        [[1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]],
        [[1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1]]
    ]
    data['vessel']['technician'] = [12, 12, 12, 12, 12]
    # data['vessel']['cost'] = [225, 300, 250, 280]
    data['vessel']['cost'] = [280, 300, 270, 300 ,290]
    # speed (/knots) randomly generated
    # data['vessel']['speed'] = [25, 17, 23, 19]
    data['vessel']['speed'] = [20, 17, 22, 18, 17]
    # knots to km/h
    data['vessel']['speed'] = [i * 1.852 for i in data['vessel']['speed']]
    data['vessel']['trans_time'] = [0.25, 0.25, 0.25, 0.25, 0.25]
    data['vessel']['availability'] = [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]
    # random generated coordinates
    data['base']['coordinate'] = [[20, 30], [60, 63]]
    data['wind_farm']['coordinate'] = [
        [[3, 7], [1, 7], [1, 0], [1, 1], [4, 1], [5, 2], [6, 1], [0, 5]],
        [[26, 33], [20, 31], [24, 30], [25, 35], [22, 32], [23, 36], [25, 35], [21, 35]],
        [[51, 32], [50, 35], [52, 35], [50, 34], [48, 35], [52, 38], [45, 37], [49, 34]]
    ]
    return data

def main():
    data = create_data_model()
    # T1 = time.perf_counter()
    opt_GB = GB2(data)
    # opt_GB.createDistanceMatrix()
    opt_GB.Solver()
    # [cost, routes] = opt_GB.Solver()
    # T2 = time.perf_counter()
    # print ('time: ', round(T2-T1, 2))


if __name__ == "__main__":
    main()
