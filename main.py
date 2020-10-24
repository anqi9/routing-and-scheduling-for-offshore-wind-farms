# -*- coding: utf-8 -*-
import time
from GB import GB

def main():
    data = create_data_model()
    opt_GB = GB2(data)
    opt_GB.Solver()


if __name__ == "__main__":
    main()
