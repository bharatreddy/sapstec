import pandas
import datetime
import math
import os
import numpy
from scipy import signal, ndimage, stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import Normalize
from matplotlib import ticker
import seaborn as sns
import trough_detect
from aacgmv2 import convert_mlt

if __name__ == "__main__":
    # Base directory where all the files are stored
    baseDir = "/home/bharat/Documents/code/new-vel-data/veldata/"
    # read data from a file
    for root, dirs, files in os.walk(baseDir):
        for fNum, fName in enumerate(files):
            currInpLosFile = root + fName    
            print "working with--->", currInpLosFile, fNum,"/",len(files)
            trghObj = trough_detect.TroughBnd( currInpLosFile )
            allTimesList = trghObj.get_all_uniq_times()
            fileNameFltrdBnds = baseDir + "/trghBnds/tec-fltrd-vals-" +\
                                 allTimesList[0].year + "--" +\
                                 allTimesList[0].month + "--" +\
                                 allTimesList[0].day + ".txt"
            fileNameBndCoeffs = baseDir + "/trghCoeffs/bnd-coeffs" +\
                                 allTimesList[0].year + "--" +\
                                 allTimesList[0].month + "--" +\
                                 allTimesList[0].day + ".txt"
            for inpDT in allTimesList:
                trghBndDF = trghObj.find_trough_loc(inpDT)
                if trghBndDF is not None:
                    fltrdTrghBndDF = trghObj.filter_trough_loc( trghBndDF )
                if fltrdTrghBndDF is not None:
                    bndCoeffsDF = trghObj.fit_boundaries( fltrdTrghBndDF )
                with open(fileNameFltrdBnds, 'a') as ftB:
                    fltrdTrghBndDF.to_csv(ftB, header=False,\
                                      index=False, sep=' ' )
                    with open(fileNameBndCoeffs, 'a') as bcF:
                    bndCoeffsDF.to_csv(bcF, header=False,\
                                      index=False, sep=' ' )