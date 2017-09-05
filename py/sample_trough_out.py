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

    inpTecFile = "/home/bharat/Documents/code/data/fullday-tecMF/tec-medFilt-20110126.txt"
    trghObj = trough_detect.TroughBnd( inpTecFile )
    allTimesList = trghObj.get_all_uniq_times()
    dtDay = datetime.datetime.utcfromtimestamp(\
            allTimesList[0].astype(int) * 1e-9)
    fileNameFltrdBnds = "/home/bharat/Documents/code/sapstec/sample-bnds.txt"
    fileNameBndCoeffs = "/home/bharat/Documents/code/sapstec/sample-coeffs.txt"
    for inpDT in allTimesList:
        currDT = datetime.datetime.utcfromtimestamp(\
            inpDT.astype(int) * 1e-9)
        trghBndDF = trghObj.find_trough_loc(currDT)
        fltrdTrghBndDF = None
        if trghBndDF is not None:
            fltrdTrghBndDF = trghObj.filter_trough_loc( trghBndDF )
        if fltrdTrghBndDF is not None:
            bndCoeffsDF = trghObj.fit_boundaries( fltrdTrghBndDF )
            with open(fileNameFltrdBnds, 'a') as ftB:
                fltrdTrghBndDF.to_csv(ftB, header=True,\
                                  index=False, sep=' ' )
            with open(fileNameBndCoeffs, 'a') as bcF:
                bndCoeffsDF.to_csv(bcF, header=True,\
                                  index=False, sep=' ' )