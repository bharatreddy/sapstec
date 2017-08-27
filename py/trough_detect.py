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

if __name__ == "__main__":
    trghObj = trough_detect.TroughBnd( \
        "../data/tec-medFilt-20110409.txt" )
    allTimesList = trghObj.get_all_uniq_times()
    # print allTimesList
    inpDT = datetime.datetime( 2011, 4, 9, 12, 0 )
    trghBndDF = trghObj.find_trough_loc(inpDT)
    print trghBndDF
    fltrdTrghBndDF = trghObj.filter_trough_loc( trghBndDF )
    print fltrdTrghBndDF


class TroughBnd(object):
    """
    A class to identify the location of TEC trough
    given a time instance. We use Evan's median filtered
    TEC data!!!
    """
    def __init__(self, tecFileName):
        import datetime
        import numpy
        import pandas
        # Some constants
        self.equTrghCutoffMLat = 40.
        self.polTrghCutoffMLat = 70.
        # We'll only look at NA continent, So we 
        # choose the corresponding MLONS - 260, 20
        self.mlonNAList = range( -100, 21, 2 )#range( -180, 181, 2 )
        # Some cutoffs to verify goodness of fit
        self.cutoffKSPval = 0.7
        self.cutoffKSDstat = 0.25
        self.cutOffTrghMinUpper = 65.
        self.cutOffTrghMinLower = 45.
        self.cutOffPrcntErrorFit = 0.2 # 20%
        self.cutOffLatCnt = 20
        self.cutOffMinTECVal = 5.
        self.cutoffUnqMlonCnt = 15.
        # set variables for trough location filtering
        self.trghLocMlonbinSize = 20
        self.trghLocMlonbinCntCutoff = 0.2
        self.cutOffGoodMlonBins = 2
        # choose columns we need to store from the tec file
        tecFileselCols = [ "dateStr", "timeStr", "Mlat",\
              "Mlon", "med_tec", "dlat", "dlon" ]
        # read the tec file for the given datetime
        self.medFltrdTecDF = pandas.read_csv(tecFileName, delim_whitespace=True,\
                                    header=None, names=tecFileselCols)
        self.medFltrdTecDF["date"] = self.medFltrdTecDF.apply(\
         self.convert_to_datetime, axis=1 )# read the tec file for the given datetime
        self.medFltrdTecDF = pandas.read_csv(tecFileName, delim_whitespace=True,\
                                    header=None, names=tecFileselCols)
        self.medFltrdTecDF["date"] = self.medFltrdTecDF.apply(\
         self.convert_to_datetime, axis=1 )


    def convert_to_datetime(self,row):
        """
        Convert date and time strings to datetime objects
        """
        import datetime
        currDateStr = str( int( row["dateStr"] ) )
        if row["timeStr"] < 10:
            currTimeStr = "000" + str( int( row["timeStr"] ) )
        elif row["timeStr"] < 100:
            currTimeStr = "00" + str( int( row["timeStr"] ) )
        elif row["timeStr"] < 1000:
            currTimeStr = "0" + str( int( row["timeStr"] ) )
        else:
            currTimeStr = str( int( row["timeStr"] ) )
        return datetime.datetime.strptime( currDateStr\
                        + ":" + currTimeStr, "%Y%m%d:%H%M" )

    def get_all_uniq_times(self):
        """
        Get all the unique times (dt objs) present
        in the inp TEC file!!!
        """
        return self.medFltrdTecDF["date"].unique()

    def gauss_function(self, x, a, x0, sigma):
        # Fit a gaussian curve to find trough
        # across a longitude.
        return a*numpy.exp(-(x-x0)**2/(2*sigma**2))

    def find_trough_loc(self, selDT):
        """
        From the median filtered tec dataframe
        identify the location of trough
        """
        BndMlonArr = []
        BndEquMlatArr = []
        BndPolMlatArr = []
        minTecMlatArr = []
        minTecValArr = []
        BndEquTecValArr = []
        BndPolTecValArr = []
        currTimeArr = []
        selDTDF = self.medFltrdTecDF[ \
                    self.medFltrdTecDF["date"] ==selDT\
                     ].reset_index(drop=True)
        for ind, sMlon in enumerate(self.mlonNAList):
            if sMlon < 0.:
                sMlon = sMlon + 360.
            selMlonDF = selDTDF[ (selDTDF["Mlon"] == sMlon) &\
                                  (selDTDF["Mlat"] >= self.equTrghCutoffMLat) &\
                                  (selDTDF["Mlat"] <= self.polTrghCutoffMLat)]
            # If no significant number of values are found discard
            if len( selMlonDF["med_tec"].values ) < 5.:
                continue
            tecGaussFitArr = numpy.max( selMlonDF["med_tec"].values ) - selMlonDF["med_tec"].values
            mlatPltArr = numpy.arange(self.equTrghCutoffMLat-20, self.polTrghCutoffMLat+20)
            try:
                popt, pcov = curve_fit(self.gauss_function, selMlonDF["Mlat"].values,\
                                       tecGaussFitArr, p0 = [2, 52., 1.])
            except:
                print "failed fit at MLON-->", sMlon
                continue
            fwhmEqu = popt[1] - popt[2]*2.355/2.
            fwhmPol = popt[1] + popt[2]*2.355/2.
            tecValFwhmEqu = self.gauss_function(fwhmEqu, *popt)
            tecValFwhmPol = self.gauss_function(fwhmPol, *popt)
            # Test goodness of fit
            ksTestTecArr = numpy.array( [ self.gauss_function(l, *popt)\
                                        for l in selMlonDF["Mlat"].values.tolist() ] )
            ksDStat, ksPVal = stats.ks_2samp( tecGaussFitArr, ksTestTecArr )
            # we setup a few conditions to discard bad fits
            # 1) location of trough min should be between 45 nad 65 MLAT
            if ((popt[1] > self.cutOffTrghMinUpper)\
                    | (popt[1] < self.cutOffTrghMinLower)):
                continue
            # 2) if percent error in any of fit parameters is more 
            # than 10 % (cutoff) then skip
            if ( ( pcov[0,0]**0.5/popt[0] > self.cutOffPrcntErrorFit  ) \
                | ( pcov[1,1]**0.5/popt[1] > self.cutOffPrcntErrorFit  ) \
                | ( pcov[2,2]**0.5/popt[2] > self.cutOffPrcntErrorFit  ) ):
                continue
            # 3) Number of latitudes should be greater than 20.
            if len( selMlonDF["Mlat"].values ) < self.cutOffLatCnt:
                continue
            # 4) If p-val from the KS TEST is low discard
            # or if KS Stat is high discard
            if ksPVal < self.cutoffKSPval:
                continue
            if ksDStat > self.cutoffKSDstat:
                continue
            # 5) if tec value at min trough is greater than 4.
            # discard
            # GET Trough min loc and tec val
            minTrghLoc = min(list(selMlonDF["Mlat"].values), key=lambda x:abs(x-popt[1]))
            minTrghTecVal = selMlonDF[ selMlonDF["Mlat"] == minTrghLoc ]["med_tec"].values[0]
            # GET Trough equ loc and tec val
            equTrghLoc = min(list(selMlonDF["Mlat"].values), key=lambda x:abs(x-fwhmEqu))
            equTrghTecVal = selMlonDF[ selMlonDF["Mlat"] == equTrghLoc ]["med_tec"].values[0]
            # GET Trough pol loc and tec val
            polTrghLoc = min(list(selMlonDF["Mlat"].values), key=lambda x:abs(x-fwhmPol))
            polTrghTecVal = selMlonDF[ selMlonDF["Mlat"] == polTrghLoc ]["med_tec"].values[0]
            if minTrghTecVal > self.cutOffMinTECVal:
                continue
            # Now we have good fits. Get appropriate boundary locations
            # and TEC values.
            # Append the values to arrays
            BndMlonArr.append( sMlon )
            BndEquMlatArr.append( equTrghLoc )
            BndPolMlatArr.append( polTrghLoc )
            minTecMlatArr.append( minTrghLoc )
            minTecValArr.append( minTrghTecVal )
            BndEquTecValArr.append( equTrghTecVal )
            BndPolTecValArr.append( polTrghTecVal )
            currTimeArr.append( selDT )
        # Check the number of unique mag lon measurements
        # we have. IF they are above cutoff, discard them
        if len( BndMlonArr ) > self.cutoffUnqMlonCnt:
            # Store data in a DF
            trghBndDF = pandas.DataFrame({
                        "BndMlon" : BndMlonArr,
                        "BndEquMlat" : BndEquMlatArr,
                        "BndPolMlat" : BndPolMlatArr,
                        "minTecMlat" : minTecMlatArr,
                        "minTecVal" : minTecValArr,
                        "BndEquTecVal" : BndEquTecValArr,
                        "BndPolTecVal" : BndPolTecValArr,
                        "date" : currTimeArr
                        })
            return trghBndDF
        return None

    def filter_trough_loc(self,trghLocDF):
        """
        Once a trough location has been calculated, we'll find 
        a few unwated locations in there. Here we'll fitler them out.
        """
        # Now we need to filter out the trough locs 
        # which are present in odd locations. We do so
        # by binning the Mlon arr in a groups 10. and count
        # how many fits we have in each bin. If we have greater
        # than 50% (actual count=5) values in that bin, we keep
        # those bins and remove the rest.
        # check if longitude goes -180 to 180 or 0 to 360.
        minEdge = self.mlonNAList[0]
        maxEdge = self.mlonNAList[-1]
        binList = [ b for b in numpy.arange(minEdge,\
                maxEdge,self.trghLocMlonbinSize) ]
        # Need to adjust the negative Mlons
        trghLocDF["adjstMlons"] = [ b - 360. if b > 180. else\
                 b for b in trghLocDF["BndMlon"] ]
        mlonFreq, mlonBins = numpy.histogram(trghLocDF["adjstMlons"].values,\
                             bins=binList)
        goodMlonValues = numpy.where( mlonFreq >= \
            self.trghLocMlonbinCntCutoff*self.trghLocMlonbinSize )
        # if there aren't more than 2 mlon bins with good fits
        # discard the datapoints.
        print goodMlonValues, len(goodMlonValues)
        print mlonFreq, mlonBins
        if len(goodMlonValues[0]) >= self.cutOffGoodMlonBins:
            return trghLocDF
        else:
            return None