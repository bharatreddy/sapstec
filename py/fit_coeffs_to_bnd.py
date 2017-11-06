import pandas
import datetime
import time
import urllib
import bs4
import os
import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.colors import Normalize
from matplotlib import ticker
import seaborn as sns
from aacgmv2 import convert_mlt
import feather

# To estimate the boundary we
# fit first order harmonics!
def trough_fit_harmonic(mlt, a0, c1, s1, phiC1, phiS1):
    phiC = (2*numpy.pi/24.) * mlt + phiC1
    phiS = (2*numpy.pi/24.) * mlt + phiS1
    cosTerm = c1 * numpy.cos(phiC)
    sinTerm = s1 * numpy.sin(phiS)
    return a0 + cosTerm + sinTerm

tecMinCutoff = 10.
delTecCutoff = 0.
delMlatCutoff = 0.

# get dst index vals from wdc kyoto website
# create a list of dates with monthly freq
date_dst_arr = []
dst_val = []
dst_time_del = datetime.timedelta(hours = 1)
start_date = datetime.datetime(2011,1,1)
end_date = datetime.datetime(2014,12,31)
daterange = pandas.date_range(start_date, end_date, freq="M")
for dt in daterange:
    if dt.month <= 9:
            monthStr = "0" + str(dt.month)
    else:
        monthStr = str(dt.month)
    if dt.year >= 2016:
        # create the url
        currUrl = "http://wdc.kugi.kyoto-u.ac.jp/" + "dst_realtime" + \
            "/" + str(dt.year) + monthStr + "/index.html"
    elif ( (dt.year >= 2014) and (dt.year <= 2015) ):
        # create the url
        currUrl = "http://wdc.kugi.kyoto-u.ac.jp/" + "dst_provisional" + \
            "/" + str(dt.year) + monthStr + "/index.html"
    else:
        # create the url
        currUrl = "http://wdc.kugi.kyoto-u.ac.jp/" + "dst_final" + \
            "/" + str(dt.year) + monthStr + "/index.html"
    conn = urllib.urlopen(currUrl)
    htmlSource = conn.read()
    soup = bs4.BeautifulSoup(htmlSource, 'html.parser')
    dataResObj = soup.find("pre", { "class" : "data" })
    # get the data as a list of strings after removing white space
    lines = dataResObj.text.strip().splitlines()
    for line in lines[6:]:
        columns = line.split()
        if len( columns ) > 0. :
            date_dst_arr.append( datetime.datetime( \
                dt.year, dt.month, int(columns[0]), 1 ) )
            for cols in range( len( columns[1:] ) ) :
                try:
                    inNumberFloatTest = float(columns[cols + 1])
                except:
                    # split these cols as well and work on them!
                    try:
                        missedCols = columns[cols + 1].split("-")[1:]
                        if len(missedCols) >= 1:
                            for mcols in missedCols:
                                dst_val.append( -1*float( mcols ) )
                                # now since we added the date earlier we need to be
                                # careful about appending date values
                                if ( len(date_dst_arr) != len(dst_val) ):
                                    date_dst_arr.append ( date_dst_arr[-1] + dst_time_del )
                    except:
                        print "something wrong with messed up vals!-->", columns[cols + 1]
                        continue
                    continue
                # I have to do this because of the messed up way Kyoto puts up the latest dst value..
                # mixed with 9999 (fillers) like if latest dst is 1 then Kyoto puts it as 199999.....
                if len( columns[ cols + 1 ] ) < 5 :
                    dst_val.append( float( columns[ cols + 1 ] ) )
                elif ( len( columns[ cols + 1 ] ) > 5 and columns[ cols + 1 ][0:3] != '999' ) :
                    mixed_messed_dst = ''
                    for jj in range(5) :
                        if columns[ cols + 1 ][jj] != '9' :
                            mixed_messed_dst = mixed_messed_dst + columns[ cols + 1 ][jj]

                    if mixed_messed_dst != '-' :
                        dst_val.append( float( mixed_messed_dst ) )
                    else :
                        dst_val.append( float( 'nan' ) )
                else :
                    dst_val.append( float( 'nan' ) )
                if cols > 0 :
                    date_dst_arr.append ( date_dst_arr[-1] + dst_time_del )
# convert dst data to a dataframe
dstDF = pandas.DataFrame(
    {'dst_date': date_dst_arr,
     'dst_index': dst_val
    })
dstDF["dateStr"] = dstDF["dst_date"].map(lambda x: x.strftime('%Y%m%d'))
dstDF["hour"] = dstDF["dst_date"].map(lambda x: x.strftime('%H'))

baseDir = "/home/bharat/Documents/code/data/trghBnds/"
fitDir = "/home/bharat/Documents/code/data/trghCoeffs/"
colNames = [ "mlatEqu", "tecEqu", "mlon",\
            "mlatPol", "tecPol", "date",\
            "mlatMin", "tecMin", "mlt", "mlonAdjst" ]
frames = []
# cnt = 0
for root, dirs, files in os.walk(baseDir):
    for fNum, fName in enumerate(files):
        currInpLosFile = root + fName
        bndDF = pandas.read_csv(currInpLosFile, delim_whitespace=True,\
                                    header=None, names=colNames,\
                                infer_datetime_format=True,\
                                parse_dates=["date"])
        frames.append( bndDF )

finBndDF = pandas.concat( frames )
finBndDF["normMLT"] = [x-24 if x >= 12 else x\
                         for x in finBndDF['mlt']]
finBndDF["delTecEqu"] = finBndDF["tecEqu"] - finBndDF["tecMin"]
finBndDF["delTecPol"] = finBndDF["tecPol"] - finBndDF["tecMin"]
finBndDF["delMlat"] = finBndDF["mlatPol"] - finBndDF["mlatEqu"]
finBndDF["timeStr"] = finBndDF["date"].dt.strftime('%H%M').astype(int)
# # discard dates where delTecEqu and delTecPol are -ve
# finBndDF["dateStr"] = finBndDF["date"].dt.strftime('%Y%m%d')
discrdDatesDelTec = finBndDF[ (finBndDF["delTecEqu"] < delTecCutoff) |\
                   (finBndDF["delTecPol"] < delTecCutoff) ]["date"].values
finBndDF = finBndDF[ ~finBndDF["date"].isin(discrdDatesDelTec) ].reset_index(drop=True)
# Discard those dates where tecMin is greater than 10.
discrdDatestecMin = finBndDF[ (finBndDF["tecMin"] > tecMinCutoff) ]["date"].values
finBndDF = finBndDF[ ~finBndDF["date"].isin(discrdDatestecMin) ].reset_index(drop=True)
# Discard locations where delMlat < 0.
finBndDF = finBndDF[ finBndDF["delMlat"] > 0. ].reset_index(drop=True)

# Read the fitting coefficients to estimate location
# of the trough from them.

coeffCols = [ "a0EquBnd", "a0MinTrgh", "a0PolBnd",\
             "c1EquBnd", "c1MinTrgh", "c1PolBnd",\
             "phiC1EquBnd", "phiC1MinTrgh", "phiC1PolBnd",\
             "phiS1EquBnd", "phiS1MinTrgh", "phiS1PolBnd",\
             "s1EquBnd", "s1MinTrgh", "s1PolBnd", "trghPredTime" ]

coeffFrames = []

# Arrays to store data
mlatEquArr = []
tecEquArr = []
mlonArr = []
mlatPolArr = []
tecPolArr = []
dateArr = []
mlatMinArr = []
tecMinArr = []
mltArr = []
mlonAdjstArr = []
normMLTArr = []
timeStrArr = []


dateStrArr = []

goodDatesList = finBndDF["date"].unique()
for nd, gd in enumerate(goodDatesList):
    ts = pandas.to_datetime(str(gd)) 
    yr = ts.strftime('%Y')
    mt = ts.strftime('%m')
    dt = ts.strftime('%d')
    if int(mt) < 10:
        mt = mt[-1]
    if int(dt) < 10:
        dt = dt[-1]
    cDtStr = ts.strftime('%y%m%d')
    print "current date--->", gd
    # if cDtStr not in dateStrArr:
    #     print cDtStr
        # dateStrArr.append(cDtStr)
    fName = fitDir + "bnd-coeffs" + yr + "-" +\
                mt + "-" + dt + ".txt"
    coeffDF = pandas.read_csv(fName, delim_whitespace=True,\
                                    header=None, names=coeffCols,\
                                infer_datetime_format=True,\
                                parse_dates=["trghPredTime"])
    # need to estimate location of trough between the MLT range where
    # we could calculate trough bnds. Get that range first!!!
    currBndDF = finBndDF[ finBndDF["date"] == gd ]
#     print currBndDF["mlon"].values
#     nMlonStart = numpy.min( currBndDF["mlon"].values )
#     nMlonEnd = numpy.max( currBndDF["mlon"].values )
    for cMlon in currBndDF["mlon"].values:
        cpMlt = round( convert_mlt( cMlon,\
                                    ts , m2a=False ) )
        if cpMlt >= 12 :
            nMlt = cpMlt - 24.
        else:
            nMlt = cpMlt
        selCoeffDF = coeffDF[ coeffDF["trghPredTime"] == ts ]
        selBndDF = currBndDF[ currBndDF["mlon"] == cMlon ]
        minTrghParams = selCoeffDF[ [ "a0MinTrgh", "c1MinTrgh",\
                                 "s1MinTrgh", "phiC1MinTrgh",\
                                 "phiS1MinTrgh" ] ].values[0]
        eqBndParams = selCoeffDF[ [ "a0EquBnd", "c1EquBnd",\
                               "s1EquBnd", "phiC1EquBnd",\
                               "phiS1EquBnd" ] ].values[0]
        polBndParams = selCoeffDF[ [ "a0PolBnd", "c1PolBnd",\
                                "s1PolBnd", "phiC1PolBnd",\
                                "phiS1PolBnd" ] ].values[0]

        minTrghMlat = trough_fit_harmonic(cpMlt, *minTrghParams)
        eqBndMlat = trough_fit_harmonic(cpMlt, *eqBndParams)
        polBndMlat = trough_fit_harmonic(cpMlt, *polBndParams)
        
        # store the data into arrays
        mlatEquArr.append( eqBndMlat ) 
        tecEquArr.append( selBndDF["tecEqu"].values[0] ) 
        mlonArr.append( cMlon ) 
        mlatPolArr.append( polBndMlat ) 
        tecPolArr.append( selBndDF["tecPol"].values[0] ) 
        dateArr.append( gd ) 
        mlatMinArr.append( minTrghMlat ) 
        tecMinArr.append( selBndDF["tecMin"].values[0] ) 
        mltArr.append( cpMlt ) 
        mlonAdjstArr.append( selBndDF["mlonAdjst"].values[0] ) 
        normMLTArr.append( nMlt ) 
        timeStrArr.append( int( ts.strftime('%H%M') ) )

trghFitDF = pandas.DataFrame(
    {'mlatEqu': mlatEquArr,
     'tecEqu': tecEquArr,
     'mlon': mlonArr,
     'mlatPol': mlatPolArr,
     'tecPol': tecPolArr,
     'date': dateArr,
     'mlatMin': mlatMinArr,
     'tecMin': tecMinArr,
     'mlt': mltArr,
     'mlonAdjst': mlonAdjstArr,
     'normMLT': normMLTArr,
     'timeStr': timeStrArr
    })

trghFitDF["delTecEqu"] = trghFitDF["tecEqu"] - trghFitDF["tecMin"]
trghFitDF["delTecPol"] = trghFitDF["tecPol"] - trghFitDF["tecMin"]
trghFitDF["delMlat"] = trghFitDF["mlatPol"] - trghFitDF["mlatEqu"]

# trghFitDF.shape
# asyDF = pandas.read_csv( "../data/Asy_processed.txt", sep=' ' )
# asyDF["date"] = pandas.to_datetime(asyDF["datetimeStr"], format='%Y%m%d-%H-%M')
trghFitDF["dateStr"] = trghFitDF["date"].map(lambda x: x.strftime('%Y%m%d'))
trghFitDF["hour"] = trghFitDF["date"].map(lambda x: x.strftime('%H'))
trghFitDF = pandas.merge( trghFitDF, dstDF, on=["dateStr", "hour"] )
dstBins = [ -150, -50, -25, -10, 10 ]
trghFitDF = pandas.concat( [ trghFitDF, \
                    pandas.cut( trghFitDF["dst_index"], \
                               bins=dstBins ) ], axis=1 )
trghFitDF.columns = [ "date", "mlatEqu", "mlatMin", "mlatPol", "mlon",\
                     "mlonAdjst", "mlt", "normMLT", "tecEqu", "tecMin",\
                     "tecPol", "timeStr", "delTecEqu", "delTecPol", "delMlat",\
                     "dateStr", "hour", "dst_date", "dst_index", "dst_bin" ]
print trghFitDF.head()
feather.write_dataframe(trghFitDF, '../data/trghBndDst-fits.feather')