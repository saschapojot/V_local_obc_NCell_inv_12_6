import sys
import glob
import re
import json
from decimal import Decimal
import pandas as pd
import numpy as np
import subprocess


#this script loads previous data
numArgErr=4
valErr=5
if (len(sys.argv)!=3):
    print("wrong number of arguments.")
    exit(numArgErr)


jsonDataFromConf =json.loads(sys.argv[1])
jsonFromSummary=json.loads(sys.argv[2])

potential_function_name=jsonDataFromConf["potential_function_name"]
U_dist_dataDir=jsonFromSummary["U_dist_dataDir"]
startingFileInd=jsonFromSummary["startingFileInd"]
startingVecPosition=jsonFromSummary["startingVecPosition"]
N=int(jsonDataFromConf["unitCellNum"])

if N<=0:
    print("N="+str(N)+"<=0")
    exit(valErr)
#search and read U_dist files
#give arbitrary values to L, d0Vec, d1Vec without reading data
UInit=6

# y0Init=1
# z0Init=1
# y1Init=1
xVec=np.array(list(range(1,2*N+1)))*0.77
# xA=coordsAll[0::2]
# xB=coordsAll[1::2]

sweepLastFile=-1

#search csv files
csvFileList=[]
sweepEndAll=[]
for file in glob.glob(U_dist_dataDir+"/*.csv"):
    csvFileList.append(file)
    matchEnd=re.search(r"sweepEnd(\d+)",file)
    if matchEnd:
        sweepEndAll.append(int(matchEnd.group(1)))

# print(csvFileList)
def create_loadedJsonData(UVal,xVec,sweepLastFileVal):
    """

    :param UVal:
    :param xVec:
    :param sweepLastFileVal:
    :return:
    """
    initDataDict={
        "U":str(UVal),
        "xVec":list(xVec),
        "sweepLastFile":str(sweepLastFileVal)
    }
    # print(initDataDict)
    return json.dumps(initDataDict)

#if no data found, return the arbitrary values

if len(csvFileList)==0:
    loadedJsonDataStr=create_loadedJsonData(UInit,xVec,sweepLastFile)
    loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
    print(loadedJsonData_stdout)
    exit(0)


#if found csv data
sortedEndInds=np.argsort(sweepEndAll)
sortedsweepEnd=[sweepEndAll[ind] for ind in sortedEndInds]
sortedCsvFileNames=[csvFileList[ind] for ind in sortedEndInds]
sweepLastFile=sortedsweepEnd[-1]

lastFileName=sortedCsvFileNames[-1]

def get_last_row(csv_file):
    result = subprocess.run(['tail', '-n', '1', csv_file], stdout=subprocess.PIPE)
    last_row = result.stdout.decode('utf-8').strip()
    return last_row


csvLastRowStr=get_last_row(lastFileName)
valsInLastRow = [float(value) for value in csvLastRowStr.split(',')]

UInit=valsInLastRow[0]

xVec=valsInLastRow[1:]


loadedJsonDataStr=create_loadedJsonData(UInit,xVec,sweepLastFile)
loadedJsonData_stdout="loadedJsonData="+loadedJsonDataStr
print(loadedJsonData_stdout)
exit(0)