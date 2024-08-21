import numpy as np

import glob
from decimal import Decimal
import pickle
import re
import matplotlib.pyplot as plt
import pandas as pd
import sys

from mk_dir import dataRoot


#This script loads  and plots the data in the last few files
def format_using_decimal(value):
    # Convert the float to a Decimal
    decimal_value = Decimal(value)
    # Remove trailing zeros and ensure fixed-point notation
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)

if (len(sys.argv)!=3):
    print("wrong number of arguments")
    exit()

T=float(sys.argv[1])
unitCellNum=int(sys.argv[2])

TStr=format_using_decimal(T)

inParamFileName="./V_inv_12_6Params.csv"
rowNum=0
inDf=pd.read_csv(inParamFileName)
oneRow=inDf.iloc[rowNum,:]
a1=float(oneRow.loc["a1"])
b1=float(oneRow.loc["b1"])
a2=float(oneRow.loc["a2"])
b2=float(oneRow.loc["b2"])

def V1(r):
    return a1*r**(-12)-b1*r**(-6)

def sort_data_files_by_swEnd(oneDataFolder):
    """

    :param oneDataFolder:
    :return:
    """


    dataFolderName=oneDataFolder
    dataFilesAll=[]
    sweepEndAll=[]

    for oneDataFile in glob.glob(dataFolderName+"/*.pkl"):
        dataFilesAll.append(oneDataFile)
        matchEnd=re.search(r"sweepEnd(\d+)",oneDataFile)
        if matchEnd:
            sweepEndAll.append(int(matchEnd.group(1)))


    endInds=np.argsort(sweepEndAll)
    # sweepStartSorted=[sweepStartAll[i] for i in startInds]
    sortedDataFiles=[dataFilesAll[i] for i in endInds]

    return sortedDataFiles

dataPath=dataRoot+"/dataAllUnitCell"+str(unitCellNum)+"/row0/T"+TStr+"/U_dist_dataFiles/"


plt_nameU="U"

inUPath=dataPath+"/"+plt_nameU+"/"

sorted_inUFiles=sort_data_files_by_swEnd(inUPath)

lastFilesNum=10

files2Plot=sorted_inUFiles[-lastFilesNum:]

arrU=np.array([])
for pkl_file in files2Plot:
    with open(pkl_file,"rb") as fptr:
        arrIn=pickle.load(fptr)
        arrU=np.append(arrU,arrIn)

avg_arrU=arrU/unitCellNum
plt.figure(figsize=(120,20))
plt.scatter(range(0,len(avg_arrU)),avg_arrU,s=1)
plt.title("T="+str(TStr)+", N="+str(unitCellNum)+", avg U in last "+str(lastFilesNum)+" files")
plt.savefig("T"+TStr+"N"+str(unitCellNum)+"lastFilesU.png")
plt.close()
# peaksU, _ = find_peaks(arrU)
# periods = np.diff(peaksU)
# average_period = np.mean(periods)
#
# # Output the results
# print(f"Periods between peaks: {periods}")
# print(f"Average period: {average_period}")
#
# print(np.sqrt(np.var(periods,ddof=1)))

ind=5

plt_name_x1="xA"+str(ind)
plt_name_x2="xB"+str(ind)

in_x1Path=dataPath+"/"+plt_name_x1+"/"
in_x2Path=dataPath+"/"+plt_name_x2+"/"


sorted_in_x1_Files=sort_data_files_by_swEnd(in_x1Path)
sorted_in_x2_Files=sort_data_files_by_swEnd(in_x2Path)

files_x1_to_plot=sorted_in_x1_Files[-lastFilesNum:]
files_x2_to_plot=sorted_in_x2_Files[-lastFilesNum:]


arr_x1=np.array([])
for pkl_file in files_x1_to_plot:
    with open(pkl_file,"rb") as fptr:
        arrIn=pickle.load(fptr)
        arr_x1=np.append(arr_x1,arrIn)


arr_x2=np.array([])
for pkl_file in files_x2_to_plot:
    with open(pkl_file,"rb") as fptr:
        arrIn=pickle.load(fptr)
        arr_x2=np.append(arr_x2,arrIn)

arr_diff=arr_x2-arr_x1

plt.figure()
plt.scatter(range(0,len(arr_diff)),arr_diff,s=1)
plt.title(plt_name_x2+" - "+plt_name_x1+" in last "+str(lastFilesNum)+" files")
plt.savefig("lastFilesDist.png")
plt.close()

def autocorrelation(x, lag):
    x_mean = np.mean(x)
    x_var = np.var(x)
    N = len(x)

    # Compute autocorrelation for given lag
    acf = np.sum((x[:N-lag] - x_mean) * (x[lag:] - x_mean)) / ((N-lag) * x_var)

    return acf

# autcU=autocorrelation(arrU,100000)
# print(autcU)
