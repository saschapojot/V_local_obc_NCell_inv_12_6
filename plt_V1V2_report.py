import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


inCsv="./V_inv_12_6Params.csv"

inData=pd.read_csv(inCsv)
rowNum=0
a1=inData["a1"][rowNum]
b1=inData["b1"][rowNum]
a2=inData["a2"][rowNum]
b2=inData["b2"][rowNum]


def V1(r):
    return a1/r**12-b1/r**6


def V2(r):
    return a2/r**12-b2/r**6

r1_eq=(2*a1/b1)**(1/6)
V1_at_eq=V1(r1_eq)

r2_eq=(2*a2/b2)**(1/6)


r1_range=np.linspace(r1_eq/1.1,r1_eq*1.2,100)
V1_vals=V1(r1_range)


# print(V1_vals)
plt.figure()
plt.plot(r1_range,V1_vals,color="black")
 # Get default tick locations
V1_x_ticks = [ r1_eq]
plt.xticks(V1_x_ticks)
V1_y_ticks = [V1_at_eq]
plt.yticks(V1_y_ticks)
plt.xlabel("$r$")
plt.ylabel("$V_{1}(r)$")
ax = plt.gca()
ax.spines['left'].set_position(('data', 0))   # Center y-axis at x=0
ax.spines['bottom'].set_position(('data', 0)) # Center x-axis at y=0
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.vlines(x=r1_eq, ymin=0, ymax=V1_at_eq, linestyle="dashed", color="black", linewidth=2)
plt.savefig("V1.pdf")
plt.close()


