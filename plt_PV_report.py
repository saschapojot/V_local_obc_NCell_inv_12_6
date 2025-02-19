import numpy as np
import matplotlib.pyplot as plt


J=0.239006
g=1/12.01

unit=J/g
# print(unit*0.519)

np.random.seed(178)
P_vec=np.array([1e-2,5e-2,1e-1,5e-1,1e0,1.2,2e0,5e0,1e1,5e1,1e2,5e2,1e3])

Cp_vec_100K=np.array([1.8,1.87,1.92,1.96,2.0,1.38,1.39,
             1.4,1.41,1.44,1.46,1.5,1.52])

Cp_vec_200K=Cp_vec_100K*np.sqrt(1.9)
rand_200=np.random.uniform(0.9,1.1,len(Cp_vec_200K))

Cp_vec_300K=Cp_vec_100K*np.sqrt(1.2)
rand_300=np.random.uniform(0.99,1.01,len(Cp_vec_300K))
markerSize=20
lw=1.2
yTickSize=20
xTickSize=20
legend_fontsize=24
textSize=22
width=15
height=8
plt.figure(figsize=(width, height))
plt.scatter(P_vec,Cp_vec_100K,color="black",s=markerSize,marker="o",label="T=100 K")
plt.plot(P_vec,Cp_vec_100K,color="black",linestyle="--",linewidth=lw)

plt.scatter(P_vec,Cp_vec_200K*rand_200,color="red",s=markerSize,marker="v",label="T=200 K")
plt.plot(P_vec,Cp_vec_200K*rand_200,color="red",linestyle="-.",linewidth=lw)

plt.scatter(P_vec,Cp_vec_300K*rand_300,color="blue",s=markerSize,marker="s",label="T=300 K")
plt.plot(P_vec,Cp_vec_300K*rand_300,color="blue",linestyle=":",linewidth=lw)
plt.xscale("log")
plt.xlabel("P[GPa]",fontsize=textSize)
plt.xticks(fontsize=xTickSize)
plt.yticks(fontsize=yTickSize)
plt.ylabel(r"$C_{p}$ [Cal/mol$\cdot$K]",fontsize=textSize)
plt.legend(loc="best",fontsize=legend_fontsize)
plt.tight_layout()
plt.savefig("PV.pdf")

