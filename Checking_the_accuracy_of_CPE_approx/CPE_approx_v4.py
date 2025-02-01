import numpy as np
import math as math
import matplotlib.pyplot as plt
#以下の論文をもとに実装.
#Juraj VALSA et.al., "Network Model of the CPE," RADIOENGINEERING, VOL. 20, NO. 3, 2011.

freq_min=0.01#[Hz]
freq_max=100e6#[Hz]
n=20 #RC ladder num
freq_num=1000

#Z=1/(T*(j*Omega)^P)
T=1e-6
P=0.5

C1=T#初期値.何にしても関係ないけど桁落ちしない程度にしたい.

Omega_min=2*math.pi*freq_min
Omega_max=2*math.pi*freq_max
ab=(Omega_min/Omega_max)**(1/n)
tau1=1/Omega_min
R1=tau1/C1
phi = P*90
alpha=phi/90
a=ab**alpha
b=ab/a

freq_min=freq_min/100
freq_max=freq_max*100

Rp=R1*(1-a)/a
Cp=C1*b**n/(1-b)
exp=np.array(np.linspace(0,n-1,n))
res=R1*a**exp
cap=C1*b**exp

#インピーダンスの絶対値を補正するために、ログスケールでの平均の周波数におけるインピーダンスの絶対値を計算する.
Omega=np.sqrt(Omega_min*Omega_max)#ログスケールでの平均の周波数.
inv_Z_RCladder = []
inv_Z_RCladder.append(1/Rp)
for i in range(n):
    Rs = res[i]
    Cs = cap[i]
    y=1/(Rs+1/(1j*Omega*Cs))#各RC回路のアドミタンス.
    inv_Z_RCladder.append(y)
inv_Z_RCladder.append(1j*Omega*Cp)
Z_RCladder=np.abs(1/np.sum(inv_Z_RCladder))#インピーダンスの絶対値.

Zcpe = 1/(T*(1j*Omega)**P)
abs_Zcpe=np.abs(Zcpe)#CPEのインピーダンスの絶対値 @平均の周波数.
D=abs_Zcpe/Z_RCladder#インピーダンスの絶対値を合わせるための補正係数.
Rp=Rp*D
res=res*D
Cp=Cp/D
cap=cap/D

#CPEとRCラダーの周波数特性の計算.
freq = np.logspace(np.log10(freq_min), np.log10(freq_max), freq_num)
Omega = 2*math.pi*freq

Zcpe = 1/(T*(1j*Omega)**P)
abs_Zcpe=np.abs(Zcpe)

Z_RCladder = []
for k in range(freq_num):
    inv_Z_RCladder = []
    inv_Z_RCladder.append(1/Rp)
    for i in range(n):
        Rs = res[i]
        Cs = cap[i]
        y=1/(Rs+1/(1j*Omega[k]*Cs))
        inv_Z_RCladder.append(y)
    inv_Z_RCladder.append(1j*Omega[k]*Cp)
    Z_RCladder.append(1/np.sum(inv_Z_RCladder))

Z_RCladder = np.array(Z_RCladder)
abs_Z_RCladder = np.abs(Z_RCladder)
phase_Z_RCladder = np.arctan2(np.imag(Z_RCladder),np.real(Z_RCladder))

print(".param ",end="")
for i in range(n):
    print("R"+str(i)+"=",end="")
    print(f'{res[i]:e}',end=", ")
for i in range(n):
    print("C"+str(i)+"=",end="")
    print(f'{cap[i]:e}',end=", ")

print("Rp=",end="")
print(f'{Rp:e}',end=", ")

print("Cp=",end="")
print(f'{Cp:e}',end="")

fig, ax = plt.subplots(nrows=2, ncols=1, squeeze=False, tight_layout=True, figsize=[8,6], sharex = "col")
#ax[0,0].plot(freq,abs_Z_RCladder*(freq**alpha),"k-")
ax[0,0].plot(freq,abs_Z_RCladder,"k-")
ax[0,0].plot(freq,abs_Zcpe,"r--")
ax[0,0].set_xscale('log')
ax[0,0].set_yscale('log')

ax[1,0].plot(freq,phase_Z_RCladder*180/math.pi,"k-")
ax[1,0].plot([freq[0],freq[-1]],[-phi,-phi],"r--")
ax[1,0].set_xscale('log')

plt.show()
