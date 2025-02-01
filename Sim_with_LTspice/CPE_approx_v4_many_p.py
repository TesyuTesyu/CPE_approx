import numpy as np
import math as math
import matplotlib.pyplot as plt

freq_min=1e-3#[Hz]
freq_max=10e9#[Hz]
n=30 #RC ladder num
freq_num=1000

#Z=1/(T*(j*Omega)^P)
T=1e-6
p_num=100

C1=T#初期値.何にしても関係ないけど桁落ちしない程度にしたい.

matp=np.linspace(0.01,0.99,p_num)
mat_res=[]
mat_cap=[]
for P in matp:
    Omega_min=2*math.pi*freq_min
    Omega_max=2*math.pi*freq_max
    ab=(Omega_min/Omega_max)**(1/n)
    tau1=1/Omega_min
    R1=tau1/C1
    phi = P*90
    alpha=phi/90
    a=ab**alpha
    b=ab/a

    Rp=R1*(1-a)/a
    Cp=C1*b**n/(1-b)
    exp=np.array(np.linspace(0,n-1,n))
    res=R1*a**exp
    cap=C1*b**exp

    k=int(n/2)
    Omega=np.sqrt(Omega_min*Omega_max)
    inv_Z_RCladder = []
    inv_Z_RCladder.append(1/Rp)
    for i in range(n):
        Rs = res[i]
        Cs = cap[i]
        y=1/(Rs+1/(1j*Omega*Cs))
        inv_Z_RCladder.append(y)
    inv_Z_RCladder.append(1j*Omega*Cp)
    Z_RCladder=np.abs(1/np.sum(inv_Z_RCladder))

    Zcpe = 1/(T*(1j*Omega)**P)
    abs_Zcpe=np.abs(Zcpe)
    D=abs_Zcpe/Z_RCladder#インピーダンスの絶対値を指定のTのそれと合わせる.
    Rp=Rp*D
    res=res*D
    Cp=Cp/D
    cap=cap/D

    mat_res.append(np.concatenate([res,[Rp]]))
    mat_cap.append(np.concatenate([cap,[Cp]]))


fig, ax = plt.subplots(nrows=2, ncols=1, squeeze=False, tight_layout=True, figsize=[8,6], sharex = "col")
cm_name = 'cool'# B->G->R
cm = plt.get_cmap(cm_name, n)#cm(i)[0:3]でRGBを取得できる.

for k in range(n+1):
    Rk=[]
    Ck=[]
    for i in range(p_num):
        Rk.append(mat_res[i][k])
        Ck.append(mat_cap[i][k])

    ax[0,0].plot(matp,Rk, color = cm(k)[0:3], alpha=1)
    ax[1,0].plot(matp,Ck, color = cm(k)[0:3], alpha=1)

    deg=5#多項式の次数.

    Rk2=np.array(Rk)
    coefs=np.polyfit(matp, np.log10(Rk), deg)
    Rk_fitted=np.poly1d(coefs)(matp)
    ax[0,0].plot(matp,10**Rk_fitted, "--", color = cm(k)[0:3], alpha=1)

    print(".param "+"R"+str(k)+"=10**(",end="")
    for i in range(deg-1):
        print(f'{coefs[i]:e}',end="")
        print('*{p}**',end="")
        print(str(deg-i),end="")
        if coefs[i+1]>0:
            print("+",end="")
    print(f'{coefs[i+1]:e}',end="")
    print('*{p}',end="")
    if coefs[i+2]>0:
        print("+",end="")
    print(f'{coefs[i+2]:e}',end="")
    print(")")

    Ck2=np.array(Ck)
    coefs=np.polyfit(matp, np.log10(Ck), deg)
    Ck_fitted=np.poly1d(coefs)(matp)
    ax[1,0].plot(matp,10**Ck_fitted, "--", color = cm(k)[0:3], alpha=1)

    print(".param "+"C"+str(k)+"=10**(",end="")
    for i in range(deg-1):
        print(f'{coefs[i]:e}',end="")
        print('*{p}**',end="")
        print(str(deg-i),end="")
        if coefs[i+1]>0:
            print("+",end="")
    print(f'{coefs[i+1]:e}',end="")
    print('*{p}',end="")
    if coefs[i+2]>0:
        print("+",end="")
    print(f'{coefs[i+2]:e}',end="")
    print(")")

ax[0,0].set_yscale('log')
ax[1,0].set_yscale('log')
plt.show()