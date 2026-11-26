# crystal.py -- crystallization temperature of the Ising arrays of regular Hadamard matrices.
# The crystal has free energy -1 - log(c_H)/(N beta) per row wire, the liquid -(2 log 2 + Psi(beta))/beta,
# with Psi(beta) = (sqrt(1+4 beta^2) - 1 - log((1+sqrt(1+4 beta^2))/2))/2, and they cross when
#     beta + log(c_H)/N = 2 log 2 + Psi(beta).
# We solve this for the 3 regular classes of order 16, compare with the maxima of the exact specific heat (from
# glowlaw.json), also computed for had16.3, check that c_H for the Walsh matrix W_16 is the number 896 of bent functions in 4 variables,
# and count the self-dual ones, with W_16 b = 4b,
# and use the known numbers of bent functions in 6 and 8 variables (Preneel, Langevin-Leander), and the bounds on
# their asymptotic number, for the Walsh matrices of order 64, 256 and N -> infinity.
import numpy as np, json, itertools
from math import log, sqrt
Psi=lambda b: 0.5*(sqrt(1+4*b*b)-1-log((1+sqrt(1+4*b*b))/2))
def Tc(s):
    lo,hi=1e-6,1e3                      # f(b) = 2log2 + Psi(b) - b - s is decreasing in b
    f=lambda b: 2*log(2)+Psi(b)-b-s
    for _ in range(200):
        mid=(lo+hi)/2
        if f(mid)>0: lo=mid
        else: hi=mid
    return 1/lo
L=json.load(open('glowlaw.json'))
print("order 16: crystallization temperature T_c, against the maximum of the exact specific heat")
for name in ['had16_0','had16_1','had16_2']:
    law=L[name]; w=np.array([int(k) for k in law],float); n=np.array([float(v) for v in law.values()])
    N=16; cH=int(law['64']); x=w/np.sqrt(N); lw=np.log(n)
    best=(0,0)
    for T in np.arange(0.2,0.7,0.0005):
        b=1/T; e=lw+b*x; p=np.exp(e-e.max()); p/=p.sum()
        c=b*b*((p*x*x).sum()-(p*x).sum()**2)/N
        if c>best[0]: best=(c,T)
    print(f"  {name}: c_H = {cH:4d}, log(c_H)/N = {log(cH)/N:.4f}, T_c = {Tc(log(cH)/N):.4f}, "
          f"specific heat maximum {best[0]:.4f} at T = {best[1]:.4f}")
law=L['had16_3']; w=np.array([int(k) for k in law],float); n=np.array([float(v) for v in law.values()])
x=w/4; lw=np.log(n); best=(0,0)
for T in np.arange(0.2,0.7,0.0005):
    b=1/T; e=lw+b*x; p=np.exp(e-e.max()); p/=p.sum(); c=b*b*((p*x*x).sum()-(p*x).sum()**2)/16
    if c>best[0]: best=(c,T)
print(f"  had16_3: no crystalline state, broad maximum of the specific heat {best[0]:.4f} at T = {best[1]:.4f}")
# W_16 and the bent functions in 4 variables
W=np.array([[(-1)**bin(j&k).count('1') for k in range(16)] for j in range(16)])
cnt=0
for bits in itertools.product([1,-1],repeat=16):
    if np.all(np.abs(W@np.array(bits))==4): cnt+=1
sd=sum(1 for bits in itertools.product([1,-1],repeat=16) if np.array_equal(W@np.array(bits),4*np.array(bits)))
print(f"\nW_16: number of b with |(Wb)_j| = 4 for all j: {cnt} (bent functions in 4 variables: 896), "
      f"self-dual ones, with Wb = 4b: {sd}")
# Walsh matrices of order 4^s: c_H = number of bent functions in 2s variables
bent={16:896,64:5425430528,256:99270589265934370305785861242880}
print("\nWalsh matrices: c_H = number of bent functions")
for N,c in bent.items():
    print(f"  N = {N:3d}: log2(c_H) = {log(c,2):8.3f}, log(c_H)/N = {log(c)/N:.4f}, T_c = {Tc(log(c)/N):.4f}")
s_up=11/32*log(2)
print(f"\nN -> infinity: log(c_H)/N between 0 and (11/32) log 2 = {s_up:.4f} (Potapov), so T_c is between "
      f"{Tc(0):.4f} and {Tc(s_up):.4f}")
print(f"for comparison: entropy crisis 1/sqrt(240) = {1/sqrt(240):.4f}")
# order 36: the 9 regular arrays of make36.py, from the exact densities of states computed by dosall.py (dos36.json)
import os
if os.path.exists('dos36.json'):
    L36=json.load(open('dos36.json')); N=36
    print("\norder 36: crystallization temperature T_c, against the maximum of the exact specific heat")
    for name,law in sorted(L36.items()):
        w=np.array([int(k) for k in law],float); n=np.array([float(v) for v in law.values()])
        x=w/np.sqrt(N); lw=np.log(n); best=(0,0)
        for T in np.arange(0.15,0.6,0.0005):
            b=1/T; e=lw+b*x; p=np.exp(e-e.max()); p/=p.sum()
            c=b*b*((p*x*x).sum()-(p*x).sum()**2)/N
            if c>best[0]: best=(c,T)
        cH=int(law['216'])
        print(f"  {name}: c_H = {cH:5d}, T_c = {Tc(log(cH)/N):.4f}, specific heat maximum {best[0]:.3f} at T = {best[1]:.4f}")
