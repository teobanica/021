# hypercubedos.py -- exact densities of states of the fully frustrated hypercube, for D <= 6, viewed as the Ising
# array of the weighing matrix A of order N = 2^(D-1) and weight D of Proposition 2.4 (hypercube.py), with energy
# E = -Omega/sqrt D, Omega = <Ab,a>.  Uses the C program dos.c (compiled here as ./dos), which groups the vectors b by
# the multiset of the |(Ab)_j|, all of the parity of D, followed by an exact convolution over a.  Prints the ground
# states, with m = Omega/(N sqrt D) <= 1, and the maximum of the specific heat per row wire, and for D = 4, where the
# bound m = 1 is reached, the crystallization temperature T_c of Section 7.
import numpy as np, subprocess, os, math, json
from collections import Counter
from math import sqrt, log
if not os.path.exists('./dos'): subprocess.run(['cc','-O3','-o','dos','dos.c'],check=True)
os.makedirs('hypercube',exist_ok=True)
def weighing(D):
    n=2**D; ev=[x for x in range(n) if bin(x).count('1')%2==0]; od=[x for x in range(n) if bin(x).count('1')%2==1]
    J={}
    for x in range(n):
        for nu in range(D):
            if not (x>>nu)&1:
                s=(-1)**bin(x&((1<<nu)-1)).count('1'); J[(x,x|(1<<nu))]=J[(x|(1<<nu),x)]=s
    A=np.array([[J.get((x,y),0) for y in od] for x in ev]); assert np.array_equal(A@A.T,D*np.eye(len(ev),dtype=int))
    return A
Psi=lambda b: 0.5*(sqrt(1+4*b*b)-1-log((1+sqrt(1+4*b*b))/2))
def Tc(s):
    lo,hi=1e-6,1e3
    for _ in range(200):
        mid=(lo+hi)/2
        if 2*log(2)+Psi(mid)-mid-s>0: lo=mid
        else: hi=mid
    return 1/lo
out={}
for D in range(2,7):
    A=weighing(D); N=len(A); np.savetxt(f'hypercube/A{D}.txt',A,fmt='%d')
    raw=subprocess.run(['./dos',f'hypercube/A{D}.txt'],capture_output=True,text=True,check=True).stdout
    tot=Counter()
    for line in raw.split('\n'):
        if not line.strip(): continue
        c,*h=map(int,line.split()); poly={0:1}
        for v,m in enumerate(h):
            r=2*v+(D%2)
            for _ in range(m):
                q=Counter()
                for w,x in poly.items(): q[w+r]+=x; q[w-r]+=x
                poly=q
        for w,x in poly.items(): tot[w]+=2*c*x
    assert sum(tot.values())==4**N
    out[D]={str(w):tot[w] for w in sorted(tot)}
    top=max(tot); w=np.array(sorted(tot),float); n=np.array([tot[int(x)] for x in w],float)
    x=w/sqrt(D); lw=np.log(n); best=(0,0)
    for T in np.arange(0.05,2.0,0.0005):
        b=1/T; e=lw+b*x; p=np.exp(e-e.max()); p/=p.sum(); cv=b*b*((p*x*x).sum()-(p*x).sum()**2)/N
        if cv>best[0]: best=(cv,T)
    line=(f"D = {D}: N = {N}, {len(tot)} levels, ground states at m = {top/(N*sqrt(D)):.4f}, {tot[top]} of them, "
          f"specific heat maximum {best[0]:.3f} at T = {best[1]:.4f}")
    if abs(top-N*sqrt(D))<1e-9: line+=f", crystallization temperature T_c = {Tc(log(tot[top])/N):.4f}"
    print(line,flush=True)
json.dump(out,open('hypercube.json','w'))
