# dosall.py -- exact densities of states of the Ising arrays of all the Hadamard matrices in a folder (of the
# same order N), computed with the C program dos.c (compiled here as ./dos), which groups the vectors b by the
# multiset of the |(Hb)_j|, followed by an exact convolution over a.  Usage: python3 dosall.py matrices24 dos24.json,
# the outputs of dos.c being kept in matrices24_raw/, or python3 dosall.py matrices28 dos28.json gaethje_output, the
# outputs of dos.c being then read from the given folders, here those of the laptop script gaethje.py.
# Then checks: total 4^N, the universality of the moments up to order 14, the order 16 formula
# E(Omega^16) - E'(Omega^16) = (2^8 16!/4!) (P4(H) - P4(H')), and groups the matrices with the same density of states.
import numpy as np, json, sys, os, glob, subprocess, math, itertools, time
from collections import Counter
from fractions import Fraction
folder,outfile=sys.argv[1],sys.argv[2]; raws=sys.argv[3:] or [f'{folder}_raw']
if not os.path.exists('./dos'): subprocess.run(['cc','-O3','-o','dos','dos.c'],check=True)
files=sorted(glob.glob(f'{folder}/*.txt'))
laws={}; t0=time.time()
for i,f in enumerate(files):
    N=len(np.loadtxt(f))
    found=[r for r in (f"{d}/{os.path.basename(f)}" for d in raws) if os.path.exists(r) and os.path.getsize(r)>0]
    if found: out=open(found[0]).read()
    else:
        out=subprocess.run(['./dos',f],capture_output=True,text=True,check=True).stdout
        os.makedirs(raws[0],exist_ok=True); open(f"{raws[0]}/{os.path.basename(f)}",'w').write(out)
    tot=Counter()
    for line in out.split('\n'):
        if not line.strip(): continue
        c,*h=map(int,line.split())
        # law over a of sum_j a_j R_j: product over v of (x^r + x^-r)^m, r = 2v + (N mod 2), as exact int64 arrays
        # indexed by the exponent, shifted by S = sum of r m (the coefficients are at most 2^N, fine in int64 for N <= 62)
        S=sum((2*v+(N%2))*m for v,m in enumerate(h)); poly=np.zeros(1,dtype=np.int64); poly[0]=1; off=0
        for v,m in enumerate(h):
            r=2*v+(N%2)
            if m==0: continue
            if r==0: poly=poly*(2**m); continue
            fac=np.zeros(2*r*m+1,dtype=np.int64)
            for k in range(m+1): fac[2*r*k]+=math.comb(m,k)
            poly=np.convolve(poly,fac); off+=r*m
        for e in np.nonzero(poly)[0]: tot[int(e)-off]+=2*c*int(poly[e])
    assert sum(tot.values())==4**N
    laws[os.path.basename(f)[:-4]]={str(w):tot[w] for w in sorted(tot)}
    print(f"{i+1}/{len(files)} {os.path.basename(f)}: {len(out.splitlines())} multisets, top level "
          f"{max(tot)} with {tot[max(tot)]} states, {time.time()-t0:.0f}s",flush=True)
json.dump(laws,open(outfile,'w'))
names=list(laws); N=len(np.loadtxt(files[0]))
def moments(law): return [Fraction(sum(int(w)**(2*p)*v for w,v in law.items()),4**N) for p in range(10)]
def P4(H):
    I=np.array(list(itertools.combinations(range(N),4)))
    p=(H[I[:,0]]*H[I[:,1]]*H[I[:,2]]*H[I[:,3]]).sum(axis=1).astype(object)
    return int(sum(x**4 for x in p))
M={n:moments(laws[n]) for n in names}
ref=names[0]
same14=all(M[n][:8]==M[ref][:8] for n in names)
print(f"\nmoments up to E(Omega^14) identical for all {len(names)} matrices: {same14}")
P={n:P4(np.loadtxt(f'{folder}/{n}.txt',dtype=int)) for n in names}
c=Fraction(2**8*math.factorial(16),24)
ok16=all(M[n][8]-M[ref][8]==c*(P[n]-P[ref]) for n in names)
print(f"order 16 formula E(Omega^16) - E'(Omega^16) = (2^8 16!/4!)(P4 - P4') holds for all: {ok16}")
print(f"distinct values of P4: {len(set(P.values()))}")
groups={}
for n in names: groups.setdefault(json.dumps(laws[n],sort_keys=True),[]).append(n)
print(f"distinct densities of states: {len(groups)}")
for g in groups.values():
    if len(g)>1: print("   same density of states:",' '.join(g))
W=sorted(set(int(w) for n in names for w in laws[n]))
X=np.array([[laws[n].get(str(w),0) for w in W] for n in names],dtype=float)
D=X-X[0]; rk=np.linalg.matrix_rank(D/np.abs(D).max() if np.abs(D).max()>0 else D)
print(f"rank of the differences of the densities of states: {rk}")
top=Counter((max(int(w) for w in laws[n]),laws[n][str(max(int(w) for w in laws[n]))]) for n in names)
print("ground state level Omega and degeneracy, with the number of matrices:")
for (w,d),k in sorted(top.items(),reverse=True): print(f"   Omega = {w} = {w/N**1.5:.4f} N^(3/2), {d} ground states: {k} matrices")
