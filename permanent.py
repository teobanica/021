# permanent.py -- the permanents of the Hadamard matrices in a folder, computed by the C program perm.c (compiled here
# as ./perm, by the formula of Glynn), compared with their densities of states, from the file written by dosall.py.
# Prints the number of distinct values of |per(H)|, the number of distinct pairs (density of states, |per(H)|), and the
# values of |per(H)|/g, with g their gcd, and writes the permanents in <folder>_perm.json.  4 matrices at a time.
# Usage: python3 permanent.py matrices24 dos24.json
import json, sys, os, glob, subprocess, time
from math import gcd
from concurrent.futures import ThreadPoolExecutor
folder,dosfile=sys.argv[1],sys.argv[2]
if not os.path.exists('./perm'): subprocess.run(['cc','-O3','-o','perm','perm.c'],check=True)
L=json.load(open(dosfile)); out=f'{folder}_perm.json'
P=json.load(open(out)) if os.path.exists(out) else {}
todo=[n for n in sorted(L) if n not in P]; t0=time.time()
def work(n): return n,int(subprocess.run(['./perm',f'{folder}/{n}.txt'],capture_output=True,text=True,check=True).stdout)
with ThreadPoolExecutor(4) as ex:
    for k,(n,p) in enumerate(ex.map(work,todo)):
        P[n]=p; json.dump(P,open(out,'w'))
        if (k+1)%20==0 or k+1==len(todo): print(f"{k+1}/{len(todo)} permanents, {time.time()-t0:.0f}s",flush=True)
dens=lambda n: json.dumps(L[n],sort_keys=True)
A={n:abs(P[n]) for n in L}; g=0
for v in A.values(): g=gcd(g,v)
nd=len(set(dens(n) for n in L)); na=len(set(A.values())); npair=len(set((dens(n),A[n]) for n in L))
print(f"{len(L)} matrices: {nd} densities of states, {na} values of |per|, {npair} pairs (density, |per|)")
print(f"same density of states <=> same |per|: {nd==na==npair}")
print(f"gcd g = {g}, values of |per|/g: {sorted(set(v//g for v in A.values()))}")
