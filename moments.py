# moments.py -- exact moments E(Omega^2p) of the real glow, from glowlaw.json, and the 4-row profile of each matrix.
# Checks that the moments of the matrices of the same order agree up to p = 7, and that the first difference, at
# p = 8, is given by the profile: E(Omega^16) - E'(Omega^16) = (2^8 16!/4!) * (P4(H) - P4(H')), where
# P4(H) = sum over the 4-subsets I of rows of (sum_j prod_{i in I} H_ij)^4.  Writes moments.json.
import numpy as np, json, itertools, math
from fractions import Fraction
L=json.load(open('glowlaw.json'))
def profile(H):
    N=len(H); P4=0; vals=[]
    for I in itertools.combinations(range(N),4):
        s=int(np.prod(H[list(I)],axis=0).sum()); P4+=s**4; vals.append(abs(s))
    return P4,dict(sorted(Counter(vals).items()))
from collections import Counter
out={}
mats={}
for name,law in L.items():
    H=np.loadtxt(f'matrices/{name}.txt').astype(int); N=len(H)
    m=[Fraction(sum(int(w)**(2*p)*v for w,v in law.items()),4**N) for p in range(11)]
    P4,prof=profile(H) if N>=4 else (0,{})
    out[name]={"N":N,"moments":[str(x) for x in m],"P4":P4,"profile":{str(k):v for k,v in prof.items()}}
    mats.setdefault(N,[]).append(name)
    print(f"{name:9s} E(Omega^2)={m[1]}  E(Omega^4)={m[2]}  4-profile {prof}",flush=True)
print()
for N,names in mats.items():
    if len(names)<2: continue
    ref=names[0]; mr=[Fraction(x) for x in out[ref]["moments"]]
    for nm in names[1:]:
        mm=[Fraction(x) for x in out[nm]["moments"]]
        first=next((p for p in range(11) if mm[p]!=mr[p]),None)
        pred=Fraction(2**8*math.factorial(16),24)*(out[nm]["P4"]-out[ref]["P4"])
        fd=f"E(Omega^{2*first})" if first is not None else "none up to E(Omega^20)"
        print(f"N={N}: {nm} vs {ref}: first different moment {fd}; "
              f"difference at 16: {mm[8]-mr[8]}, predicted by the profile: {pred}, equal: {mm[8]-mr[8]==pred}")
for N in sorted(mats):
    n=mats[N][0]; m=[Fraction(x) for x in out[n]["moments"]]
    print(f"N={N}: E(Omega^4) = {m[2]}, 3N^4-6N^3+4N^2 = {3*N**4-6*N**3+4*N**2}")
json.dump(out,open('moments.json','w'),indent=1)
