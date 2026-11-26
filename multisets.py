# multisets.py -- a finer invariant than the density of states: for a Hadamard matrix H of order N, the distribution
# of the multiset of the numbers |(Hb)_j|, over the 2^(N-1) vectors b with b_1 = 1.  The density of states only
# depends on this distribution, so if 2 matrices have the same distribution, they have the same density of states.
# We compare the 3 classes of order 20, and also print the local field profile counts.
import numpy as np, sys
from collections import Counter
def profile(H):
    N=len(H); counts=Counter(); chunk=1<<14; nb=1<<(N-1)
    for start in range(0,nb,chunk):
        k=np.arange(start,min(nb,start+chunk))
        B=1-2*((k[:,None]>>np.arange(N-1)[None,:])&1)
        B=np.hstack([np.ones((len(k),1),int),B])
        R=np.sort(np.abs(B@H.T),axis=1)
        u,c=np.unique(R,axis=0,return_counts=True)
        for row,m in zip(map(tuple,u),c): counts[row]+=int(m)
    return counts
names=sys.argv[1:] or ['had20_0','had20_1','had20_2']
P={n:profile(np.loadtxt(f'matrices/{n}.txt',dtype=int)) for n in names}
for n in names: print(f"{n}: {len(P[n])} distinct multisets", flush=True)
ref=P[names[0]]
for n in names[1:]:
    keys=set(ref)|set(P[n]); diff={k:(ref.get(k,0),P[n].get(k,0)) for k in keys if ref.get(k,0)!=P[n].get(k,0)}
    print(f"{names[0]} vs {n}: {'identical' if not diff else str(len(diff))+' multisets differ'}")
    for k,(x,y) in sorted(diff.items(),key=lambda t:-abs(t[1][0]-t[1][1]))[:12]:
        prof=Counter(k); print("   ", dict(sorted(prof.items())), x, y)
