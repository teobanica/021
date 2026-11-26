# glowlaw.py -- exact density of states of the Ising array of a real Hadamard matrix H, that is, the exact law of
# Omega = sum_ij a_i b_j H_ij over the 4^N pairs of sign vectors (a,b).  For fixed b the law of sum_i a_i R_i, with
# R = Hb, only depends on the multiset of the |R_i|, so we group the vectors b by this multiset, and convolve.
# Writes glowlaw.json: for each matrix, the counts n(w) of pairs (a,b) with Omega = w, as exact integers.
import numpy as np, json, sys, time, glob, os
from collections import Counter
def law(H):
    N=len(H); counts=Counter(); chunk=1<<14
    nb=1<<(N-1)                                    # b_1 = 1; the vectors b and -b give the same law
    for start in range(0,nb,chunk):
        k=np.arange(start,min(nb,start+chunk))
        B=1-2*((k[:,None]>>np.arange(N-1)[None,:])&1)      # entries b_2..b_N
        B=np.hstack([np.ones((len(k),1),int),B])
        R=np.sort(np.abs(B@H.T),axis=1)
        u,c=np.unique(R,axis=0,return_counts=True)
        for row,m in zip(map(tuple,u),c): counts[row]+=int(m)
    tot=Counter()
    for row,m in counts.items():
        p={0:1}
        for r in row:
            q=Counter()
            for w,v in p.items(): q[w+r]+=v; q[w-r]+=v
            p=q
        for w,v in p.items(): tot[w]+=2*m*v
    assert sum(tot.values())==4**N
    return dict(sorted(tot.items())),len(counts)
if __name__=="__main__":
    out={}
    files=sorted(glob.glob('matrices/*.txt'),key=lambda f:(len(np.loadtxt(f)),f))
    for f in files:
        t0=time.time(); H=np.loadtxt(f).astype(int); name=os.path.basename(f)[:-4]
        L,nm=law(H); out[name]={str(w):v for w,v in L.items()}
        top=max(L); print(f"{name:9s} N={len(H):2d}: {len(L)} energy levels, {nm} multisets of |Hb|, "
              f"max Omega = {top} with {L[top]} pairs (a,b)  [{time.time()-t0:.1f}s]",flush=True)
    json.dump(out,open('glowlaw.json','w'),indent=0)
    same=lambda a,b: out[a]==out[b]
    print("laws at N=20 coincide:",same('had20_0','had20_1') and same('had20_1','had20_2'))
