# groundstates.py -- the ground states of the Ising arrays: the pairs of sign vectors (a,b) maximizing
# Omega = <a, H b>.  For each matrix we give the maximal value, the number of ground states, the profiles |Hb| of the
# ground states, and the number of free spins, the rows with (Hb)_i = 0, whose spin a_i does not change the energy.
# For had16.3 we also sort the ground states into orbits of the symmetry group (needs aut.npy, from autgroup.py run
# on matrices/had16_3.txt; skipped if absent).
import numpy as np, itertools, os
from collections import Counter
def ground(H):
    N=len(H)
    B=1-2*((np.arange(2**(N-1))[:,None]>>np.arange(N-1)[None,:])&1); B=np.hstack([np.ones((len(B),1),int),B])
    best=-1; keep=[]
    for s in range(0,len(B),1<<16):
        R=B[s:s+(1<<16)]@H.T; S=np.abs(R).sum(1); m=S.max()
        if m>best: best=m; keep=[]
        if m==best: keep+= [(B[s+i],R[i]) for i in np.nonzero(S==m)[0]]
    states=[]
    for b,r in keep:
        z=np.nonzero(r==0)[0]
        for sg in itertools.product((1,-1),repeat=len(z)):
            a=np.sign(r); a[z]=sg
            states+= [(tuple(a),tuple(b)),(tuple(-a),tuple(-b))]
    prof=Counter(tuple(sorted(np.abs(r))) for b,r in keep)
    return best,states,prof
for name in ['walsh4','walsh8','had12','had16_0','had16_1','had16_2','had16_3','had16_4','had20_0','had20_1','had20_2']:
    H=np.loadtxt(f'matrices/{name}.txt').astype(int); N=len(H)
    best,states,prof=ground(H)
    print(f"{name:8s} N={N:2d}: max Omega = {best} ({best/N**1.5:.4f} N^(3/2)), ground states (a,b): {len(states)}")
    for p,c in sorted(prof.items()):
        z=p.count(0); print(f"     |Hb| = {dict(sorted(Counter(p).items()))}: {2*c} vectors b, {z} free spins each")
    if name=='had16_3' and os.path.exists('aut.npy'):
        Qs=np.load('aut.npy')
        Ps=np.rint(np.einsum('ij,kjl,ml->kim',H,np.transpose(Qs,(0,2,1)),H)/N).astype(int)   # P = H Q^t H^t / N
        rest=set(states); orbits=[]
        while rest:
            a0,b0=map(np.array,next(iter(rest)))
            orb={(tuple(P.T@a0),tuple(Q@b0)) for P,Q in zip(Ps,Qs)}
            assert orb<=rest
            orbits.append((len(orb),int((H@b0==0).sum()))); rest-=orb
        print("     orbits of the symmetry group on the ground states (size, free spins):",orbits)
