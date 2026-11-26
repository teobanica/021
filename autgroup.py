# autgroup.py -- the symmetry group of a real Hadamard matrix H: all pairs (P,Q) of signed permutation matrices with
# (usage: autgroup.py matrix.txt [output.npy]; writes the matrices Q, as aut.npy by default, about 180 MB for had16.3,
# not kept in the folder, rerun to regenerate)
# P H Q = H.  Backtracking on the images of the columns: after choosing the signed images of columns 0..t-1, the
# multiset of rows of the chosen submatrix, each normalized up to sign, must equal that of the original submatrix.
import numpy as np, sys, time
from collections import Counter
H=np.loadtxt(sys.argv[1]).astype(int); N=len(H)
def rowsig(M):
    # multiset of rows up to sign: normalize each row so that its first entry is +1
    s=M*M[:,:1]; return Counter(map(bytes,(s>0).astype(np.uint8)))
orig=[None]+[rowsig(H[:,:t]) for t in range(1,N+1)]
sols=[]
def bt(img,sg):
    t=len(img)
    if t==N:
        Q=np.zeros((N,N),int)
        for j,(c,s) in enumerate(zip(img,sg)): Q[c,j]=s          # H Q has column j = s * H[:,c]
        HQ=H@Q
        # find P with P (H Q) = H : rows of HQ, up to sign, are a permutation of rows of H
        P=np.zeros((N,N),int)
        for i in range(N):
            k=[r for r in range(N) if np.array_equal(HQ[r],H[i]) or np.array_equal(HQ[r],-H[i])][0]
            P[i,k]=1 if np.array_equal(HQ[k],H[i]) else -1
        assert np.array_equal(P@H@Q,H); sols.append((P,Q)); return
    for c in range(N):
        if c in img: continue
        for s in (1,-1):
            M=H[:,img+[c]]*np.array(sg+[s])
            if rowsig(M)==orig[t+1]: bt(img+[c],sg+[s])
t0=time.time(); bt([],[])
print(f"|Aut(H)| = {len(sols)} pairs (P,Q), including (-1,-1); found in {time.time()-t0:.1f}s")
np.save(sys.argv[2] if len(sys.argv)>2 else 'aut.npy',np.array([Q for P,Q in sols]))
