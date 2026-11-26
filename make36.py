# make36.py -- 9 Hadamard matrices of order 36: 8 regular ones, H = J - 2A with A the adjacency matrix of the Latin
# square graph of a random Latin square of order 6 (a strongly regular graph with parameters (36,15,6,6), so that
# H H^t = 36), keeping one matrix for each pair of 4-row profiles of H and H^t found, and the Paley II matrix from
# q = 17.  Writes matrices36/.  The densities of states are then computed by frazier.py, on a laptop.
import numpy as np, random, itertools, os
from collections import Counter
random.seed(7)
def latin(n):
    while True:
        L=[[None]*n for _ in range(n)]
        def bt(pos):
            if pos==n*n: return True
            r,c=divmod(pos,n)
            opts=[s for s in range(n) if all(L[r][j]!=s for j in range(c)) and all(L[i][c]!=s for i in range(r))]
            random.shuffle(opts)
            for s in opts:
                L[r][c]=s
                if bt(pos+1): return True
            L[r][c]=None
            return False
        if bt(0): return L
def hadamard_from_latin(L):
    n=len(L); N=n*n; A=np.zeros((N,N),int)
    cells=[(r,c) for r in range(n) for c in range(n)]
    for i,(r,c) in enumerate(cells):
        for j,(s,d) in enumerate(cells):
            if i!=j and (r==s or c==d or L[r][c]==L[s][d]): A[i,j]=1
    H=1-2*A
    assert np.array_equal(H@H.T,N*np.eye(N,dtype=int)); return H
def paley2(q):
    chi=lambda x: 0 if x%q==0 else (1 if pow(x,(q-1)//2,q)==1 else -1)
    Q=np.array([[chi(j-i) for j in range(q)] for i in range(q)])
    C=np.zeros((q+1,q+1),int); C[0,1:]=1; C[1:,0]=1; C[1:,1:]=Q
    H=np.kron(C,np.array([[1,1],[1,-1]]))+np.kron(np.eye(q+1,dtype=int),np.array([[1,-1],[-1,-1]]))
    assert np.array_equal(H@H.T,(2*q+2)*np.eye(2*q+2,dtype=int)); return H
N=36; I=np.array(list(itertools.combinations(range(N),4)))
def prof(H):
    p=np.abs((H[I[:,0]]*H[I[:,1]]*H[I[:,2]]*H[I[:,3]]).sum(axis=1)); return tuple(sorted(Counter(p.tolist()).items()))
seen={}
for t in range(300):
    H=hadamard_from_latin(latin(6)); k=(prof(H),prof(H.T))
    if k not in seen: seen[k]=H
print("distinct 4-profile pairs among Latin square matrices:",len(seen))
os.makedirs('matrices36',exist_ok=True)
for i,H in enumerate(seen.values()): np.savetxt(f'matrices36/had36_latin{i:02d}.txt',H,fmt='%d')
P=paley2(17); np.savetxt('matrices36/had36_paley.txt',P,fmt='%d')
print("Paley II row sums:",sorted(set(np.abs(P.sum(axis=1)).tolist()))[:5], "profile new:",(prof(P),prof(P.T)) not in seen)
