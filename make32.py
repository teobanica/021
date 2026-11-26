# make32.py -- Hadamard matrices of order 32, for the computations of golovkin.py and valentina.py (laptop scripts):
# the Walsh matrix W_32, the Kronecker products had_16.k x W_2, the Paley matrix from q = 31, and doublings
# [[A, B], [A, -B]] with A one of the 5 matrices of order 16 and B a random signed permutation of one of them, which
# are Hadamard for any such A, B.  The matrices are kept only when their pair of 4-row profiles (of H and H^t) is new,
# so they are certainly inequivalent.  Writes matrices32/ and matrices32.log, with the 2-rank of the core of each.
import numpy as np, itertools, os, glob
from collections import Counter
rng=np.random.default_rng(32)
had=[np.loadtxt(f'matrices/had16_{k}.txt',dtype=int) for k in range(5)]
W2=np.array([[1,1],[1,-1]])
def walsh(n): return np.array([[(-1)**bin(i&j).count('1') for j in range(n)] for i in range(n)])
def paley(q):
    chi=lambda x: 0 if x%q==0 else (1 if pow(x,(q-1)//2,q)==1 else -1)
    Q=np.array([[chi(j-i) for j in range(q)] for i in range(q)])
    S=np.zeros((q+1,q+1),int); S[0,1:]=1; S[1:,0]=-1; S[1:,1:]=Q
    return S+np.eye(q+1,dtype=int)
def signperm(H):
    n=len(H); P=rng.permutation(n); Q=rng.permutation(n)
    return (rng.choice([-1,1],n)[:,None]*H[P][:,Q]*rng.choice([-1,1],n)[None,:])
I4=np.array(list(itertools.combinations(range(32),4)))
def prof(H):
    p=np.abs((H[I4[:,0]]*H[I4[:,1]]*H[I4[:,2]]*H[I4[:,3]]).sum(axis=1)); return tuple(sorted(Counter(p.tolist()).items()))
def rank2(H):
    H=H*H[0][None,:]; H=H*H[:,0][:,None]; M=((1-H[1:,1:])//2)%2; M=M.copy(); r=0
    for c in range(M.shape[1]):
        piv=next((i for i in range(r,M.shape[0]) if M[i,c]),None)
        if piv is None: continue
        M[[r,piv]]=M[[piv,r]]
        for i in range(M.shape[0]):
            if i!=r and M[i,c]: M[i]^=M[r]
        r+=1
    return r
os.makedirs('matrices32',exist_ok=True); seen=set(); kept=[]
def keep(name,H):
    assert np.array_equal(H@H.T,32*np.eye(32,dtype=int))
    k=(prof(H),prof(H.T))
    if k in seen: return
    seen.add(k); kept.append(name); np.savetxt(f'matrices32/{name}.txt',H,fmt='%d')
keep('had32_walsh',walsh(32))
for k in range(1,5): keep(f'had32_kron{k}',np.kron(had[k],W2))
keep('had32_paley',paley(31))
pairs=[(i,j) for i in range(5) for j in range(i,5)]
t=0
while len(kept)<360 and t<20000:
    i,j=pairs[t%len(pairs)]; t+=1
    A=had[i]; B=signperm(had[j]); keep(f'had32_d{i}{j}_{t:05d}',np.block([[A,B],[A,-B]]))
with open('matrices32.log','w') as f:
    for n in kept:
        line=f"{n}: 2-rank of the core {rank2(np.loadtxt(f'matrices32/{n}.txt',dtype=int))}"; print(line); f.write(line+'\n')
print(f"{len(kept)} inequivalent matrices of order 32, from {t} doublings tried")
