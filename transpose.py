# transpose.py -- the arrays of H and of its transpose H^t have the same density of states, since Omega = <Hb,a> =
# <H^t a,b>.  For a given order, we compare the distributions of the multisets of the |(Hb)_j| (rows, in folder_raw)
# and of the |(H^t a)_k| (columns, folderT_raw), computed by dos.c, and check that 2 matrices have the same density
# of states precisely when they have the same unordered pair of such distributions.  For the pairs of matrices with
# the same density of states, we also compare a finer invariant, the multiset over the rows of the distributions of
# the numbers |sum_k prod_{i in I} H_ik| over the sets I of 4 rows containing the given row, for H and for H^t:
# when the row invariants of H, H' differ, but those of H' and H^t agree, the pair is consistent with H' = H^t.
# Usage: python3 transpose.py matrices24 dos24.json, or python3 transpose.py matrices28 dos28.json gaethje_output, with
# the folder of the outputs of dos.c for the rows, by default matrices24_raw for matrices24, as for dosall.py.
import json, sys, itertools, numpy as np
from collections import Counter
folder,dosfile=sys.argv[1],sys.argv[2]; raw=sys.argv[3] if len(sys.argv)>3 else f'{folder}_raw'
L=json.load(open(dosfile)); names=sorted(L)
N=len(np.loadtxt(f'{folder}/{names[0]}.txt'))
def ms(path): return frozenset(tuple(map(int,l.split())) for l in open(path) if l.strip())
Mr={n:ms(f'{raw}/{n}.txt') for n in names}
Mc={n:ms(f'{folder}T_raw/{n}.txt') for n in names}
key=lambda n: json.dumps(L[n],sort_keys=True)
groups={}
for n in names: groups.setdefault(key(n),[]).append(n)
pairdist=lambda n: frozenset([Mr[n],Mc[n]])
ok=all(len({pairdist(n) for n in g})==1 for g in groups.values()) and len(groups)==len(set(pairdist(n) for n in names))
print(f"N = {N}: {len(names)} matrices, {len(groups)} densities of states, {len(set(Mr.values()))} distinct row "
      f"distributions, {sum(Mr[n]==Mc[n] for n in names)} matrices with equal row and column distributions")
print(f"same density of states <=> same pair of row and column distributions: {ok}")
I=np.array(list(itertools.combinations(range(N),4)))
def local(H):
    p=np.abs((H[I[:,0]]*H[I[:,1]]*H[I[:,2]]*H[I[:,3]]).sum(axis=1))
    per=[Counter() for _ in range(N)]
    for idx,v in zip(I,p):
        for j in idx: per[j][int(v)]+=1
    return tuple(sorted(tuple(sorted(c.items())) for c in per))
R={};C={}
for n in names:
    H=np.loadtxt(f'{folder}/{n}.txt',dtype=int); R[n]=local(H); C[n]=local(H.T)
cnt=Counter()
for g in groups.values():
    for a,b in itertools.combinations(g,2):
        if R[a]!=R[b] and R[a]==C[b] and C[a]==R[b]: cnt['consistent with transposition only']+=1
        elif R[a]==R[b] and C[a]==C[b] and not (R[a]==C[b]): cnt['not transposition']+=1
        elif R[a]==R[b]: cnt['undecided']+=1
        else: cnt['other']+=1
print("pairs of matrices with the same density of states, by the finer invariant:",dict(cnt))
