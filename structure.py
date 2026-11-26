# structure.py -- the structure of the densities of states of a given order, from the file written by dosall.py.
# Computes the exact rank of the differences n_H - n_H' (over the rationals), the order of the first differing moment
# for all the pairs of distinct densities of states, and the rank of the differences of the moments of order 16, ..., 22.
# For N = 24 it checks the decomposition n_H = A + n_H(112) V + P4(H) U, with A, U fixed and
# V the coefficients of (x^4 - x^-4)^20 (x^8 - x^-8)^4, where 112 is the ground state level and P4 the 4-profile.
# For N = 28 it checks that n_H is an affine function of P4(H) and of the numbers of vectors b, with b_1 = 1, such that
# the multiset of the |(Hb)_j| is {2^7, 6^21}, {4^21, 8^7}, or {2^9, 6^18, 10}, read from the outputs of dos.c, and that
# these numbers A, B, C give n_H(140) = 2(A+B) and n_H(136) = 2(C+7A).
# Usage: python3 structure.py dos24.json matrices24, or python3 structure.py dos28.json matrices28 gaethje_output, with
# the folder of the outputs of dos.c, by default matrices24_raw for matrices24, as for dosall.py.
import json, sys, itertools, numpy as np
from fractions import Fraction
L=json.load(open(sys.argv[1])); folder=sys.argv[2]; names=list(L); raw=sys.argv[3] if len(sys.argv)>3 else f'{folder}_raw'
N=len(np.loadtxt(f'{folder}/{names[0]}.txt'))
W=sorted(set(int(w) for n in names for w in L[n])); W=list(range(min(W),max(W)+1,4))
vec=lambda n:[L[n].get(str(w),0) for w in W]
def rank(rows):
    rows=[[Fraction(x) for x in r] for r in rows if any(r)]; rk=0
    if not rows: return 0
    for c in range(len(rows[0])):
        piv=next((i for i in range(rk,len(rows)) if rows[i][c]!=0),None)
        if piv is None: continue
        rows[rk],rows[piv]=rows[piv],rows[rk]
        for i in range(len(rows)):
            if i!=rk and rows[i][c]!=0:
                f=rows[i][c]/rows[rk][c]; rows[i]=[a-f*b for a,b in zip(rows[i],rows[rk])]
        rk+=1
    return rk
v0=vec(names[0]); D={n:[a-b for a,b in zip(vec(n),v0)] for n in names}
print(f"N = {N}: {len(names)} matrices, {len(set(json.dumps(L[n],sort_keys=True) for n in names))} distinct densities of states")
print(f"rank of the differences of the densities of states: {rank(list(D.values()))}")
# first differing moment, over all the pairs of distinct densities of states: the partition of the distinct densities
# is refined by the successive even moments, and we count the pairs separated at each order
dist={}
for n in names: dist.setdefault(json.dumps(L[n],sort_keys=True),n)
reps=list(dist.values()); blocks=[reps]; sep={}; p=0
while any(len(B)>1 for B in blocks):
    p+=1; nb=[]
    for B in blocks:
        g={}
        for n in B: g.setdefault(sum(int(w)**(2*p)*v for w,v in L[n].items()),[]).append(n)
        parts=list(g.values())
        if len(parts)>1:
            tot=len(B)*(len(B)-1)//2-sum(len(P)*(len(P)-1)//2 for P in parts); sep[2*p]=sep.get(2*p,0)+tot
        nb+=parts
    blocks=nb
print(f"first differing moment, over the {len(reps)*(len(reps)-1)//2} pairs of distinct densities of states: {sep}")
mom=lambda n,p: Fraction(sum(int(w)**(2*p)*v for w,v in L[n].items()),4**N)
M=[[mom(n,p)-mom(names[0],p) for p in range(8,12)] for n in names]
print(f"rank of the differences of the moments of order 16, 18, 20, 22: {rank(M)}")
if N==24:
    P4={n:sum(int(np.prod(H[list(I)],axis=0).sum())**4 for H in [np.loadtxt(f'{folder}/{n}.txt',dtype=int)]
              for I in itertools.combinations(range(N),4)) for n in names}
    def mul(p,q):
        r={}
        for a,x in p.items():
            for b,y in q.items(): r[a+b]=r.get(a+b,0)+x*y
        return r
    V={0:1}
    for _ in range(20): V=mul(V,{4:1,-4:-1})
    for _ in range(4): V=mul(V,{8:1,-8:-1})
    Vv=[V.get(w,0) for w in W]; ref=names[0]
    R={n:[d-(L[n].get('112',0)-L[ref].get('112',0))*v for d,v in zip(D[n],Vv)] for n in names}
    U=None; ok=True
    for n in names:
        dp=P4[n]-P4[ref]
        if dp==0: ok&=not any(R[n]); continue
        u=[Fraction(x,dp) for x in R[n]]
        if U is None: U=u
        ok&=(u==U)
    print(f"N = 24: n_H = A + n_H(112) V + P4(H) U, with V from (x^4 - x^-4)^20 (x^8 - x^-8)^4, holds for all: {ok}")
    print(f"   {len(set(P4.values()))} values of P4, {len(set(L[n].get('112',0) for n in names))} values of n_H(112), "
          f"{len(set((P4[n],L[n].get('112',0)) for n in names))} pairs")
    print("   U(w) for w = 0, 4, ..., 108:",[str(U[W.index(w)]) for w in range(0,109,4)])
if N==28:
    I=np.array(list(itertools.combinations(range(N),4)))
    def p4(n):
        H=np.loadtxt(f'{folder}/{n}.txt',dtype=int)
        p=(H[I[:,0]]*H[I[:,1]]*H[I[:,2]]*H[I[:,3]]).sum(axis=1).astype(object); return int(sum(x**4 for x in p))
    def ms(n):
        d={}
        for l in open(f'{raw}/{n}.txt'):
            if l.strip():
                c,*h=map(int,l.split()); d[tuple(h)]=c
        return d
    def key(d):
        h=[0]*(N//2+1)
        for r,m in d.items(): h[r//2]=m
        return tuple(h)
    K=[key({2:7,6:21}),key({4:21,8:7}),key({2:9,6:18,10:1})]
    prm={}
    for n in names:
        m=ms(n); prm[n]=[p4(n)]+[m.get(k,0) for k in K]
    ref=names[0]
    X={n:[Fraction(a-b) for a,b in zip(prm[n],prm[ref])] for n in names}
    sel=[]
    for n in names:
        if rank([X[m] for m in sel]+[X[n]])>len(sel): sel.append(n)
        if len(sel)==4: break
    A=[X[n][:]+[Fraction(x) for x in D[n]] for n in sel]
    for c in range(4):
        p=next(i for i in range(c,4) if A[i][c]!=0); A[c],A[p]=A[p],A[c]
        A[c]=[x/A[c][c] for x in A[c]]
        for i in range(4):
            if i!=c and A[i][c]!=0:
                f=A[i][c]; A[i]=[a-f*b for a,b in zip(A[i],A[c])]
    C=[row[4:] for row in A]
    ok=all([sum(X[n][i]*C[i][j] for i in range(4)) for j in range(len(W))]==D[n] for n in names)
    print(f"N = 28: n_H is an affine function of P4(H) and of the numbers of b with |Hb| = (2^7,6^21), (4^21,8^7), "
          f"(2^9,6^18,10): {ok}")
    ok2=all(L[n].get('140',0)==2*(prm[n][1]+prm[n][2]) and L[n].get('136',0)==2*(prm[n][3]+7*prm[n][1]) for n in names)
    print(f"N = 28: with A, B, C these 3 numbers, n_H(140) = 2(A+B) and n_H(136) = 2(C+7A), for all: {ok2}")
    def F(n): return {int(w)//4:v for w,v in L[n].items()}
    def mulp(p,q):
        r={}
        for a,x in p.items():
            for b,y in q.items(): r[a+b]=r.get(a+b,0)+x*y
        return r
    if 'had28_242' in L and 'had28_313' in L:
        G={0:16}
        for _ in range(28): G=mulp(G,{1:1,-1:-1})
        for _ in range(7): G=mulp(G,{1:1,-1:1})
        d={k:F('had28_242').get(k,0)-F('had28_313').get(k,0) for k in set(F('had28_242'))|set(F('had28_313'))}
        print("N = 28: F_242 - F_313 = 16 (y - 1/y)^28 (y + 1/y)^7, with F_H(y) = sum n_H(w) y^(w/4):",
              all(d.get(k,0)==G.get(k,0) for k in set(d)|set(G)))
    print(f"   {len(set(tuple(v) for v in prm.values()))} distinct parameter vectors; ground state levels: "
          f"{sorted(set(max(int(w) for w in L[n]) for n in names))}")
if N==36:
    def mul(p,q):
        r={}
        for a,x in p.items():
            for b,y in q.items(): r[a+b]=r.get(a+b,0)+x*y
        return r
    C={0:1}
    for _ in range(36): C=mul(C,{6:1,-6:-1})
    Cv=[C.get(w,0) for w in W]
    print(f"N = 36: rank with the crystal vector (x^6 - x^-6)^36 added: {rank(list(D.values())+[Cv])}")
