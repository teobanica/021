# residues.py -- the density of states modulo 16, at all the orders where we know it exactly: N = 4, ..., 20 from
# glowlaw.json, N = 24, 28 from dos24.json, dos28.json (all the classes), N = 32 from dos32.json (the 360 matrices of
# make32.py, computed by golovkin.py and valentina.py), and N = 36 from dos36.json (the 9 regular matrices of make36.py).  For each order we list the 2-ranks r of the cores of the dephased binary matrices, and the
# values of p_4 = P(Omega = 4(16)), p_12 and p_0 - p_8, and we check the formulas p_0 - p_8 = (-1)^((N-4)/8) 2^(-N/2-1)
# for N = 4(8), p_0 - p_8 = -3 2^(-N/2-1) for N = 8(16), and p_0 - p_8 = 2^(-N/2) + 2^(-2-r) for N = 0(16).  For N <= 20 we
# also check the formula p_0 - p_8 = i^n E((-1)^(n(alpha_1+beta_1)) i^<alpha,H beta>), with n = N/4 and H dephased, by
# summing over alpha in closed form: E_alpha = prod_j (1 + s_j i^((H beta)_j))/2, with s_1 = (-1)^n and s_j = 1 otherwise.
# Writes residues.dat, with the points (N, -log_2 |p_0 - p_8|) of the figure, and the sign of p_0 - p_8.
import numpy as np, json
from fractions import Fraction
from collections import defaultdict
def rank2(A):
    A=A.copy()%2; r=0; rows,cols=A.shape
    for c in range(cols):
        piv=[i for i in range(r,rows) if A[i,c]]
        if not piv: continue
        A[[r,piv[0]]]=A[[piv[0],r]]
        for i in range(rows):
            if i!=r and A[i,c]: A[i]^=A[r]
        r+=1
    return r
def core_rank(H):
    D=H*H[:,[0]]; D=D*D[[0],:]; return rank2(((1-D)//2)[1:,1:].astype(np.int64))
def formula(N,r):
    if N%8==4: return Fraction((-1)**((N-4)//8),2**(N//2+1))
    if N%16==8: return Fraction(-3,2**(N//2+1))
    if N%16==0: return Fraction(1,2**(N//2))+Fraction(1,2**(2+r))
sources=[('glowlaw.json','matrices')]+[(f'dos{n}.json',f'matrices{n}') for n in (24,28,32,36)]
table=defaultdict(lambda: defaultdict(set)); count=defaultdict(int)
for fn,folder in sources:
    for name,law in json.load(open(fn)).items():
        H=np.loadtxt(f'{folder}/{name}.txt').astype(int); N=len(H)
        if N<4: continue
        r=core_rank(H); tot=sum(law.values())
        P=[Fraction(sum(v for w,v in law.items() if int(w)%16==k),tot) for k in range(16)]
        assert P[4]==P[12] and sum(P[k] for k in range(16) if k%4)==0
        table[N][r].add((P[4],P[0]-P[8])); count[N]+=1
        assert P[0]-P[8]==formula(N,r), (name,P[0]-P[8])
        if N<=20:
            n=N//4; D=H*H[:,[0]]; D=D*D[[0],:]
            B=((np.arange(2**N)[:,None]>>np.arange(N)[None,:])&1); S=(B@D.T)%4
            s=np.ones(N); s[0]=(-1)**n
            z=(1j**n*((-1.0)**(n*B[:,0]))*np.prod((1+s*1j**S)/2,axis=1)).mean()
            assert abs(z-float(P[0]-P[8]))<1e-12 and abs(z.imag)<1e-12, (name,z)
for N in sorted(table):
    for r in sorted(table[N]):
        for p4,d in sorted(table[N][r]):
            print(f"N = {N:2d}: {count[N]:3d} matrices, r = {r:2d}, p_4 = p_12 = {p4}, p_0 - p_8 = {d}, formula {formula(N,r)}")
print("all the formulas hold, and for N <= 20 the imaginary temperature formula holds too")
from math import log2
with open('residues.dat','w') as f:
    for N in sorted(table):
        for r in sorted(table[N]):
            for p4,d in sorted(table[N][r]): f.write(f"{N} {r} {-log2(abs(d)):.3f} {'+' if d>0 else '-'}\n")
