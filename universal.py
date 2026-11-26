# universal.py -- the universal moments of the real glow.  For p <= 7, E(Omega^2p) is a polynomial U_p(N) of degree 2p,
# the same for all real Hadamard matrices of order N, N = 1, 2 included.  We interpolate U_p from the orders
# N = 1, 2, 4, 8, 12, 16, 20 (moments.json, together with the order 1), check it on the remaining orders when there
# are some, and give the cumulants kappa_2p, whose leading terms are (-1)^(p-1) (2p-1)! Cat_(p-1) N^(p+1).
import json
from fractions import Fraction as F
from math import comb, factorial
M=json.load(open('moments.json'))
data={1:[F(1)]*11}
for name in ['walsh2','walsh4','walsh8','had12','had16_0','had20_0']:
    data[M[name]['N']]=[F(x) for x in M[name]['moments']]
Ns=sorted(data)
def interp(pts):
    """Lagrange interpolation with exact fractions; returns the coefficients, constant term first"""
    n=len(pts); coef=[F(0)]*n
    for i,(xi,yi) in enumerate(pts):
        basis=[F(1)]
        den=F(1)
        for j,(xj,_) in enumerate(pts):
            if j==i: continue
            basis=[F(0)]+basis
            for k in range(len(basis)-1): basis[k]-=xj*basis[k+1]
            den*=(xi-xj)
        for k in range(n): coef[k]+=yi*basis[k]/den
    return coef
def show(c):
    terms=[]
    for k in range(len(c)-1,-1,-1):
        if c[k]==0: continue
        s=f"{c[k]}"+(f"N^{k}" if k>1 else ("N" if k==1 else ""))
        terms.append(s)
    return " + ".join(terms).replace("+ -","- ")
def cumul(m):
    mm={}; ka={}
    for n in range(1,len(m)*2):
        mm[n]=m[n//2] if n%2==0 else F(0)
    for n in range(1,len(m)*2):
        ka[n]=mm[n]-sum(comb(n-1,j-1)*ka[j]*mm[n-j] for j in range(1,n))
    return ka
for p in range(1,6):
    deg=2*p; pts=[(N,data[N][p]) for N in Ns][:deg+1]
    if len(pts)<deg+1: break
    c=interp(pts)+[F(0)]*0
    c=c[:deg+1]
    ok=all(sum(ck*N**k for k,ck in enumerate(c))==data[N][p] for N in Ns)
    print(f"E(Omega^{2*p}) = {show(c)}   (checked on all {len(Ns)} orders: {ok}, using {len(pts)} of them)")
cum={N:cumul(data[N][:8]) for N in Ns}
for p in range(1,5):
    deg=p+1; pts=[(N,cum[N][2*p]) for N in Ns][:deg+1]
    c=interp(pts)[:deg+1]
    ok=all(sum(ck*N**k for k,ck in enumerate(c))==cum[N][2*p] for N in Ns)
    lead=(-1)**(p-1)*factorial(2*p-1)*comb(2*p-2,p-1)//p
    print(f"kappa_{2*p} = {show(c)}   (checked on all orders: {ok}; predicted leading coefficient {lead})")
