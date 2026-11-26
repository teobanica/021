# twoadic.py -- the density of states modulo 16.  From the exact laws (glowlaw.json) we compute the probabilities of the
# residues of Omega modulo 16, the 2-rank r of the core of the dephased binary matrix M, with H = (-1)^M, and we compare
# P(0) - P(8) with 2^-8 + 2^-(2+r) at N = 16, and with the complex temperature formula: for N divisible by 16 and H
# dephased, P(0) - P(8) = E(i^(alpha^t H beta)), over alpha, beta in {0,1}^N, which is computed independently here.
import numpy as np, json
from fractions import Fraction
L=json.load(open('glowlaw.json'))
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
for name,law in L.items():
    H=np.loadtxt(f'matrices/{name}.txt').astype(int); N=len(H)
    if N<8: continue
    D=H*H[:,[0]]; D=D*D[[0],:]; r=rank2(((1-D)//2)[1:,1:].astype(np.int64))
    P=[Fraction(sum(v for w,v in law.items() if int(w)%16==k),4**N) for k in range(16)]
    line=f"{name:9s} N={N:2d} r={r:2d}  P(4)={P[4]}, P(12)={P[12]}, P(0)-P(8)={P[0]-P[8]}"
    if N==16: line+=f",  2^-8+2^-(2+r)={Fraction(1,256)+Fraction(1,2**(2+r))}"
    if N%16==0:
        B=((np.arange(2**N)[:,None]>>np.arange(N)[None,:])&1)
        S=(B@D.T)%4
        z=np.prod((1+1j**S)/2,axis=1).mean()
        line+=f",  E(i^(a^t H b)) = {z.real:.10f}"
    print(line)
