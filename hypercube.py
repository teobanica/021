# hypercube.py -- the fully frustrated hypercube of Marinari, Parisi and Ritort, as an Ising array.
# Signs on the edges of {0,1}^D, in the gauge J(x, x+e_nu) = (-1)^(x_1+...+x_{nu-1}); we check that all the square
# faces are frustrated, and that the matrix W of couplings between the even and the odd vertices, of size 2^(D-1),
# satisfies W W^t = D, i.e. is a weighing matrix of weight D.  For D <= 5 we also compute the ground state energy
# per spin, E/2^D = -max_b sum_j |(Wb)_j| / (sqrt(D) 2^D), to be compared with the bound -1/2, reached for D = 1, 4.
import numpy as np, itertools
from math import sqrt
def coupling(D):
    J={}
    for x in range(2**D):
        for nu in range(D):
            if not (x>>nu)&1:
                s=(-1)**bin(x&((1<<nu)-1)).count('1')
                J[(x,x|(1<<nu))]=J[(x|(1<<nu),x)]=s
    return J
for D in range(1,9):
    J=coupling(D); n=2**D
    faces=frus=0
    for x in range(n):
        for mu in range(D):
            for nu in range(mu+1,D):
                if not (x>>mu)&1 and not (x>>nu)&1:
                    y,z,w=x|(1<<mu),x|(1<<nu),x|(1<<mu)|(1<<nu)
                    faces+=1; frus+=(J[(x,y)]*J[(y,w)]*J[(w,z)]*J[(z,x)]==-1)
    ev=[x for x in range(n) if bin(x).count('1')%2==0]; od=[x for x in range(n) if bin(x).count('1')%2==1]
    W=np.array([[J.get((x,y),0) for y in od] for x in ev])
    ok=np.array_equal(W@W.T,D*np.eye(len(ev),dtype=int))
    line=f"D = {D}: {frus}/{faces} faces frustrated, W of size {len(ev)} with W W^t = D: {ok}"
    if D<=5:
        best=max(np.abs(W@np.array((1,)+b)).sum() for b in itertools.product([1,-1],repeat=len(od)-1))
        line+=f", ground state energy per spin {-best/(sqrt(D)*n):.4f}"
    print(line, flush=True)
