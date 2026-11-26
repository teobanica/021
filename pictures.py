# pictures.py -- the data of the pictures and tables of Sections 2 to 7.
# (1) The fully frustrated cube, D = 3, with the signs J(x,x+e_v) = (-1)^(x_1+...+x_(v-1)): its edges, the frustration
#     of its 6 faces, and its 4 x 4 weighing matrix A, between the even and the odd vertices, with A A^t = 3.
# (2) The full density of states of had16.3, from glowlaw.json.
# (3) A crystalline state of W_16, built from the bent function x1 x2 + x3 x4, and a ground state of had16.3, of the type
#     with |Hb| given by 6 on 6 rows and 2 on 10 rows, as sign matrices D_a H D_b, with rows and columns sorted by sums.
# (4) The entropy per row wire s(T) of the large N limit, and of the exact arrays of had16.0, had16.3 and of the Paley
#     matrix of order 36, from glowlaw.json and dos36.json.
# (5) The crystallization temperatures of the Walsh matrices W_4, W_16, W_64, W_256, from the numbers of bent functions.
# (6) The specific heat per row wire of the arrays of order 32 of make32.py, from dos32.json: W_32, the Kronecker
#     products, the Paley matrix, and the envelope of the curves of the 355 doublings, with the ground states.
# (7) The numbers of ground states of the 355 doublings of order 32, all at Omega = 168, against the 2-rank r of the
#     core, read from matrices32.log: for each r, the number of matrices, and the minimum, mean and maximum.
# (8) The difference of the densities of states of the matrices had28_242 and had28_313 of order 28, whose moments
#     coincide up to order 26, from dos28.json, in units of 10^6 configurations.
# (9) For the 27 densities of states of order 24, from dos24.json: the ground state level and degeneracy, and the
#     temperatures and heights of the maxima of the specific heat, on a grid of step 0.001.
# Writes pictures.dat, and prints a summary.
import numpy as np, json, itertools
from math import log, sqrt
out=open('pictures.dat','w')
def w(s): out.write(s+"\n")
# (1) the cube
D=3; V=list(itertools.product([0,1],repeat=D))
J={}
for x in V:
    for v in range(D):
        if x[v]==0:
            y=list(x); y[v]=1; J[(x,tuple(y))]=(-1)**sum(x[:v])
w("% cube: edges x -> y with sign"); [w(f"{''.join(map(str,x))} {''.join(map(str,y))} {s:+d}") for (x,y),s in J.items()]
def sgn(x,y): return J.get((x,y),J.get((y,x),0))
faces=0
for v1,v2 in itertools.combinations(range(D),2):
    for x in V:
        if x[v1]==0 and x[v2]==0:
            a=list(x); a[v1]=1; a=tuple(a); b=list(x); b[v2]=1; b=tuple(b); c=list(a); c[v2]=1; c=tuple(c)
            faces+=sgn(x,a)*sgn(a,c)*sgn(c,b)*sgn(b,x)==-1
even=[x for x in V if sum(x)%2==0]; odd=[x for x in V if sum(x)%2==1]
A=np.array([[sgn(x,y) for y in odd] for x in even])
w("% cube: even vertices "+" ".join(''.join(map(str,x)) for x in even)+", odd vertices "+" ".join(''.join(map(str,y)) for y in odd))
w("% cube: weighing matrix A\n"+"\n".join(" ".join(f"{e:+d}" for e in r) for r in A))
print(f"cube: {faces} frustrated faces out of 6, A A^t = 3: {np.array_equal(A@A.T,3*np.eye(4,dtype=int))}")
# (2) the density of states of had16.3
L=json.load(open('glowlaw.json')); law=L['had16_3']
w("% density of states of had16.3, Omega n(Omega)")
for o in range(-64,65,4): w(f"{o} {law.get(str(o),0)}")
print(f"had16.3: {sum(law.values())} configurations, at Omega = 0: {law['0']}")
# (3) the crystal and the frustrated ground state
def walsh(n): return np.array([[(-1)**bin(j&k).count('1') for k in range(n)] for j in range(n)])
def picture(H,b,name):
    a=np.sign(H@b); a[a==0]=1; X=(a[:,None]*H)*b[None,:]
    X=X[np.argsort(-X.sum(1),kind='stable')]; X=X[:,np.argsort(-X.sum(0),kind='stable')]
    w(f"% {name}: rows of D_a H D_b, 1 for +1 and 0 for -1, then row sums and column sums")
    for r in X: w("{"+",".join('1' if e>0 else '0' for e in r)+"}")
    w("row sums "+" ".join(map(str,X.sum(1)))); w("column sums "+" ".join(map(str,X.sum(0))))
    return X
W=walsh(16); b=np.array([(-1)**(((k>>3)&1)*((k>>2)&1)+((k>>1)&1)*(k&1)) for k in range(16)])
X=picture(W,b,'crystal of W16'); print(f"W16 crystal: row sums {set(X.sum(1).tolist())}, column sums {set(X.sum(0).tolist())}")
H=np.loadtxt('matrices/had16_3.txt',dtype=int)
for bits in itertools.product([1,-1],repeat=15):
    b=np.array((1,)+bits); s=sorted(np.abs(H@b))
    if s==[2]*10+[6]*6: break
X=picture(H,b,'ground state of had16.3'); print(f"had16.3 ground state: sum {X.sum()}, row sums {X.sum(1).tolist()}, column sums {X.sum(0).tolist()}")
# the 3 types of ground states of had16.3, by the multisets of the row sums and of the column sums of D_a H D_b
B=np.array(list(itertools.product([1,-1],repeat=16))); HB=B@H.T; S=np.abs(HB).sum(1); types={}
for i in np.nonzero(S==S.max())[0]:
    hb=HB[i]; z=np.nonzero(hb==0)[0]
    for sg in itertools.product([1,-1],repeat=len(z)):
        a=np.sign(hb); a[z]=sg; X=a[:,None]*H*B[i][None,:]
        key=(tuple(sorted(X.sum(1).tolist(),reverse=True)),tuple(sorted(X.sum(0).tolist(),reverse=True)))
        types[key]=types.get(key,0)+1
for (r,c),v in types.items(): print(f"had16.3: {v} ground states with row sums {r} and column sums {c}")
# (4) entropies
def entropy(law,T):
    x=np.array([int(k) for k in law],float); n=np.array([float(v) for v in law.values()])
    N=int(round(log(n.sum())/log(4))); x/=sqrt(N); lw=np.log(n); res=[]
    for t in T:
        e=lw+x/t; M=e.max(); p=np.exp(e-M); Z=p.sum(); p/=Z; res.append((M+log(Z)-(p*x).sum()/t)/N)
    return res
T=np.round(np.concatenate([np.arange(0.03,0.6,0.02),np.arange(0.6,1.21,0.1)]),2)
slim=[2*log(2)-0.5*log((1+sqrt(1+4/t**2))/2) for t in T]
D36=json.load(open('dos36.json'))
for name,lw in [('large N',None),('had16_0',L['had16_0']),('had16_3',L['had16_3']),('had36_paley',D36['had36_paley'])]:
    s=slim if lw is None else entropy(lw,T)
    w(f"% entropy s(T) {name}\n"+" ".join(f"({t:g},{v:.4f})" for t,v in zip(T,s)))
print(f"entropy at T = 0.03: had16.0 {entropy(L['had16_0'],[0.03])[0]:.4f} (log 896/16 = {log(896)/16:.4f}), "
      f"had16.3 {entropy(L['had16_3'],[0.03])[0]:.4f} (log 43008/16 = {log(43008)/16:.4f}), "
      f"Paley 36 {entropy(D36['had36_paley'],[0.03])[0]:.4f} (log 816/36 = {log(816)/36:.4f})")
# (5) the Walsh matrices
Psi=lambda b: 0.5*(sqrt(1+4*b*b)-1-log((1+sqrt(1+4*b*b))/2))
def Tc(s):
    lo,hi=1e-6,1e3
    for _ in range(200):
        mid=(lo+hi)/2
        if 2*log(2)+Psi(mid)-mid-s>0: lo=mid
        else: hi=mid
    return 1/lo
bent={4:8,16:896,64:5425430528,256:99270589265934370305785861242880}
for N,c in bent.items():
    line=f"W{N}: c_H = {c}, log2(c_H) = {log(c,2):.3f}, log(c_H)/N = {log(c)/N:.4f}, T_c = {Tc(log(c)/N):.4f}"
    w("% "+line); print(line)
# (6) the arrays of order 32
L32=json.load(open('dos32.json'))
def heat(law,Ts):
    x=np.array([int(k) for k in law],float); n=np.array([float(v) for v in law.values()])
    N=int(round(log(n.sum())/log(4))); x/=sqrt(N); lw=np.log(n); res=[]
    for t in Ts:
        e=lw+x/t; p=np.exp(e-e.max()); p/=p.sum(); res.append(((p*x*x).sum()-(p*x).sum()**2)/(t*t*N))
    return np.array(res)
T=np.round(np.arange(0.06,1.001,0.02),2)
for name in ['had32_walsh','had32_kron1','had32_kron2','had32_kron3','had32_paley']:
    c=heat(L32[name],T); top=max(int(k) for k in L32[name]); k=int(np.argmax(c))
    w(f"% specific heat {name}\n"+" ".join(f"({t:g},{v:.4f})" for t,v in zip(T,c)))
    print(f"{name}: ground states at Omega = {top}, {L32[name][str(top)]} of them, specific heat maximum {c[k]:.3f} at T = {T[k]:g}")
D=[n for n in L32 if n.startswith('had32_d')]; C=np.array([heat(L32[n],T) for n in D])
w("% specific heat of the doublings, upper envelope\n"+" ".join(f"({t:g},{v:.4f})" for t,v in zip(T,C.max(0))))
w("% specific heat of the doublings, lower envelope\n"+" ".join(f"({t:g},{v:.4f})" for t,v in zip(T,C.min(0))))
tops=[max(int(k) for k in L32[n]) for n in D]; gs=[L32[n][str(t)] for n,t in zip(D,tops)]
km=C.argmax(1)
print(f"{len(D)} doublings: ground states at Omega = {set(tops)}, from {min(gs)} to {max(gs)} of them, specific heat maxima "
      f"from {C.max(1).min():.3f} to {C.max(1).max():.3f}, at T from {T[km].min():g} to {T[km].max():g}")
# (7) ground states of the doublings against the 2-rank
rk={l.split(':')[0]:int(l.split()[-1]) for l in open('matrices32.log')}
by={}
for n,t,g in zip(D,tops,gs): by.setdefault(rk[n],[]).append(g)
w("% ground states of the doublings of order 32: r, number of matrices, minimum, mean, maximum")
for r in sorted(by):
    v=by[r]; w(f"{r} {len(v)} {min(v)} {round(sum(v)/len(v))} {max(v)}")
    print(f"doublings with r = {r}: {len(v)} matrices, from {min(v)} to {max(v)} ground states, mean {round(sum(v)/len(v))}")
# (8) the pair of order 28 first differing at order 28
L28=json.load(open('dos28.json')); a,b=L28['had28_242'],L28['had28_313']
d={o:a.get(str(o),0)-b.get(str(o),0) for o in range(-140,141,4)}
w("% difference of the densities of states of had28_242 and had28_313, Omega and difference/10^6")
w(",".join(f"{o}/{d[o]/1e6:.3f}" for o in range(-140,141,4) if d[o]))
mom=lambda p: sum(o**p*v for o,v in d.items())
print(f"had28_242 - had28_313: {sum(1 for v in d.values() if v)} nonzero levels, first nonzero moment of order "
      f"{min(p for p in range(2,60,2) if mom(p))}, largest difference {max(abs(v) for v in d.values())}")
# (9) the specific heat of the arrays of order 24
L24=json.load(open('dos24.json')); seen={}
for n in sorted(L24): seen.setdefault(json.dumps(L24[n],sort_keys=True),[]).append(n)
T=np.arange(0.08,0.6,0.001)
w("% order 24: matrices with this density of states, ground level, ground states, maxima (T, c) of the specific heat")
for k,ns in seen.items():
    law=L24[ns[0]]; top=max(int(x) for x in law); c=heat(law,T)
    loc=[i for i in range(1,len(T)-1) if c[i]>c[i-1] and c[i]>=c[i+1]]
    line=f"{len(ns)} {top} {law[str(top)]} "+" ".join(f"{T[i]:.3f}/{c[i]:.3f}" for i in loc)
    w(line); print("order 24:",line)
out.close()
