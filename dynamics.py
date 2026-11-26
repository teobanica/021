# dynamics.py -- Metropolis dynamics of the Ising arrays, with the C program mc.c (compiled here as ./mc).
# Matrices: the Walsh matrices W_64, W_256, which have crystalline states, the Paley matrix of order 252 (q = 251),
# which has none since 252 is not a square, and a random sign matrix of order 256, for comparison (its couplings are
# not orthogonal, so its curve is not comparable with the paramagnetic one).
# Protocols: slow cooling from a random configuration, from T = 1 down to T = 0.02, several seeds, followed by a
# quench at T = 0, and for the Walsh matrices, slow heating starting from a crystalline state, built from the bent
# function (-1)^{x.y} of Maiorana and McFarland.  Writes dynamics.dat (the curves) and prints a summary.
import numpy as np, subprocess, os, sys, time
from math import sqrt, log
if not os.path.exists('./mc'): subprocess.run(['cc','-O3','-o','mc','mc.c','-lm'],check=True)
os.makedirs('dyn',exist_ok=True)
def walsh(n): return np.array([[(-1)**bin(i&j).count('1') for j in range(n)] for i in range(n)])
def paley(q):
    chi=lambda x: 0 if x%q==0 else (1 if pow(x,(q-1)//2,q)==1 else -1)
    Q=np.array([[chi(j-i) for j in range(q)] for i in range(q)])
    S=np.zeros((q+1,q+1),int); S[0,1:]=1; S[1:,0]=-1; S[1:,1:]=Q
    return S+np.eye(q+1,dtype=int)
rng=np.random.default_rng(1)
mats={'W64':walsh(64),'W256':walsh(256),'Paley252':paley(251),'random256':rng.choice([-1,1],size=(256,256))}
for k,H in mats.items():
    N=len(H)
    if k!='random256': assert np.array_equal(H@H.T,N*np.eye(N,dtype=int)),k
    np.savetxt(f'dyn/{k}.txt',H,fmt='%d')
def crystal(N):
    s=int(round(np.log2(N)))//2; m=1<<s
    b=np.array([(-1)**bin((i>>s)&(i&(m-1))).count('1') for i in range(N)])
    R=walsh(N)@b; assert np.all(np.abs(R)==int(sqrt(N))); a=np.sign(R)
    return a,b
for k in ('W64','W256'):
    a,b=crystal(len(mats[k])); np.savetxt(f'dyn/{k}_crystal.txt',np.concatenate([a,b]),fmt='%d')
Tcool=[round(x,3) for x in list(np.arange(1.0,0.049,-0.025))+[0.04,0.03,0.02]]
Theat=Tcool[::-1]
sweeps=int(sys.argv[1]) if len(sys.argv)>1 else 10000
seeds=range(1,9)
def run(k,seed,temps,start=None):
    cmd=['./mc',f'dyn/{k}.txt',str(seed),str(sweeps)]+[str(t) for t in temps]+(['-start',start] if start else [])
    out=subprocess.run(cmd,capture_output=True,text=True,check=True).stdout.split('\n')
    rows=[list(map(float,l.split())) for l in out if l and not l.startswith('quench')]
    q=float([l for l in out if l.startswith('quench')][0].split()[1])
    return np.array(rows),q
mpara=lambda T:(sqrt(1+4/T**2)-1)/(2/T)
dat=open('dynamics.dat','w'); t0=time.time()
for k in mats:
    N=len(mats[k]); C=[];Q=[]
    for s in seeds:
        r,q=run(k,s,Tcool); C.append(r); Q.append(q)
    C=np.array(C)
    print(f"{k}: cooling, {len(seeds)} runs of {sweeps} sweeps per temperature, {time.time()-t0:.0f}s",flush=True)
    print("   T      m(cooling)   sd     m(paramagnet)  specific heat")
    for i,T in enumerate(Tcool):
        m=C[:,i,1].mean(); sd=C[:,i,1].std(); c=(C[:,i,2]/T**2).mean()
        print(f"   {T:.2f}   {m:.4f}   {sd:.4f}   {mpara(T):.4f}   {c:.3f}")
        dat.write(f"{k} cool {T:.3f} {m:.5f} {sd:.5f} {mpara(T):.5f} {c:.4f}\n")
    print(f"   quench at T=0: m = {np.mean(Q):.4f} (min {min(Q):.4f}, max {max(Q):.4f}); best m seen {C[:,:,3].max():.4f}")
    dat.write(f"{k} quench {np.mean(Q):.5f} {min(Q):.5f} {max(Q):.5f} {C[:,:,3].max():.5f}\n")
    if k in ('W64','W256'):
        Hh=[run(k,s,Theat,f'dyn/{k}_crystal.txt')[0] for s in list(seeds)[:4]]; Hh=np.array(Hh)
        print(f"{k}: heating from a crystalline state, {time.time()-t0:.0f}s")
        for i,T in enumerate(Theat):
            m=Hh[:,i,1].mean(); print(f"   {T:.2f}   {m:.4f}"); dat.write(f"{k} heat {T:.3f} {m:.5f}\n")
    dat.flush()
