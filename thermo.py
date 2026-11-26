# thermo.py -- exact thermodynamics of the Ising arrays of order 16 and 20, from the densities of states in
# glowlaw.json.  The energy is E = -Omega/sqrt(N), and per row wire we compute the normalized mean energy
# m = <Omega>/N^(3/2), the specific heat c = beta^2 Var(Omega/sqrt N)/N and the entropy s = (log Z - beta <Omega>/sqrt N)/N,
# and compare with the large N predictions m(beta) = (sqrt(1+4 beta^2)-1)/(2 beta) and
# s(beta) = 2 log 2 - log((1+sqrt(1+4 beta^2))/2)/2.  Writes thermo.json and the tables of the figures, thermo.dat.
import numpy as np, json
L=json.load(open('glowlaw.json'))
def thermo(name,beta):
    law=L[name]; w=np.array([int(k) for k in law],float); n=np.array([float(v) for v in law.values()])
    N=round(np.sqrt(np.sqrt(n.sum()))**1) ; N=int(round(np.log(n.sum())/np.log(4)))
    x=w/np.sqrt(N); lw=np.log(n)
    out=[]
    for b in beta:
        e=lw+b*x; M=e.max(); p=np.exp(e-M); Z=p.sum(); p/=Z
        mean=(p*x).sum(); var=(p*x*x).sum()-mean**2
        logZ=M+np.log(Z)
        out.append((mean/N,b*b*var/N,(logZ-b*mean)/N))
    return np.array(out),N
mlim=lambda b:(np.sqrt(1+4*b*b)-1)/(2*b) if b>0 else 0.0
slim=lambda b:2*np.log(2)-0.5*np.log((1+np.sqrt(1+4*b*b))/2)
beta=np.round(np.arange(0,8.01,0.25),2)
res={}
dat=open('thermo.dat','w')
for name in ['had16_0','had16_1','had16_2','had16_3','had20_0']:
    T,N=thermo(name,beta); res[name]=T.tolist()
    dat.write(f"# {name}: beta, m, c, s\n")
    for b,(m,c,s) in zip(beta,T): dat.write(f"{b:5.2f} {m:.4f} {c:.4f} {s:.4f}\n")
dat.write("# large N: beta, m, s\n")
for b in beta: dat.write(f"{b:5.2f} {mlim(b):.4f} {slim(b):.4f}\n")
json.dump(res,open('thermo.json','w'))
print(" beta   m(16.0)  m(16.3)  m(20)   m(limit)   c(16.0)  c(16.3)   s(16.0)  s(16.3)  s(limit)")
for i,b in enumerate(beta):
    if b in (0.25,0.5,1,1.5,2,3,4,6,8):
        print(f"{b:5.2f}  {res['had16_0'][i][0]:.4f}   {res['had16_3'][i][0]:.4f}   {res['had20_0'][i][0]:.4f}   {mlim(b):.4f}"
              f"     {res['had16_0'][i][1]:.4f}   {res['had16_3'][i][1]:.4f}    {res['had16_0'][i][2]:.4f}   {res['had16_3'][i][2]:.4f}   {slim(b):.4f}")
bc=np.sqrt(240); print(f"\nlarge N: the entropy s(beta) vanishes at beta = sqrt(240) = {bc:.4f}, T = {1/bc:.5f}, where m = {mlim(bc):.5f}")
for name in ['had16_0','had16_1','had16_2','had16_3','had20_0']:
    law=L[name]; top=max(int(k) for k in law); N=int(round(np.log(sum(law.values()))/np.log(4)))
    print(f"{name}: ground state Omega = {top} = {top/N**1.5:.4f} N^(3/2), degeneracy {law[str(top)]}, entropy per row wire log(deg)/N = {np.log(law[str(top)])/N:.4f}")
