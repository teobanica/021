# figures.py -- the coordinates of the plots of the paper, computed from the exact densities of states (glowlaw.json):
# the normalized energy m = <Omega>/N^(3/2) as a function of beta, for the 4 spectra of order 16 and the large N curve,
# and the specific heat per row wire c = beta^2 Var(Omega/sqrt N)/N as a function of the temperature T = 1/beta.
# Writes figures.dat, with one TikZ coordinate list per curve.
import numpy as np, json
L=json.load(open('glowlaw.json'))
def curves(name,betas):
    law=L[name]; w=np.array([int(k) for k in law],float); n=np.array([float(v) for v in law.values()])
    N=int(round(np.log(n.sum())/np.log(4))); x=w/np.sqrt(N); out=[]
    for b in betas:
        e=np.log(n)+b*x; p=np.exp(e-e.max()); p/=p.sum()
        mean=(p*x).sum(); var=(p*x*x).sum()-mean**2
        out.append((mean/N,b*b*var/N))
    return out
f=open('figures.dat','w')
betas=np.arange(0,8.001,0.25)
for name in ['had16_0','had16_1','had16_2','had16_3','had20_0']:
    m=curves(name,betas)
    f.write(f"% m(beta) {name}\n"+" ".join(f"({b:g},{v[0]:.4f})" for b,v in zip(betas,m))+"\n")
f.write("% m(beta) large N\n"+" ".join(f"({b:g},{((np.sqrt(1+4*b*b)-1)/(2*b) if b>0 else 0):.4f})" for b in betas)+"\n")
Ts=np.arange(0.05,1.5001,0.025)
for name in ['had16_0','had16_1','had16_2','had16_3','had20_0']:
    c=curves(name,1/Ts)
    f.write(f"% c(T) {name}\n"+" ".join(f"({T:.3f},{v[1]:.4f})" for T,v in zip(Ts,c))+"\n")
# large N: c = beta^2 dm/dbeta with m = (sqrt(1+4b^2)-1)/(2b): dm/db = (1 - 1/sqrt(1+4b^2))/(2 b^2)... computed numerically
mm=lambda b:(np.sqrt(1+4*b*b)-1)/(2*b)
f.write("% c(T) large N\n"+" ".join(f"({T:.3f},{(1/T)**2*(mm(1/T+1e-6)-mm(1/T-1e-6))/2e-6:.4f})" for T in Ts)+"\n")
f.close(); print(open('figures.dat').read()[:3000])
for name in ['had16_0','had16_1','had16_2','had16_3','had20_0']:
    c=curves(name,1/Ts); k=int(np.argmax([v[1] for v in c])); print(name,"specific heat maximum",round(c[k][1],4),"at T =",round(Ts[k],3))
