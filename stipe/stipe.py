# stipe.py -- Metropolis dynamics of the Ising arrays of the Walsh matrices W_1024, W_4096 and of the Paley matrix of
# order 1020 (q = 1019), with the C program mc.c: slow cooling from random configurations, slower cooling in the
# window 0.45 >= T >= 0.1, and heating from a crystalline state built from the bent function (-1)^{<x,y>}.
# The matrices are written in stipe_data/, and 4 runs are done at a time, each run writing stipe_output/<run>.txt (one
# line per temperature: T, mean m, N times the variance of m, largest m seen, then the energy after a final quench).
# The progress is printed every minute and kept in stipe_output/stipe.log.  Runs already finished are skipped, so the
# script can be stopped and restarted.
import numpy as np, subprocess, os, time
os.chdir(os.path.dirname(os.path.abspath(__file__)))
subprocess.run(['cc','-O3','-o','mc','mc.c','-lm'],check=True)
os.makedirs('stipe_data',exist_ok=True); os.makedirs('stipe_output',exist_ok=True)
log=open('stipe_output/stipe.log','a'); t0=time.time()
def say(msg): print(msg,flush=True); log.write(msg+'\n'); log.flush()
def walsh(n):
    pc=np.array([bin(i).count('1') for i in range(n)])
    i=np.arange(n); return (1-2*(pc[np.bitwise_and.outer(i,i)]&1)).astype(np.int8)
def paley(q):
    sq=np.zeros(q,int); sq[(np.arange(1,q)**2)%q]=1
    chi=np.where(sq==1,1,-1); chi[0]=0
    i=np.arange(q); Q=chi[(i[None,:]-i[:,None])%q]
    S=np.zeros((q+1,q+1),int); S[0,1:]=1; S[1:,0]=-1; S[1:,1:]=Q
    return (S+np.eye(q+1,dtype=int)).astype(np.int8)
def crystal(n):
    s=int(round(np.log2(n)))//2; m=1<<s; i=np.arange(n)
    pc=np.array([bin(x).count('1') for x in range(m)])
    b=1-2*(pc[(i>>s)&(i&(m-1))]&1); R=walsh(n).astype(np.int64)@b
    assert np.all(np.abs(R)==int(np.sqrt(n))); return np.sign(R),b
for name,H in (('W1024',None),('W4096',None),('Paley1020',None)):
    if not os.path.exists(f'stipe_data/{name}.txt'):
        say(f"writing the matrix {name} in stipe_data/")
        H=walsh(int(name[1:])) if name[0]=='W' else paley(1019)
        if name[0]=='P': assert np.array_equal(H.astype(np.int64)@H.T.astype(np.int64),1020*np.eye(1020,dtype=np.int64))
        np.savetxt(f'stipe_data/{name}.txt',H,fmt='%d')
    if name[0]=='W' and not os.path.exists(f'stipe_data/{name}_crystal.txt'):
        a,b=crystal(int(name[1:])); np.savetxt(f'stipe_data/{name}_crystal.txt',np.concatenate([a,b]),fmt='%d')
Tcool=[round(x,3) for x in list(np.arange(1.0,0.049,-0.025))+[0.04,0.03,0.02]]
Theat=Tcool[::-1]
Tslow=[round(x,3) for x in np.arange(0.45,0.099,-0.025)]
runs=[('W1024_slow1','W1024',11,100000,Tslow,None),('W1024_slow2','W1024',12,100000,Tslow,None),
      ('W4096_cool1','W4096',21,3000,Tcool,None),('W4096_heat1','W4096',22,3000,Theat,'stipe_data/W4096_crystal.txt')]
runs+=[(f'W1024_cool{s}','W1024',s,10000,Tcool,None) for s in (1,2,3,4)]
runs+=[(f'W1024_heat{s}','W1024',30+s,10000,Theat,'stipe_data/W1024_crystal.txt') for s in (1,2)]
runs+=[(f'Paley1020_cool{s}','Paley1020',40+s,10000,Tcool,None) for s in (1,2)]
finished=lambda r: os.path.exists(f'stipe_output/{r}.txt') and 'quench' in open(f'stipe_output/{r}.txt').read()
todo=[r for r in runs if not finished(r[0])]; running={}
last=0
while todo or running:
    while todo and len(running)<4:
        name,mat,seed,sw,temps,start=todo.pop(0)
        cmd=['./mc',f'stipe_data/{mat}.txt',str(seed),str(sw)]+[str(t) for t in temps]+(['-start',start] if start else [])
        running[name]=(subprocess.Popen(cmd,stdout=open(f'stipe_output/{name}.txt','w')),len(temps))
    time.sleep(5)
    for name,(p,nt) in list(running.items()):
        if p.poll() is not None:
            say(f"{(time.time()-t0)/60:.1f} min: {name} finished"+("" if p.returncode==0 else f" with error {p.returncode}")); del running[name]
    if time.time()-last>60 and running:
        last=time.time()
        st=', '.join(f"{n} {sum(1 for _ in open(f'stipe_output/{n}.txt'))}/{nt} temperatures" for n,(p,nt) in running.items())
        say(f"{(time.time()-t0)/60:.1f} min: {sum(finished(r[0]) for r in runs)}/{len(runs)} runs finished; running: {st}")
say(f"{(time.time()-t0)/60:.1f} min: {sum(finished(r[0]) for r in runs)}/{len(runs)} runs finished"); log.close()
