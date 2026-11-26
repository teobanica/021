# valentina.py -- exact densities of states of the Ising arrays of Hadamard matrices of order 32, the odd-numbered
# half of those made by make32.py (W_32, Kronecker products, Paley, and doublings [[A,B],[A,-B]] of matrices of order
# 16), the other half being done by golovkin.py. The C program dos.c enumerates the 2^31 vectors b with b_1 = 1 and
# counts them according to the multiset of the |(Hb)_j|, writing one file per matrix in valentina_output/, as used by
# dosall.py, the progress of each run being kept in valentina_progress/.  4 matrices at a time, about 1-2 minutes each
# on a laptop; the progress is printed every minute and kept in valentina_output/valentina.log, and the matrices
# already done are skipped, so the script can be stopped and restarted.
import subprocess, os, glob, time
os.chdir(os.path.dirname(os.path.abspath(__file__)))
subprocess.run(['cc','-O3','-o','dos','dos.c'],check=True)
out,prog='valentina_output','valentina_progress'
os.makedirs(out,exist_ok=True); os.makedirs(prog,exist_ok=True)
files=sorted(glob.glob('matrices32/*.txt'))[1::2]
done=lambda f: os.path.exists(f'{out}/{os.path.basename(f)}') and os.path.getsize(f'{out}/{os.path.basename(f)}')>0
todo=[f for f in files if not done(f)]; running={}; t0=time.time(); last=0; log=open(f'{out}/valentina.log','a')
def say(msg): print(msg,flush=True); log.write(msg+'\n'); log.flush()
def pct(n):
    try:
        l=open(f'{prog}/{n}').read().strip().split('\n')[-1].split()[1].split('/'); return 100*int(l[0])/int(l[1])
    except Exception: return 0.0
while todo or running:
    while todo and len(running)<4:
        f=todo.pop(0); n=os.path.basename(f)
        running[n]=subprocess.Popen(['./dos',f],stdout=open(f'{out}/{n}.tmp','w'),stderr=open(f'{prog}/{n}','w'))
    time.sleep(5)
    for n,p in list(running.items()):
        if p.poll() is not None:
            if p.returncode==0: os.rename(f'{out}/{n}.tmp',f'{out}/{n}')
            say(f"{(time.time()-t0)/60:.1f} min: {n} "+("done" if p.returncode==0 else f"error {p.returncode}")); del running[n]
    if time.time()-last>60 and running:
        last=time.time()
        say(f"{(time.time()-t0)/60:.1f} min: {sum(done(f) for f in files)}/{len(files)} done; running: "+', '.join(f"{n} {pct(n):.0f}%" for n in running))
say(f"{(time.time()-t0)/60:.1f} min: {sum(done(f) for f in files)}/{len(files)} done"); log.close()
