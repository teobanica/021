# gaethje.py -- exact densities of states of the Ising arrays of the 487 Hadamard matrices of order 28. Compiles dos.c,
# then runs it on the matrices in matrices28/, 4 at a time, writing one file per matrix in gaethje_output/, as used by
# dosall.py, with the progress printed and kept in gaethje_output/gaethje.log.  Can be stopped and restarted, the
# matrices already done being skipped.
import subprocess, os, glob, time
from concurrent.futures import ThreadPoolExecutor
os.chdir(os.path.dirname(os.path.abspath(__file__)))
subprocess.run(['cc','-O3','-o','dos','dos.c'],check=True)
out='gaethje_output'; os.makedirs(out,exist_ok=True)
files=sorted(glob.glob('matrices28/*.txt'))
todo=[f for f in files if not (os.path.exists(f'{out}/{os.path.basename(f)}') and os.path.getsize(f'{out}/{os.path.basename(f)}')>0)]
log=open(f'{out}/gaethje.log','a'); t0=time.time(); done=len(files)-len(todo)
def say(msg): print(msg,flush=True); log.write(msg+'\n'); log.flush()
def work(f):
    res=subprocess.run(['./dos',f],capture_output=True,text=True,check=True).stdout
    tmp=f'{out}/{os.path.basename(f)}.tmp'; open(tmp,'w').write(res); os.rename(tmp,f'{out}/{os.path.basename(f)}')
    return f,len(res.splitlines())
say(f"{done}/{len(files)} done")
with ThreadPoolExecutor(4) as ex:
    for f,k in ex.map(work,todo):
        done+=1; say(f"{done}/{len(files)} {os.path.basename(f)}: {k} multisets, {time.time()-t0:.0f}s elapsed")
log.close()
