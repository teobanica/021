# transposes.py -- writes the transposes of the matrices of a folder into folderT/, and runs dos.c on them, writing
# the outputs into folderT_raw/, as used by transpose.py.  Usage: python3 transposes.py matrices24
import numpy as np, glob, os, sys, subprocess
src=sys.argv[1]
os.makedirs(src+'T',exist_ok=True); os.makedirs(src+'T_raw',exist_ok=True)
files=sorted(glob.glob(src+'/*.txt'))
for i,f in enumerate(files):
    n=os.path.basename(f)
    np.savetxt(f'{src}T/{n}',np.loadtxt(f,dtype=int).T,fmt='%d')
    if not (os.path.exists(f'{src}T_raw/{n}') and os.path.getsize(f'{src}T_raw/{n}')>0):
        open(f'{src}T_raw/{n}','w').write(subprocess.run(['./dos',f'{src}T/{n}'],capture_output=True,text=True,check=True).stdout)
    print(f"{i+1}/{len(files)} {n}",flush=True)
