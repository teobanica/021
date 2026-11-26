# mckay.py -- converts McKay's files of Hadamard matrices (one matrix per line, each row a hexadecimal number,
# whose binary digits give the row, 1 -> +1 and 0 -> -1) into text files of +-1 entries, and checks H H^t = N.
# Usage: python3 mckay.py had24.txt 24 matrices24   (writes matrices24/had24_00.txt, had24_01.txt, ...)
import numpy as np, sys, os
src,N,out=sys.argv[1],int(sys.argv[2]),sys.argv[3]
os.makedirs(out,exist_ok=True)
lines=[l.split() for l in open(src) if l.strip()]
for i,l in enumerate(lines):
    assert len(l)==N
    H=np.array([[1 if (int(h,16)>>(N-1-k))&1 else -1 for k in range(N)] for h in l],dtype=int)
    assert np.array_equal(H@H.T,N*np.eye(N,dtype=int)), (i,'not Hadamard')
    np.savetxt(f"{out}/had{N}_{i:02d}.txt" if len(lines)<100 else f"{out}/had{N}_{i:03d}.txt",H,fmt='%d')
print(f"{len(lines)} Hadamard matrices of order {N} converted and checked, in {out}/")
