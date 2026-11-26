# dynlarge.py -- summary of the runs of stipe.py (laptop script, in stipe/), for W_1024, W_4096 and the Paley matrix
# of order 1020: for each run, the energy after the final quench, the temperature where the cooling curve leaves the
# paramagnetic curve by more than 0.01, and for the heating runs, the melting temperature, where m drops below 0.9.
import glob, os
from math import sqrt
mpara=lambda T:(sqrt(1+4/T**2)-1)/(2/T)
for f in sorted(glob.glob('stipe/stipe_output/*.txt')):
    rows=[l.split() for l in open(f) if l.strip()]
    q=float([r for r in rows if r[0]=='quench'][0][1])
    data=[(float(r[0]),float(r[1])) for r in rows if r[0]!='quench']
    name=os.path.basename(f)[:-4]
    if 'heat' in name:
        melt=next(T for T,m in data if m<0.9)
        print(f"{name}: heating from a crystal, m = {dict(data)[0.35]:.4f} at T = 0.35, melted (m < 0.9) at T = {melt}")
    else:
        dep=next((T for T,m in data if abs(m-mpara(T))>0.01),None)
        print(f"{name}: cooling, leaves the paramagnetic curve at T = {dep}, final quench m = {q:.4f}")
