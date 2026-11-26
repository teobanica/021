# spectrum16.py -- the densities of states of the 5 Ising arrays of order 16 form a one-parameter family:
# n_H(64 - 8k) = n(64 - 8k) + (-1)^k binom(16,k) n_H(64),  n_H(w) = n(w) for w = 4 mod 8,
# where n is the density of states of had16.3 and n_H(64) = 2^(14-r) - 128 is the number of crystal ground states,
# r being the 2-rank.  Equivalently Z_H(beta) = Z_163(beta) + n_H(64) (2 sinh beta)^16.  Checked exactly here.
import json
from math import comb
L=json.load(open('glowlaw.json'))
ref={int(k):v for k,v in L['had16_3'].items()}
rank={'had16_0':4,'had16_1':5,'had16_2':6,'had16_3':7,'had16_4':7}
for name,r in rank.items():
    law={int(k):v for k,v in L[name].items()}
    c=law.get(64,0)
    ok=all(law.get(w,0)==ref.get(w,0)+((-1)**((64-w)//8)*comb(16,(64-w)//8)*c if (64-w)%8==0 else 0) for w in range(-64,65,4))
    print(f"{name}: n(64) = {c} = 2^(14-r) - 128 with r = {r}: {c==2**(14-r)-128};  one-parameter formula holds: {ok}")
print("density of states of had16.3, n(w) for w = 0, 4, ..., 64:")
print({w:ref.get(w,0) for w in range(0,65,4)})
