"""
CY-piece of V_lambda via Kloosterman sheaf:
  Kl(t) = sum_{x in F_p^*} e((x+1/x)t/p),  Sym^4 trace S4(t) = Kl^4 - 3p Kl^2 + p^2.
  c_p(lambda) := -(1/p) * sum_{t != 0} S4(t) e(-lambda t/p)     (integer, weight 3 expected)
  also s2 := (1/p)* sum_{t!=0} (Kl^2 - p) e(-lambda t/p)        (weight-1 piece)
"""
import numpy as np, sys
from pauli_cy3_counts import primes
def kl(p):
    x=np.arange(1,p); inv=np.array([pow(int(a),p-2,p) for a in x]); v=(x+inv)%p
    t=np.arange(p)[:,None]
    return np.real(np.exp(2j*np.pi*t*v[None,:]/p).sum(1))
def pieces(p,lams):
    K=kl(p); t=np.arange(1,p); K=K[1:]
    S4=K**4-3*p*K**2+p**2; S2=K**2-p
    out=[]
    for l in lams:
        e=np.exp(-2j*np.pi*l*t/p)
        out.append((int(round(-np.real((S4*e).sum())/p)), int(round(np.real((S2*e).sum())/p))))
    return out
if __name__=="__main__":
    P=int(sys.argv[1]); lams=[int(x) for x in sys.argv[2].split(',')]
    for p in primes(P):
        if p<5: continue
        r=pieces(p,lams)
        print(p, '  '.join(f"{l}:{c:>7} [{s:>4}] r={c/p**1.5:+.2f}" for l,(c,s) in zip(lams,r)))
