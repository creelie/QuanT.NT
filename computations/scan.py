import sys
from eta import FORMS4, coeffs
from sym4 import pieces
from pauli_cy3_counts import primes
P=int(sys.argv[1]); lams=list(range(int(sys.argv[2]),int(sys.argv[3])))
F4={k:coeffs(v,P+1) for k,v in FORMS4.items()}; F4={k:v for k,v in F4.items() if v is not None}
ps=[p for p in primes(P) if p>=5]
def leg(a,p): a%=p; return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
C={p:dict(zip(lams,[x[0] for x in pieces(p,lams)])) for p in ps}
tws=[1,-1,2,-2,3,-3,5,-5,6,-6,7,-7,10,-10,15,-15]
for lam in lams:
    D=lam*(lam*lam-16)*(lam*lam-64)
    good=[p for p in ps if (D==0 and p>3) or (D!=0 and D%p and p>3)]
    good=[p for p in good if 6 % p]  # keep p>=5
    hits=[f"a_p({n})+{s}*({t}|p)*p" for n,f in F4.items() for t in tws for s in (1,-1)
          if all(C[p][lam]-s*leg(t,p)*p==f[p] for p in good)]
    hits+= [f"a_p({n})" for n,f in F4.items() if all(C[p][lam]==f[p] for p in good)]
    mx=max(abs(C[p][lam])/p**1.5 for p in good)
    print(f"lam={lam:3d} #p={len(good)} max|c|/p^1.5={mx:.2f}", hits)
