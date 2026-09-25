import itertools, sys
from eta import eta_q
from sym4 import pieces
from pauli_cy3_counts import primes
lam=int(sys.argv[1]); P=150
ps=[p for p in primes(P) if p>=5]
def leg(a,p): a%=p; return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
C={p:pieces(p,[lam])[0][0] for p in ps}
targets={}
for t in [1,-1,2,-2,3,-3,6,-6]:
    for s in (1,-1,0):
        targets[(t,s)]={p:C[p]-s*leg(t,p)*p for p in ps}
found=set()
for N in [12,16,18,24,32,36,48,64]:
    divs=[d for d in range(1,N+1) if N%d==0]
    for combo in itertools.product(range(-4,9), repeat=len(divs)) if len(divs)<=6 else []:
        pass
    # restricted search: up to 4 nonzero exponents
    for k in range(1,5):
        for ds in itertools.combinations(divs,k):
            for ex in itertools.product([e for e in range(-8,17) if e],repeat=k):
                if sum(ex)!=8 or sum(d*e for d,e in zip(ds,ex))!=24: continue
                spec=dict(zip(ds,ex))
                try: f=eta_q(spec,P+1)
                except AssertionError: continue
                if f[1]!=1: continue
                for (t,s),T in targets.items():
                    if all(T[p]==f[p] for p in ps):
                        key=(tuple(sorted(spec.items())),t,s)
                        if key not in found:
                            found.add(key); print(f"lam={lam}: c_p = a_p(eta {spec}) + {s}*({t}|p)*p")
print("done", len(found))
