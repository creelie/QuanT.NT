from pauli_cy3_counts import counts, primes
import sys
P=int(sys.argv[1]); lams=[int(x) for x in sys.argv[2].split(',')]
print('p', *lams)
for p in primes(P):
    if p<5: continue
    N=counts(p); base=p**3-4*p**2+6*p-4
    print(p, *[int(N[l%p])-base for l in lams])
