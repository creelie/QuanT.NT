"""
Periods of V_lambda: constant terms W_{2n} = CT[(sum_i (x_i+1/x_i))^{2n}] = # closed 2n-step walks on Z^4.
Test truncated-period congruences  S_p(t) = sum_{n=0}^{p-1} W_{2n} t^{2n}  mod p^r  against a_p(f).
"""
import sys
from math import comb
from eta import eta_q
from pauli_cy3_counts import primes
M = 2*int(sys.argv[1]) if len(sys.argv) > 1 else 800
# 1-D walks: c1(2k) = C(2k,k). Z^4 walks via EGF product: W_{2n} = sum over 2k1+..+2k4 = 2n of (2n)!/prod(2ki)! prod C(2ki,ki)
def walks(nmax):
    # a_k = C(2k,k)/(2k)!  -> use integer convolution of multinomials instead
    # 2-D: W2(2n)=C(2n,n)^2 ; Z^4 = convolution of two Z^2 via multinomial
    W2 = [comb(2*n, n)**2 for n in range(nmax+1)]
    return [sum(comb(2*n, 2*k)*W2[k]*W2[n-k] for k in range(n+1)) for n in range(nmax+1)]
W = walks(M//2)
print("W_2n:", W[:8])
F = {'f6': eta_q({1:2,2:2,3:2,6:2}, M+2), 'f8': eta_q({2:4,4:4}, M+2)}
for lam in (8, 4, 2, 6, 12, 16):
    for name, f in F.items():
        for r in (1, 2, 3):
            ok = True; checked = 0
            for p in primes(M//2):
                if p < 5 or lam % p == 0: continue
                inv = pow(lam*lam, -1, p**r)
                S = sum(W[n]*pow(inv, n, p**r) for n in range(p)) % p**r
                a = int(f[p])
                if (S - a) % p**r and (S + a) % p**r: ok = False; break
                checked += 1
            if ok: print(f"lambda={lam}: sum_(n<p) W_2n/lambda^(2n) == ±a_p({name}) mod p^{r}  [{checked} primes]")
