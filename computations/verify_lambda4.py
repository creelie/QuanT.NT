"""
Exact checks for the lambda = +-4 fibre and the mod-p link lemma (paper/main.tex, Theorems D, E, Lemma 3.3).

  (L)  N_p(lam) == -7 - sum_{k<p} C(2k,k) D_k lam^{-2k}   (mod p)   for all lam in F_p^*   [proved; checked]
  (D)  N_p(+-4) = p^3 - 4p^2 + 4p(1 - (-1/p)) - 7 - a_p(f_12)
  (E)  sum_{k=0}^{p-1} C(2k,k) D_k / 16^k == a_p(f_12)   (mod p^3)

f_12 = q + 3q^3 - 18q^5 + 8q^7 + 9q^9 + 36q^11 - ...  is the unique newform in S_4(Gamma_0(12)) (LMFDB 12.4.a.a).
Requires numpy and cypari2 (pip install cypari2).  Usage: python3 verify_lambda4.py [bound] [bound for (L)]
"""
import sys
from math import comb
import cypari2
from pauli_cy3_counts import counts, primes

B = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
BL = int(sys.argv[2]) if len(sys.argv) > 2 else 200
pari = cypari2.Pari()
f = pari(f"L=mfinit([12,4],0); mfcoefs(mfeigenbasis(L)[1],{B})")
a = [int(f[n]) for n in range(B + 1)]
assert a[:8] == [0, 1, 0, 3, 0, -18, 0, 8]
D = [sum(comb(n, k)**2 * comb(2*k, k) * comb(2*n - 2*k, n - k) for k in range(n + 1)) for n in range(B)]
W = [comb(2*n, n) * D[n] for n in range(B)]

ps = [p for p in primes(B) if p >= 5]
failD, failE, modp4 = [], [], []
for p in ps:
    N = counts(p)
    rhs = p**3 - 4*p**2 + 4*p*(1 - (1 if p % 4 == 1 else -1)) - 7 - a[p]
    if int(N[4 % p]) != rhs or int(N[-4 % p]) != rhs: failD.append(p)
    M = p**4; z = pow(16, -1, M)
    s = sum(W[k] * pow(z, k, M) for k in range(p)) % M
    if (s - a[p]) % p**3: failE.append(p)
    if (s - a[p]) % p**4 == 0: modp4.append(p)
failL = 0
for p in [q for q in primes(BL) if q >= 3]:
    N = counts(p)
    for lam in range(1, p):
        if (int(N[lam]) + 7 + sum(W[k] * pow(lam, -2*k, p) for k in range(p))) % p: failL += 1
print(f"(D),(E): {len(ps)} primes 5<=p<={B}; failures D={failD} E={failE}; (E) also holds mod p^4 only at {modp4}")
print(f"(L): all lam in F_p^*, primes 3<=p<={BL}; violations: {failL}")
