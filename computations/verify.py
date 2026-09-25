"""
Exact-integer verification of every identity stated in paper/main.tex.

  (A)  N_p(0) = p^3 - 2p^2 - 7 - a_p(f_8)
  (A') sum_{a in F_p^*} (a/p) K(a)^4 = -3p^2 - p a_p(f_8)          [twisted 4th Kloosterman moment]
  (B)  N_p(+-8) = p^3 - 4p^2 + 3p - 7 - 3p a_p(g_24) - a_p(f_6)
  (C)  sum_{k=0}^{p-1} C(2k,k) D_k / 64^k == a_p(f_6)   (mod p^3)   [D_k = Domb numbers]

N_p(lam) = #{x in (F_p^*)^4 : sum_i (x_i + 1/x_i) = lam},   K(a) = sum_{x in F_p^*} e((x + a/x)/p),
  f_8  = eta(2z)^4 eta(4z)^4,  f_6 = (eta(z)eta(2z)eta(3z)eta(6z))^2,  g_24 = eta(2z)eta(4z)eta(6z)eta(12z).
Usage: python3 verify.py [bound for A,B,C] [bound for A']
"""
import sys
from math import comb
import numpy as np
from pauli_cy3_counts import counts, primes, legendre_table
from eta import eta_q

B = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
BK = int(sys.argv[2]) if len(sys.argv) > 2 else 400
f8 = eta_q({2: 4, 4: 4}, B + 1)
f6 = eta_q({1: 2, 2: 2, 3: 2, 6: 2}, B + 1)
g24 = eta_q({2: 1, 4: 1, 6: 1, 12: 1}, B + 1)
D = [sum(comb(n, k)**2 * comb(2*k, k) * comb(2*n - 2*k, n - k) for k in range(n + 1)) for n in range(B)]

def kloosterman_real(p):
    x = np.arange(1, p); inv = np.array([pow(int(a), p - 2, p) for a in x])
    a = np.arange(1, p)[:, None]
    return np.cos(2*np.pi*((x[None, :] + a*inv[None, :]) % p)/p).sum(1)   # K(a) is real

ps = [p for p in primes(B) if p >= 5]
fail = {k: [] for k in "ABC"}; failK = []
for p in ps:
    N = counts(p)
    if int(N[0]) != p**3 - 2*p**2 - 7 - int(f8[p]): fail["A"].append(p)
    rhs = p**3 - 4*p**2 + 3*p - 7 - 3*p*int(g24[p]) - int(f6[p])
    if int(N[8 % p]) != rhs or int(N[-8 % p]) != rhs: fail["B"].append(p)
    m = p**3; inv = pow(64, -1, m)
    if (sum(comb(2*k, k)*D[k]*pow(inv, k, m) for k in range(p)) - int(f6[p])) % m: fail["C"].append(p)
    if p <= BK:
        K = kloosterman_real(p); chi = legendre_table(p)[1:]
        if round(float((chi*K**4).sum())) != -3*p*p - p*int(f8[p]): failK.append(p)
print(f"(A),(B),(C): {len(ps)} primes 5<=p<={B}; failures: {fail}")
print(f"(A'): {len([p for p in ps if p <= BK])} primes 5<=p<={BK}; failures: {failK}")
