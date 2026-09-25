"""
Point counts for the 'Pauli-expectation' Calabi-Yau family
    V_lambda :  sum_{i=1}^4 (x_i + 1/x_i) = lambda      on the torus (F_p^*)^4,
the affine chart of the (2,2,2,2) hypersurface in (P^1)^4 (4-qubit product states).

N_p(lambda) = #{x in (F_p^*)^4 : sum f(x_i) = lambda},  f(x) = x + 1/x.
Value distribution of f: g(v) = #{x : x+1/x = v} = 1 + (v^2-4 | p).
N_p = g*g*g*g (additive convolution mod p).
"""
import numpy as np, sys

def primes(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**.5)+1):
        if s[i]: s[i*i::i] = False
    return [int(i) for i in np.nonzero(s)[0]]

def legendre_table(p):
    t = np.full(p, -1, dtype=np.int64); t[0] = 0
    t[(np.arange(1, p)**2) % p] = 1
    return t

def conv(a, b, p):
    # cyclic convolution over Z/p, exact integers
    out = np.zeros(p, dtype=np.int64)
    for s in range(p):
        out += a[s] * np.roll(b, s)
    return out

def counts(p):
    chi = legendre_table(p)
    v = np.arange(p)
    g = 1 + chi[(v*v - 4) % p]
    g2 = conv(g, g, p)
    return conv(g2, g2, p)     # N_p(lambda) for lambda = 0..p-1

if __name__ == "__main__":
    P = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    for p in primes(P):
        if p < 5: continue
        N = counts(p)
        print(p, {lam: int(N[lam % p]) for lam in (0, 1, 2, 3, 4, 5, 6, 8, -8, 12, 16)})
