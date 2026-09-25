"""Blind scan: for Y qubits, is S_p = sum_{k<p} W^{(Y)}_{2k} m^{-k} (W = C(2k,k) A_Y(k), kind W) or sum A_Y(k) m^{-k} (kind A)
anomalously small mod p^2 / p^3?  Prints (Y, kind, m, mean log_p|S_p mod p^2|, mean log_p|S_p mod p^3|); random ~1.6 / ~2.7.
Genuine mod-p^3 signals in Y=4: m=16 (Theorem E) and m=64 (Theorem C). Slow (~30 min)."""
import math, sys
from math import comb, log
def primes(n):
    s=bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
N=260
A={1:[1]*N}
for Y in range(2,8):
    A[Y]=[sum(comb(n,k)**2*A[Y-1][k] for k in range(n+1)) for n in range(N)]
C2=[comb(2*n,n) for n in range(N)]
PS=[p for p in primes(N) if p>=60]
res=[]
for Y in range(2,8):
  for kind in ("W","A"):
    seq=[C2[n]*A[Y][n] if kind=="W" else A[Y][n] for n in range(N)]
    for m in list(range(-128,129)):
        if m in (0,): continue
        es={2:[],3:[]}
        for p in PS:
            if m%p==0: continue
            for r in (2,3):
                M=p**r; z=pow(m,-1,M); s=0; zk=1
                for k in range(p):
                    s=(s+seq[k]*zk)%M; zk=zk*z%M
                if s>M//2: s-=M
                es[r].append(log(abs(s)+1)/log(p))
        e2=sum(es[2])/len(es[2]); e3=sum(es[3])/len(es[3])
        if e2<1.6 or e3<2.4:
            res.append((Y,kind,m,round(e2,2),round(e3,2))); print(Y,kind,m,round(e2,2),round(e3,2),flush=True)
