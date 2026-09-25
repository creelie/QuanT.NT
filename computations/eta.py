"""q-expansions of multiplicative eta-quotients of weight 4 (Martin's list) and weight 2."""
import numpy as np
def eta_prod(spec, N):
    # spec: {m: exponent}; returns coeffs of prod eta(m z)^e without the q^{sum m e/24} shift
    c=np.zeros(N,dtype=object); c[0]=1
    for m,e in spec.items():
        for _ in range(abs(e)):
            for n in range(1,N):
                k=n*m
                if k>=N: break
                if e>0:   # multiply by (1-q^k)
                    c[k:]=c[k:]-c[:N-k]
                else:     # divide by (1-q^k)
                    for i in range(k,N): c[i]+=c[i-k]
    return c
def eta_q(spec, N):
    shift=sum(m*e for m,e in spec.items())
    assert shift%24==0
    s=shift//24; c=eta_prod(spec,N)
    out=np.zeros(N,dtype=object); out[s:]=c[:N-s]; return out
FORMS4={'5':{1:4,5:4},'6':{1:2,2:2,3:2,6:2},'8':{2:4,4:4},'9':{3:8},
        '32':{4:16,2:-4,8:-4},'12':{2:4,6:4}}
FORMS2={'11':{1:2,11:2},'14':{1:1,2:1,7:1,14:1},'15':{1:1,3:1,5:1,15:1},'20':{2:2,10:2},
        '24':{2:1,4:1,6:1,12:1},'27':{3:2,9:2},'32':{4:2,8:2},'36':{6:4}}
def coeffs(spec,N):
    try: return eta_q(spec,N)
    except AssertionError: return None
