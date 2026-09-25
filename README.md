# QuanT.NT — arithmetic of Calabi–Yau sections of the four-qubit QISM

A draft math.NT paper plus reproducible code. It connects the QuanTY / Quantum Inner State
Manifold (QISM) framework to number theory.

**Core idea.** The Y-qubit QISM is `M_Y = (CP^1)^Y` (product of Bloch spheres). As a variety
over Z, its anticanonical sections are Calabi–Yau. For Y = 4 they are exactly the quadric
sections of the separable-state (Segre) variety. The *Pauli family*

    V_λ :  Σ_i q(ψ_i) Π_{j≠i} m(ψ_j) = λ Π_j m(ψ_j),   q = ψᵀψ,  m = ½ ψᵀσ_xψ
           (torus chart:  Σ_{i=1}^4 (x_i + 1/x_i) = λ)

has point counts given by twisted Kloosterman moments and periods given by closed walks on Z⁴
(`C(2n,n)·Domb_n`).

## Results (exact integer checks for all primes 5 ≤ p ≤ 1000)

| | Statement | Status |
|---|---|---|
| Prop. 1 | No projective-bundle QISM total space is Calabi–Yau (c₁ restricted to a fibre is N·h ≠ 0) | proved |
| Prop. 2 | Calabi–Yau sections of `M_4` are exactly quadric sections of the 4-qubit separable-state variety | proved |
| A | `#V_0°(F_p) = p³ − 2p² − 7 − a_p(η(2τ)⁴η(4τ)⁴)` | verified, conjectural |
| A′ | `Σ_a (a/p) K(a)⁴ = −3p² − p·a_p(f_8)` (twisted 4th Kloosterman moment) | verified for p ≤ 400, conjectural |
| B | `#V_{±8}°(F_p) = p³ − 4p² + 3p − 7 − 3p·a_p(24a) − a_p((η₁η₂η₃η₆)²)` | verified, conjectural |
| C | `Σ_{k<p} C(2k,k) D_k / 64^k ≡ a_p((η₁η₂η₃η₆)²) (mod p³)`, fails mod p⁴ | verified, conjectural |

## Layout

- `paper/main.tex`: the arXiv draft (amsart). Not compiled in this environment, so compile locally.
- `computations/verify.py`: checks A, A′, B, C in exact arithmetic (`python3 verify.py 1000 400`, ~30 s).
- `computations/pauli_cy3_counts.py`: fast point counts via convolution of `x + 1/x` value distributions.
- `computations/eta.py`: q-expansions of eta quotients.
- `computations/supercong.py`: truncated-period congruence search.
- `computations/sym4.py`, `scan.py`, `etasearch.py`, `motive.py`: exploratory scripts (Sym⁴ Kloosterman splitting, λ-scans).

Requires `numpy`.
