# Erdős–Straus class-depth census

An independent, from-scratch implementation of the deficit/class framework for the
Erdős–Straus conjecture, and a census of the minimum class `C_k` for primes up to 10^12.

## Result

- **No prime below 10^12 requires `k > 27`.** This is established in two parts: the
  36 residue classes mod 9240 are swept exhaustively below 10^12 (`c27`), and every
  prime outside those classes is already in `C_3` by Ionascu & Wilson, Theorem 2.2
  (see *Method*). No prime in the swept set needed δ > 107.
- The deepest known class remains **k = 27**, held by `n = 8,803,369` (δ = 107).
- For `10^10 ≤ n < 10^12`, the maximum is **k = 15** (δ = 59, at n = 514,597,369,441).
- Prior published census: 4,146,894,049 (Ionascu & Wilson, 2010).

## The mathematics

`4/n = 1/a + 1/b + 1/c`, for integers `n ≥ 2`.

Even `n` reduces by a three-term identity; `n ≡ 3 (mod 4)` falls to a closed form;
multiplicativity reduces the general case to primes. What remains is `n ≡ 1 (mod 4)`,
in practice primes.

Write `a = (n + δ)/4`, so `δ = 4a − n ≡ 3 (mod 4)` is the **deficit**. Define

    C_k = { n : a solution exists with a ≤ (n + 4k − 1)/4 },   i.e. δ ≤ 4k − 1.

**Exact criterion.** For `n ≡ 1 (mod 4)` and `δ ≡ 3 (mod 4)`, put `g = gcd(δ, na)`,
`p = δ/g`, `q = na/g`. The δ-family solves `4/n` **iff** `p = 1`, or `q²` has a
divisor `D` with `D ≡ −q (mod p)`.

No search is required: the test is arithmetic.

## Method

By Ionascu & Wilson, Theorem 2.2, every prime outside 36 residues mod 9240 is already
in `C_3` (or `C_8`/`C_40` for the two special residues 1201 and 6001, which are
included explicitly). Therefore **no prime outside those 36 classes can have k > 27**,
and the census only needs to sweep them.

Each residue class is an arithmetic progression. The sieve runs inside the progression:
for a prime `p`, the multiples in the class form a single progression, so the step in
`k = (n−r)/9240` is exactly `p`. Memory stays at ~1 MB per thread regardless of range.

## Validation

Three checks, described below with their limits:

1. **Self-test** — the engine reproduces `8,803,369 → δ = 107` and confirms no smaller δ works.
2. **Cross-check** — on `[10^9, 10^10)` `c27` returns max δ = 83, matching `es_big`
   run exhaustively over all `n ≡ 1 (mod 4)`. These two programs share the `works()`
   criterion but factor and sieve differently, so this check catches implementation
   errors, not a wrong criterion. The `es_big` run also reports
   `cofactor-survived-trial-division`, which is `0` across `[2, 10^10)`.
3. **Prime count** — 696,585,317 primes found, versus ~6.97×10^8 predicted by the
   prime number theorem in arithmetic progressions (agreement to 0.02%).

**Known gap:** none of the three checks validates the exact criterion itself against
first principles. A brute-force enumeration of `(a, b, c)` for small `n` would close
this, and is not yet implemented.

## Build and run

    make
    ./c27 10000000000 1000000000000

`es_big` is the exhaustive deficit census over all `n ≡ 1 (mod 4)` (used up to 10^10).

## Scope

Primes only (as in the prior published census). Composites are not classified.
This is a census, not a proof of the conjecture.

## References

- L. J. Mordell, *Diophantine Equations* (1967) — mod-840 exceptions.
- E. Ionascu, A. Wilson, *On the Erdős–Straus conjecture*, arXiv:1001.1100 (2010).
- S. E. Salez, arXiv:1406.6307 (2014) — verification to 10^17.

## Attribution

The mathematics is due to Mordell and to Ionascu & Wilson. This repository is an
independent re-derivation and implementation, written with an AI assistant and
validated against the published record and an independent exhaustive census.
