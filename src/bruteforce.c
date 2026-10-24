/*
 * bruteforce.c — independent validation of the exact criterion used by
 * c27.c / es_big.c / bench.c.
 *
 * The other programs in this repo share a single `works()` predicate: given
 * n and a deficit delta, decide whether the delta-family solves 4/n by
 * enumerating reachable residues mod p = delta/gcd(delta,na).
 *
 * That predicate is arithmetic and never enumerated directly. The three
 * checks documented in the README (self-test, c27-vs-es_big cross-check,
 * prime count) all assume the predicate is correct: they catch coding
 * errors, not a wrong criterion.
 *
 * This program closes that gap. For small n it searches for solutions
 * 4/n = 1/a + 1/b + 1/c by direct enumeration — no residue sets, no
 * modular arithmetic beyond integer division — and reports the minimum
 * deficit actually observed. Compare its output against c27 on the same
 * range.
 *
 * METHOD
 *
 *   Fix a first. From 4/n = 1/a + 1/b + 1/c the remaining two terms are
 *
 *       1/b + 1/c = 4/n - 1/a = (4a - n) / (n a)
 *
 *   Write delta = 4a - n and na = n*a. Reduce by g = gcd(delta, na) to get
 *   coprime  r = delta/g,  B = na/g,  so  1/b + 1/c = r/B.
 *
 *   For fixed r and B, the solutions are exactly the divisors of... rather,
 *   for each candidate b the required c is
 *
 *       c = B b / (b r - B)
 *
 *   which is a positive integer exactly when  b r > B  and  (b r - B) | B b.
 *   The search window is  B/r < b <= 2B/r, since b <= c forces
 *   1/b >= r/(2B).
 *
 *   NOTE: this window is independent of sqrt(B). An earlier version bounded
 *   the scan by sqrt(B) and reported false negatives for n = 409 and
 *   n = 577, whose solutions have b >> sqrt(B). Do not reintroduce that
 *   bound.
 *
 *   a is scanned upward from floor(n/4)+1 so the first hit is the minimum
 *   delta, which is what C_k classifies.
 *
 * SCOPE
 *
 *   Deliberately slow and simple. Intended for n up to ~1e5, where the
 *   inner loop stays tractable. The point is to be obviously correct, not
 *   fast.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static inline uint64_t gcd64(uint64_t a, uint64_t b) {
    while (b) { uint64_t t = a % b; a = b; b = t; }
    return a;
}

/*
 * Smallest delta <= dmax admitting a solution, found by direct enumeration.
 * Returns 0 if none found within the cap.
 */
static uint32_t brute_min_delta(uint64_t n, uint32_t dmax) {
    for (uint32_t delta = 3; delta <= dmax; delta += 4) {
        if ((n + delta) % 4 != 0) continue;
        uint64_t a  = (n + delta) / 4;
        uint64_t na = n * a;
        uint64_t g  = gcd64(delta, na);
        uint64_t r  = delta / g;   /* reduced numerator   */
        uint64_t B  = na / g;      /* reduced denominator */

        if (r == 1) return delta;  /* identity solution */

        uint64_t b_lo = B / r + 1;
        uint64_t b_hi = (2 * B) / r;

        for (uint64_t b = b_lo; b <= b_hi; b++) {
            uint64_t den = b * r - B;
            if (den == 0) continue;
            uint64_t num = B * b;
            if (num % den) continue;
            uint64_t c = num / den;
            if (c >= b) return delta;
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    uint64_t N   = argc > 1 ? strtoull(argv[1], 0, 10) : 100000ULL;
    uint32_t cap = argc > 2 ? (uint32_t)strtoull(argv[2], 0, 10) : 107;

    printf("brute force over n < %llu, delta cap %u\n",
           (unsigned long long)N, cap);

    uint64_t checked = 0, unsolved = 0;
    uint32_t gmax = 0;
    uint64_t gmaxn = 0;

    /* n ≡ 1 (mod 4) only: the other classes reduce by the identities
     * described in the README, so they carry no depth information. */
    for (uint64_t n = 5; n < N; n += 4) {
        uint32_t d = brute_min_delta(n, cap);
        checked++;
        if (!d) { unsolved++; continue; }
        if (d > gmax) { gmax = d; gmaxn = n; }
    }

    printf("n checked: %llu | unsolved at cap: %llu\n",
           (unsigned long long)checked, (unsigned long long)unsolved);
    printf("MAX delta=%u at n=%llu\n", gmax, (unsigned long long)gmaxn);
    return 0;
}
