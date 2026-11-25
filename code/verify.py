#!/usr/bin/env python3
"""Finite checks for "Commuting permutations with few joint orbits"
(equation and theorem numbers as in manuscript.pdf).

Suites (each prints PASS/FAIL; exit status 1 on any failure; no assert, so the
script also fails correctly under python -O):

 1. enumeration   Direct enumeration of commuting l-tuples in S_n (n <= 5, l <= 4)
                  against the coefficient formula (2.2); the first rows (Table 3).
 2. identities    Exact rational identities through n = 40, l <= 6: the two
                  Dirichlet-convolution descriptions (2.1)/(2.4), integrality,
                  and the cubic lower bound of Section 6.1.
 3. asymptotic    Exact integers A(l,n,k) (exponential-formula recursion) for
                  l = 3, 4, k = 1..5, n <= 1000, compared with the main term of
                  Theorem 2.2 (Euler product over primes < 10^4, 40 digits);
                  the explicit bound (2.14) is checked at every sampled point.
 4. constants     Exact rational checks of the constants in Sections 7-8.
 5. thresholds    Smallest n (as log10) at which the finite criterion
                  64 k E_{l,k+1}(n) <= 1 of Theorem 2.6 holds, for a few (l,k).

Numerical items (suite 3, 5) use mpmath with 40 digits; they are diagnostics,
not proofs.  Exact items use Python integers and fractions.
"""
from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import permutations
from math import comb, factorial, isqrt
import json
import sys
import time

import mpmath as mp

mp.mp.dps = 40
FAILURES = []
OUT = {}


def check(name, ok):
    if not ok:
        FAILURES.append(name)
    print(('PASS: ' if ok else 'FAIL: ') + name, flush=True)


def sublattice_counts(ell, N):
    """B_l(n) = (id^0 *_D id^1 *_D ... *_D id^{l-1})(n), n = 0..N (integers)."""
    prev = [0] + [1] * N
    for power in range(1, ell):
        cur = [0] * (N + 1)
        for d in range(1, N + 1):
            dp = d ** power
            for m in range(1, N // d + 1):
                cur[d * m] += dp * prev[m]
        prev = cur
    return prev


# ---------------------------------------------------------------- suite 1
def direct_counts(n, max_rank=4):
    perms = list(permutations(range(n)))
    size = len(perms)
    commuting = []
    for a in perms:
        mask = 0
        for j, b in enumerate(perms):
            if all(a[b[i]] == b[a[i]] for i in range(n)):
                mask |= 1 << j
        commuting.append(mask)

    @lru_cache(None)
    def joined(partition, j):
        parent = list(range(n))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for i in range(n):
            x, y = find(i), find(partition[i])
            if x != y:
                parent[y] = x
        for i, t in enumerate(perms[j]):
            x, y = find(i), find(t)
            if x != y:
                parent[y] = x
        groups = {}
        for i in range(n):
            groups.setdefault(find(i), []).append(i)
        res = [0] * n
        for g in groups.values():
            for i in g:
                res[i] = min(g)
        return tuple(res)
    counts = [defaultdict(int) for _ in range(max_rank + 1)]

    def visit(depth, allowed, partition):
        counts[depth][len(set(partition))] += 1
        if depth == max_rank:
            return
        bits = allowed
        while bits:
            bit = bits & -bits
            j = bit.bit_length() - 1
            visit(depth + 1, allowed & commuting[j], joined(partition, j))
            bits -= bit
    visit(0, (1 << size) - 1, tuple(range(n)))
    return counts


def A_table(ell, N, K):
    """Exact A(l,n,k), 0<=n<=N, 0<=k<=K, from the exponential formula (2.3):
    A(n,k) = sum_{m=1}^{n} C(n-1,m-1) (m-1)! B(m) A(n-m,k-1)."""
    B = sublattice_counts(ell, N)
    A = [[0] * (K + 1) for _ in range(N + 1)]
    A[0][0] = 1
    fact = [factorial(i) for i in range(N + 1)]
    for n in range(1, N + 1):
        for k in range(1, min(n, K) + 1):
            s = 0
            for m in range(1, n - k + 2):
                s += comb(n - 1, m - 1) * fact[m - 1] * B[m] * A[n - m][k - 1]
            A[n][k] = s
    return A


def suite_enumeration():
    rows = []
    ok = True
    for n in range(1, 6):
        counts = direct_counts(n)
        for ell in range(1, 5):
            A = A_table(ell, n, n)
            row = [A[n][k] for k in range(1, n + 1)]
            direct = [counts[ell].get(k, 0) for k in range(1, n + 1)]
            ok = ok and row == direct
            rows.append({'n': n, 'ell': ell, 'A_k_1_to_n': row})
    check('enumeration: formula (2.2)/(2.3) equals direct count for n<=5, l<=4 (60 values)', ok)
    table = {(r['ell'], r['n']): r['A_k_1_to_n'] for r in rows}
    expected = {(2, 3): [8, 9, 1], (2, 4): [42, 59, 18, 1], (2, 5): [144, 450, 215, 30, 1],
                (3, 3): [26, 21, 1], (3, 4): [210, 251, 42, 1], (3, 5): [744, 2870, 995, 70, 1],
                (4, 3): [80, 45, 1], (4, 4): [930, 995, 90, 1], (4, 5): [3744, 16650, 4175, 150, 1]}
    check('enumeration: Table 2 rows reproduced', all(table[key] == val for key, val in expected.items()))
    OUT['enumeration'] = rows


# ---------------------------------------------------------------- suite 2
def suite_identities():
    N = 40
    ok_int, ok_dir, ok_low = True, True, True
    for ell in range(1, 7):
        B = sublattice_counts(ell, N)
        b = [Fraction(0)] + [Fraction(B[n], n) for n in range(1, N + 1)]
        A = A_table(ell, N, N)
        # integrality and agreement of (2.2) with the recursion
        powers = [[Fraction(int(n == 0)) for n in range(N + 1)]]
        for k in range(1, 9):
            prev = powers[-1]
            powers.append([sum((prev[n - j] * b[j] for j in range(1, n + 1)), Fraction(0)) for n in range(N + 1)])
        for n in range(1, 21):
            for k in range(1, min(n, 8) + 1):
                val = powers[k][n] * Fraction(factorial(n), factorial(k))
                ok_int = ok_int and val.denominator == 1 and val.numerator == A[n][k]
        if ell >= 2:
            a = ell - 2
            g = [Fraction(0)] + [Fraction(1)] * N
            for j in range(1, ell):
                h = [Fraction(0)] + [Fraction(1, n ** j) for n in range(1, N + 1)]
                new = [Fraction(0)] * (N + 1)
                for d in range(1, N + 1):
                    for m in range(1, N // d + 1):
                        new[d * m] += g[d] * h[m]
                g = new
            ok_dir = ok_dir and all(b[n] == n ** a * g[n] for n in range(1, N + 1))
            s3 = powers[3]
            ok_low = ok_low and all(s3[n] >= factorial(a) ** 3 * comb(n - 1, 3 * a + 2) for n in range(3 * a + 3, N + 1))
    check('identities: (2.2) is an integer and equals the exponential-formula count (n<=20, k<=8, l<=6)', ok_int)
    check('identities: b_l(n) = n^{l-2} g_l(n), (2.4), through n=40, 2<=l<=6', ok_dir)
    check('identities: cubic lower bound of Section 6.1 through n=40, 2<=l<=6', ok_low)


# ---------------------------------------------------------------- suite 3
def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


PR = primes_upto(10000)


def local_L(ell, p, r, v):
    pp = mp.mpf(p)
    T0 = mp.fprod(1 / (1 - pp ** (-j)) for j in range(2, ell + 1))
    M = v + 2
    w = [mp.mpf(1)] + [mp.mpf(0)] * M
    for j in range(2, ell + 1):
        rr = pp ** (-j)
        for e in range(1, M + 1):
            w[e] += rr * w[e - 1]
    q, T = [], T0
    for m in range(M + 1):
        q.append(T / T0)
        T -= w[m]
    return mp.fsum(pp ** m * (q[m] ** r - q[m + 1] ** r) for m in range(v + 1))


def singular(ell, r, n):
    """S_{l,r}(n) = C_l^r prod_p L_{p,r}(v_p(n)), primes < 10^4 (truncation lowers it by < 1e-7 relative)."""
    C = mp.fprod(mp.zeta(j) for j in range(2, ell + 1))
    val = C ** r
    for p in PR:
        v, m = 0, n
        while m % p == 0:
            m //= p
            v += 1
        val *= local_L(ell, p, r, v)
    return val


def E_bound(ell, r, n):
    L = ell - 1
    D = mp.fprod(mp.zeta(j) for j in range(2, ell))
    C = mp.fprod(mp.zeta(j) for j in range(2, ell + 1))
    # U_l = sum_d h_l(d)^2 computed numerically (Euler product of squares).
    U = mp.mpf(1)
    for p in PR:
        # h(p^e): complete homogeneous polynomial of degree e in p^{-1},...,p^{-(l-1)}
        hs = [mp.mpf(1)] + [mp.mpf(0)] * 39
        for j in range(1, ell):
            x = mp.mpf(p) ** (-j)
            for e in range(1, 40):
                hs[e] += x * hs[e - 1]
        U *= mp.fsum(h * h for h in hs)
    H = 1 + mp.log(n)
    beta_inv = mp.gamma(r * L) / mp.gamma(L) ** r
    return (4 * mp.e * D ** r * H ** r * beta_inv + 2 * r * D * C ** (r - 1) * U) / n, U


def suite_asymptotic():
    rows = []
    ok_bound = True
    worst = {}
    for ell in (3, 4):
        L = ell - 1
        A = A_table(ell, 1000, 5)
        for k in (2, 3, 4):
            beta = mp.gamma(L) ** k / mp.gamma(k * L)
            for n in (100, 200, 500, 997, 1000):
                main = mp.mpf(factorial(n)) / factorial(k) * beta * mp.mpf(n) ** (k * L - 1) * singular(ell, k, n)
                rel = mp.mpf(A[n][k]) / main - 1
                E, U = E_bound(ell, k, n)
                ok_bound = ok_bound and abs(rel) <= E
                rows.append({'ell': ell, 'k': k, 'n': n, 'relative_error': mp.nstr(rel, 8),
                             'scaled_n_rel_over_(1+log n)^k': mp.nstr(n * abs(rel) / (1 + mp.log(n)) ** k, 6),
                             'E_bound_(2.13)': mp.nstr(E, 6)})
            for n in (1000,):
                R = mp.mpf(A[n][k]) ** 2 / (mp.mpf(A[n][k - 1]) * A[n][k + 1])
                worst[(ell, k)] = mp.nstr(R, 10)
            print('  asymptotic l=%d k=%d done' % (ell, k), flush=True)
    check('asymptotic: |A/main - 1| <= E_{l,k}(n) of (2.13) at all 30 sampled points', ok_bound)
    small = all(abs(mp.mpf(r['relative_error'])) < mp.mpf('0.02') for r in rows if r['n'] >= 997)
    check('asymptotic: relative error below 2% at n>=997 (l=3,4; k=2,3,4)', small)
    OUT['asymptotic'] = rows
    OUT['exact_R_at_n_1000'] = {'l=%d,k=%d' % key: val for key, val in worst.items()}


# ---------------------------------------------------------------- suite 4
def suite_constants():
    F = Fraction
    # (7.4), (7.5)
    check('Section 7.2 (first sum): (27/8)(3^-6 + 3^-5/5) = 1/135', F(27, 8) * (F(1, 3 ** 6) + F(1, 5 * 3 ** 5)) == F(1, 135))
    check('Section 7.2 (second sum): (3/2)(3^-3 + 3^-2/2) = 5/36', F(3, 2) * (F(1, 27) + F(1, 18)) == F(5, 36))
    check('(7.4) first bound at t=3: (1/135)(216/215) = 8/1075', F(1, 135) * F(216, 215) == F(8, 1075))
    check('(7.4) second bound at r=2: (5/36)(12/11) = 5/33', F(5, 36) * F(12, 11) == F(5, 33))
    # Section 7.4
    G23 = F(4, 3) * F(factorial(3) * factorial(7), factorial(5) ** 2)
    check('G_{2,3} = 14/5', G23 == F(14, 5))
    P3, Q2, Q4 = F(7, 8), F(3, 2), F(15, 14)
    check('(2.9) k=3 prefactor: G_{2,3} P3^2/(Q2 Q4) = 2401/1800', G23 * P3 ** 2 / (Q2 * Q4) == F(2401, 1800))
    s = F(16, 1075) + F(5, 33) + F(5, 1188)
    check('Section 7.4: 16/1075 + 5/33 + 5/1188 < 7/40', s < F(7, 40))
    check('Section 7.4: (2401/1800)(33/40) = 79233/72000 > 11/10', F(2401, 1800) * F(33, 40) == F(79233, 72000) > F(11, 10))
    check('Section 7.5: 4/15 + 16/1075 < 3/10', F(4, 15) + F(16, 1075) < F(3, 10))
    check('Section 7.5: 1/3 + 5/33 < 1/2', F(1, 3) + F(5, 33) < F(1, 2))
    check('Section 7.5: 1/15 + 5/1188 < 3/40', F(1, 15) + F(5, 1188) < F(3, 40))
    check('Section 7.5: 3/k - 7/2^k > 0 for 4<=k<=200', all(F(3, k) > F(7, 2 ** k) for k in range(4, 201)))
    check('Section 7.3: ((k+1)/k)^4 >= G_{2,k} for 3<=k<=200',
          all(F(k + 1, k) ** 4 >= F((k + 1) * (2 * k + 1), (k - 1) * (2 * k - 1)) for k in range(3, 201)))
    check('Section 7.6: (2^-6 + 2^-5/5)(64/63)*2 = (128/45) 4^-3',
          (F(1, 64) + F(1, 5 * 32)) * F(64, 63) * 2 == F(896, 315) * F(1, 64))
    check('Section 7.6: (2^-3 + 2^-2/2)(8/7) = (32/7) 4^-2', (F(1, 8) + F(1, 8)) * F(8, 7) == F(32, 7) * F(1, 16))
    # (8.1) margins
    check('(8.1): 3/k - 7/2^k >= 1/k for 4<=k<=200', all(F(3, k) - F(7, 2 ** k) >= F(1, k) for k in range(4, 201)))
    check('(8.1): log(11/10) > 1/12', mp.log(mp.mpf(11) / 10) > mp.mpf(1) / 12)
    check('(8.1): 2/k - 25/4^k >= 1/(4k) for 3<=k<=200', all(F(2, k) - F(25, 4 ** k) >= F(1, 4 * k) for k in range(3, 201)))
    # (8.4) multinomial gamma bound, exact
    ok = all(Fraction(factorial(r * L - 1), factorial(L - 1) ** r) <= r ** (r * L - 1) * L ** (r - 1)
             for L in range(1, 9) for r in range(2, 21))
    check('Section 8.1 multinomial bound: Gamma(rL)/Gamma(L)^r <= r^{rL-1} L^{r-1}, 1<=L<=8, 2<=r<=20', ok)
    # (9.12) alpha_{2,2} = 30/7
    val = 3 * mp.zeta(2) * mp.zeta(6) / mp.zeta(4) ** 2
    check('Section 9.3: 3 zeta(2) zeta(6)/zeta(4)^2 = 30/7', abs(val - mp.mpf(30) / 7) < mp.mpf(10) ** -35)
    # Section 8.1 integral of (1+log x)^2/x^2 on [1,inf) = 5
    check('Section 8.1: integral of (1+log x)^2/x^2 over [1,inf) = 5',
          abs(mp.quad(lambda x: (1 + mp.log(x)) ** 2 / x ** 2, [1, mp.inf]) - 5) < mp.mpf(10) ** -30)


# ---------------------------------------------------------------- suite 5
def suite_thresholds():
    res = {}
    for ell, k in ((2, 3), (3, 3), (3, 4), (4, 3)):
        L = ell - 1
        D = mp.fprod(mp.zeta(j) for j in range(2, ell))
        C = mp.fprod(mp.zeta(j) for j in range(2, ell + 1))
        _, U = E_bound(ell, 2, 10)
        r = k + 1
        beta_inv = mp.gamma(r * L) / mp.gamma(L) ** r

        def crit(logn):
            H = 1 + logn
            E = (4 * mp.e * D ** r * H ** r * beta_inv + 2 * r * D * C ** (r - 1) * U) * mp.exp(-logn)
            return 64 * k * E - 1
        lo, hi = mp.mpf(1), mp.mpf(400)
        for _ in range(200):
            mid = (lo + hi) / 2
            if crit(mid) > 0:
                lo = mid
            else:
                hi = mid
        res['l=%d,k=%d' % (ell, k)] = mp.nstr(hi / mp.log(10), 5)
    OUT['finite_criterion_log10_threshold'] = res
    print('  thresholds (log10 n):', res, flush=True)
    check('thresholds computed', len(res) == 4)


def main():
    t0 = time.time()
    for name, fn in (('enumeration', suite_enumeration), ('identities', suite_identities),
                     ('asymptotic', suite_asymptotic), ('constants', suite_constants),
                     ('thresholds', suite_thresholds)):
        print('== suite', name, flush=True)
        fn()
    OUT['elapsed_seconds'] = round(time.time() - t0, 1)
    OUT['failures'] = FAILURES
    path = sys.argv[sys.argv.index('--output') + 1] if '--output' in sys.argv else 'results/verification.json'
    with open(path, 'w') as fh:
        json.dump(OUT, fh, indent=2)
        fh.write('\n')
    print('ALL CHECKS PASSED' if not FAILURES else 'FAILURES: %s' % FAILURES)
    sys.exit(1 if FAILURES else 0)


if __name__ == '__main__':
    main()
