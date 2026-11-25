#!/usr/bin/env python3
"""Computations for Sections 10 and 11 of the manuscript (numbering as in manuscript.pdf).

Part A (Theorem 2.5, Table 2): the exact limit inferior and limit
superior of the log-concavity ratio R_{l,k}(n) for fixed k >= 3 equal
G_{L,k} * prod_p inf_v K_{p,k}(v) and G_{L,k} * prod_p sup_v K_{p,k}(v).
The table evaluates these Euler products numerically (mpmath, 40 digits,
primes below PRIME_CUTOFF, valuations 0..VMAX plus the v -> infinity limit).
These values are numerical evaluations, not interval certificates.  At l = 2
they are compared with Starr's closed formula (Corollary 2.4 of arXiv
2601.18599, multiplied by (k+1)/k to pass from k!A/n! to A).

Part B (Proposition 2.9, Section 11): rigorous interval-arithmetic certificates that
R_{l,2}(n) < 1 for explicit integers n at ranks l = 2, 3, 4.  The certificate
uses only Theorem 2.6 (with the bounds D_l, C_l <= e and U_l <= 6e^2 of
Section 8.1), the exact one-orbit identity, exact rational local factors for
primes up to PRIME_CERT, and the bound prod_{p > P} K_{p,2}(0) <= exp(2/P)
proved in Section 11.  Every floating step is carried out in mpmath interval
arithmetic (iv, 60 digits), so the printed upper bounds are rigorous.

The script exits with status 1 if any check fails (no use of assert, so it
also fails correctly under python -O).
"""
from fractions import Fraction
from math import factorial, isqrt
import json
import sys
import time

import mpmath as mp
from mpmath import iv

mp.mp.dps = 40
iv.dps = 60
PRIME_CUTOFF = 100000
PRIME_CERT = 100000
VMAX = 60
FAILURES = []


def check(name, ok):
    if not ok:
        FAILURES.append(name)
    print(('PASS: ' if ok else 'FAIL: ') + name, flush=True)


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, isqrt(n) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


PRIMES = primes_upto(max(PRIME_CUTOFF, PRIME_CERT))


# ---------- exact local data ------------------------------------------------
def weights_exact(ell, p, M):
    """w(e) = p^{-e} h_l(p^e), e = 0..M-1: coefficients of prod_{j=2}^{l} (1-p^{-j} z)^{-1}."""
    w = [Fraction(1)] + [Fraction(0)] * (M - 1)
    for j in range(2, ell + 1):
        r = Fraction(1, p ** j)
        for e in range(1, M):
            w[e] += r * w[e - 1]
    return w


def q_exact(ell, p, M):
    """Tail probabilities q_p(m), m = 0..M, as exact fractions."""
    T0 = Fraction(1)
    for j in range(2, ell + 1):
        T0 /= (1 - Fraction(1, p ** j))
    w = weights_exact(ell, p, M)
    q, T = [], T0
    for m in range(M + 1):
        q.append(T / T0)
        if m < M:
            T -= w[m]
    return q


def L_exact(ell, p, r, v):
    q = q_exact(ell, p, v + 1)
    return sum(Fraction(p) ** m * (q[m] ** r - q[m + 1] ** r) for m in range(v + 1))


def K_exact(ell, k, p, v):
    return L_exact(ell, p, k, v) ** 2 / (L_exact(ell, p, k - 1, v) * L_exact(ell, p, k + 1, v))


def G_exact(ell, k):
    L = ell - 1
    return Fraction(k + 1, k) * Fraction(factorial((k - 1) * L - 1) * factorial((k + 1) * L - 1),
                                         factorial(k * L - 1) ** 2)


# ---------- Part A: limit inferior and limit superior -----------------------
def q_float(ell, p, M):
    pp = mp.mpf(p)
    T0 = mp.fprod(1 / (1 - pp ** (-j)) for j in range(2, ell + 1))
    w = [mp.mpf(1)] + [mp.mpf(0)] * (M - 1)
    for j in range(2, ell + 1):
        r = pp ** (-j)
        for e in range(1, M):
            w[e] += r * w[e - 1]
    q, T = [], T0
    for m in range(M + 1):
        q.append(T / T0)
        if m < M:
            T -= w[m]
    return q


def K_profile(ell, k, p, vmax):
    """K_{p,k}(v) for v = 0..vmax, and the limit as v -> infinity."""
    M = vmax + (80 if p < 50 else (20 if p < 1000 else 8))
    q = q_float(ell, p, M)
    pp = mp.mpf(p)
    out = []
    Ls = {r: mp.mpf(0) for r in (k - 1, k, k + 1)}
    for m in range(M):
        for r in Ls:
            Ls[r] += pp ** m * (q[m] ** r - q[m + 1] ** r)
        if m <= vmax:
            out.append(Ls[k] ** 2 / (Ls[k - 1] * Ls[k + 1]))
    limit = Ls[k] ** 2 / (Ls[k - 1] * Ls[k + 1])
    return out, limit


def liminf_table():
    rows = []
    for ell in (2, 3, 4, 5, 6):
        for k in (3, 4, 5):
            G = G_exact(ell, k)
            loginf = mp.mpf(0)
            logsup = mp.mpf(0)
            nonmonotone = []
            for p in PRIMES:
                if p > PRIME_CUTOFF:
                    break
                vmax = VMAX if p < 50 else (12 if p < 1000 else 3)
                prof, lim = K_profile(ell, k, p, vmax)
                vals = prof + [lim]
                loginf += mp.log(min(vals))
                logsup += mp.log(max(vals))
                if any(vals[i + 1] > vals[i] * (1 + mp.mpf(10) ** -30) for i in range(len(vals) - 1)):
                    nonmonotone.append(p)
            liminf = mp.mpf(G.numerator) / G.denominator * mp.exp(loginf)
            limsup = mp.mpf(G.numerator) / G.denominator * mp.exp(logsup)
            row = {'ell': ell, 'k': k, 'G': str(G), 'G_float': mp.nstr(mp.mpf(G.numerator) / G.denominator, 12),
                   'liminf': mp.nstr(liminf, 12), 'limsup': mp.nstr(limsup, 12),
                   'primes_where_K_not_nonincreasing_in_v': nonmonotone}
            if ell == 2:
                z = mp.zeta
                starr = (mp.mpf(k) / (k - 1)) * z(2 * k + 2) * z(2 * k - 2) / z(2 * k) ** 2 \
                    * z(2 * k - 1) ** 2 / (z(2 * k - 3) * z(2 * k + 1))
                starr_R = starr * mp.mpf(k + 1) / k
                row['starr_closed_form_times_(k+1)/k'] = mp.nstr(starr_R, 12)
                rel = abs(liminf / starr_R - 1)
                row['relative_difference_from_starr'] = mp.nstr(rel, 3)
                check('l=2, k=%d: Euler-product lim inf equals Starr closed form (rel. diff < 1e-8)' % k, rel < mp.mpf(10) ** -8)
            rows.append(row)
            print('  l=%d k=%d  G=%s  liminf=%s  limsup=%s' % (ell, k, row['G_float'], row['liminf'], row['limsup']), flush=True)
    return rows


# ---------- Part B: certified two-orbit counterexamples ----------------------
def lo(x):
    return mp.mpf(x._mpi_[0])


def hi(x):
    return mp.mpf(x._mpi_[1])


def ivfrac(x):
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)


def certify(ell, exps):
    """Rigorous upper bound for R_{l,2}(n), n = prod p^e over exps (dict)."""
    L = ell - 1
    n = 1
    for p, e in exps.items():
        n *= p ** e
    G = G_exact(ell, 2)
    prod = ivfrac(G)
    for p in PRIMES:
        if p > PRIME_CERT:
            break
        v = exps.get(p, 0)
        if v == 0:
            q1 = 1 - Fraction(1)
            prodj = Fraction(1)
            for j in range(2, ell + 1):
                prodj *= 1 - Fraction(1, p ** j)
            q1 = 1 - prodj
            K0 = (1 - q1 ** 2) ** 2 / ((1 - q1) * (1 - q1 ** 3))
            prod *= ivfrac(K0)
        else:
            prod *= ivfrac(K_exact(ell, 2, p, v))
    # primes above PRIME_CERT do not divide n (checked below); their factor is <= exp(2/P).
    if max(exps) > PRIME_CERT:
        return None
    lead_upper = prod * iv.exp(iv.mpf(2) / PRIME_CERT)
    lead_lower = prod
    # Theorem 2.6 error bounds with D, C <= e and U <= 6 e^2 (Section 8.1).
    logn = iv.log(iv.mpf(n))
    H = 1 + logn
    e = iv.e

    def Ebound(r):
        ratio = Fraction(factorial(r * L - 1), factorial(L - 1) ** r)
        return (4 * e * e ** r * H ** r * ivfrac(ratio) + 2 * r * e * e ** (r - 1) * 6 * e ** 2) / iv.mpf(n)
    E2, E3 = Ebound(2), Ebound(3)
    R_upper = lead_upper * (1 + E2) ** 2 / (1 - E3)
    return {'ell': ell, 'log10_n': mp.nstr(mp.log10(n), 8), 'digits_of_n': len(str(n)),
            'factorisation': {str(p): e for p, e in sorted(exps.items())},
            'leading_ratio_interval': [mp.nstr(lo(lead_lower), 15), mp.nstr(hi(lead_upper), 15)],
            'E_l2_upper': mp.nstr(hi(E2), 5), 'E_l3_upper': mp.nstr(hi(E3), 5),
            'R_upper_bound': mp.nstr(hi(R_upper), 15), 'n': str(n) if len(str(n)) < 120 else 'see factorisation', 'certified_R_below_one': bool(hi(R_upper) < 1)}


def greedy_exponents(ell, max_steps=20000):
    """Greedy search (heuristic, not part of the certificate) for valuations driving
    the leading two-orbit ratio below one with small log n."""
    from math import log
    v = {}
    Kcache = {}

    def logK(p, e):
        if (p, e) not in Kcache:
            Kcache[(p, e)] = float(mp.log(K_profile_point(ell, p, e)))
        return Kcache[(p, e)]
    # leading log ratio at n = 1 (all v = 0)
    total = float(mp.log(mp.mpf(G_exact(ell, 2).numerator) / G_exact(ell, 2).denominator))
    for p in PRIMES[:20000]:
        total += logK(p, 0)
    steps = 0
    while total >= -1e-9 and steps < max_steps:
        best = None
        for p in PRIMES:
            e = v.get(p, 0)
            d = (logK(p, e + 1) - logK(p, e)) / log(p)
            if best is None or d < best[0]:
                best = (d, p)
            if e == 0:
                break
        d, p = best
        total += d * log(p)
        v[p] = v.get(p, 0) + 1
        steps += 1
    return v


def K_profile_point(ell, p, e):
    M = e + (80 if p < 50 else 10)
    q = q_float(ell, p, M)
    pp = mp.mpf(p)
    Ls = {}
    for r in (1, 2, 3):
        Ls[r] = mp.fsum(pp ** m * (q[m] ** r - q[m + 1] ** r) for m in range(e + 1))
    return Ls[2] ** 2 / (Ls[1] * Ls[3])


def main():
    t0 = time.time()
    out = {'scope': 'Part A: numerical evaluation (mpmath %d digits, primes < %d). Part B: interval-arithmetic '
                    'certificates (iv, %d digits) resting on Theorem 2.4.' % (mp.mp.dps, PRIME_CUTOFF, iv.dps)}
    # sanity: exact K agrees with float K
    for ell in (2, 3, 4):
        for p in (2, 3, 7):
            for v in (0, 1, 5):
                a = K_exact(ell, 3, p, v)
                b = K_profile(ell, 3, p, v)[0][v]
                check('exact vs float local factor l=%d p=%d v=%d' % (ell, p, v),
                      abs(mp.mpf(a.numerator) / a.denominator - b) < mp.mpf(10) ** -30)
    # (9.9): K_{p,k}(0) > 1
    for ell in (2, 3, 4):
        for k in (2, 3, 4):
            check('K_{p,%d}(0) > 1 for l=%d, p<100' % (k, ell),
                  all(K_exact(ell, k, p, 0) > 1 for p in PRIMES[:25]))
    print('Part A: limit inferior / superior table', flush=True)
    out['liminf_table'] = liminf_table()
    for row in out['liminf_table']:
        check('lim inf > 1 for l=%d, k=%d' % (row['ell'], row['k']), mp.mpf(row['liminf']) > 1)
    print('Part B: certified two-orbit counterexamples', flush=True)
    certs = []
    # rank 2: Charlton-Heim-Stumpenhusen first counterexample (consistency check of the method).
    lam = {2: 6, 3: 3, 5: 2, 7: 2, 11: 1, 13: 1, 17: 1, 19: 1, 23: 1, 29: 1}
    c2 = certify(2, lam)
    c2['note'] = 'n = 65 214 507 758 400 (first rank-two counterexample, Charlton et al. 2026a)'
    certs.append(c2)
    # Heuristic scale check (Section 11): where the greedy search first drives the
    # rank-two leading ratio below one (compare the true first failure, about 6.5e13).
    g2 = greedy_exponents(2)
    n2 = 1
    for p, e in g2.items():
        n2 *= p ** e
    out['greedy_rank_two_log10_n'] = mp.nstr(mp.log10(n2), 6)
    print('  greedy search, rank two: log10 n = %s' % out['greedy_rank_two_log10_n'], flush=True)
    for ell in (3, 4):
        exps = greedy_exponents(ell)
        c = certify(ell, exps)
        c['note'] = 'valuations found by a greedy heuristic search; not claimed to be minimal'
        certs.append(c)
    for c in certs:
        print('  l=%d log10 n=%s  R <= %s  certified=%s' % (c['ell'], c['log10_n'], c['R_upper_bound'],
                                                          c['certified_R_below_one']), flush=True)
        check('certified R_{%d,2}(n) < 1 at log10 n = %s' % (c['ell'], c['log10_n']), c['certified_R_below_one'])
    out['two_orbit_certificates'] = certs
    out['elapsed_seconds'] = round(time.time() - t0, 1)
    out['failures'] = FAILURES
    path = sys.argv[sys.argv.index('--output') + 1] if '--output' in sys.argv else 'results/extensions.json'
    with open(path, 'w') as fh:
        json.dump(out, fh, indent=2)
        fh.write('\n')
    print('ALL CHECKS PASSED' if not FAILURES else 'FAILURES: %s' % FAILURES)
    sys.exit(1 if FAILURES else 0)


if __name__ == '__main__':
    main()
