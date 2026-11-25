# AI index: Commuting permutations with few joint orbits (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1 is
a plain introduction with the question, prior work, the results, the mechanism and a roadmap).
DOI 10.5281/zenodo.24541319 (concept 10.5281/zenodo.24541318). Status: unrefereed candidate; internal AI-assisted
review and computer checks only.

## Normalisation

A(l,n,k) = number of ordered l-tuples of pairwise commuting permutations of [n] with exactly k joint orbits (orbits of
the generated group). l = rank, L = l - 1. A(1,n,k) = unsigned Stirling numbers of the first kind; A(2,n,k) =
normalised d'Arcais coefficients. B_l(n) = number of subgroups of index n in Z^l; A(l,n,1) = (n-1)! B_l(n).
R_{l,k}(n) = A(l,n,k)^2 / (A(l,n,k-1) A(l,n,k+1)); the conjecture (Abdesselam, Conjecture 1.1) is R_{l,k}(n) >= 1 for
2 <= k <= n-1. Starr's normalisation alpha(n,k) = k! A(2,n,k)/n! has ratio (k/(k+1)) R_{2,k}(n) (Remark 2.1).
S_{l,k}(n) = singular series (eq. 2.5), bounded between 1 and C_{l,k}; beta_{l,k} = Gamma(L)^k / Gamma(kL).

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| A | A(l,n,k) = (n!/k!) beta_{l,k} n^{kL-1} [S_{l,k}(n) + O((1+log n)^k / n)], fixed l >= 2, k >= 2, uniform in the factorisation of n | Theorem 2.2 | written proof (Sections 3-4); qualitative leading term follows the Ingham-Halberstam route (Lemke Oliver-Shrestha-Thorne; Starr at rank two) |
| B | Two-orbit failure in every rank: R_{l,2}(m!) -> 0; for every eps > 0, R_{l,2}(n) < eps on all large multiples of some Q | Theorem 2.3 | written proof (Section 6); rank two due to Starr |
| C | Eventual strict log-concavity at every fixed k >= 3; explicit uniform bounds (> 1.10 at k = 3 for l >= 3) | Theorem 2.4 | written proof (Section 7); rank two due to Starr |
| D | Exact liminf and limsup of R_{l,k}(n), k >= 3, as G_{L,k} prod_p inf/sup_v K_{p,k}(v); rank two gives Starr Corollary 2.4 times (k+1)/k | Theorem 2.5; Table 2 | written proof (Section 10); table values are numerical evaluations |
| E | Explicit error bound E_{l,r}(n) and finite criterion 64k E_{l,k+1}(n) <= 1 => R_{l,k}(n) >= exp(1/(8k)) | Theorem 2.6 | written proof (Section 8) |
| F | Strict log-concavity for all 3 <= k <= c log n / log log n, c < 1/l, n large | Corollary 2.7 | written proof (Section 8) |
| G | Continuous limiting distribution of R_{l,k}(n) for fixed l, k; density of two-orbit failures exists, 0 < delta_{l,2}(1) < 1 | Theorem 2.8 | written proof (Section 9); applies the Erdos-Wintner theorem; existence new also at rank two (Charlton-Heim-Stumpenhusen gave bounds) |
| H | Certified failures: R_{3,2}(n_3) < 0.9926, n_3 = 2^9 3^5 5^3 7^2 11^2 13^2 prod_{17<=p<=109} p ~ 1.45e53 (54 digits); R_{4,2}(n_4) < 0.99961, n_4 ~ 2.8e844 (845 digits, largest prime 1913); confirms R_{2,2} < 1 at 65 214 507 758 400 | Proposition 2.9; Section 11 | rigorous interval certificate resting on Theorem 2.6; not claimed minimal; independent exact recheck R_{3,2}(n_3) <= 0.9925002 by the proof auditor |
| I | Every rank: the reduced polynomial P_{l,n}(X)/X has a root in the open right half-plane for n in a set of positive lower density; the negative-second-Hurwitz-determinant set has a density in (0,1) | Section 9.4 | written argument; rank two due to Charlton-Heim-Neuhauser-Stumpenhusen-Troger |

## Evidence

| Check | Script | Output | Arithmetic |
|---|---|---|---|
| Suite 1: brute-force enumeration of commuting l-tuples (n <= 5, l <= 4, 60 values) equals formula (2.2); Table 3 rows | [code/verify.py](code/verify.py) | [results/verification.json](results/verification.json) | exact integers |
| Suite 2: identities (2.1)/(2.4) through n = 40, l <= 6; integrality; cubic lower bound of Section 6.1 | [code/verify.py](code/verify.py) | [results/verification.json](results/verification.json) | exact fractions |
| Suite 3: exact A(l,n,k) vs the main term of Theorem 2.2 for l = 3, 4, k = 2, 3, 4, n <= 1000, and the explicit bound of Theorem 2.6 at all 30 sampled points (largest relative error 1.19% at n = 1000, l = 3, k = 4) | [code/verify.py](code/verify.py) | [results/verification.json](results/verification.json) | exact integers; mpmath 40 digits for the main term |
| Suite 4: every rational constant of Sections 7-8 (2401/1800, 79233/72000, 30/7, ...) | [code/verify.py](code/verify.py) | [results/verification.json](results/verification.json) | exact fractions |
| Suite 5: finite-criterion thresholds log10 n = 9.5427 (2,3), 13.969 (3,3), 18.254 (3,4), 17.35 (4,3) | [code/verify.py](code/verify.py) | [results/verification.json](results/verification.json) | mpmath 40 digits |
| Table 2 (15 rows), rank-two agreement with Starr's closed form (relative difference <= 4.2e-12), monotonicity of K_{p,k}(v) in every computed case | [code/extensions.py](code/extensions.py) | [results/extensions.json](results/extensions.json) | mpmath 40 digits, primes < 1e5, valuations to 60 plus the limit |
| Certificates at ranks 2, 3, 4 (upper bounds 0.99612942, 0.99251008, 0.99960597); greedy rank-two stop at log10 n = 14.46 | [code/extensions.py](code/extensions.py) | [results/extensions.json](results/extensions.json) | exact rational local factors, 60-digit interval arithmetic |

Numbers quoted on the release page and in its audio: 65 214 507 758 400 (Charlton-Heim-Stumpenhusen); a 54-digit n,
about 1.45 x 10^53 (rank three); an 845-digit n (rank four); about 2.52 against 1.10 (l = 3, k = 3); density bounds
2.47 x 10^-16 and 6.8 x 10^-4 (Charlton-Heim-Stumpenhusen, rank two). All are proved in the paper, reproduced by the
scripts or taken from the cited sources.

## Limits

Fixed numbers of orbits only; the range between about log n / log log n and the typical orbit counts is open
(Section 12, open problem 1). The explicit failures are not claimed to be the smallest; higher-rank failures are
astronomically rare (heuristic estimate in Section 11). Monotonicity of K_{p,k}(v) is observed, not proved. Whether the
limit laws are singular is undecided. Table 2 values are numerical evaluations; the certificates are rigorous.

## Replay and records

Replay commands and the comparison policy (elapsed_seconds ignored; integers, rationals, text and verdicts identical;
mpmath decimals to 10 significant digits): [README.md](README.md), [ENVIRONMENT.txt](ENVIRONMENT.txt),
[REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI [.github/workflows/replay.yml](.github/workflows/replay.yml) with
[.github/ci/compare_json.py](.github/ci/compare_json.py). Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Prior art and
limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json). Bibliography: [refs.bib](refs.bib). Licences:
[LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff). Hashes: [MANIFEST.sha256](MANIFEST.sha256).
