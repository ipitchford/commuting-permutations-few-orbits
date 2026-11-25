# Commuting permutations with few joint orbits: asymptotics, log-concavity and the density of failures in every rank

Version 0.1.0-candidate, 25 November 2026. DOI [10.5281/zenodo.24541319](https://doi.org/10.5281/zenodo.24541319)
(all versions: [10.5281/zenodo.24541318](https://doi.org/10.5281/zenodo.24541318)). Author: Anonymous.
Release page: https://evidencepress.org/releases/commuting-permutations-few-orbits/

**Status: unrefereed research candidate.** The results have written proofs that were checked by internal AI-assisted
review, with computer checks. No human referee has assessed them.

## What this is about

Take a few shuffles (permutations) of n cards that all commute with one another: doing them in any order gives the same
result. Together they split the cards into groups, called *joint orbits*, that the shuffles never mix. Count how many
families of ℓ such shuffles split the cards into exactly k groups; call the count A(ℓ, n, k). For a single shuffle
(ℓ = 1) the groups are its cycles and the counts are the classical Stirling numbers. The same counts describe coverings
of an ℓ-dimensional torus and, for ℓ = 2, the coefficients of the d'Arcais polynomials behind powers of the partition
generating function.

Stirling numbers rise and then fall smoothly: they are *log-concave*, meaning each one squared is at least the product
of its two neighbours. Abdesselam conjectured the same for every ℓ and n. In 2026 Starr showed that for ℓ = 2 it fails
at k = 2 (two groups) for some huge n, though it holds at every fixed larger k once n is large, and he named the
higher-rank case the most important open problem. Charlton, Heim and Stumpenhusen found the first rank-two failure,
n = 65 214 507 758 400, and bounds on how often failures occur, but did not show that their frequency is well defined.

## What the paper shows

For every fixed number of shuffles ℓ ≥ 2 and every fixed number of groups k:

1. **The conjecture fails in every rank, and in the fixed-k regime only at two groups.** The log-concavity ratio at
   k = 2 tends to zero along n = m! and is below any ε > 0 on a set of n of positive density. At every fixed k ≥ 3
   the conjecture holds once n is large (Theorems 2.3 and 2.4).
2. **Exact limiting worst and best values.** At fixed k ≥ 3 the limit inferior and limit superior of the ratio are
   products over primes (Theorem 2.5). For example, for three shuffles and three groups the limit inferior is about
   2.52, where the elementary argument gives only the bound 1.10. At rank two the formula reproduces Starr's.
3. **The frequency of failures exists.** The ratio has a continuous limiting distribution, so the proportion of n at
   which two-group log-concavity fails exists and lies strictly between 0 and 1 in every rank (Theorem 2.8). This is
   new also at rank two.
4. **Explicit certified failures.** A 54-digit n (about 1.45 × 10^53) for three commuting shuffles and an 845-digit n
   for four (Proposition 2.9), certified by rigorous interval arithmetic.
5. **An explicit error bound** for the asymptotic formula for A(ℓ, n, k), uniform in the prime factorisation of n
   (Theorems 2.2 and 2.6), and a growing window 3 ≤ k ≤ c log n / log log n (c < 1/ℓ) of strict log-concavity.

**Why two groups are special.** One group means the shuffles act transitively, and that count carries the factor
σ(n)/n (the sum of the divisors of n, divided by n), which is large when n has many small prime factors, as for
factorials. Counts with two or more groups carry only bounded arithmetic factors. The two-group ratio has the one-group
count in its denominator, so it collapses exactly when σ(n)/n is large: a close cousin of the abundant numbers.

**Limits.** Only a fixed number of groups is covered; the typical group counts were treated by Abdesselam, and the
range in between is open beyond k up to about log n / log log n. The qualitative leading asymptotic for ℓ ≥ 3 follows
the classical Ingham–Halberstam method, as Starr noted at rank two, and the density result applies the classical
Erdős–Wintner theorem. The explicit failures are not claimed to be the smallest, and failures in ranks three and four
are astronomically rare. The values in the table of limits are numerical evaluations; the certificates are rigorous
interval computations.

## Contents

| File | What it is |
|---|---|
| [manuscript.pdf](manuscript.pdf), [manuscript.tex](manuscript.tex) | The paper (20 pages); Section 1 is a plain introduction with a roadmap |
| [code/verify.py](code/verify.py) | Five suites: brute-force enumeration of commuting tuples (n ≤ 5, ℓ ≤ 4) against the coefficient formula; exact Dirichlet-convolution identities (n ≤ 40, ℓ ≤ 6); exact A(ℓ,n,k) against the asymptotic formula and the explicit bound for ℓ = 3, 4, n ≤ 1000; every rational constant of Sections 7-8; thresholds of the finite criterion |
| [code/extensions.py](code/extensions.py) | Table 2 (limit inferior and superior, 15 rows, with the rank-two comparison with Starr's closed form) and the interval certificates of Proposition 2.9 at ranks two, three and four |
| [results/verification.json](results/verification.json), [results/extensions.json](results/extensions.json) | Recorded outputs of the two scripts |
| [AI_INDEX.md](AI_INDEX.md) | Machine-oriented map of claims, evidence and limits |
| [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) | What internal checking was done and how objections were handled |
| [RESEARCH_GATES.json](RESEARCH_GATES.json) | Prior-art and limitation record |
| [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [ENVIRONMENT.txt](ENVIRONMENT.txt) | Replay runs and environment |
| [refs.bib](refs.bib) | Bibliography of the paper |
| [CITATION.cff](CITATION.cff), [LICENSE.md](LICENSE.md), [MANIFEST.sha256](MANIFEST.sha256) | Citation, licences, hashes |

## Replay (about three minutes)

From the package root, with [uv](https://docs.astral.sh/uv/):

```sh
uv run --python 3.12 --with mpmath==1.3.0 python3 -I -B -O code/verify.py --output results/verification.json
uv run --python 3.12 --with mpmath==1.3.0 python3 -I -B -O code/extensions.py --output results/extensions.json
```

Each prints `ALL CHECKS PASSED` and exits 0; any failed check prints `FAIL` and exits 1, also under `python -O`.
Comparison policy (declared before tagging; also in CI, ENVIRONMENT.txt and REPLAY_RECEIPT.md): the key
`elapsed_seconds` (run time) is ignored; everything else must be identical, except that numbers computed with mpmath
(the 40-digit diagnostics of verify.py suites 3 and 5, Table 2 and the printed certificate intervals) may differ across
platforms in the last places and must agree to 10 significant digits. Integers, booleans, exact rationals, the
enumeration, the identities and the certificate verdicts must be identical. GitHub Actions runs both scripts on Linux,
compares under this policy (`.github/ci/compare_json.py`) and runs two negative controls.

## Assurance

Internal AI-assisted review only: a panel of five Claude reviewers, one revision round with two extensions, a
re-review, a proof audit of every new or rewritten proof with an independent exact recheck of the rank-three
certificate, and clean replays with negative controls. No Astra/Codex review took place. This is not peer review,
independent reproduction or formal verification.

## Licence and citation

Text, paper, recorded outputs and records: CC0-1.0. Code: MIT. See [LICENSE.md](LICENSE.md). Cite as: Anonymous,
*Commuting permutations with few joint orbits: asymptotics, log-concavity and the density of failures in every rank*,
Evidence Press, version 0.1.0-candidate (2026), https://doi.org/10.5281/zenodo.24541319.
