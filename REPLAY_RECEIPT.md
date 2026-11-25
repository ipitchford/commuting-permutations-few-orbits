# Replay receipt (v0.1.0-candidate, 25 November 2026)

Release replay by the publisher from fresh copies of the approved candidate (code and requirements only; results
regenerated), working directory = package root, macOS arm64. Interpreter as printed inside the runs: CPython 3.12.11,
mpmath 1.3.0 (see [ENVIRONMENT.txt](ENVIRONMENT.txt)).

| Run | Command | Exit | Output vs recorded |
|---|---|---|---|
| verify -O | `uv run --python 3.12 --with mpmath==1.3.0 python3 -I -B -O code/verify.py --output results/verification.json` | 0 | identical apart from elapsed_seconds (290 values); 29 PASS lines, "ALL CHECKS PASSED" |
| extensions -O | `uv run --python 3.12 --with mpmath==1.3.0 python3 -I -B -O code/extensions.py --output results/extensions.json` | 0 | identical apart from elapsed_seconds (463 values); 57 PASS lines, "ALL CHECKS PASSED" |
| verify, extensions without -O | as above without -O | 0, 0 | identical apart from elapsed_seconds |
| Negative control 1 | verify.py with the constant of (2.9) altered (`F(2401, 1800)` replaced by `F(2400, 1800)` in the check on line 284), under -O | 1 | `FAIL: (2.9) k=3 prefactor: G_{2,3} P3^2/(Q2 Q4) = 2401/1800` |
| Negative control 2 | extensions.py with the certificate threshold tightened (`hi(R_upper) < 1` replaced by `< 0.99`), under -O | 1 | `FAIL: certified R_{2,2}(n) < 1`, and the same for ranks three and four |

stderr was empty for every passing run. The lead's earlier clean replay (24-25 November 2026, also CPython 3.12.11)
produced the same printed output. After the replay only the title-page DOI line of manuscript.tex, the rebuilt PDF and
the release records changed; code/ and results/ are byte-identical to the replayed candidate.

Comparison policy for other platforms (declared before tagging; also in CI, README.md and ENVIRONMENT.txt):
elapsed_seconds ignored; integers, booleans, exact rationals, text and certificate verdicts identical; mpmath decimal
values agree to 10 significant digits.
