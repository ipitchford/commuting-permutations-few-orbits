# Review summary (internal, AI-assisted)

This is a public summary of the internal checking of version 0.1.0-candidate, without scores. It is not peer review,
independent reproduction or formal verification. No human referee has assessed the paper.

## What was done (24-25 November 2026)

1. **Panel of five Claude reviewers** (editor, methodology and proof checking, domain and literature, perspective, and
   devil's advocate) on the first draft, with the cited sources read. The methodology reviewer recomputed all constants
   and found no mathematical error; no reviewer raised a critical issue. The major items concerned positioning against
   the literature, a buried headline, readability, and the strength of the explicit bounds.
2. **One revision round with two extensions.** (a) A new theorem giving the exact limit inferior and superior of the
   log-concavity ratio as Euler products (Theorem 2.5, Table 2), which replaces the weak uniform bound 1.10 by about
   2.52 at three shuffles and three orbits and reduces to Starr's formula at rank two. (b) Certified explicit failures
   in ranks three and four (Proposition 2.9), with rank-two consistency at the first known counterexample. Also: a
   context-first introduction and roadmap; credit stating that the leading asymptotic follows the classical
   Ingham–Halberstam route and that the density result applies the Erdős–Wintner theorem; six added references; a full
   proof of Lemma 3.1; thresholds of the finite criterion stated; new scripts without assert, with brute-force
   enumeration and an exact-versus-asymptotic comparison.
3. **Re-review** of the revision: every earlier major item addressed or declined with a statement in the paper; no
   mathematical error found in the new material; editorial items applied.
4. **Proof audit** of every new or rewritten proof (Lemma 3.1, Theorem 2.5, Section 11 and Proposition 2.9, the
   Hurwitz density paragraph, Corollary 2.7): all verified, with an independent exact-rational recheck of the rank-three
   certificate (R_{3,2}(n_3) ≤ 0.9925002).
5. **Citation of Gadiyar, Murty and Padma (2014)** reinstated after the paper was obtained and read page by page; their
   Theorem 10 is stated as in the paper and Theorem 13(a) identified as the case ℓ = 2, k = 2.
6. **Clean replays** of both scripts from fresh copies with negative controls, by the lead and again by the publisher
   (see [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md)).

No Astra/Codex second-model review took place in this batch.

## Main objections and how they were handled

| Objection | Disposition |
|---|---|
| "The 1.10 lower bound at k = 3 is weak; the exact limit inferior is computable." | New Theorem 2.5: exact limit inferior and superior as Euler products, evaluated in Table 2. |
| "The scale of failures is hidden; are there explicit higher-rank counterexamples?" | Proposition 2.9: rigorous interval certificates at a 54-digit n (rank three) and an 845-digit n (rank four); the heuristic estimate of the first failure is labelled as heuristic. |
| "The leading asymptotic for ℓ ≥ 3 is routine Ingham–Halberstam." | Agreed and stated in Section 1.5, with credit to Lemke Oliver–Shrestha–Thorne and Starr's own remark; the paper's analytic contribution is the explicit uniform error and the local-factor inequalities. |
| "The density theorem is a routine application of Erdős–Wintner." | Framed as an application, with no new probabilistic input claimed. |
| "The finite criterion is far from sharp." | Thresholds stated (about 10^9.5 to 10^18.3); a sharper relative form is open problem 5. |
| "Normalisation differs from Starr's." | Remark 2.1 gives the exact factor (k+1)/k. |
| Not addressed | The range between about log n / log log n and the typical orbit counts, monotonicity of the local factors, the first rank-three failure and the numerical density, and singularity of the limit law are listed as open problems. |

## Scope of assurance

Internal AI-assisted review only: written proofs checked by AI reviewers, computer checks replayed, sources read. No
human specialist review, no independent reproduction by an unaffiliated party, no proof-assistant formalisation.
Literature searches were of limited scope (see [RESEARCH_GATES.json](RESEARCH_GATES.json)).
