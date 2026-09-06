# Stage C1 frozen-v2 execution implementation concretization

**Date:** 2026-09-06  
**Status:** **LOCKED BEFORE ANY STAGE C1 BETA-VALUE BIOLOGICAL RESULT**

This note records implementation details needed to execute the already-approved 2026-08-30 Stage C1 v2 scientific contract. It does not change the scientific question, thresholds, features, null logic, multiplicity families, promotion ladder, or claim ceiling.

## Scientific firewall

At the time this implementation was frozen:

- internal P0 predictive development was already closed;
- no Stage C1 biological beta-value result had been calculated or inspected;
- P0/D4/P3 results were not used to change C1 science;
- historical `CV/2` is not used as biological chi;
- no damping, exceptional-point, temporal-inheritance, causal, clinical, optimum, or `chi=1` claim is admitted.

Canonical C1 scientific contract: `STAGE_C1_ANALYSIS_CONTRACT_DRAFT_V2_20260830.md`, frozen by the preregistration as source blob SHA-1 `a39414c5990c234dc513569bd405b0237117d434`.

## Frozen input identities

- methylation source SHA-256: `5934c497882fbe8178d128a3a7f71e765480af6bbd460e0398de3428cd075b77`
- Stage A Hallmark-union RNA cache SHA-256: `e65f6788aa6037fef407169794f29d63322de2769343bb6e594fe469dfeb8e63`
- Hallmark v2026.1.Hs membership SHA-256: `bc6a9a33d7421dc407d33a66859760ba25e47b2f398e1a43c9156f80c71b3900`
- Patient_Split SHA-256: `12b85d67d06e57c6d1be914444c65aa526f2b821e117725dd036142cc6b0a825`
- Sample_Eligibility SHA-256: `4a906d2d80ab36af0cc1dc8f9d8dcb834cfc717021a355c78c089f5cc665e3c2`
- ABSOLUTE purity SHA-256: `f430a975433d82e0098d7405619d4f12a0c765fcd97e7d63cc9b1de7f2d763cd`
- leukocyte-fraction SHA-256: `5a8268caedbf8dc98a75be0528d583238d7355761d9fc746e42002f223a982d9`
- exact C1A annotation export SHA-256: `a7f83233f97c3933752d74b8042e967de88df20eba2cf477a536136631a8da17`
- exact C1A Chen export SHA-256: `078e95716af2b20c3515f59d09310d20c43deb5ed8ba8d1b70885810acde2179`
- C1A source summary SHA-256: `b19297c4855ec12bc1fe4dd473742a72bef23db0b1b446258aa6a51b3679205b`

The exact recovered C1A annotation package is Bioconductor 3.8 `IlluminaHumanMethylation450kanno.ilmn12.hg19` **v0.6.0**.

## Mechanical/statistical concretizations

These details were not chosen from C1 outcomes.

1. The identical patient draw across primary/masked/RNA views uses seed tokens `track=SHARED`, `stratum=ALL`, `null_type=resample`, `replicate=0` inside the frozen SHA-256 seed token order.
2. The singular construction, patient, and label permutation specified for each observed resample is implemented as one deterministic null replicate per resample. The 100 resamples provide the stability distribution; resamples are not treated as independent inferential units.
3. Exact-zero ties in the pan-cancer sign test are excluded from Binomial `n`; no continuity correction or pseudocount is used.
4. The five secondary regulatory strata form one ten-test BH family per probe track: five strata x `{patient-null, label-null}`. These tests are contextual and cannot rescue `PROMOTER_CORE`.
5. Adjusted patient-alignment nulls permute the already residualized methylation representation relative to residualized RNA, preserving the frozen B1 composition projection while destroying cross-assay patient identity.
6. The secondary topology label null reuses the same deterministic Hallmark-label permutation used for `H3b` in that cancer/track/resample.
7. Sample-space eigenvector signs are stored reproducibly by orienting each retained top-five eigenvector so its largest-absolute sample loading is positive. This changes no eigenspace, eigenvalue, contribution, CKA, or principal angle.
8. Probe contribution formulas are evaluated in float64, checked to sum to one before storage, then stored as float32 to keep the required full top-five vectors tractable. Matching ordered probe-ID sidecars and SHA-256 hashes are retained.
9. The exact Stage A RNA Hallmark rule retains its minimum 15 eligible mapped genes; finite-coverage/z-score/standardized-zero-imputation behavior remains the Stage A rule.
10. Reconstructed C1A gzip bytes are allowed to differ from the historical compressed files only when exact source exports/historical code reproduce the frozen semantic counts. Historical byte identity is never falsely claimed.

## Required output completeness

The execution retains/reports:

- every sample eligibility record;
- every cancer x probe eligibility decision, finite count, and imputation median;
- imputation burden by cancer and probe track;
- all 29 observed and construction-null normalized spectral positions;
- top-five observed sample-space eigenvectors with resample membership mapping;
- full top-five probe contribution vectors locally plus ordered probe sidecars and hashes;
- top-five methylation/RNA principal angles;
- raw and composition-adjusted signed Hallmark correlations;
- primary and secondary null effects;
- cancer-level medians and exact sign tests;
- separate raw and adjusted BH families;
- explicit C1-1/C1-2/C1-3 promotion decision and frozen claim ceiling.

## Runtime/recovery controls

- 5 GB methylation projection is disk-backed and resumable.
- RNA expression member is extracted once and memory-mapped.
- probe rules are written atomically.
- each cancer checkpoint is hash-sidecar verified and keyed to frozen contract/input identities.
- every biological resample checkpoints after both tracks finish.
- worker crashes are automatically restarted up to three times from verified checkpoints.
- large immutable input SHA verification is cached only when path, size, and modification timestamp are unchanged.
- result files are atomically written and SHA-manifested.

## Frozen runtime source hashes

- `src/run_stage_c1.py`: `9afe3be8dfd0236d002afc8c45f4605bab227f7f031b36e3b03c8c40a5687408`
- `src/stage_c1_core.py`: `06911532db113198167e4aa32c35632c9cc4e0debe2bb8b30e3c1235622fe1f2`
- `src/c1a_historical_rebuild.py`: `ccb2a1108643f46cba81ef8e0ea00182f691d6ddb13c13fa764e384220660a8f`
- `src/windows_launcher.py`: `f27963fa33997efcd78c55a6a203cd8dddfa8be3354347905772a02270d253a9`

Regression/contract suite before packaging: **17/17 PASS**, including a one-resample synthetic end-to-end checkpoint/output-contract test.
