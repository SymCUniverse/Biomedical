# P0 P3-v2 FINAL_HOLDOUT independent audit

**Date:** 2026-09-06  
**Status:** **CLOSED — FINAL_HOLDOUT EVALUATION COMPLETE**

## 1. Provenance and firewall

This evaluation uses the D4-informed P3-v2 amendment that was frozen before FINAL_HOLDOUT molecular values were opened. It does not rewrite the original categorical P3-v1 result; P3-v1 remains `NOT_EVALUABLE_CATEGORICAL_SELECTIVE_TEST`.

Uploaded checkpoint integrity:

- `Final_Methylation_Scores.pkl` SHA-256: `6cd6fba96656edf3c176d7fb203d403f59c88d00009fda2c41e44d5e01fee33e`
- methylation checkpoint key: `ab99a87da30476caba4c8eebc1f9b64c50bf43c6d2bd1325484cee02443388d9`
- `Final_RNA_Scores.pkl` SHA-256: `78c918c60745d40fdcd61462efdd8379e3414df2c13615c644825026ecf4b189`
- RNA checkpoint key: `6ec465353feb03e7d2ad18ce236491989fc8168cf1dff45901915f66b4496fdd`

Both checkpoint keys exactly match those recomputed from the frozen P3-v2 implementation/config identities and the previously audited source identities. Structural validators pass:

- methylation score rows: **151,740**
- RNA score rows: **84,300**
- FINAL_HOLDOUT eligible participants: **1,686**
- composition-complete participants: **1,534**
- primary evaluation cancers: **18**
- PAAD remains pre-value excluded from primary evaluation because its frozen composition-complete count is 26, below the inherited n=30 rule.

No model, transform, composition residualizer, participant partition, covariate value, Hallmark definition, or prediction threshold was refit from FINAL_HOLDOUT.

## 2. Mechanical evaluation repair

The checkpoint-generating run reached the frozen projected-score checkpoints but did not produce its return ZIP. Re-execution exposed two implementation issues after the scientific quantities had already been frozen:

1. pandas Copy-on-Write produced a read-only array inside the rank-correlation helper;
2. the literal 99,999-permutation implementation was needlessly slow despite an algebraically equivalent rank-vector calculation.

The evaluator was repaired mechanically by materializing writable ranked vectors and evaluating the same deterministic permutation stream with equivalent dot-product formulas. A 1,000-permutation equivalence test against the original algorithm with only the Copy-on-Write compatibility fix produced identical Spearman coefficients and p-values for all three association tests. AURC differed only by `5.55e-17` floating-point roundoff and had the identical permutation p-value. No scientific rule or RNG stream changed.

## 3. Final P1 direct prediction

Primary comparison: cancer-median normalized MSE for unchanged `ALL_METHYLATION_RIDGE` versus unchanged `COVARIATE_ONLY`.

- all-methylation better: **17 / 18 cancers**
- exact one-sided sign test: **p = 7.2479248046875e-05**
- sole loss: **PCPG**

Median across the 18 primary cancers:

| Model | Median normalized MSE | Median held-out R2 |
|---|---:|---:|
| ALL_METHYLATION_RIDGE | **0.4440** | **0.5074** |
| COVARIATE_ONLY | 0.7018 | 0.1878 |
| SAME_HALLMARK_ONLY | 0.6357 | 0.2493 |
| MEAN_ONLY | 0.8871 | 0.0000 |

The median cancer-level improvement `COVARIATE_ONLY - ALL_METHYLATION_RIDGE` is approximately **0.2842 normalized-MSE units**.

This independently repeats the central D4 internal-replication result in a second untouched TCGA partition. It remains internal TCGA evidence, and the frozen 24-cancer promotion floor is not met.

## 4. P3-v2 continuous confidence test

Frozen D4-informed selector:

`GLOBAL_GEOMETRY_CONFIDENCE = min(raw PRIMARY_PUBLICATION Delta_CKA, raw MASKED_TECHNICAL Delta_CKA)`

Primary FINAL_HOLDOUT association:

- Spearman rho with final prediction risk: **-0.56037**
- one-sided permutation p: **0.00795**
- permutations: **99,999**

Thus stronger discovery global cross-layer geometry prospectively ranks lower FINAL_HOLDOUT prediction risk under the D4-informed pre-FINAL formulation.

Development-set stability comparator:

- D4 replication risk versus FINAL_HOLDOUT risk: rho **+0.78328**, p **0.00011**

Simple sample-size comparator in the prespecified favorable direction:

- DISCOVERY_N versus FINAL_HOLDOUT risk: rho **+0.39112**
- one-sided negative-direction p **0.94619**

The confidence result therefore is not explained by the trivial proposition that larger discovery cancers have lower risk.

## 5. Coverage-risk

Frozen geometry-confidence AURC:

- AURC: **0.4026622**
- one-sided permutation p: **0.00526**
- permutations: **99,999**

At 50% coverage, selecting the 9/18 cancers with highest frozen geometry confidence gives mean normalized-MSE risk **0.41372**, versus **0.44061** over all 18 cancers, about a **6.1% reduction**.

This is a modest but prospectively supported ranking effect, not a calibrated clinical acceptance threshold.

## 6. PCPG stress case

PCPG again is the only primary cancer in which the global all-methylation model loses to the covariate-only baseline:

- ALL_METHYLATION_RIDGE nMSE: **0.41656**
- COVARIATE_ONLY nMSE: **0.39993**
- difference `covariate - all`: **-0.01663**
- SAME_HALLMARK_ONLY nMSE: **0.36438**

The global loss is smaller than in D4 but recurs in the same direction. At target level, ALL_METHYLATION_RIDGE beats COVARIATE_ONLY on **17/50** targets and loses on **33/50**.

The recurrent pattern supports treating PCPG as a transport/representation stress case. It does not justify exclusion or post-result model retuning. The fact that SAME_HALLMARK_ONLY again outperforms the global model is a concrete lead for a future routing/subtype-conditioned investigation on new data.

## 7. Post-FINAL component diagnostics

These are explicitly exploratory and cannot replace the frozen P3-v2 result.

Discovery-to-final Spearman correlations show that the predictive-confidence signal is concentrated in global CKA geometry:

- adjusted masked Delta_CKA: about **-0.606**
- adjusted primary Delta_CKA: about **-0.602**
- raw masked Delta_CKA: about **-0.577**
- raw primary Delta_CKA: **-0.560** (the prospectively frozen P3-v2 axis)

By contrast, Hallmark-label-specific Delta_A_label is weak or wrong-signed in raw analyses and near zero after adjustment. Patient-alignment and same-Hallmark strength retain moderate association but less than the global geometry signal.

This suggests a useful separation of roles: global geometry is a candidate predictive-confidence coordinate; semantic label-null machinery is better retained as an interpretation firewall than as a predictive selector. Adjusted global geometry is a promising **future external-validation hypothesis**, not a result promoted from this FINAL_HOLDOUT.

## 8. Claim ceiling

What is established internally:

1. discovery-trained methylation Hallmark architecture adds substantial predictive information for RNA Hallmark states beyond purity/leukocyte covariates in two untouched TCGA partitions;
2. a D4-informed continuous discovery global-geometry confidence coordinate prospectively stratifies FINAL_HOLDOUT prediction risk and coverage-risk;
3. the original categorical ACCEPT/CAUTION/REFUSE selector is not validated by P3-v1 and remains non-evaluable;
4. Hallmark-label-specific semantic sharing remains unsupported as a general promoted branch;
5. PCPG is a reproducible stress case for the global transport model.

Not established: external generalization; clinical utility; subtype-specific mechanism; causal methylation-to-RNA regulation; temporal inheritance; treatment response; damping, exceptional-point behavior, or biological chi.

## 9. Artifact identity

`Final_Holdout_Result.zip` SHA-256:

`b5e54917f0030ba6fc0841d79bcd38671922c0c0d37e620005d667f50c7e305f`

**Gate decision: P3-v2 CLOSED. Proceed to component retention and the already-frozen Stage C1 biological execution without changing Stage C1 science.**
