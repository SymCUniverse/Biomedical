# P0 component-retention decision after D4 and P3-v2

**Date:** 2026-09-06  
**Status:** **CLOSED — INTERNAL TOOL ARCHITECTURE NARROWED BEFORE EXTERNAL VALIDATION**

The retention question is not whether each component is mathematically interesting. It is whether the component earned a role through prediction, prospective risk ranking, robustness, or scientifically useful traceability.

## Retain

### 1. ALL_METHYLATION_RIDGE predictive core — RETAIN

Reason: improved over COVARIATE_ONLY in 18/19 untouched D4 cancers and 17/18 primary FINAL_HOLDOUT cancers. This is the strongest internally replicated tool-level capability.

Role: cross-layer prediction benchmark/core, not causal mechanism.

### 2. Continuous global cross-layer geometry confidence — RETAIN AS CANDIDATE CONFIDENCE LAYER

Reason: the D4-informed, pre-FINAL frozen geometry-confidence coordinate prospectively correlated with lower FINAL_HOLDOUT prediction risk (rho=-0.5604, one-sided permutation p=0.00795) and favorable AURC (p=0.00526).

Role: continuous confidence/risk-ranking signal. No categorical threshold is licensed from the current data.

### 3. Technical-mask track — RETAIN AS ROBUSTNESS / TRACEABILITY GUARDRAIL

Reason: it detects technical-track instability and prevented overpromotion in prior stages. Current evidence does not establish unique incremental predictive gain from the mask itself, so it is retained for robustness and provenance rather than marketed as a predictive feature.

### 4. Composition attack — RETAIN AS ROBUSTNESS / INTERPRETATION GUARDRAIL

Reason: COAD and other cases show that composition adjustment can materially change interpretation. Post-FINAL exploratory analysis suggests adjusted global geometry may even be useful for future confidence modeling, but that is a new hypothesis for external data.

### 5. Modal / within-layer organization — RETAIN FOR DIAGNOSIS AND TRACEABILITY

Reason: it distinguishes internal-layer organization from cross-layer transport and helped characterize PCPG's failure mode. It is not yet independently validated as a selector.

### 6. Hallmark-label null — RETAIN AS CLAIM FIREWALL

Reason: the label-null repeatedly blocks unsupported semantic identity claims. Its predictive-risk value is weak, but its falsification value is strong. This component earns its place as an epistemic guardrail rather than a prediction feature.

## Demote or retire

### 7. Categorical ACCEPT / CAUTION / REFUSE as the primary predictive selector — DEMOTE / RETIRE FOR CURRENT TOOL VERSION

Reason: all 19 discovery cancers collapsed to CAUTION, making P3-v1 non-identifiable. The categorical selector did not earn prospective risk-stratification utility.

The labels may remain descriptive reason codes, but they should not be advertised as a validated selective-prediction system.

### 8. Hallmark-label-specific semantic sharing as a promoted biological claim — RETIRE

Reason: failed broadly in discovery/replication and does not show useful predictive-risk behavior in FINAL_HOLDOUT diagnostics. Same-Hallmark/local relationships may still be useful per cancer or target, especially PCPG, but no general label-specific sharing claim survives.

### 9. Regulatory autonomy scalar — REMAINS DROPPED

Reason: F4 already classified it `AUTONOMY_REDUNDANT_NARROW`; no later prospective evidence restores it.

### 10. Biological chi — NOT ADMITTED

No same-coordinate dynamical derivation exists in this static TCGA program. No scalar is renamed or promoted to chi.

## New leads, not retained claims

1. **Composition-adjusted global Delta_CKA** may be a stronger confidence coordinate than raw geometry. This was learned after FINAL_HOLDOUT and must be tested externally/future, not promoted here.
2. **PCPG local-vs-global routing:** repeated superiority of SAME_HALLMARK_ONLY over the global all-methylation model suggests a heterogeneous transport regime. A subtype/genotype/local-route model is a candidate external protocol, not a current rescue.
3. **Continuous rather than categorical uncertainty:** the evidence favors preserving graded geometry information instead of collapsing it prematurely into ACCEPT/CAUTION/REFUSE.

## Frozen next architecture for external validation

The next external tool version should be built around:

- unchanged or independently retrained all-methylation cross-layer prediction;
- continuous geometry-based confidence;
- technical and composition attacks as guardrails;
- semantic label null as an interpretation firewall;
- modal/conglomeration outputs for traceability;
- explicit failure routing for cancers/domains such as PCPG;
- no biological chi and no causal/clinical claim without separate evidence.

External cohort compatibility and feature mapping must be frozen from assay compatibility before external outcomes are inspected.

## Immediate program consequence

Internal P0 predictive development is now closed. The next two tracks can proceed in parallel without changing one another:

1. execute the already-frozen Stage C1 biological architecture case study exactly as preregistered;
2. prepare an independent external validation protocol for the narrowed predictive/confidence architecture.
