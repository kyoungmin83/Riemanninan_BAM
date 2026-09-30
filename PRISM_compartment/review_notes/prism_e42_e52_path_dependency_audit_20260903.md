# PRISM e42/e52 frozen path-dependency audit

Date: 2026-09-03

## Purpose

This validation-only audit determines which fitted pathway actually carries
module recovery and donor ADNC ordering before changing the target of automatic
pathology-rank learning. It compares the frozen rank-8 epoch-42 checkpoint with
the frozen compartmental/path-rank epoch-52 checkpoint on the same sampled
validation nuclei.

The official nine-donor test split is not opened. No model parameter is updated.

## Paired inference variants

For every sampled nucleus, the model is run once. Each OFF result is obtained by
subtracting one explicitly reconstructed score contribution before converting
ordinal scores to probabilities.

- full fitted model;
- z-derived named pathology decoder OFF;
- common global pathology term OFF;
- zero-centred common cell-type pathology delta OFF;
- whole common pathology term OFF;
- observed-pathology interaction OFF;
- common plus interaction OFF;
- all observed-pathology-driven terms OFF, including personal response;
- all observed-pathology-driven terms plus the z-derived named route OFF;
- rank-2 personal baseline OFF;
- personal pathology response OFF;
- module-local term OFF;
- optional module-local nonlinear subpath OFF;
- all personal terms OFF.

The common global and cell-type-delta pieces are recomputed from the checkpoint
parameters and checked numerically to sum back to the stored common coefficient
on every batch. The z-derived named route is removed after its learned state
mixer gate, which makes its score subtraction exact.

## Readouts

- validation ordinal NLL;
- donor-centred general-module recovery;
- same-donor and donor-disjoint disease-module direction recovery;
- five-axis disease-module recovery;
- module-to-ADNC leave-one-donor-out Spearman correlation;
- donor variance, AUROC, and average precision of the z-derived named pathology
  coordinates.

The historical module-to-ADNC implementation is reproduced for parity. A second,
preferred value fits feature centring and scaling inside each leave-one-donor-out
training fold.

## Interpretation boundary

This is a conditional-reliance audit of a frozen fitted model. A deterioration
after removing a path shows that the current checkpoint depends on that path. It
does not prove how a model retrained without that path would behave and does not
establish causal biological importance.

## Pre-registered design decision

- Rank the common cell-type delta only if removing that delta materially and
  reproducibly harms held-out module/disease recovery beyond removing the global
  common term alone.
- If the global common term carries most of the signal, factorising only the
  cell-type delta is the wrong target.
- Do not retain or rank the z-derived named route merely because it exists. It
  must show donor-varying pathology readout and useful OFF-minus-full behavior.
- Treat module-to-ADNC as a safety/dependency readout, not a clinical-performance
  headline, especially if it collapses when observed-pathology-driven terms are
  removed.

## Artifacts

- Runner: `scripts/analysis/run_sv6_e42_e52_path_dependency_audit_20260903.sh`
- Extractor: `scripts/analysis/extract_prism_path_ablation_fingerprints.py`
- Summarizer: `scripts/analysis/summarize_prism_path_ablation.py`
- Remote output: `/gstorage_data/kmlee/project/riemann_bam/analysis_outputs/prism_common_rank_dependency_audit_20260903`
- Local copy: `analysis_outputs/prism_common_rank_dependency_audit_20260903`

## Results

The paired audit completed on 14,986 validation nuclei from 16 donors and 24
cell types, with at most 40 nuclei per donor and cell type. There were 404
evaluable modules. The summarizer used 100 donor bootstrap resamples and 100
donor-disjoint splits. These resampling counts are sufficient for the design
decision below, but publication-grade interval estimation should use more
resamples.

### Rank-8 epoch 42

| Removed path | OFF minus full NLL | General centred | AD-module direction | Donor-disjoint AD | Module-to-ADNC |
|---|---:|---:|---:|---:|---:|
| None (full) | 0.000000 | 0.630 | 0.849 | 0.538 | 0.793 |
| Named z-derived route | 0.000000 | 0.630 | 0.849 | 0.538 | 0.793 |
| Common global | +0.006298 | 0.624 | 0.597 | 0.297 | 0.714 |
| Common cell-type delta | +0.004540 | 0.581 | 0.853 | 0.499 | 0.743 |
| Whole common pathology | +0.010907 | 0.599 | 0.551 | 0.242 | 0.673 |
| All observed-pathology terms | +0.012091 | 0.536 | 0.643 | 0.247 | -0.194 |
| Personal pathology response | +0.002055 | 0.543 | 0.822 | 0.550 | 0.816 |

The named z-derived route is identically inactive in this checkpoint. Both the
global and cell-type-delta pieces of the common label-driven pathway help NLL.
The global piece carries most of the same-donor AD direction and much of the
donor-disjoint AD signal, while the cell-type delta is most visible in general
centred recovery.

### Compartmental epoch 52

| Removed path | OFF minus full NLL | General centred | AD-module direction | Donor-disjoint AD | Module-to-ADNC |
|---|---:|---:|---:|---:|---:|
| None (full) | 0.000000 | 0.841 | 0.854 | 0.402 | 0.289 |
| Named z-derived route | +0.000236 | 0.841 | 0.850 | 0.399 | 0.289 |
| Common global | -0.009645 | 0.791 | 0.803 | 0.336 | 0.016 |
| Common cell-type delta | +0.001360 | 0.844 | 0.844 | 0.343 | 0.374 |
| Whole common pathology | -0.007576 | 0.757 | 0.795 | 0.267 | -0.022 |
| All observed-pathology terms | -0.007513 | 0.710 | 0.782 | 0.274 | 0.057 |
| Personal pathology response | +0.000224 | 0.815 | 0.833 | 0.415 | 0.340 |
| Module-local | +0.000777 | 0.732 | 0.845 | 0.481 | 0.701 |

The named z-derived route has a substantial score RMS of 0.257, but its removal
has negligible effects on every primary readout. Its NLL interval includes zero
(-0.000199 to +0.000766), and its recovery changes are approximately zero. It is
therefore large but functionally redundant in the fitted model.

Removing the common global term improves NLL but reduces biological recovery
and collapses module-to-ADNC from 0.289 to 0.016. Removing the cell-type delta
barely changes same-donor recovery, but donor-disjoint AD recovery falls from
0.402 to 0.343. Thus global and cell-type-specific common effects contribute in
different ways, and NLL alone would select the wrong biological decomposition.

The module-local path shows a second trade-off. Removing it lowers general
centred recovery from 0.841 to 0.732 and five-axis recovery from 0.784 to 0.611,
but raises donor-disjoint AD recovery from 0.402 to 0.481 and module-to-ADNC from
0.289 to 0.701. This is evidence of objective interference, not proof that the
path should simply be deleted; a retrained ON/OFF arm is required.

### Named-coordinate audit

Only 2.1% to 8.4% of the e52 named-coordinate variance is between donors. Its
donor-level AUROCs are 0.476 for Braak, 0.587 for Thal, 0.650 for LATE, and
0.714 for Lewy. With only 16 validation donors, this is at most weak evidence of
donor-varying pathology readout and does not justify making this redundant path
the target of automatic rank learning.

### Module-to-ADNC interpretation

Turning off all terms directly driven by observed pathology reduces
module-to-ADNC from 0.793 to -0.194 at e42 and from 0.289 to 0.057 at e52.
Module-to-ADNC is therefore mainly a check that the reconstructed modules retain
the supplied pathology ordering. It is useful as a dependency/safety readout,
but it is not independent clinical prediction evidence.

## Design decision

1. Do not learn pathology rank on the z-derived named route. It is inactive at
   e42 and large but redundant at e52.
2. Do not factorise only the common cell-type delta. The global common component
   is the dominant ADNC and same-donor AD carrier, whereas the delta contributes
   to general and donor-disjoint recovery.
3. If automatic pathology rank is retained, apply it to the observed-label
   common pathway. Keep separate rank gates for the global matrix and the
   centred cell-type-delta matrix so that the global signal cannot hide or prune
   the transferable cell-type-specific signal.
4. Fade the z-derived named route out during Phase II, as in the rank-8 lineage,
   and keep its pathology head only as a detached monitoring probe unless a
   retrained ablation demonstrates unique utility.
5. Test module-local OFF versus a constrained module-local redesign in a
   retrained proof-of-life arm. Select with a multi-objective gate containing
   NLL, centred module recovery, donor-disjoint disease recovery, module-to-ADNC,
   and leakage—not NLL alone.

These decisions are based on frozen-checkpoint conditional removals. They define
the next proof-of-life experiment; they do not by themselves authorize or
replace a retrained ablation.
