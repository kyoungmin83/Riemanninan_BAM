# Integrated PRISM Cell-type/AD-module Curriculum v2

Status: design started; training not launched  
Date: 2026-08-24  
Decision: one integrated architecture, one continuous run, one checkpoint history

Persistent user-locked contract: the next primary run must execute this entire document as one merged design. The canonical CT64 Phase-I rank8 biological warmup and the later cell-type/general-module, cell-type/AD-module, PHU, and automatic-generator curriculum are cumulative requirements, not alternative experiments. This decision is also recorded in the workspace `AGENTS.md`.

## 1. Objective

The next primary run must improve both of the following without removing intentional cell-type content:

1. cell-type-specific recovery of the general donor-varying module program;
2. cell-type-specific recovery of AD/pathology-associated module directions.

The model is compared against the frozen CLS PRISM personal-rank2 epoch-20 baseline with seed 42, the identical donor split, and an identical effective optimization budget. The 9 official test donors remain sealed throughout training and checkpoint selection.

## 2. Evidence driving the redesign

| Evidence | Current integrated result | Design consequence |
|---|---:|---|
| SV6 guarded best | epoch 36 | Use module-aware selection, not loss-only selection |
| e36 AD-module median / pooled | 0.787 / 0.651 | Median is near baseline, but pooled recovery needs direct protection |
| Frozen rank2 e20 median / pooled | 0.803 / 0.703 | Replacement requires recovery of the broad pooled program |
| e36 clean mini-batch z participation | 2.54 / 32 | Scratch representation warmup was insufficient |
| e49 global z participation | 5.59 / 32 | Increasing z rank alone does not guarantee biological alignment |
| Oligodendrocyte e31 / e36 / e49 | 0.432 / 0.399 / 0.295 | Low-amplitude non-neuronal signals are progressively lost |
| Oligodendrocyte train cells per donor, median | about 3,045 | The problem is not cell scarcity |
| Oligodendrocyte robust target scale | 0.331, lowest of 24 cell types | Use larger pseudobulk aggregates and reliability-aware losses |
| Current non-neuronal rescue | legacy 2-cell blocks | Replace with an 8-cell aggregate pseudobulk loss |
| Current rescue state-readout updates | disabled | Keep disabled in the primary CLS model; the existing allowlist is valid only for persistent AGP pooling |
| Current unique-gene protection during pruning | disabled | Protect module/gene coverage and per-cell-type recovery |
| Current integrated source config | Phase-II rank2 PHU-relative config | The run was a Phase-II scratch curriculum, not a Phase-I/II union |
| Current integrated first-order pathology decoder / auxiliary head | disabled / disabled from epoch 1 | Restore the exact CT64 Phase-I rank8 biological warmup inside the one-run model |

The strong frozen rank2 lineage used an initial CT64 regime with learning rate 0.0001 and gradient accumulation 1. The scratch integrated run instead used learning rate 0.00002 and gradient accumulation 2 from the beginning. The current “warmup” therefore had one fifth of the learning rate and about half as many optimizer updates per data pass.

The current SV6 run is not accepted as a faithful integrated equivalent of the frozen lineage. Its generated manifest names the already-final `rank2 PHU-relative` configuration as its source, sets `use_pathology_decoder=false`, disables `pathology_aux`, and starts from scratch. It therefore omitted the CT64 Phase-I rank8 pathology-conditioned representation learning that shaped the frozen rank2 lineage. Its results remain useful as an ablation of that omission, but not as evidence against a correctly integrated Phase-I/II design.

Persistent AGP is not selected for the new primary design. It raised global z participation to about 13 / 32, but its epoch-44 AD-module median and pooled recovery were only 0.632 and 0.522. The evidence favors a well-trained CLS representation over a pooling change.

## 3. Fixed architecture

- CLS pooling, latent width 32, explicit 64-dimensional cell-type embedding.
- Cell type remains biological content. It is not projected out and is not an adversarial removal target.
- Personal rank remains 2.
- The 414-module registry and full 414-generator initialization remain fixed.
- Reproduce the CT64 Phase-I first-order pathology decoder exactly: four named axes, pathology rank 8, no pairwise pathology term, no free pathology interactions, and no disease gate.
- Reproduce the CT64 Phase-I pathology auxiliary head: Braak, Thal, LATE, and Lewy supervision with its canonical weights and pathology-direction tie.
- First-order pathology, common, personal, response, PHU, module-rescue, and generator-gate components are all instantiated before epoch 1. Curriculum gates change their contribution, not the architecture.
- Pairwise and higher-order pathology interaction terms remain disabled in the primary run unless a separate predeclared ablation is authorized.
- No warm-started model weights are used. The run starts from scratch but reproduces the useful optimization regimes of the old lineage inside one checkpoint history.
- Use disjoint deterministic initialization streams so dormant Phase-II modules cannot change the canonical Phase-I initialization under seed 42.

## 4. Core optimization strategy

### 4.1 Rebuild the strong representation before decomposition

- Use learning rate 0.0001 and gradient accumulation 1 during the representation stage.
- Keep PRISM decomposition, PHU relative weighting, module rescue, and generator pruning at zero initially.
- For epochs 1–12, reproduce the selected CT64 Phase-I contract: reconstruction, ordinal prediction, decoder thresholds, reference anchoring, sex/technology nuisance controls, projection-aware whitening, the active rank8 first-order pathology decoder, and donor-balanced pathology auxiliary supervision.
- Use the canonical Phase-I values `lambda_path_aux=0.3`, `lambda_axis_specificity=0.05`, and `lambda_pathology_tie=0.1`; persist the donor-by-cell-type pathology EMA bank in checkpoints instead of rewarming it after resume.
- During epochs 13–16, smoothly cross-fade the direct Phase-I pathology-decoder contribution and auxiliary loss toward the exact Phase-II setting while common/personal PRISM ramps up. The rank8 parameters remain instantiated and checkpointed; they are not deleted or reinitialized.
- After the cross-fade, preserve the learned shared encoder/decoder state exactly as a warm-started Phase II would, but continue inside the same optimizer, scheduler, run directory, and checkpoint history.
- Do not add automatic pathology-rank selection or a new centering constraint in the primary parity run. Those are scientifically reasonable but must be separate ablations; the primary run first restores the omitted canonical rank8 Phase-I contract.
- Provenance audit: 73 SV6 configurations and the available git history use fixed `pathology_rank=8`; no `learnable_pathology_rank` implementation or executed pathology-rank search was found. The historical learnable-rank run targeted Lie-generator rank with maximum rank 6 and did not alter pathology rank.
- Treat rank 8 as maximum capacity, not as a claim that eight directions are used. At every epoch and stage boundary, log the per-axis singular-value spectrum, participation rank, and 95%-energy rank of the learned cell-type pathology-correction matrix. These are diagnostics only and do not prune or change the primary parity model.
- Do not directly maximize participation ratio. z must gain dimensions because they predict stable biological aggregates, not because an isotropic rank penalty fills unused axes with noise.

### 4.2 Dual module rescue

The module objective has two complementary views.

1. General centered recovery
   - Match train-only donor-by-cell-type module pseudobulk after removing the cell-type module mean.
   - Retain Huber and CCC terms over all 414 modules.
   - Optimize the full decoder view and the explicit PRISM branch view.

2. Pathology-axis recovery
   - Build train-only donor-by-cell-type-by-region module pseudobulk targets.
   - Estimate module coefficient vectors for Braak, Thal, CERAD, LATE, and Lewy while controlling the other pathology axes, sex, age, and region.
   - Use fixed train-donor cross-fitting and split-half reliability masks. Validation donors are never used to construct the target.
   - Align predicted and observed coefficient vectors with correlation/CCC plus a small magnitude-calibration term.
   - Weight unreliable cell-type/module/axis combinations down rather than forcing the model to fit noise.

### 4.3 Cell-type tail protection

- Give neuronal and non-neuronal lineages equal total module-rescue mass.
- Start equal within each lineage.
- Maintain an EMA of each cell type's ceiling-normalized module deficit.
- Apply a bounded deficit multiplier between 0.5 and 2.0 and renormalize within lineage.
- Optimize a worst-tail term over the weakest 25% of cell types, but exclude target combinations that fail the train-only reliability mask.
- Oligodendrocyte is not assigned an arbitrary permanent weight. It receives additional weight only while a reliable deficit persists.

### 4.4 Non-neuronal pseudobulk correction

- Set `non_neuronal_pseudobulk_mode` to `aggregate`.
- Combine the four deterministic 2-cell draws into one 8-cell donor-by-cell-type pseudobulk loss.
- Balance DLPFC and MTG contributions within donor where both regions are available; retain an explicit missing-region mask otherwise.
- Keep state-readout rescue updates disabled in the primary CLS run. The existing state-readout rescue allowlist targets `persistent_attention_pool` and may be enabled only in a separately declared persistent-AGP arm.

### 4.5 PHU and relative weighting

- Bootstrap and persist the PHU donor-by-cell-type bank before applying PHU loss.
- Start alignment only after decoder reconstruction and ordinal discrimination are stable.
- Start relative reconstruction weighting later than alignment.
- Normalize PHU weights to mean one inside each cell type so a difficult lineage is not globally suppressed.
- Use bounded weights, nominally 0.80 to 1.25, and require effective sample-size fraction of at least 0.90.

### 4.6 Generator-cardinality learning

- Do not wait until the final curriculum stage to learn generator importance. Begin a protected search when structured module rescue starts, so the gate receives both reconstruction and biological module gradients.
- Use epochs 25–28 as shadow search: update importance, dual, and constraint EMA state from counterfactual gate views while the actual forward path keeps all 414 generators on.
- Use epochs 29–35 as high-temperature soft gating: let the decoder adapt to continuous gate attenuation while retaining a paired all-on reference path.
- Begin reversible hard-mask learning at epoch 36. Continue automatic cardinality search through epoch 55, but physically delete no generator before the final model is frozen.
- Do not specify a target generator count or a user-chosen minimum count. Optimize the normalized expected active fraction directly, with primal-dual multipliers enforcing the biological and reconstruction constraints.
- Enable `protect_unique_gene_coverage` and retain singleton protection.
- Evaluate paired all-on and gated full reconstruction, generator-isolated reconstruction, and general and AD-module recovery for every cell type on fixed train-only architecture folds.
- Increase sparsity pressure only while all confidence-adjusted constraints remain feasible. Constraint violations automatically raise their dual multipliers and oppose further pruning.
- Accept a new hard mask only after three consecutive train-only architecture audits pass. Otherwise freeze or roll back to the last safe mask.
- Never use validation donors to update gate logits, dual variables, constraint EMA state, or the safe mask. Validation remains checkpoint-selection evidence only.
- Generator count is a secondary efficiency result. Among checkpoints that are biologically indistinguishable within uncertainty and pass every gate, the smaller feasible mask is the tie-breaker.

This gives 31 epochs of generator-importance learning, 27 epochs of soft adaptation, and 20 epochs of protected hard-mask adaptation. The final count may be below, near, or above 97; the frozen rank2 count is a comparator, not a target. The effective lower bound emerges from unique-gene coverage and reconstruction/module non-inferiority rather than a preset cardinality.

## 5. Nominal 55-epoch continuous curriculum

The 55-epoch budget mirrors the selected frozen model lineage: 12 representation epochs, 12 PRISM preparation epochs, 11 module-rescue epochs, and 20 final epochs. The builder must match the baseline's exact optimizer-update count, not only the nominal epoch count.

| Stage | Nominal epochs | Active components | Optimizer regime | Exit evidence |
|---|---:|---|---|---|
| A. CT64 Phase-I bootstrap | 1–3 | Canonical reconstruction/ordinal/reference/nuisance contract plus active first-order pathology rank8 decoder and pathology auxiliary head | learning rate 0.0001; accumulation 1 | finite logits, moving rank8 pathology parameters, and stable decoder thresholds |
| B. CT64 Phase-I biological warmup | 4–12 | Continue the exact canonical Phase-I contract and optimizer regime | learning rate 0.0001; accumulation 1 | epoch-12 Phase-I parity audit and pathology/z diagnostics pass |
| C. Continuous Phase-I to Phase-II cross-fade | 13–24 | Cross-fade the Phase-I direct pathology route at 13–16 while common, personal rank2, and response terms ramp; all parameters remain in one model | decay to learning rate 0.00002; accumulation 2 | Phase-II branch prediction is reliable without losing the Phase-I representation |
| D. Structured module rescue and soft gate search | 25–35 | General centered rescue, AD-axis rescue, 8-cell non-neuronal aggregation, reference EMA; shadow gate search at 25–28 and soft gating at 29–35 | learning rate 0.00002; accumulation 2 | module metrics improve and gate ranking stabilizes without a hard mask |
| E. PHU-relative consolidation and protected hardening | 36–43 | PHU alignment and bounded relative weighting; rescue remains active; reversible hard-mask candidates with no target count | learning rate 0.00002; accumulation 2 | PHU ESS, reconstruction, module, and gate-mask gates pass |
| F. Automatic cardinality consolidation | 44–55 | Continue primal-dual search for the smallest feasible mask, then adapt to the last safe mask | learning rate 0.00002; accumulation 2 | no further feasible reduction, stable safe mask, and post-pruning adaptation complete |

Stage boundaries are controller states stored in every checkpoint. A transition requires both the nominal minimum step and its validation gates. A delayed transition consumes the fixed total budget by shortening the later pruning stage; total training exposure is not increased opportunistically.

## 6. z-specific policy

- Use global balanced-extract participation ratio at stage boundaries; use mini-batch participation only as a live trend.
- Provisional Stage B exit floors are mini-batch clean participation of at least 6 / 32 and global balanced participation of at least 10 / 32.
- The target is not maximal rank. A practical healthy range is global participation around 10–20 / 32 with stable module recovery.
- Track the gap between raw, clean, and sampled participation. A large sampled-to-mean gap is a collapse alarm; a large raw-to-clean loss is a nuisance-refill alarm.
- Require module improvement alongside z expansion. z rank never overrides a module-recovery regression.

## 7. Validation and checkpoint gates

### 7.1 Safety gates

- Reconstruction and branch NLL remain non-inferior to the frozen baseline under one unified evaluator.
- Ordinal nonzero accuracy and balanced recall do not regress beyond predeclared bootstrap uncertainty.
- Generator full and isolated safety violations remain non-positive.
- PHU effective sample-size fraction remains at least 0.90 once relative weighting is active.
- Sex and technology probes remain bounded; cell-type predictability is reported as content retention and is not an exclusion gate.
- No official test donor is loaded.

### 7.2 Module replacement gates

The frozen rank2 e20 values define the validation reference.

- AD-module median must be at least 0.783.
- Pooled AD-module recovery must be at least 0.683.
- Blind-centered recovery must be at least the baseline value minus 0.02 under the same evaluator.
- At least 12 of 24 cell types must be within 0.02 of baseline.
- No reliable cell type may fall more than 0.05 below baseline.
- Oligodendrocyte must reach at least 0.484, corresponding to the baseline value minus 0.05.
- Neuronal and non-neuronal lineage medians must independently pass their baseline-minus-0.02 floors.

Cross-disjoint donor recovery is reported with bootstrap uncertainty. Because validation has only 16 donors, it is a confirmation metric rather than the sole early-stopping target.

The following values freeze the per-cell-type AD-module comparator before the new run. `Provisional floor` is the rank2 value minus 0.05. It is enforced only for cell types passing the train-only split-half reliability contract; reliability is never decided from validation performance. The preferred target for every cell type is at least its rank2 value, not merely its floor.

| Cell type | Lineage | Rank2 e20 AD-module | Provisional floor |
|---|---|---:|---:|
| Astrocyte | non-neuronal | 0.496 | 0.446 |
| Chandelier | neuronal | 0.795 | 0.745 |
| Endothelial | non-neuronal | 0.790 | 0.740 |
| L2/3 IT | neuronal | 0.872 | 0.822 |
| L4 IT | neuronal | 0.909 | 0.859 |
| L5 ET | neuronal | 0.916 | 0.866 |
| L5 IT | neuronal | 0.912 | 0.862 |
| L5/6 NP | neuronal | 0.765 | 0.715 |
| L6 CT | neuronal | 0.632 | 0.582 |
| L6 IT | neuronal | 0.786 | 0.736 |
| L6 IT Car3 | neuronal | 0.904 | 0.854 |
| L6b | neuronal | 0.664 | 0.614 |
| Lamp5 | neuronal | 0.836 | 0.786 |
| Lamp5 Lhx6 | neuronal | 0.812 | 0.762 |
| Microglia-PVM | non-neuronal | 0.723 | 0.673 |
| OPC | non-neuronal | 0.628 | 0.578 |
| Oligodendrocyte | non-neuronal | 0.534 | 0.484 |
| Pax6 | neuronal | 0.826 | 0.776 |
| Pvalb | neuronal | 0.815 | 0.765 |
| Sncg | neuronal | 0.718 | 0.668 |
| Sst | neuronal | 0.857 | 0.807 |
| Sst Chodl | neuronal | 0.605 | 0.555 |
| VLMC | non-neuronal | 0.826 | 0.776 |
| Vip | neuronal | 0.845 | 0.795 |

These are non-inferiority guards, not 24 independent tuning targets. Optimization uses lineage-balanced and tail-robust aggregate losses so validation noise in one cell type cannot drive the run.

### 7.3 Checkpoint ranking

First discard checkpoints that fail any safety or module floor. Rank the remaining checkpoints with this validation-only score:

- 30% pooled AD-module recovery;
- 25% cell-type median AD-module recovery;
- 20% blind-centered recovery;
- 15% ceiling-normalized worst-six-cell-type recovery;
- 10% cross-disjoint recovery.

Reconstruction, generator count, and z participation are gates or diagnostics, not substitutes for the biological ranking.

## 8. Evaluation cadence

- Save every epoch and persist optimizer, scheduler, curriculum stage, PHU state, pathology aggregate bank, module-rescue bank, generator dual/EMA state, and gate mask.
- Run lightweight cap-10 validation module screening every 2 epochs after Stage D begins.
- Run cap-40 full module recovery only at stage boundaries and for final shortlisted checkpoints.
- Fit probes on train and evaluate on validation. Never fit a probe on validation labels.
- Keep the 9 test donors sealed until a final model is frozen independently of test performance.

## 9. Implementation map

### Existing mechanisms to reuse

- non-neuronal aggregate pseudobulk mode;
- state-readout rescue allowlist for persistent-AGP arms only;
- PHU state persistence and bounded relative weighting;
- subgroup CVaR infrastructure;
- learned generator-count gate and unique-gene protection;
- current validation module evaluator and cap-10/cap-40 sampling contracts.
- canonical CT64 Phase-I rank8 pathology decoder and donor-balanced pathology auxiliary head.

### Required implementation changes

1. Build the integrated configuration from the union of the canonical CT64 Phase-I and frozen Phase-II configurations; a Phase-II-only source config is forbidden.
2. Add a checkpointed curriculum controller with step-based transitions, dynamic learning rate, dynamic gradient accumulation, and a persisted Phase-I/II cross-fade scale.
3. Persist donor-by-cell-type pathology aggregate state rather than rewarming it after resume.
4. Add train-only region-balanced general module targets and pathology-axis coefficient targets.
5. Add the pathology-axis module-rescue loss and its split-half reliability mask.
6. Add ceiling-normalized cell-type deficit EMA and module-tail CVaR weighting.
7. Add shadow, soft, and protected-hard gate modes with checkpointed mode transitions.
8. Replace fixed generator-count floors with a train-only primal-dual smallest-feasible-cardinality controller.
9. Add per-cell-type general/AD-module constraints to generator-cardinality learning.
10. Add automatic safe-mask rollback and persist the last safe generator mask.
11. Replace loss-only best-checkpoint selection with the gated module-ranking policy.

## 10. Verification before launch

- Unit tests for every curriculum transition and exact resume equivalence.
- A Phase-I parity test proving that, under seed 42 and the same first-12-epoch update budget, the integrated epoch-12 shared parameters and metrics match the canonical CT64 Phase-I contract within declared numerical tolerance.
- A launch assertion requiring `use_pathology_decoder=true`, `pathology_rank=8`, `pathology_pairwise=false`, `pathology_aux.enabled=true`, and the canonical pathology-loss weights before epoch 1.
- A state-dict test proving rank8 pathology parameters exist at epoch 1, receive gradients during Phase I, remain checkpointed after cross-fade, and are never silently dropped on resume.
- Tests for deterministic per-axis pathology-correction singular-value, participation-rank, and 95%-energy-rank diagnostics under DDP and resume.
- DDP tests for aggregate pseudobulk gradients, cell-type weights, and EMA synchronization.
- Tests proving that validation and test donors cannot enter train-side banks or coefficient targets.
- Tests that an 8-cell aggregate is the exact mean of its four non-overlapping 2-cell draws.
- Tests that CLS configurations reject state-readout rescue updates and persistent-AGP configurations alone may enable them after their declared stage.
- Tests that generator pruning rolls back on a synthetic cell-type module regression.
- Tests that no runtime field imposes a target or minimum generator count and that a synthetic feasible optimum is recovered from constraints.
- Tests proving validation batches cannot update gate logits, dual variables, constraint EMA state, or the safe mask.
- A compressed smoke run that traverses all stages and resumes once; smoke metrics are functional checks only.
- A dry-run manifest fixing seed, split hashes, registry hash, baseline checkpoint hash, step budget, and sealed-test contract.

## 11. Launch status

This document authorizes design work only. No new training run has been launched. Implementation, smoke verification, and a final configuration audit must complete before requesting launch approval.
