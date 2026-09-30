# PRISM v5: online joint generator-count learning

Date: 2026-08-05

## Decision

The generator count is learned inside the ordinary full training run.  The
W44/R10/C10 fixed-​K grid is not used by this experiment.  All 64 training
donors remain in the optimizer dataset; official validation controls early
stopping and the official test split is not opened during training.

PRISM retains its 414 named module coefficients for the explicit
region/age/common/personal/response decomposition.  The learned count applies
to the decoder's 414 Lie-generator routes.  A gate closes the complete route,
including both its linear Lie action and affine translation; it does not
silently delete named PRISM coefficient axes.

## Online architecture variable

Each generator (j) has one fp32 hard Binary-Concrete logit
​\(\alpha_j\).  A training forward samples an exact binary switch
​\(g_j \in \{0,1\}\), while the backward pass uses a straight-through
relaxation.  Validation uses a deterministic binary mask.  The learned count
is

$$K = \sum_{j=1}^{414}\mathbf{1}(g_j = 1).$$

Epochs 1–2 are an exact all-on warm-up.  From epoch 3, model weights, PRISM
weights and gate logits are optimized in the same run.  Temperature anneals
from 2.0 to 0.3 through epoch 18.  The run permits at most 36 epochs and uses
validation early stopping after epoch 18 with patience 8.

## Minimum-cardinality objective

For every ordinary minibatch, the gated candidate and an all-on
counterfactual use the same cells, latent state and current weights.  Two
paired guards are computed:

1. full ordinal NLL, including the normal full decoder;
2. generator-isolated ordinal NLL, with Direct, PRISM, mixer, pathology and
   interaction bypass routes excluded.

The optimization target is

$$L_{\mathrm{count}} = \frac{\mathbb{E}[K]}{414} + \sum_m \lambda_m[v_m]_+ + \frac{\rho}{2}\sum_m[v_m]_+^2,$$

where each violation equals the gated NLL minus its paired all-on NLL and its
predeclared non-inferiority margin.  The full margin equals 1% of paired NLL;
the isolated margin equals 1.5%.  Dual variables are averaged across DDP ranks
and updated once per real optimizer step.  The weighted term added to the
ordinary objective equals (0.05L_{\mathrm{count}}).

This is not magnitude pruning.  The forward gate is exactly zero or one, so a
route cannot evade sparsity by shrinking its gate and inflating its weights.

## Safety contracts

- The gate and its probability calculations remain fp32 under bf16 model
  autocast.
- The exact same sampled gate is used throughout one ordinary forward.
- The 10 singleton generators are always active.  Unique-coverage modules are
  not blanket-protected: doing so would pre-fix 326 of 414 routes and defeat
  count learning.  Their removal is instead governed by the paired full and
  isolated non-inferiority constraints.
- At least 32 generators must remain active.
- PRISM module-rescue uses the epoch's deterministic exact-hard mask.  Its
  straight-through gradient may update only the gate logits plus the existing
  rescue allowlist, so module recovery can keep necessary routes alive without
  unfreezing the encoder, nuisance projector, baselines or shared mixer.
- CP-BAM and metric-contrast decoder re-forwards are fail-closed while this
  mode is enabled because those paths do not yet accept the sampled gate.
- External `generator_budget_search` must remain disabled.
- The old R10 frozen-signature/leakage backend is irrelevant to this online
  experiment and is not silently substituted.

## Checkpoint and audit

Gate logits, protected mask, dual state and current gate epoch are persistent
model state and therefore travel with every checkpoint.  Each epoch also
writes `joint_generator_count_latest.json`, containing deterministic (K),
expected (K), temperature, selected generator IDs and a mask SHA256.

The best checkpoint is selected by the existing validation composite augmented
with gated full NLL, gated isolated NLL and expected active count.  An epoch is
not allowed to replace `checkpoint_best.pt` unless both paired validation
non-inferiority violations satisfy \(v_{\mathrm{full}} \le 0\) and
\(v_{\mathrm{isolated}} \le 0\).  If no pruned epoch is feasible, selection
fails closed to the last best feasible all-on/larger architecture rather than
releasing an underperforming small \(K\).  The official test split remains
unopened.  A final fixed-mask consolidation run may be done after this
discovery run, but it must consume the discovered mask and may not run a
candidate-​(K) search.

## Immutable starting assets

- CT64 canonical checkpoint SHA256:
  `ef3f6241052b4e691132b9b3fd098fb848f82cfb79a43e9473f30bd230fab44a`
- Previous PRISM v2a epoch-11 checkpoint SHA256:
  `9bee5aba99e38453499d6cd2333831509b811945a33e77c2b381eb3ea655c555`

The new run uses a fresh output directory and cannot overwrite either asset.
