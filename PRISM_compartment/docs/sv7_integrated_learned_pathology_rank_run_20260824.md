# SV7 integrated PRISM learned-pathology-rank run (2026-08-24)

## User-authorized design

- One integrated architecture, one continuous Phase I to Phase II curriculum,
  one optimizer lineage, and one checkpoint history.
- Keep the SV6 integrated-primary design except that the pathology-conditioned
  decoder correction rank is learned instead of fixed at rank 8.
- The user explicitly authorized the initial SV7 launch and unattended restart
  within this same design. A new design or destructive data change still
  requires a new decision.
- The nine official test donors remain sealed.

## Learned pathology rank

- `decoder.pathology_rank = 0` is an AUTO sentinel, not a rank-zero model.
- The tensor allocation is the algebraic capacity:

$$R_{\mathrm{capacity}} = \min(n_{\mathrm{celltype}}n_{\mathrm{axis}}, n_{\mathrm{gene}}) = \min(24\times4, 13498) = 96.$$

- All 96 hard candidates are open at initialization with identical gate
  initialization. There is no target rank, user-chosen maximum rank, or
  minimum rank; rank zero is permitted if the task finds no useful
  cell-type-specific correction.
- Epochs 1--2: task gradients learn all soft gates; no cardinality pressure.
- Epochs 3--12: task gradients plus a target-free cardinality penalty.
- Epochs 13--16: exact hard forward with straight-through adaptation during
  the pathology cross-fade.
- Epoch 17 onward: freeze the learned mask. This prevents meaningless collapse
  after the pathology route reaches zero contribution.
- Rank-gate learning rate is `3e-5`; cardinality loss fraction is `0.0025`.
- Gate parameters, current epoch, final mask, optimizer state, and curriculum
  state are checkpointed.

## SV6-equivalent global budgets on four SV7 GPUs

$$B_{\mathrm{Phase\ I}} = 4\times96 = 384.$$

$$B_{\mathrm{Phase\ II}} = 4\times48\times2 = 384.$$

- Both phases use 5,255 optimizer updates per epoch, exactly matching the SV6
  six-GPU geometry.
- Module rescue uses 24 steps per rank on four GPUs, preserving the SV6 global
  draw contract:

$$N_{\mathrm{rescue}} = 24\times4 = 16\times6 = 96.$$

- Generator cardinality remains SV6-equivalent: shadow epochs 13--16, soft
  epochs 17--24, hard pruning from epoch 25, and no manual minimum generator
  count.
- BAM pathology-aware uncertainty, PHU alignment, within-cell-type relative
  weighting, and the minimum relative ESS fraction of 0.90 remain enabled.
- PRISM personal rank remains 2; this is distinct from the learned pathology
  correction rank.

## Immutable runtime and provenance

- Runtime: `/home/kmlee/project_local/kmlee_bam_integrated_pathranklearn_20260824`
- Base config:
  `configs/generated/train_config_prism_integrated_pathranklearn_s42_sv7_20260824.json`
- Base-config SHA256:
  `97abe29ed657d9189bd4c5c2d706b26275be7eaf84bddb6cd30690143b8d1f36`
- Train-only module-rescue artifact SHA256:
  `b32a5156d835bf2049dbe61f838e686e460aaa1cd3130c3cf99b8ca1d7de857c`
- Output:
  `/gstorage_data/kmlee/project/riemann_bam/outputs/SEA_AD_OUTPUTS/train_runs/kmlee_bam_prism_integrated_pathranklearn_s42_sv7_20260824`
- Preflight report:
  `artifacts/preflight_prism_integrated_pathranklearn_sv7_20260824.json`

## Unattended recovery

- User service: `prism-integrated-pathranklearn-sv7.service`
- The service is enabled under the lingering user manager, so it starts again
  after host reboot.
- Transient process, CUDA, NCCL, or I/O failure causes automatic retry without
  additional user approval.
- Once an epoch-boundary checkpoint exists, recovery restores the system,
  optimizer, scheduler, AMP scaler, PHU/EMA banks, generator/rank gates, module
  rescue state, history, and integrated curriculum state in the same output
  lineage.
- Before the first epoch checkpoint, a failure restarts epoch 1 from scratch;
  exact mid-epoch sampler/RNG state is intentionally not claimed.

## Launch audit and first live evidence

- The first attempted launch exposed a six-GPU module-rescue step count in the
  four-GPU runtime (`16` configured versus `24` required). It failed before any
  optimizer step or checkpoint and was corrected to preserve 96 global draws.
- A subsequent live check exposed the 1.5-fold update-budget mismatch that
  would result from retaining per-GPU batch 64 on four GPUs. It was stopped
  before a checkpoint and corrected to the SV6-equivalent 96/48 geometry.
- The final run started on 2026-08-24 and reached optimizer step 500 in epoch 1.
  At that point: loss `1.997`, gradient norm `0.58`, projected-z effective
  dimension `16.0/32`, PHU ESS `0.99`, and pathology decoder magnitude
  `0.0142`. Four GPUs used approximately 43.3--43.6 GB of 49.1 GB each with no
  OOM.
