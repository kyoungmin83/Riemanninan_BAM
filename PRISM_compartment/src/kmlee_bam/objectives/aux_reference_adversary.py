from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional

import numpy as np
import torch
from torch.utils.data import DataLoader, Subset
from torch.utils.data._utils.collate import default_collate


@dataclass(frozen=True)
class AuxReferenceAdversaryConfig:
    """
    Configuration for the v15 *auxiliary reference batch* that drives the
    reference-only celltype adversary (:class:`ReferenceCelltypeAdversary`).

    The v14 failure mode: the adversary only ever saw the ~3 reference cells that
    happened to land in a 32-cell main batch, so its celltype CLASSIFIER could not
    learn and its gradient-reversal signal to the encoder was useless. v15 gives
    the adversary a DEDICATED, celltype-balanced reference batch of
    ``aux_batch_size`` cells, forwarded FRESH every ``aux_every`` steps (the
    encoder changes every step, so a cached / stale z_perp would give wrong
    gradients).

    OFF by default: ``enabled=False`` or ``aux_batch_size<=0`` ⇒ no aux runner is
    constructed and the training loop is byte-identical to the v14 path.
    """

    enabled: bool = False
    aux_batch_size: int = 0
    aux_every: int = 1
    classifier_pretrain_epochs: int = 0
    grl_pretrain_strength: float = 0.0
    # GRL-STRENGTH ramp (NOT loss scaling — Adam normalizes a constant loss
    # scale away). After the classifier-pretrain window the encoder pressure
    # would otherwise jump to the full configured ``grl_strength`` abruptly;
    # instead ramp the GRL strength linearly from 0 up to ``grl_strength`` over
    # this many adversarial epochs. ``<=0`` ⇒ no ramp (full strength on the
    # first post-pretrain epoch).
    grl_ramp_epochs: int = 3
    celltype_balanced: bool = True
    num_workers: int = 0
    seed: int = 42
    min_cells_per_celltype: int = 1


@dataclass(frozen=True)
class AuxReferenceAdversaryStepStats:
    loss: float
    balanced_accuracy: float
    accuracy: float
    n_reference: float
    pretrain: bool
    grl_strength: float


class ReferenceIndexBank:
    """
    Precomputed bank of dataset-LOCAL indices of reference/control cells,
    organized per cell type for celltype-balanced sampling.

    Reuses the same dataset attributes the :class:`ReferenceAnchorBatchSampler`
    relies on (``row_idx`` maps dataset-local → global zarr rows;
    ``celltype_ids`` and ``is_reference_origin`` are indexed by global rows).
    """

    def __init__(self, dataset: Any, *, min_cells_per_celltype: int = 1) -> None:
        row_idx = np.asarray(getattr(dataset, "row_idx"), dtype=np.int64)
        celltype_ids = np.asarray(getattr(dataset, "celltype_ids"), dtype=np.int64)
        is_reference = np.asarray(
            getattr(dataset, "is_reference_origin"), dtype=bool
        )

        # Global-row reference mask projected onto dataset-local positions.
        ref_global = is_reference[row_idx]
        local_ref = np.nonzero(ref_global)[0].astype(np.int64)
        ct_of_local_ref = celltype_ids[row_idx[local_ref]].astype(np.int64)

        by_ct: Dict[int, np.ndarray] = {}
        for ct in np.unique(ct_of_local_ref):
            members = local_ref[ct_of_local_ref == int(ct)]
            if members.size >= int(min_cells_per_celltype):
                by_ct[int(ct)] = members

        self.by_celltype: Dict[int, np.ndarray] = by_ct
        self.celltypes: np.ndarray = np.asarray(sorted(by_ct.keys()), dtype=np.int64)
        self.all_reference_local: np.ndarray = local_ref
        self.n_reference: int = int(local_ref.size)

    def has_reference(self) -> bool:
        return self.n_reference > 0 and self.celltypes.size > 0

    def sample(
        self,
        n: int,
        rng: np.random.Generator,
        *,
        celltype_balanced: bool = True,
    ) -> np.ndarray:
        """Return ``n`` dataset-local reference indices (sampled with replacement).

        ``celltype_balanced`` round-robins cell types so rare types are not
        swamped by abundant ones (the balanced-accuracy metric the adversary
        optimizes is per-class). Otherwise samples uniformly over all reference
        cells.
        """
        if n <= 0 or not self.has_reference():
            return np.empty(0, dtype=np.int64)

        if not celltype_balanced:
            return rng.choice(self.all_reference_local, size=int(n), replace=True)

        cts = self.celltypes
        # Round-robin cell types in a shuffled order; draw one cell per slot.
        order = rng.permutation(cts.size)
        out = np.empty(int(n), dtype=np.int64)
        for i in range(int(n)):
            ct = int(cts[order[i % cts.size]])
            pool = self.by_celltype[ct]
            out[i] = int(pool[rng.integers(0, pool.size)])
        return out


class AuxReferenceAdversaryRunner:
    """
    Drives the reference-only celltype adversary on a dedicated reference batch.

    The runner only PRODUCES the adversary loss + diagnostics from a FRESH
    encoder forward on ``aux_batch_size`` reference cells; the owning trainer
    performs the ``zero_grad → backward → optimizer.step`` (so the aux update
    reuses the trainer's grad-clip / optimizer machinery and stays fully
    decoupled from the main batch's autograd graph).

    DDP-safety (the v15 crash fix): the aux forward runs on the UNWRAPPED module
    (``system.module`` if DDP-wrapped, else ``system``) via ``_unwrap_system``,
    NOT through the DDP wrapper. A second DDP-wrapped forward-backward per
    iteration would leave DDP's reducer "unfinished" (the aux loss reduces only
    encoder+adversary params, not the decoder the DDP forward registered),
    crashing the next main forward with "Expected to have finished reduction in
    the prior iteration". Running unwrapped keeps the main forward the only DDP
    forward-backward; the trainer then MANUALLY all-reduces the local aux grads
    across ranks (see ``V3Trainer._run_aux_celltype_adv_step``) so the replicas
    stay in sync.

    Encoder-only: the forward is called with ``compute_decoder=False``, so the
    gene-level ordinal decoder (and, in :class:`V3OrdinalBAMSystem`, the
    in-forward tech/celltype adversary attachments) are SKIPPED. The runner only
    needs ``encoder_out.z_s`` / ``mu_p`` / ``sigma_p`` to form the σ-detached
    z_perp and feed the adversary head; the full decoder / reconstruction stack
    would be pure waste every aux step.

    Sampling mirrors :class:`ReferencePriorAnchor`: build a ``Subset`` +
    ``DataLoader`` of reference cells, default-collate one batch, move it to
    device, and run the (unwrapped, encoder-only) forward to obtain ``z_s`` /
    ``mu_p`` / ``sigma_p``.
    """

    def __init__(
        self,
        *,
        dataset: Any,
        config: AuxReferenceAdversaryConfig,
    ) -> None:
        self.config = config
        self.dataset = dataset
        self.bank = ReferenceIndexBank(
            dataset, min_cells_per_celltype=int(config.min_cells_per_celltype)
        )
        self._call_count = 0

    def is_active(self) -> bool:
        return (
            bool(self.config.enabled)
            and int(self.config.aux_batch_size) > 0
            and self.bank.has_reference()
        )

    def is_pretrain_epoch(self, epoch_index: Optional[int]) -> bool:
        """Classifier pre-train window (1-based epoch). During pre-train the
        adversary trains ONLY its classifier (GRL strength forced to
        ``grl_pretrain_strength``, default 0 ⇒ no encoder reversal), so it is a
        strong detective before the encoder pressure ramps in."""
        if epoch_index is None:
            return False
        return int(epoch_index) <= int(self.config.classifier_pretrain_epochs)

    def should_run_this_step(self, step_idx: int) -> bool:
        every = max(1, int(self.config.aux_every))
        return (int(step_idx) % every) == 0

    def grl_ramp_fraction(self, epoch_index: Optional[int]) -> float:
        """Fraction in ``[0, 1]`` of the configured ``grl_strength`` to apply on
        an ADVERSARIAL (non-pretrain) epoch.

        Linearly ramps from 0 (right after the classifier-pretrain window) up to
        1 over ``grl_ramp_epochs`` epochs, so the encoder pressure eases in
        instead of jumping to full strength. ``grl_ramp_epochs<=0`` ⇒ always 1
        (no ramp). Pretrain epochs are handled separately (GRL forced to
        ``grl_pretrain_strength``) and never call this.

        With 1-based ``epoch_index`` and ``classifier_pretrain_epochs=P``:
        the first adversarial epoch is ``P+1`` ⇒ fraction
        ``clamp((epoch - P) / grl_ramp_epochs, 0, 1)`` ⇒ ``1/grl_ramp_epochs``
        on that first epoch, reaching 1 at epoch ``P + grl_ramp_epochs``.
        """
        ramp = int(self.config.grl_ramp_epochs)
        if ramp <= 0 or epoch_index is None:
            return 1.0
        epochs_after_pretrain = int(epoch_index) - int(self.config.classifier_pretrain_epochs)
        frac = float(epochs_after_pretrain) / float(ramp)
        return max(0.0, min(1.0, frac))

    def _sample_batch(self, epoch_index: Optional[int], step_idx: int) -> Dict[str, torch.Tensor]:
        # Deterministic per (epoch, step, call) so all DDP ranks stay schedule-
        # aligned and reruns are reproducible. Different ranks draw different
        # reference cells (rank offset) so the global aux batch is more diverse.
        rank = 0
        if torch.distributed.is_available() and torch.distributed.is_initialized():
            rank = torch.distributed.get_rank()
        seed = (
            int(self.config.seed)
            + 1_000_003 * int(epoch_index or 0)
            + 9_973 * int(step_idx)
            + 131 * int(rank)
        )
        rng = np.random.default_rng(seed)
        local_idx = self.bank.sample(
            int(self.config.aux_batch_size),
            rng,
            celltype_balanced=bool(self.config.celltype_balanced),
        )
        subset = Subset(self.dataset, local_idx.tolist())
        loader = DataLoader(
            subset,
            batch_size=int(local_idx.size),
            shuffle=False,
            num_workers=max(0, int(self.config.num_workers)),
            pin_memory=False,
            drop_last=False,
            collate_fn=default_collate,
        )
        return next(iter(loader))

    def aux_loss_and_stats(
        self,
        trainer: Any,
        *,
        epoch_index: Optional[int],
        step_idx: int,
        sample_latent: bool = True,
    ) -> Optional[tuple[torch.Tensor, AuxReferenceAdversaryStepStats]]:
        """Fresh encoder-only forward on a dedicated reference batch + adversary
        loss.

        Returns ``(loss, stats)`` or ``None`` when inactive / no reference cells.
        The caller backward()s ``loss``, manually all-reduces the aux grads across
        ranks (DDP sync), then steps the optimizer.

        The forward runs on the UNWRAPPED module (``_unwrap_system``) with
        ``compute_decoder=False`` so it (a) never engages DDP's reducer — the v15
        crash fix — and (b) skips the gene-level decoder it does not need.
        """
        if not self.is_active():
            return None

        self._call_count += 1
        batch = self._sample_batch(epoch_index, step_idx)
        batch = trainer._move_batch_to_device(batch)

        system = trainer.system
        # Run the aux forward on the underlying module, NOT the DDP wrapper, so
        # this second-per-iteration forward-backward never engages DDP's reducer
        # (which would crash the next main forward with an unfinished reduction).
        # The replicas are re-synced afterwards by a MANUAL all-reduce of the aux
        # grads in V3Trainer._run_aux_celltype_adv_step.
        unwrapped_system = _unwrap_system(system)
        celltype_adversary = _get_celltype_adversary(system)
        if celltype_adversary is None:
            return None

        pretrain = self.is_pretrain_epoch(epoch_index)
        original_strength = float(celltype_adversary.grl.strength)
        # GRL strength for THIS aux step:
        #   - pretrain      ⇒ grl_pretrain_strength (default 0: no encoder pull,
        #                     classifier-only warmup);
        #   - adversarial   ⇒ configured grl_strength × ramp_fraction (C), easing
        #                     the encoder pressure in instead of an abrupt jump.
        ramp_fraction = 1.0
        if pretrain:
            effective_strength = float(self.config.grl_pretrain_strength)
        else:
            ramp_fraction = self.grl_ramp_fraction(epoch_index)
            effective_strength = original_strength * ramp_fraction
        celltype_adversary.grl.strength = effective_strength

        try:
            with trainer._autocast_context():
                # UNWRAPPED + encoder-only forward (compute_decoder=False): no DDP
                # reducer, no decoder. V3 also skips its in-forward adversary
                # attachments under compute_decoder=False, so model_out carries no
                # celltype_adversary_out — we recompute the adversary on the
                # σ-detached residual below regardless.
                model_out = unwrapped_system(
                    batch,
                    sample_latent=bool(sample_latent),
                    return_all_hidden_states=False,
                    return_attn_diagnostics=False,
                    compute_decoder=False,
                )
                # (A) σ-DETACH: a z_perp = (z_s − μ_c)/σ_c that keeps σ_c attached
                # flows the adversary's REVERSED gradient into σ_c
                # (raw_sigma_embedding is trainable), letting the adversary "cheat"
                # by inflating σ to shrink z_perp instead of cleaning the encoder.
                # Recompute z_perp here with σ DETACHED so the reversed gradient
                # reaches only the encoder (z_s), NOT the prior scale. (All aux
                # cells are reference cells, so the CE is over the full N-cell
                # batch.)
                z_perp_detached = _residualize_sigma_detached(unwrapped_system, model_out)
                is_reference = batch.get("is_reference")
                if is_reference is None:
                    is_reference = torch.ones(
                        z_perp_detached.shape[0],
                        dtype=torch.bool,
                        device=z_perp_detached.device,
                    )
                adv_out = celltype_adversary(
                    z_perp_detached,
                    batch["celltype_id"],
                    is_reference,
                )
                loss = adv_out.loss
        finally:
            celltype_adversary.grl.strength = original_strength

        stats = AuxReferenceAdversaryStepStats(
            loss=float(loss.detach().cpu()),
            balanced_accuracy=float(adv_out.balanced_accuracy.detach().cpu()),
            accuracy=float(adv_out.accuracy.detach().cpu()),
            n_reference=float(adv_out.n_reference.detach().cpu()),
            pretrain=bool(pretrain),
            grl_strength=float(effective_strength),
        )
        return loss, stats


def _unwrap_system(system: Any):
    """Unwrap a possibly DDP-wrapped system to the underlying module."""
    target = system
    if hasattr(target, "module"):
        target = target.module
    return target


def _residualize_sigma_detached(system: Any, model_out: Any) -> torch.Tensor:
    """Recompute z_perp = (z_s − μ_c) / (σ_c.detach() + eps), clamped — the SAME
    transform as :meth:`CellTypePrior.residualize` but with the prior scale σ_c
    DETACHED.

    Rationale (refinement A): μ_c (``mu_embedding``) is frozen, but σ_c
    (``raw_sigma_embedding``) is trainable, so the adversary's gradient-reversal
    would otherwise flow into σ_c and let it cheat by inflating σ to shrink
    z_perp. Detaching σ in the AUX forward routes the reversed gradient to the
    encoder (z_s) only. The MAIN forward's z_perp is untouched.

    Uses the model_out's already-computed ``mu_p`` / ``sigma_p`` (so no second
    prior embedding lookup) and the prior module's own ``residual_eps`` /
    ``residual_clamp`` constants when available (falling back to the
    :class:`CellTypePrior` defaults of 1e-5 / 10.0).
    """
    z_s = model_out.encoder_out.z_s
    mu_p = model_out.mu_p
    sigma_p = model_out.sigma_p

    prior = getattr(_unwrap_system(system), "prior", None)
    residual_eps = float(getattr(prior, "residual_eps", 1e-5))
    residual_clamp = float(getattr(prior, "residual_clamp", 10.0))

    z_perp = (z_s - mu_p) / (sigma_p.detach() + residual_eps)
    return z_perp.clamp(-residual_clamp, residual_clamp)


def _get_celltype_adversary(system: Any):
    """Unwrap a possibly DDP-wrapped system and return its celltype adversary."""
    return getattr(_unwrap_system(system), "celltype_adversary", None)
