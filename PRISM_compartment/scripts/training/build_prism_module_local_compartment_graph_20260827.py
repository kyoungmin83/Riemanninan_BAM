#!/usr/bin/env python3
"""Build the immutable sparse module graph used by nonlinear module-local PRISM.

Edges depend only on overlap between gene-module memberships in the supplied
registry.  No cell, donor, pathology, validation, or test outcome is opened.
For every module we retain the strongest directed overlaps and row-normalize
them.  Isolated singleton/residual modules retain an all-zero row and therefore
remain independent compartments at runtime.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile

import numpy as np

from kmlee_bam.data.module_local_reliability import (
    MODULE_LOCAL_COMPARTMENT_GRAPH_SCHEMA,
    sha256_file,
)
from kmlee_bam.gene_module.gene_module_registry import GeneModuleRegistry


def build_adjacency(
    membership: np.ndarray,
    *,
    topk: int,
    minimum_jaccard: float,
) -> np.ndarray:
    if membership.ndim != 2:
        raise ValueError("membership must have shape [module,gene]")
    if topk <= 0:
        raise ValueError("topk must be positive")
    if not 0.0 <= minimum_jaccard <= 1.0:
        raise ValueError("minimum_jaccard must lie in [0,1]")
    binary = np.asarray(membership > 0, dtype=np.float64)
    intersection = binary @ binary.T
    size = binary.sum(axis=1)
    union = size[:, None] + size[None, :] - intersection
    jaccard = np.divide(
        intersection,
        np.maximum(union, 1.0),
        out=np.zeros_like(intersection),
        where=union > 0.0,
    )
    np.fill_diagonal(jaccard, 0.0)
    adjacency = np.zeros_like(jaccard, dtype=np.float32)
    for row in range(jaccard.shape[0]):
        eligible = np.flatnonzero(jaccard[row] >= float(minimum_jaccard))
        if eligible.size == 0:
            continue
        # Stable tie-breaking: greater Jaccard first, then lower module index.
        order = np.lexsort((eligible, -jaccard[row, eligible]))
        chosen = eligible[order[: int(topk)]]
        weights = jaccard[row, chosen]
        total = float(weights.sum())
        if total > 0.0:
            adjacency[row, chosen] = (weights / total).astype(np.float32)
    return adjacency


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--topk", type=int, default=4)
    parser.add_argument("--minimum-jaccard", type=float, default=0.05)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    registry_path = Path(args.registry).expanduser().resolve()
    registry = GeneModuleRegistry.load(registry_path)
    adjacency = build_adjacency(
        registry.membership_binary.detach().cpu().numpy(),
        topk=int(args.topk),
        minimum_jaccard=float(args.minimum_jaccard),
    )
    output = Path(args.output).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        suffix=".npz", dir=output.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
    try:
        np.savez_compressed(
            temporary,
            schema_version=np.asarray(MODULE_LOCAL_COMPARTMENT_GRAPH_SCHEMA),
            adjacency=adjacency,
            module_names=np.asarray(registry.module_names),
            registry_sha256=np.asarray(sha256_file(registry_path)),
            topk=np.asarray(int(args.topk)),
            minimum_jaccard=np.asarray(float(args.minimum_jaccard)),
            validation_donors_used=np.asarray(False),
            test_donors_used=np.asarray(False),
        )
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()
    nonempty = adjacency.sum(axis=1) > 0.0
    print(
        "[module-local compartment graph] "
        f"modules={adjacency.shape[0]} nonempty={int(nonempty.sum())} "
        f"edges={int((adjacency > 0.0).sum())} "
        f"sha256={sha256_file(output)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
