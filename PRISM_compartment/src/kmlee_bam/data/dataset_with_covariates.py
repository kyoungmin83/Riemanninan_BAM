from __future__ import annotations

from typing import Optional

import numpy as np
import torch

try:
    from kmlee_bam.data.ordinal_dataset import (
        OrdinalScDataset,
        read_zarr_dataframe,
    )
except ImportError:
    from kmlee_bam.data.ordinal_dataset import OrdinalScDataset, read_zarr_dataframe


class V3OrdinalScDataset(OrdinalScDataset):
    """
    v3 dataset wrapper that adds covariates needed for calibration losses.

    The base dataset already returns `is_reference`, `celltype_id`, and the
    configured technology `batch_id`. v3 additionally needs donor ids for
    donor-balanced reference centering and optional depth covariates for
    uncertainty audits.
    """

    def __init__(
        self,
        *args,
        donor_key: str = "donor_id",
        depth_key: Optional[str] = "Number of UMIs",
        detected_genes_key: Optional[str] = "Genes detected",
        **kwargs,
    ) -> None:
        super().__init__(*args, **kwargs)

        obs = read_zarr_dataframe(self.root, "obs")
        if donor_key not in obs.columns:
            raise KeyError(f"'{donor_key}' not found in zarr obs")

        donor_values = obs[donor_key].astype(str)
        donor_vocab = np.asarray(sorted(donor_values.unique().tolist()), dtype=object)
        donor_to_id = {str(v): i for i, v in enumerate(donor_vocab.tolist())}
        self.donor_key = donor_key
        self.donor_vocab = donor_vocab
        self.donor_ids = np.asarray(
            [donor_to_id[str(v)] for v in donor_values.values],
            dtype=np.int64,
        )

        self.depth_key = depth_key if depth_key in obs.columns else None
        if self.depth_key is None:
            self.depth_values = None
        else:
            self.depth_values = (
                np.asarray(obs[self.depth_key].astype(float).values, dtype=np.float32)
            )

        self.detected_genes_key = (
            detected_genes_key if detected_genes_key in obs.columns else None
        )
        if self.detected_genes_key is None:
            self.detected_genes_values = None
        else:
            self.detected_genes_values = np.asarray(
                obs[self.detected_genes_key].astype(float).values,
                dtype=np.float32,
            )

    def __getitem__(self, idx: int):
        item = super().__getitem__(idx)
        r = int(item["row_index"].item())
        item["donor_id"] = torch.tensor(int(self.donor_ids[r]), dtype=torch.long)
        if self.depth_values is not None:
            item["depth_value"] = torch.tensor(float(self.depth_values[r]), dtype=torch.float32)
        if self.detected_genes_values is not None:
            item["detected_genes"] = torch.tensor(
                float(self.detected_genes_values[r]),
                dtype=torch.float32,
            )
        return item

