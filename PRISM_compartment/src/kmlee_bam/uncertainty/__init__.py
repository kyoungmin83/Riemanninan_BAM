"""
Pathology-Aware Hierarchical Uncertainty (PHU / model v7a).

This package implements the bank/ANCOVA machinery that gives the BAM
uncertainty head a meaningful supervision target. See:
    doc/model_v7_design_pathology_aware_uncertainty.md  (rev 3, design)
    doc/model_v7a_implementation_plan.md                (this code's spec)

Public surface:
    DifficultyProfile, DifficultyProfileConfig
    stratified_cell_sample
    compute_noise_score, normalize_rec_err
    NoiseScoreBank
    ANCOVAFit, ANCOVAConfig
    compute_u_total
"""

from kmlee_bam.uncertainty.sampler import stratified_cell_sample
from kmlee_bam.uncertainty.difficulty_profile import (
    DifficultyProfile,
    DifficultyProfileConfig,
)
from kmlee_bam.uncertainty.noise_score import (
    compute_noise_score,
    normalize_rec_err,
)
from kmlee_bam.uncertainty.bank import NoiseScoreBank
from kmlee_bam.uncertainty.ancova import ANCOVAConfig, ANCOVAFit
from kmlee_bam.uncertainty.decomposition import compute_u_total

__all__ = [
    "stratified_cell_sample",
    "DifficultyProfile",
    "DifficultyProfileConfig",
    "compute_noise_score",
    "normalize_rec_err",
    "NoiseScoreBank",
    "ANCOVAConfig",
    "ANCOVAFit",
    "compute_u_total",
]
