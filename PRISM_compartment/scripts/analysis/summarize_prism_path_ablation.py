#!/usr/bin/env python
"""Summarize paired PRISM path ablations from module fingerprints."""
from __future__ import annotations

import argparse
import csv
import json
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata


DISPLAY = {
    "full": "Full model",
    "common_pathology_off": "Common pathology OFF",
    "pathology_interaction_off": "Pathology interaction OFF",
    "all_shared_pathology_off": "All shared pathology OFF",
    "personal_rank2_baseline_off": "Rank-2 personal baseline OFF",
    "personal_pathology_response_off": "Personal pathology response OFF",
    "module_local_off": "Module-local path OFF",
    "module_local_nonlinear_off": "Compartmental nonlinearity OFF",
    "all_personal_off": "All personal paths OFF",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-npz", required=True)
    parser.add_argument("--provenance-json", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--disjoint-splits", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260902)
    return parser.parse_args()


def spearman(left: np.ndarray, right: np.ndarray) -> float:
    left = np.asarray(left, dtype=np.float64).reshape(-1)
    right = np.asarray(right, dtype=np.float64).reshape(-1)
    valid = np.isfinite(left) & np.isfinite(right)
    if valid.sum() < 3:
        return float("nan")
    left = left[valid]
    right = right[valid]
    if np.ptp(left) <= 1.0e-12 or np.ptp(right) <= 1.0e-12:
        return float("nan")
    return float(np.corrcoef(rankdata(left), rankdata(right))[0, 1])


def zscore_with_missing(value: np.ndarray) -> np.ndarray | None:
    value = np.asarray(value, dtype=np.float64)
    valid = np.isfinite(value)
    if valid.sum() == 0:
        return None
    filled = np.where(valid, value, np.nanmean(value[valid]))
    scale = float(np.std(filled))
    if not np.isfinite(scale) or scale <= 1.0e-12:
        return None
    return (filled - float(np.mean(filled))) / scale


def partial_coefficients(response: np.ndarray, predictors: list[np.ndarray]) -> np.ndarray | None:
    columns = []
    for predictor in predictors:
        standardized = zscore_with_missing(predictor)
        if standardized is not None:
            columns.append(standardized)
    if not columns:
        return None
    design = np.column_stack((np.ones(response.shape[0]), *columns))
    coefficients, *_ = np.linalg.lstsq(design, response, rcond=None)
    return coefficients[1:]


def evaluate_variant(
    predicted: np.ndarray,
    observed: np.ndarray,
    count: np.ndarray,
    covariates: np.ndarray,
    donor_indices: np.ndarray,
) -> dict[str, float]:
    raw_values = []
    centered_values = []
    ad_values = []
    axis_values = []
    for celltype in range(predicted.shape[1]):
        indices = donor_indices[count[donor_indices, celltype] >= 10]
        if indices.size < 8:
            continue
        model = predicted[indices, celltype]
        truth = observed[indices, celltype]
        raw_values.append(spearman(model, truth))
        centered_values.append(
            spearman(
                model - np.nanmean(model, axis=0, keepdims=True),
                truth - np.nanmean(truth, axis=0, keepdims=True),
            )
        )

        ad_predictors = [
            covariates[indices, 0],
            covariates[indices, 4],
            covariates[indices, 5],
            covariates[indices, 6],
            covariates[indices, 7],
        ]
        model_beta = partial_coefficients(model, ad_predictors)
        truth_beta = partial_coefficients(truth, ad_predictors)
        if model_beta is not None and truth_beta is not None:
            ad_values.append(spearman(model_beta[0], truth_beta[0]))

        five_axes = [covariates[indices, column] for column in (1, 2, 3, 4, 5)]
        model_axes = partial_coefficients(model, five_axes)
        truth_axes = partial_coefficients(truth, five_axes)
        if model_axes is not None and truth_axes is not None:
            axis_values.extend(
                spearman(model_axes[axis], truth_axes[axis])
                for axis in range(min(model_axes.shape[0], truth_axes.shape[0]))
            )

    def median(values: list[float]) -> float:
        finite = np.asarray([value for value in values if np.isfinite(value)])
        return float(np.median(finite)) if finite.size else float("nan")

    return {
        "general_raw_median": median(raw_values),
        "general_centered_median": median(centered_values),
        "ad_module_direction_median": median(ad_values),
        "five_axis_direction_median": median(axis_values),
        "n_celltypes": int(np.isfinite(centered_values).sum()),
    }


def donor_disjoint(
    predicted: np.ndarray,
    observed: np.ndarray,
    count: np.ndarray,
    covariates: np.ndarray,
    n_splits: int,
    seed: int,
) -> dict[str, float | list[float]]:
    rng = np.random.default_rng(seed)
    ad_cross = []
    axis_cross = []
    donors = np.arange(predicted.shape[0])
    for split_index in range(n_splits):
        permutation = rng.permutation(donors)
        left = permutation[: len(permutation) // 2]
        right = permutation[len(permutation) // 2 :]
        for celltype in range(predicted.shape[1]):
            left_valid = left[count[left, celltype] >= 10]
            right_valid = right[count[right, celltype] >= 10]
            if left_valid.size < 6 or right_valid.size < 6:
                continue
            model_left = predicted[left_valid, celltype]
            truth_right = observed[right_valid, celltype]
            model_ad = partial_coefficients(
                model_left,
                [
                    covariates[left_valid, 0],
                    covariates[left_valid, 4],
                    covariates[left_valid, 5],
                ],
            )
            truth_ad = partial_coefficients(
                truth_right,
                [
                    covariates[right_valid, 0],
                    covariates[right_valid, 4],
                    covariates[right_valid, 5],
                ],
            )
            if model_ad is not None and truth_ad is not None:
                ad_cross.append(spearman(model_ad[0], truth_ad[0]))

            model_axes = partial_coefficients(
                model_left,
                [covariates[left_valid, column] for column in (1, 2, 3, 4, 5)],
            )
            truth_axes = partial_coefficients(
                truth_right,
                [covariates[right_valid, column] for column in (1, 2, 3, 4, 5)],
            )
            if model_axes is not None and truth_axes is not None:
                axis_cross.extend(
                    spearman(model_axes[axis], truth_axes[axis])
                    for axis in range(min(model_axes.shape[0], truth_axes.shape[0]))
                )

    def summary(values: list[float]) -> tuple[float, list[float]]:
        finite = np.asarray([value for value in values if np.isfinite(value)])
        if finite.size == 0:
            return float("nan"), [float("nan"), float("nan")]
        return float(np.median(finite)), [
            float(np.quantile(finite, 0.025)),
            float(np.quantile(finite, 0.975)),
        ]

    ad_median, ad_interval = summary(ad_cross)
    axis_median, axis_interval = summary(axis_cross)
    return {
        "ad_module_cross_disjoint_median": ad_median,
        "ad_module_cross_disjoint_split_95_interval": ad_interval,
        "five_axis_cross_disjoint_median": axis_median,
        "five_axis_cross_disjoint_split_95_interval": axis_interval,
        "n_ad_comparisons": len(ad_cross),
        "n_axis_comparisons": len(axis_cross),
    }


def bootstrap_interval(values: np.ndarray) -> list[float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return [float("nan"), float("nan")]
    return [float(np.quantile(finite, 0.025)), float(np.quantile(finite, 0.975))]


def main() -> None:
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    data = np.load(args.input_npz, allow_pickle=False)
    variants = [str(value) for value in data["variant_names"]]
    model_fp = data["model_fp"].astype(np.float64)
    observed_fp = data["observed_fp"].astype(np.float64)
    count = data["group_count"].astype(np.int64)
    covariates = data["donor_covariates"].astype(np.float64)
    donor_nll = data["donor_nll"].astype(np.float64)
    donor_names = [str(value) for value in data["donor_names"]]
    celltype_names = [str(value) for value in data["celltype_names"]]
    module_names = [str(value) for value in data["module_names"]]
    with open(args.provenance_json, encoding="utf-8") as handle:
        provenance = json.load(handle)

    donor_indices = np.arange(len(donor_names))
    point = {
        variant: evaluate_variant(
            model_fp[index], observed_fp, count, covariates, donor_indices
        )
        for index, variant in enumerate(variants)
    }
    disjoint = {
        variant: donor_disjoint(
            model_fp[index],
            observed_fp,
            count,
            covariates,
            args.disjoint_splits,
            args.seed,
        )
        for index, variant in enumerate(variants)
    }

    rng = np.random.default_rng(args.seed)
    full_index = variants.index("full")
    metric_names = (
        "general_centered_median",
        "ad_module_direction_median",
        "five_axis_direction_median",
    )
    bootstrap_metric_delta = {
        variant: {metric: [] for metric in metric_names} for variant in variants
    }
    bootstrap_nll_delta = {variant: [] for variant in variants}
    for bootstrap_index in range(args.bootstrap):
        sampled = rng.choice(donor_indices, size=len(donor_indices), replace=True)
        full_metrics = evaluate_variant(
            model_fp[full_index], observed_fp, count, covariates, sampled
        )
        for variant_index, variant in enumerate(variants):
            variant_metrics = evaluate_variant(
                model_fp[variant_index], observed_fp, count, covariates, sampled
            )
            for metric in metric_names:
                bootstrap_metric_delta[variant][metric].append(
                    variant_metrics[metric] - full_metrics[metric]
                )
            bootstrap_nll_delta[variant].append(
                float(np.mean(donor_nll[variant_index, sampled] - donor_nll[full_index, sampled]))
            )

    rows = []
    for variant_index, variant in enumerate(variants):
        nll_delta_by_donor = donor_nll[variant_index] - donor_nll[full_index]
        row = {
            "variant": variant,
            "label": DISPLAY.get(variant, variant),
            "mean_nll": float(np.mean(donor_nll[variant_index])),
            "delta_nll_off_minus_full": float(np.mean(nll_delta_by_donor)),
            "delta_nll_bootstrap_95ci": bootstrap_interval(
                np.asarray(bootstrap_nll_delta[variant])
            ),
            **point[variant],
            **disjoint[variant],
        }
        for metric in metric_names:
            row[f"delta_{metric}_off_minus_full"] = (
                point[variant][metric] - point["full"][metric]
            )
            row[f"delta_{metric}_bootstrap_95ci"] = bootstrap_interval(
                np.asarray(bootstrap_metric_delta[variant][metric])
            )
        rows.append(row)

    report = {
        "interpretation": (
            "Frozen-checkpoint conditional reliance. Positive OFF-minus-full NLL means the "
            "removed path helped the fitted model. Recovery changes are not retrained ablations "
            "and must not be interpreted causally."
        ),
        "n_donors": len(donor_names),
        "n_celltypes": len(celltype_names),
        "n_evaluable_modules": len(module_names),
        "bootstrap_resamples": args.bootstrap,
        "donor_disjoint_splits": args.disjoint_splits,
        "source_provenance": provenance,
        "rows": rows,
    }
    with open(os.path.join(args.outdir, "path_ablation_summary.json"), "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)

    tsv_path = os.path.join(args.outdir, "path_ablation_summary.tsv")
    scalar_columns = [
        "variant",
        "label",
        "mean_nll",
        "delta_nll_off_minus_full",
        "general_raw_median",
        "general_centered_median",
        "delta_general_centered_median_off_minus_full",
        "ad_module_direction_median",
        "delta_ad_module_direction_median_off_minus_full",
        "ad_module_cross_disjoint_median",
        "five_axis_direction_median",
        "delta_five_axis_direction_median_off_minus_full",
        "five_axis_cross_disjoint_median",
    ]
    with open(tsv_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=scalar_columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    off_rows = [row for row in rows if row["variant"] != "full"]
    labels = [row["label"] for row in off_rows]
    y = np.arange(len(off_rows))[::-1]
    navy = "#183B6B"
    teal = "#159A8C"
    orange = "#E87924"
    fig, axes = plt.subplots(1, 3, figsize=(18, 7.4), constrained_layout=True)

    nll_value = np.asarray([row["delta_nll_off_minus_full"] for row in off_rows])
    nll_ci = np.asarray([row["delta_nll_bootstrap_95ci"] for row in off_rows])
    axes[0].errorbar(
        nll_value,
        y,
        xerr=np.vstack((nll_value - nll_ci[:, 0], nll_ci[:, 1] - nll_value)),
        fmt="o",
        color=navy,
        ecolor="#8DA3BF",
        capsize=3,
    )
    axes[0].axvline(0, color="#8B97A5", linestyle="--", linewidth=1)
    axes[0].set_yticks(y, labels)
    axes[0].set_xlabel("Validation NLL change after path removal\n(OFF minus full; positive = helpful path)")
    axes[0].set_title("A. Conditional reconstruction reliance", loc="left", fontweight="bold")

    centered = np.asarray([row["general_centered_median"] for row in off_rows])
    axes[1].scatter(centered, y, color=teal, s=55, zorder=3)
    axes[1].axvline(point["full"]["general_centered_median"], color=teal, linestyle="--")
    for ypos, row in zip(y, off_rows):
        axes[1].plot(
            [point["full"]["general_centered_median"], row["general_centered_median"]],
            [ypos, ypos],
            color="#B7DCD7",
            linewidth=2,
            zorder=1,
        )
    axes[1].set_yticks(y, [])
    axes[1].set_xlabel("Donor-centered general-module recovery (Spearman rho)")
    axes[1].set_title("B. General-module recovery", loc="left", fontweight="bold")

    ad = np.asarray([row["ad_module_direction_median"] for row in off_rows])
    dis = np.asarray([row["ad_module_cross_disjoint_median"] for row in off_rows])
    axes[2].scatter(ad, y + 0.12, color=orange, s=55, label="Same-donor AD-module direction")
    axes[2].scatter(dis, y - 0.12, color=navy, marker="s", s=42, label="Donor-disjoint AD-module direction")
    axes[2].axvline(point["full"]["ad_module_direction_median"], color=orange, linestyle="--")
    axes[2].axvline(disjoint["full"]["ad_module_cross_disjoint_median"], color=navy, linestyle=":")
    axes[2].set_yticks(y, [])
    axes[2].set_xlabel("AD-module direction recovery (Spearman rho)")
    axes[2].set_title("C. Disease-axis recovery", loc="left", fontweight="bold")
    axes[2].legend(frameon=False, fontsize=9, loc="lower right")

    for axis in axes:
        axis.grid(axis="x", color="#E5EAF0", linewidth=0.8)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "SV6 compartmental e52: frozen-path conditional reliance audit",
        fontsize=18,
        fontweight="bold",
        color=navy,
    )
    fig.savefig(os.path.join(args.outdir, "path_ablation_summary.png"), dpi=220, bbox_inches="tight")
    fig.savefig(os.path.join(args.outdir, "path_ablation_summary.pdf"), bbox_inches="tight")
    plt.close(fig)

    markdown = [
        "# SV6 compartmental epoch-52 frozen path-ablation audit",
        "",
        "This is a paired, validation-only conditional-reliance audit. The checkpoint is frozen,",
        "the same sampled cells are used for every variant, and one explicit score component is",
        "subtracted at inference. It is not a retrained ablation and does not establish causality.",
        "",
        f"- Donors: {len(donor_names)} validation donors",
        f"- Cell types: {len(celltype_names)}",
        f"- Evaluable modules: {len(module_names)} of 414 registered modules",
        "- Technical-zero OFF: not performed because it changes training-loss weighting, not an inference path",
        "",
        "| Path removed | NLL change | General centered | AD-module | AD donor-disjoint |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in off_rows:
        markdown.append(
            f"| {row['label']} | {row['delta_nll_off_minus_full']:+.6f} | "
            f"{row['general_centered_median']:.3f} | {row['ad_module_direction_median']:.3f} | "
            f"{row['ad_module_cross_disjoint_median']:.3f} |"
        )
    with open(os.path.join(args.outdir, "README.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(markdown) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
