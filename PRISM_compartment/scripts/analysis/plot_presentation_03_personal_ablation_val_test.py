#!/usr/bin/env python3
"""Render validation/test personal and response branch ablations."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


PRISM_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_VALIDATION = (
    PRISM_ROOT.parent
    / "kmlee_bam"
    / "analysis_outputs"
    / "prism_rank2_personal_audit_20260811"
    / "integrated_audit.json"
)
DEFAULT_TEST = (
    PRISM_ROOT
    / "analysis_outputs"
    / "presentation_results_20260819"
    / "00_locked_test_rank2_e20"
    / "rank2_e20_branch_ablation_test.json"
)
DEFAULT_OUTPUT = (
    PRISM_ROOT
    / "analysis_outputs"
    / "presentation_results_20260819"
    / "03_personal_ablation_val_test"
)

PANELS = {
    "personal": [
        (
            "full_without_personal_minus_full",
            "직접 personal 항 제거",
            "병리와 무관한 개인 기본 특성 제거",
        ),
        (
            "personal_matched_shuffle_minus_full",
            "Donor-matched\npersonal code 셔플",
            "비슷한 공변량의 다른 donor code로 교체",
        ),
    ],
    "response": [
        (
            "full_without_response_minus_full",
            "Personal×병리반응 제거",
            "같은 병리에 다르게 반응하는 항 제거",
        ),
        (
            "personal_zero_all_minus_full",
            "Personal 두 항 모두 제거",
            "기본 특성과 병리반응을 동시에 제거",
        ),
    ],
}

COLORS = {
    "navy": "#123B73",
    "blue": "#1F67B1",
    "orange": "#E67E32",
    "dark": "#172033",
    "muted": "#667085",
    "grid": "#D9E1EA",
    "connector": "#AEB9C6",
    "panel_a": "#F1F7FC",
    "panel_b": "#FFF7F0",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_style() -> None:
    regular = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
    bold = Path("/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc")
    for font_path in (regular, bold):
        if font_path.exists():
            font_manager.fontManager.addfont(font_path)
    family = font_manager.FontProperties(fname=regular).get_name() if regular.exists() else "DejaVu Sans"
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [family, "DejaVu Sans"],
            "axes.unicode_minus": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def extract_split(path: Path, split: str) -> tuple[pd.DataFrame, dict]:
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    integrated_validation = "conditions" in payload
    if integrated_validation:
        comparisons = {
            key: value["donor_balanced_delta_nll_variant_minus_full"]
            for key, value in payload["conditions"].items()
            if "donor_balanced_delta_nll_variant_minus_full" in value
        }
    else:
        comparisons = payload["donor_bootstrap_nll_comparisons"]
    rows = []
    for panel, definitions in PANELS.items():
        for order, (key, label, detail) in enumerate(definitions):
            result = comparisons[key]
            lower, upper = result["bootstrap_95ci"]
            rows.append(
                {
                    "split": split,
                    "panel": panel,
                    "order": order,
                    "comparison": key,
                    "label": label,
                    "detail": detail,
                    "n_donors": int(result["n_donors"]),
                    "delta_nll": float(result["mean_nll_left_minus_right"]),
                    "ci_low": float(lower),
                    "ci_high": float(upper),
                    "fraction_perturbation_worse": float(result["fraction_left_worse"]),
                }
            )
    if integrated_validation:
        first_condition = next(iter(payload["conditions"].values()))
        summary = {
            "selected_cells": int(payload["provenance"]["n_validation_cells"]),
            "nll": float(first_condition["full_nll_gene_cell_weighted"]),
            "accuracy": None,
            "balanced_recall": None,
            "n_donors": int(payload["provenance"]["n_validation_donors"]),
            "checkpoint_sha256": payload["provenance"].get("checkpoint_sha256"),
            "config_sha256": payload["provenance"].get("config_sha256"),
        }
    else:
        overall = payload["overall"]["full"][split]
        summary = {
            "selected_cells": int(payload["provenance"]["selected_cells"]),
            "nll": float(overall["nll"]),
            "accuracy": float(overall["accuracy"]),
            "balanced_recall": float(overall["balanced_recall"]),
            "n_donors": int(rows[0]["n_donors"]),
        }
    return pd.DataFrame(rows), summary


def panel_limits(data: pd.DataFrame) -> tuple[float, float]:
    values = data[["ci_low", "ci_high"]].to_numpy(float) * 1000.0
    low = float(np.nanmin(values))
    high = float(np.nanmax(values))
    span = max(high - low, 0.25)
    pad = 0.16 * span
    return min(low - pad, -0.05), max(high + pad, 0.05)


def draw_panel(axis: plt.Axes, data: pd.DataFrame, title: str, subtitle: str, background: str) -> None:
    axis.set_facecolor(background)
    rows = data[["order", "label", "detail"]].drop_duplicates().sort_values("order")
    y_base = np.arange(len(rows), dtype=float)
    lookup = {int(row.order): index for index, row in enumerate(rows.itertuples(index=False))}
    offsets = {"val": -0.11, "test": 0.11}
    for split, color, filled in (
        ("val", COLORS["blue"], False),
        ("test", COLORS["navy"], True),
    ):
        subset = data.loc[data["split"].eq(split)].sort_values("order")
        y = np.asarray([lookup[int(value)] for value in subset["order"]], dtype=float) + offsets[split]
        x = subset["delta_nll"].to_numpy(float) * 1000.0
        low = subset["ci_low"].to_numpy(float) * 1000.0
        high = subset["ci_high"].to_numpy(float) * 1000.0
        axis.errorbar(
            x,
            y,
            xerr=np.vstack([x - low, high - x]),
            fmt="o",
            markersize=8.5,
            markerfacecolor=color if filled else "white",
            markeredgecolor=color,
            markeredgewidth=1.8,
            ecolor=color,
            elinewidth=2.2,
            capsize=4.5,
            capthick=1.6,
            zorder=4,
        )
        for x_value, y_value in zip(x, y, strict=True):
            axis.text(
                x_value,
                y_value - 0.18 if split == "validation" else y_value + 0.18,
                f"{x_value:+.2f}",
                ha="center",
                va="center",
                fontsize=9.2,
                color=color,
                fontweight=700,
            )

    x_min, x_max = panel_limits(data)
    axis.set_xlim(x_min, x_max)
    axis.set_ylim(len(rows) - 0.55, -0.55)
    axis.axvline(0.0, color="#7B8797", linewidth=1.4, zorder=1)
    axis.grid(axis="x", color=COLORS["grid"], linewidth=0.9, zorder=0)
    axis.set_yticks(y_base)
    axis.set_yticklabels([row.label for row in rows.itertuples(index=False)], fontsize=12.0, fontweight=700, color=COLORS["dark"])
    axis.tick_params(axis="y", length=0, pad=12)
    axis.tick_params(axis="x", length=0, labelsize=10.5, colors=COLORS["muted"], pad=7)
    axis.set_xlabel("구성요소 교란 후 ΔNLL × 1,000  (양수 = 복원 악화)", fontsize=11.5, color=COLORS["muted"], labelpad=12)
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_color("#AAB5C2")
    axis.set_title(title, loc="left", fontsize=17.5, color=COLORS["dark"], fontweight=800, pad=36)
    axis.text(0.0, 1.055, subtitle, transform=axis.transAxes, ha="left", va="bottom", fontsize=10.8, color=COLORS["muted"])


def draw_figure(data: pd.DataFrame, summaries: dict[str, dict], output: Path) -> dict:
    configure_style()
    figure = plt.figure(figsize=(16, 9), facecolor="white")
    grid = figure.add_gridspec(1, 2, left=0.165, right=0.970, bottom=0.225, top=0.700, wspace=0.30)
    axis_personal = figure.add_subplot(grid[0, 0])
    axis_response = figure.add_subplot(grid[0, 1])
    draw_panel(
        axis_personal,
        data.loc[data["panel"].eq("personal")],
        "A. Personal 기본 특성",
        "병리와 무관한 개인 code가 donor 고유 정보를 더하는가?",
        COLORS["panel_a"],
    )
    draw_panel(
        axis_response,
        data.loc[data["panel"].eq("response")],
        "B. Personal×병리반응",
        "같은 병리에서도 donor마다 다른 반응을 설명하는가?",
        COLORS["panel_b"],
    )

    figure.text(
        0.045,
        0.955,
        "PRISM의 개인 성분을 제거·셔플해 validation과 locked test에서 직접 검증한다",
        ha="left",
        va="top",
        fontsize=24.5,
        fontweight=800,
        color=COLORS["navy"],
    )
    figure.text(
        0.047,
        0.904,
        "Frozen personal-rank2 · epoch 20  |  donor-bootstrap 95% CI  |  같은 네 가지 intervention을 두 split에 적용",
        ha="left",
        va="top",
        fontsize=13.0,
        color=COLORS["muted"],
    )
    figure.text(
        0.955,
        0.910,
        f"FULL NLL  validation = {summaries['validation']['nll']:.4f}  ·  test = {summaries['test']['nll']:.4f}",
        ha="right",
        va="center",
        fontsize=10.2,
        color=COLORS["navy"],
        fontweight=800,
        bbox={
            "boxstyle": "round,pad=0.45,rounding_size=0.2",
            "facecolor": "#EDF3FB",
            "edgecolor": "#B9CBE2",
            "linewidth": 1.0,
        },
    )
    legend = [
        Line2D([0], [0], marker="o", color=COLORS["blue"], markerfacecolor="white", markeredgewidth=1.8, markersize=8, linewidth=2, label=f"Validation · {summaries['validation']['n_donors']} donors"),
        Line2D([0], [0], marker="o", color=COLORS["navy"], markerfacecolor=COLORS["navy"], markeredgewidth=1.8, markersize=8, linewidth=2, label=f"Locked test · {summaries['test']['n_donors']} donors"),
    ]
    figure.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 0.825), ncol=2, frameon=False, fontsize=11.4, handletextpad=0.7, columnspacing=2.2)

    test_rows = data.loc[data["split"].eq("test")]
    supported = int(np.sum(test_rows["ci_low"].to_numpy(float) > 0))
    total = int(len(test_rows))
    banner = FancyBboxPatch(
        (0.040, 0.068),
        0.920,
        0.057,
        boxstyle="round,pad=0.006,rounding_size=0.008",
        transform=figure.transFigure,
        facecolor=COLORS["navy"],
        edgecolor=COLORS["navy"],
        linewidth=0.0,
        zorder=2,
    )
    figure.add_artist(banner)
    figure.text(
        0.500,
        0.096,
        "핵심  |  개인 기본 code는 locked test에서 지지됐지만, Personal×병리반응의 추가 설명력은 아직 확인되지 않았다.",
        ha="center",
        va="center",
        fontsize=12.3,
        color="white",
        fontweight=700,
        zorder=3,
    )
    figure.text(
        0.047,
        0.030,
        "점 = donor별 mean-cell NLL 차이의 평균, 선 = donor-bootstrap 95% CI. 양수이면 구성요소 제거·셔플 후 복원이 악화된다. "
        "Matched shuffle은 공변량과 target-excluded support 구성을 가능한 한 맞춘 다른 donor code를 사용한다. Test는 재조정에 사용하지 않았다.",
        ha="left",
        va="center",
        fontsize=9.0,
        color=COLORS["muted"],
    )

    output.mkdir(parents=True, exist_ok=True)
    base = output / "figure_03_personal_ablation_validation_test"
    figure.savefig(base.with_suffix(".png"), dpi=220, facecolor="white")
    figure.savefig(base.with_suffix(".pdf"), facecolor="white")
    figure.savefig(base.with_suffix(".svg"), facecolor="white")
    plt.close(figure)
    return {"test_supported_interventions": supported, "n_interventions": total}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation", type=Path, default=DEFAULT_VALIDATION)
    parser.add_argument("--test", type=Path, default=DEFAULT_TEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    validation_path = args.validation.resolve()
    test_path = args.test.resolve()
    output = args.output.resolve()
    validation, validation_summary = extract_split(validation_path, "val")
    test, test_summary = extract_split(test_path, "test")
    data = pd.concat([validation, test], ignore_index=True)
    summaries = {"validation": validation_summary, "test": test_summary}
    figure_summary = draw_figure(data, summaries, output)

    table_path = output / "figure_03_personal_ablation_values.tsv"
    data.to_csv(table_path, sep="\t", index=False)
    manifest = {
        "schema_version": "prism.presentation.figure03.v1",
        "model": "PRISM personal-rank2 epoch 20",
        "analysis": "personal and personal-by-pathology branch ablations",
        "validation_source": str(validation_path),
        "validation_source_sha256": sha256(validation_path),
        "test_source": str(test_path),
        "test_source_sha256": sha256(test_path),
        "summaries": summaries,
        **figure_summary,
        "outputs": {
            "png": "figure_03_personal_ablation_validation_test.png",
            "pdf": "figure_03_personal_ablation_validation_test.pdf",
            "svg": "figure_03_personal_ablation_validation_test.svg",
            "data": table_path.name,
        },
    }
    (output / "figure_03_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
