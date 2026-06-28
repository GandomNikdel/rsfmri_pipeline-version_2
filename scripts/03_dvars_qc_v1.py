cd ~/wrkdir/HBNFINAL                                                                                  cat > 03_dvars_qc_v1.py <<'PY'
#!/usr/bin/env python3
"""
Stage 2B: DVARS QC v1

Purpose:
  Compute preliminary DVARS-based artifact QC from MCFLIRT-corrected
  resting-state fMRI data.

Input:
  sub-*/derivatives/mcflirt_v2/func_mc_v2.nii.gz

Outputs:
  Per subject:
    sub-*/derivatives/dvars_qc_v1/

  Dataset level:
    logs/dvars_qc_v1/

Method:
  FSL fsl_motion_outliers --dvars --nomoco

Important:
  This stage does NOT make final exclusion decisions.
  Final decisions will be made in the next stage by combining:
    FD + DVARS + usable-data criteria.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Dict, Any, Tuple


PIPELINE_NAME = "dvars_qc_v1"
PIPELINE_VERSION = "1.0.0"

BASE = Path.home() / "wrkdir" / "HBNFINAL"

MCFLIRT_REL_PATH = Path("derivatives") / "mcflirt_v2" / "func_mc_v2.nii.gz"

SUBJECT_DERIV_DIRNAME = "dvars_qc_v1"
LOG_DIR = BASE / "logs" / "dvars_qc_v1"

BOXPLOT_IQR_MULTIPLIER = 1.5
FIXED_DVARS_THRESHOLD = 5.0
DVARS_OUTLIER_PERCENT_REVIEW = 20.0


def now_stamp() -> str:
    return dt.datetime.now().strftime("%Y%m%d_%H%M%S")


RUN_STAMP = now_stamp()


def safe_mkdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def run_command(cmd: List[str], log_path: Path | None = None) -> Tuple[int, str]:
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    output = proc.stdout or ""

    if log_path is not None:
        with open(log_path, "a", encoding="utf-8") as f:
            f.write("\n$ " + " ".join(cmd) + "\n")
            f.write(output)
            f.write(f"\n[return_code={proc.returncode}]\n")

    return proc.returncode, output


def command_exists(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def list_subjects(base: Path) -> List[Path]:
    return sorted([p for p in base.glob("sub-*") if p.is_dir()])


def fslnvols(nifti_path: Path, log_path: Path) -> int | None:
    rc, out = run_command(["fslnvols", str(nifti_path)], log_path)
    if rc != 0:
        return None
    try:
        return int(out.strip().split()[0])
    except Exception:
        return None


def read_metric_values(path: Path) -> List[float]:
    if not path.exists():
        return []

    text = path.read_text(encoding="utf-8", errors="replace").strip()
    if not text:
        return []

    vals: List[float] = []
    for tok in text.replace(",", " ").split():
        try:
            vals.append(float(tok))
        except ValueError:
            pass

    return vals


def percentile(values: List[float], p: float) -> float:
    if not values:
        return float("nan")

    xs = sorted(values)

    if len(xs) == 1:
        return xs[0]

    pos = (len(xs) - 1) * (p / 100.0)
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))

    if lo == hi:
        return xs[lo]

    frac = pos - lo
    return xs[lo] * (1.0 - frac) + xs[hi] * frac


def mean(values: List[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def median(values: List[float]) -> float:
    return percentile(values, 50.0)


def fmt(x: Any) -> str:
    if isinstance(x, float):
        if math.isnan(x):
            return "NA"
        return f"{x:.6f}"
    if x is None:
        return "NA"
    return str(x)


def align_dvars_to_volumes(values: List[float], nvol: int) -> Tuple[List[float], str]:
    """
    fsl_motion_outliers DVARS output may contain either:
      - nvol values
      - nvol-1 values

    If nvol-1 values are returned, we pad the first volume with 0.0 because
    DVARS at the first volume has no previous volume.
    """

    if len(values) == nvol:
        return values, "n_values_equals_nvol"

    if len(values) == nvol - 1:
        return [0.0] + values, "n_values_equals_nvol_minus_1_padded_first_volume_zero"

    return values, f"UNEXPECTED_LENGTH_metric_values_{len(values)}_nvol_{nvol}"


def write_vector(path: Path, values: List[Any]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for v in values:
            f.write(f"{v}\n")


def write_timeseries_tsv(
    path: Path,
    dvars_aligned: List[float],
    box_mask: List[int],
    fixed_mask: List[int],
) -> None:
    fixed_name = str(FIXED_DVARS_THRESHOLD).replace(".", "p")

    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "volume_index_0based",
            "dvars_fsl_scaled",
            "dvars_outlier_boxplot",
            f"dvars_gt_{fixed_name}",
        ])

        for i, val in enumerate(dvars_aligned):
            w.writerow([i, fmt(val), box_mask[i], fixed_mask[i]])


def write_spike_regressors(path: Path, mask: List[int], prefix: str) -> None:
    outlier_indices = [i for i, m in enumerate(mask) if int(m) == 1]

    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")

        if not outlier_indices:
            w.writerow([f"{prefix}_none"])
            for _ in mask:
                w.writerow([0])
            return

        header = [f"{prefix}_{idx:04d}" for idx in outlier_indices]
        w.writerow(header)

        for i in range(len(mask)):
            row = [1 if i == idx else 0 for idx in outlier_indices]
            w.writerow(row)


def process_subject(sub_dir: Path) -> Dict[str, Any]:
    subject = sub_dir.name

    out_dir = sub_dir / "derivatives" / SUBJECT_DERIV_DIRNAME
    safe_mkdir(out_dir)

    subj_log = out_dir / "dvars_qc_subject.log"
    subj_log.write_text(
        f"Subject: {subject}\n"
        f"Pipeline: {PIPELINE_NAME} {PIPELINE_VERSION}\n"
        f"Run stamp: {RUN_STAMP}\n"
        f"Date: {dt.datetime.now().isoformat()}\n",
        encoding="utf-8",
    )

    func_mc = sub_dir / MCFLIRT_REL_PATH

    row: Dict[str, Any] = {
        "subject": subject,
        "qc_status_preliminary": "UNKNOWN",
        "input_func_mc": str(func_mc),
        "nvol": "NA",
        "n_dvars_metric_values_raw": "NA",
        "dvars_alignment": "NA",
        "dvars_mean": "NA",
        "dvars_median": "NA",
        "dvars_p75": "NA",
        "dvars_p95": "NA",
        "dvars_max": "NA",
        "dvars_boxplot_threshold": "NA",
        "dvars_boxplot_outlier_count": "NA",
        "dvars_boxplot_outlier_pct": "NA",
        "dvars_fixed_threshold": FIXED_DVARS_THRESHOLD,
        "dvars_fixed_outlier_count": "NA",
        "dvars_fixed_outlier_pct": "NA",
        "flag_dvars_boxplot_pct_gt_20": "NA",
        "flag_dvars_fixed_pct_gt_20": "NA",
        "fsl_return_code": "NA",
        "notes": "",
    }

    if not func_mc.exists():
        row["qc_status_preliminary"] = "MISSING_MCFLIRT_OUTPUT"
        row["notes"] = "MCFLIRT-corrected func_mc_v2.nii.gz not found."
        return row

    nvol = fslnvols(func_mc, subj_log)

    if nvol is None or nvol <= 1:
        row["qc_status_preliminary"] = "FAILED_NVol"
        row["notes"] = "Could not determine valid number of volumes."
        return row

    row["nvol"] = nvol

    fsl_metric = out_dir / "dvars_fsl_metric_values.txt"
    fsl_plot = out_dir / "dvars_fsl_metric_plot.png"
    fsl_confounds = out_dir / "dvars_fsl_outlier_confounds.tsv"

    timeseries_tsv = out_dir / "dvars_qc_timeseries.tsv"
    box_mask_path = out_dir / "dvars_outlier_mask_boxplot.txt"

    fixed_name = str(FIXED_DVARS_THRESHOLD).replace(".", "p")
    fixed_mask_path = out_dir / f"dvars_outlier_mask_gt_{fixed_name}.txt"

    box_spikes_path = out_dir / "dvars_spike_regressors_boxplot.tsv"
    fixed_spikes_path = out_dir / f"dvars_spike_regressors_gt_{fixed_name}.tsv"

    summary_path = out_dir / "dvars_qc_subject_summary.tsv"
    metadata_path = out_dir / "dvars_qc_metadata.json"

    cmd = [
        "fsl_motion_outliers",
        "-i", str(func_mc),
        "-o", str(fsl_confounds),
        "--dvars",
        "--nomoco",
        "-s", str(fsl_metric),
        "-p", str(fsl_plot),
        "-v",
    ]

    rc, _ = run_command(cmd, subj_log)
    row["fsl_return_code"] = rc

    if rc != 0:
        row["qc_status_preliminary"] = "FAILED_FSL_MOTION_OUTLIERS"
        row["notes"] = "fsl_motion_outliers failed; inspect subject log."
        return row

    raw_vals = read_metric_values(fsl_metric)
    row["n_dvars_metric_values_raw"] = len(raw_vals)

    if not raw_vals:
        row["qc_status_preliminary"] = "FAILED_EMPTY_DVARS"
        row["notes"] = "No DVARS metric values were written."
        return row

    aligned_vals, alignment = align_dvars_to_volumes(raw_vals, nvol)
    row["dvars_alignment"] = alignment

    if len(aligned_vals) != nvol:
        row["qc_status_preliminary"] = "FAILED_DVARS_LENGTH_MISMATCH"
        row["notes"] = f"DVARS length mismatch after alignment: len={len(aligned_vals)}, nvol={nvol}"
        return row

    q25 = percentile(raw_vals, 25.0)
    q75 = percentile(raw_vals, 75.0)
    iqr = q75 - q25
    box_threshold = q75 + BOXPLOT_IQR_MULTIPLIER * iqr

    box_mask = [1 if v > box_threshold else 0 for v in aligned_vals]
    fixed_mask = [1 if v > FIXED_DVARS_THRESHOLD else 0 for v in aligned_vals]

    box_count = sum(box_mask)
    fixed_count = sum(fixed_mask)

    box_pct = 100.0 * box_count / nvol
    fixed_pct = 100.0 * fixed_count / nvol

    flag_box_pct = int(box_pct > DVARS_OUTLIER_PERCENT_REVIEW)
    flag_fixed_pct = int(fixed_pct > DVARS_OUTLIER_PERCENT_REVIEW)

    row.update({
        "dvars_mean": mean(raw_vals),
        "dvars_median": median(raw_vals),
        "dvars_p75": q75,
        "dvars_p95": percentile(raw_vals, 95.0),
        "dvars_max": max(raw_vals),
        "dvars_boxplot_threshold": box_threshold,
        "dvars_boxplot_outlier_count": box_count,
        "dvars_boxplot_outlier_pct": box_pct,
        "dvars_fixed_threshold": FIXED_DVARS_THRESHOLD,
        "dvars_fixed_outlier_count": fixed_count,
        "dvars_fixed_outlier_pct": fixed_pct,
        "flag_dvars_boxplot_pct_gt_20": flag_box_pct,
        "flag_dvars_fixed_pct_gt_20": flag_fixed_pct,
    })

    if flag_box_pct or flag_fixed_pct:
        row["qc_status_preliminary"] = "HIGH_DVARS_REVIEW"
    else:
        row["qc_status_preliminary"] = "PASS_PRELIMINARY_DVARS"

    row["notes"] = (
        "Preliminary DVARS QC only; final exclusion requires combined FD + DVARS + usable-data decision."
    )

    write_vector(box_mask_path, box_mask)
    write_vector(fixed_mask_path, fixed_mask)
    write_timeseries_tsv(timeseries_tsv, aligned_vals, box_mask, fixed_mask)
    write_spike_regressors(box_spikes_path, box_mask, "dvars_boxplot_spike")
    write_spike_regressors(fixed_spikes_path, fixed_mask, "dvars_gt5_spike")

    metadata = {
        "pipeline_name": PIPELINE_NAME,
        "pipeline_version": PIPELINE_VERSION,
        "run_stamp": RUN_STAMP,
        "subject": subject,
        "input_func_mc": str(func_mc),
        "method": "FSL fsl_motion_outliers --dvars --nomoco",
        "nvol": nvol,
        "dvars_alignment": alignment,
        "boxplot_threshold_rule": "P75 + 1.5*IQR",
        "fixed_dvars_threshold_auxiliary": FIXED_DVARS_THRESHOLD,
        "preliminary_review_rule": (
            f"HIGH_DVARS_REVIEW if boxplot outlier percent > {DVARS_OUTLIER_PERCENT_REVIEW}% "
            f"or fixed-threshold outlier percent > {DVARS_OUTLIER_PERCENT_REVIEW}%"
        ),
        "important_note": "This stage creates preliminary DVARS review flags only, not final exclusions.",
    }

    with open(metadata_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    with open(summary_path, "w", encoding="utf-8", newline="") as f:
        fields = list(row.keys())
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerow({k: fmt(v) for k, v in row.items()})

    return row


def write_dataset_outputs(rows: List[Dict[str, Any]]) -> None:
    safe_mkdir(LOG_DIR)

    summary_file = LOG_DIR / f"dvars_qc_v1_dataset_summary_{RUN_STAMP}.tsv"
    review_file = LOG_DIR / f"dvars_qc_v1_review_candidates_{RUN_STAMP}.tsv"
    config_file = LOG_DIR / f"dvars_qc_v1_config_{RUN_STAMP}.json"
    master_log = LOG_DIR / f"dvars_qc_v1_master_{RUN_STAMP}.log"
    stage_summary = LOG_DIR / f"dvars_qc_v1_stage_summary_{RUN_STAMP}.txt"

    fields = [
        "subject",
        "qc_status_preliminary",
        "input_func_mc",
        "nvol",
        "n_dvars_metric_values_raw",
        "dvars_alignment",
        "dvars_mean",
        "dvars_median",
        "dvars_p75",
        "dvars_p95",
        "dvars_max",
        "dvars_boxplot_threshold",
        "dvars_boxplot_outlier_count",
        "dvars_boxplot_outlier_pct",
        "dvars_fixed_threshold",
        "dvars_fixed_outlier_count",
        "dvars_fixed_outlier_pct",
        "flag_dvars_boxplot_pct_gt_20",
        "flag_dvars_fixed_pct_gt_20",
        "fsl_return_code",
        "notes",
    ]

    with open(summary_file, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: fmt(r.get(k, "NA")) for k in fields})

    review_rows = [
        r for r in rows
        if r.get("qc_status_preliminary") != "PASS_PRELIMINARY_DVARS"
    ]

    with open(review_file, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in review_rows:
            w.writerow({k: fmt(r.get(k, "NA")) for k in fields})

    status_counts: Dict[str, int] = {}
    for r in rows:
        status = str(r.get("qc_status_preliminary", "UNKNOWN"))
        status_counts[status] = status_counts.get(status, 0) + 1

    def numeric_list(key: str) -> List[float]:
        out = []
        for r in rows:
            try:
                val = float(r.get(key, "nan"))
                if not math.isnan(val):
                    out.append(val)
            except Exception:
                pass
        return out

    def write_dist(f, label: str, vals: List[float]) -> None:
        f.write(f"{label}:\n")
        if not vals:
            f.write("- N: 0\n\n")
            return
        f.write(f"- N: {len(vals)}\n")
        f.write(f"- min: {min(vals):.6f}\n")
        f.write(f"- median: {median(vals):.6f}\n")
        f.write(f"- mean: {mean(vals):.6f}\n")
        f.write(f"- p75: {percentile(vals, 75):.6f}\n")
        f.write(f"- p90: {percentile(vals, 90):.6f}\n")
        f.write(f"- max: {max(vals):.6f}\n\n")

    config = {
        "pipeline_name": PIPELINE_NAME,
        "pipeline_version": PIPELINE_VERSION,
        "run_stamp": RUN_STAMP,
        "base": str(BASE),
        "input": str(MCFLIRT_REL_PATH),
        "method": "FSL fsl_motion_outliers --dvars --nomoco",
        "boxplot_iqr_multiplier": BOXPLOT_IQR_MULTIPLIER,
        "fixed_dvars_threshold_auxiliary": FIXED_DVARS_THRESHOLD,
        "dvars_outlier_percent_review": DVARS_OUTLIER_PERCENT_REVIEW,
        "python": sys.version,
        "platform": platform.platform(),
        "commands": {
            "fsl_motion_outliers": shutil.which("fsl_motion_outliers"),
            "fslnvols": shutil.which("fslnvols"),
        },
        "important_note": "This is preliminary DVARS QC, not final exclusion.",
    }

    with open(config_file, "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    box_pct_values = numeric_list("dvars_boxplot_outlier_pct")
    fixed_pct_values = numeric_list("dvars_fixed_outlier_pct")
    mean_values = numeric_list("dvars_mean")
    max_values = numeric_list("dvars_max")

    top_box = sorted(
        [r for r in rows if str(r.get("dvars_boxplot_outlier_pct", "NA")) != "NA"],
        key=lambda r: float(r.get("dvars_boxplot_outlier_pct", -1)),
        reverse=True,
    )[:30]

    with open(stage_summary, "w", encoding="utf-8") as f:
        f.write("Stage: DVARS QC\n")
        f.write(f"Stage ID: {PIPELINE_NAME}\n")
        f.write(f"Version: {PIPELINE_VERSION}\n")
        f.write(f"Run stamp: {RUN_STAMP}\n")
        f.write(f"Base: {BASE}\n\n")

        f.write("Purpose:\n")
        f.write("Compute preliminary DVARS-based artifact QC from MCFLIRT-corrected resting-state fMRI data.\n\n")

        f.write("Method:\n")
        f.write("FSL fsl_motion_outliers --dvars --nomoco\n")
        f.write("Boxplot threshold: P75 + 1.5*IQR\n")
        f.write(f"Auxiliary fixed DVARS threshold: {FIXED_DVARS_THRESHOLD}\n")
        f.write("No final exclusions are made in this stage.\n\n")

        f.write("Processing status counts:\n")
        for k in sorted(status_counts):
            f.write(f"- {k}: {status_counts[k]}\n")
        f.write("\n")

        write_dist(f, "DVARS mean distribution", mean_values)
        write_dist(f, "DVARS max distribution", max_values)
        write_dist(f, "DVARS boxplot outlier percentage distribution", box_pct_values)
        write_dist(f, "DVARS fixed-threshold outlier percentage distribution", fixed_pct_values)

        f.write("Top 30 subjects by DVARS boxplot outlier percentage:\n")
        f.write("subject,qc_status,dvars_mean,dvars_max,dvars_boxplot_outlier_pct,dvars_fixed_outlier_pct\n")
        for r in top_box:
            f.write(",".join([
                str(r.get("subject", "NA")),
                str(r.get("qc_status_preliminary", "NA")),
                fmt(r.get("dvars_mean", "NA")),
                fmt(r.get("dvars_max", "NA")),
                fmt(r.get("dvars_boxplot_outlier_pct", "NA")),
                fmt(r.get("dvars_fixed_outlier_pct", "NA")),
            ]) + "\n")

        f.write("\nKey files:\n")
        f.write(f"- {summary_file}\n")
        f.write(f"- {review_file}\n")
        f.write(f"- {config_file}\n")
        f.write(f"- {master_log}\n")

    with open(master_log, "w", encoding="utf-8") as f:
        f.write(f"Pipeline: {PIPELINE_NAME} {PIPELINE_VERSION}\n")
        f.write(f"Run stamp: {RUN_STAMP}\n")
        f.write(f"Base: {BASE}\n")
        f.write(f"Subjects processed/listed: {len(rows)}\n")
        f.write("Status counts:\n")
        for k in sorted(status_counts):
            f.write(f"{k}\t{status_counts[k]}\n")
        f.write("\nDataset-level outputs:\n")
        f.write(f"summary\t{summary_file}\n")
        f.write(f"review_candidates\t{review_file}\n")
        f.write(f"config\t{config_file}\n")
        f.write(f"stage_summary\t{stage_summary}\n")

    print("\nDONE: DVARS QC v1")
    print(f"Dataset summary: {summary_file}")
    print(f"Review candidates: {review_file}")
    print(f"Config: {config_file}")
    print(f"Master log: {master_log}")
    print(f"Stage summary: {stage_summary}")


def main() -> int:
    if not BASE.exists():
        print(f"ERROR: BASE does not exist: {BASE}", file=sys.stderr)
        return 1

    safe_mkdir(LOG_DIR)

    missing_cmds = [
        cmd for cmd in ["fsl_motion_outliers", "fslnvols"]
        if not command_exists(cmd)
    ]

    if missing_cmds:
        print("ERROR: Missing required FSL commands:", ", ".join(missing_cmds), file=sys.stderr)
        print("Make sure FSL is loaded/available in your shell.", file=sys.stderr)
        return 1

    subjects = list_subjects(BASE)

    if not subjects:
        print(f"ERROR: no sub-* directories found in {BASE}", file=sys.stderr)
        return 1

    print(f"Pipeline: {PIPELINE_NAME} {PIPELINE_VERSION}")
    print(f"Base: {BASE}")
    print(f"Subjects found: {len(subjects)}")
    print(f"Run stamp: {RUN_STAMP}")
    print("Starting DVARS QC...")

    rows: List[Dict[str, Any]] = []

    for idx, sub_dir in enumerate(subjects, start=1):
        print(f"[{idx}/{len(subjects)}] {sub_dir.name}")
        rows.append(process_subject(sub_dir))

    write_dataset_outputs(rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
PY

chmod +x 03_dvars_qc_v1.py
