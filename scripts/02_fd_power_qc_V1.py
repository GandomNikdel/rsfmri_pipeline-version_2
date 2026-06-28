cd ~/wrkdir/HBNFINAL

cat > 02_fd_power_qc_v1.py <<'PY'
#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
import platform
import shutil
import statistics
import subprocess
import sys
from datetime import datetime
from pathlib import Path

PIPELINE_NAME = "fd_power_qc_v1"
PIPELINE_VERSION = "1.0.0"

BASE = Path.home() / "wrkdir" / "HBNFINAL"

# Power-style FD settings
RADIUS_MM = 50.0
FD_THRESHOLDS_MM = [0.2, 0.3, 0.5]
PRIMARY_FD_THRESHOLD_MM = 0.5

# Preliminary review flags only; not final exclusion.
MEAN_FD_REVIEW_MM = 0.30
FD05_PERCENT_REVIEW = 20.0
MIN_REMAINING_VOLUMES_REVIEW = 150

RUN_ID = datetime.now().strftime("%Y%m%d_%H%M%S")
LOGDIR = BASE / "logs" / PIPELINE_NAME
LOGDIR.mkdir(parents=True, exist_ok=True)

MASTERLOG = LOGDIR / f"{PIPELINE_NAME}_master_{RUN_ID}.log"
DATASET_SUMMARY = LOGDIR / f"{PIPELINE_NAME}_dataset_summary_{RUN_ID}.tsv"
REVIEW_CANDIDATES = LOGDIR / f"{PIPELINE_NAME}_review_candidates_{RUN_ID}.tsv"
CONFIG_JSON = LOGDIR / f"{PIPELINE_NAME}_config_{RUN_ID}.json"
STAGE_SUMMARY = LOGDIR / f"{PIPELINE_NAME}_stage_summary_{RUN_ID}.txt"


def log(msg: str) -> None:
    line = f"[{datetime.now().strftime('%F %T')}] {msg}"
    print(line)
    with open(MASTERLOG, "a") as f:
        f.write(line + "\n")


def run_cmd(cmd):
    try:
        return subprocess.check_output(cmd, stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return "NA"


def safe_int(x):
    try:
        return int(str(x))
    except Exception:
        return None


def safe_float(x):
    try:
        y = float(str(x))
        if math.isfinite(y):
            return y
        return None
    except Exception:
        return None


def lab(t: float) -> str:
    return str(t).replace(".", "p")


def fmt(x):
    if isinstance(x, float):
        return f"{x:.6f}"
    return x


def percentile(sorted_values, p):
    if not sorted_values:
        return "NA"
    if len(sorted_values) == 1:
        return sorted_values[0]
    k = (len(sorted_values) - 1) * p
    lo = math.floor(k)
    hi = math.ceil(k)
    if lo == hi:
        return sorted_values[int(k)]
    return sorted_values[lo] * (hi - k) + sorted_values[hi] * (k - lo)


def summary_stats(values):
    if not values:
        return {
            "mean": "NA",
            "median": "NA",
            "min": "NA",
            "max": "NA",
            "std": "NA",
            "p75": "NA",
            "p90": "NA",
            "p95": "NA",
        }
    x = sorted(values)
    return {
        "mean": statistics.mean(values),
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
        "std": statistics.pstdev(values) if len(values) > 1 else 0.0,
        "p75": percentile(x, 0.75),
        "p90": percentile(x, 0.90),
        "p95": percentile(x, 0.95),
    }


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        w.writeheader()
        for row in rows:
            clean = {}
            for k in fieldnames:
                v = row.get(k, "NA")
                if isinstance(v, float):
                    clean[k] = f"{v:.10g}"
                else:
                    clean[k] = v
            w.writerow(clean)


def write_vector(path, values):
    with open(path, "w") as f:
        for v in values:
            f.write(str(v) + "\n")


def read_mcflirt_par(par_file):
    rows = []
    with open(par_file, "r", errors="ignore") as f:
        for lineno, line in enumerate(f, 1):
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 6:
                raise ValueError(f"{par_file}: line {lineno} has fewer than 6 columns")
            vals = [float(x) for x in parts[:6]]
            if not all(math.isfinite(v) for v in vals):
                raise ValueError(f"{par_file}: line {lineno} contains NaN/Inf")
            rows.append(vals)
    return rows


def compute_power_fd(par_rows):
    out = []
    previous = None

    for i, vals in enumerate(par_rows):
        # MCFLIRT .par order:
        # columns 1-3 = rotations in radians
        # columns 4-6 = translations in mm
        rot_x, rot_y, rot_z = vals[0], vals[1], vals[2]
        trans_x, trans_y, trans_z = vals[3], vals[4], vals[5]

        if previous is None:
            drot_x = drot_y = drot_z = 0.0
            dtrans_x = dtrans_y = dtrans_z = 0.0
        else:
            drot_x = rot_x - previous[0]
            drot_y = rot_y - previous[1]
            drot_z = rot_z - previous[2]
            dtrans_x = trans_x - previous[3]
            dtrans_y = trans_y - previous[4]
            dtrans_z = trans_z - previous[5]

        trans_fd_mm = abs(dtrans_x) + abs(dtrans_y) + abs(dtrans_z)
        rot_fd_mm = RADIUS_MM * (abs(drot_x) + abs(drot_y) + abs(drot_z))
        fd_power_mm = trans_fd_mm + rot_fd_mm

        out.append({
            "volume_index": i,
            "rot_x_rad": rot_x,
            "rot_y_rad": rot_y,
            "rot_z_rad": rot_z,
            "trans_x_mm": trans_x,
            "trans_y_mm": trans_y,
            "trans_z_mm": trans_z,
            "drot_x_rad": drot_x,
            "drot_y_rad": drot_y,
            "drot_z_rad": drot_z,
            "dtrans_x_mm": dtrans_x,
            "dtrans_y_mm": dtrans_y,
            "dtrans_z_mm": dtrans_z,
            "trans_fd_mm": trans_fd_mm,
            "rot_fd_mm": rot_fd_mm,
            "fd_power_2012_mm": fd_power_mm,
        })

        previous = vals

    return out


CONFIG = {
    "pipeline_name": PIPELINE_NAME,
    "pipeline_version": PIPELINE_VERSION,
    "run_id": RUN_ID,
    "base_dir": str(BASE),
    "created": datetime.now().isoformat(),
    "fd_formula": "FD = sum(abs(diff(translations_mm))) + radius_mm * sum(abs(diff(rotations_rad)))",
    "radius_mm": RADIUS_MM,
    "mcflirt_par_column_order": {
        "columns_1_to_3": "rot_x_rad, rot_y_rad, rot_z_rad",
        "columns_4_to_6": "trans_x_mm, trans_y_mm, trans_z_mm",
    },
    "fd_thresholds_mm_saved": FD_THRESHOLDS_MM,
    "primary_fd_threshold_mm_for_preliminary_flags": PRIMARY_FD_THRESHOLD_MM,
    "mean_fd_review_mm": MEAN_FD_REVIEW_MM,
    "fd05_percent_review": FD05_PERCENT_REVIEW,
    "min_remaining_volumes_review": MIN_REMAINING_VOLUMES_REVIEW,
    "important_note": "This pipeline computes FD and preliminary review flags only. It does not make final motion-exclusion decisions.",
    "software": {
        "python": sys.version.replace("\n", " "),
        "platform": platform.platform(),
        "host": platform.node(),
        "fslversion": run_cmd(["fslversion"]),
    },
}

with open(CONFIG_JSON, "w") as f:
    json.dump(CONFIG, f, indent=2)

base_fields = [
    "subject",
    "status",
    "nvol_func",
    "nrow_par",
    "tr_sec",
    "scan_minutes",
    "mean_fd",
    "mean_fd_excluding_first",
    "median_fd",
    "min_fd",
    "max_fd",
    "std_fd",
    "p75_fd",
    "p90_fd",
    "p95_fd",
    "mean_trans_fd",
    "mean_rot_fd",
    "max_trans_fd",
    "max_rot_fd",
    "max_abs_translation_x",
    "max_abs_translation_y",
    "max_abs_translation_z",
    "max_abs_rotation_x_rad",
    "max_abs_rotation_y_rad",
    "max_abs_rotation_z_rad",
    "max_abs_rotation_x_deg",
    "max_abs_rotation_y_deg",
    "max_abs_rotation_z_deg",
]

threshold_fields = []
for t in FD_THRESHOLDS_MM:
    L = lab(t)
    threshold_fields += [
        f"fd_gt_{L}_n",
        f"fd_gt_{L}_pct",
        f"remaining_volumes_fd_le_{L}",
        f"remaining_minutes_fd_le_{L}",
    ]

flag_fields = [
    "flag_mean_fd_gt_0p30",
    "flag_fd05_pct_gt_20",
    "flag_fd05_remaining_lt_150",
    "qc_status_preliminary",
    "notes",
]

dataset_fields = base_fields + threshold_fields + flag_fields

fd_timeseries_fields = [
    "volume_index",
    "rot_x_rad",
    "rot_y_rad",
    "rot_z_rad",
    "trans_x_mm",
    "trans_y_mm",
    "trans_z_mm",
    "drot_x_rad",
    "drot_y_rad",
    "drot_z_rad",
    "dtrans_x_mm",
    "dtrans_y_mm",
    "dtrans_z_mm",
    "trans_fd_mm",
    "rot_fd_mm",
    "fd_power_2012_mm",
]

subjects = sorted([p for p in BASE.glob("sub-*") if p.is_dir()])

log(f"Starting {PIPELINE_NAME} v{PIPELINE_VERSION}")
log(f"BASE: {BASE}")
log(f"Subjects found: {len(subjects)}")
log(f"Config JSON: {CONFIG_JSON}")
log(f"FSL version: {CONFIG['software']['fslversion']}")

all_rows = []

for sub in subjects:
    subject = sub.name
    mcdir = sub / "derivatives" / "mcflirt_v2"
    func = mcdir / "func_mc_v2.nii.gz"
    par = mcdir / "func_mc_v2.par"

    outdir = sub / "derivatives" / PIPELINE_NAME

    row = {k: "NA" for k in dataset_fields}
    row["subject"] = subject
    notes = []

    try:
        if not func.is_file():
            raise FileNotFoundError(f"missing MCFLIRT corrected NIfTI: {func}")
        if not par.is_file():
            raise FileNotFoundError(f"missing MCFLIRT .par: {par}")

        if outdir.exists():
            shutil.rmtree(outdir)
        outdir.mkdir(parents=True, exist_ok=True)

        nvol_func = run_cmd(["fslnvols", str(func)])
        tr_sec = run_cmd(["fslval", str(func), "pixdim4"])

        nvol_i = safe_int(nvol_func)
        tr_f = safe_float(tr_sec)

        par_rows = read_mcflirt_par(par)
        fd_rows = compute_power_fd(par_rows)

        if nvol_i is not None and nvol_i != len(par_rows):
            notes.append(f"nvol_func_{nvol_i}_neq_nrow_par_{len(par_rows)}")

        fd_values = [float(r["fd_power_2012_mm"]) for r in fd_rows]
        fd_excluding_first = fd_values[1:] if len(fd_values) > 1 else fd_values
        trans_fd_values = [float(r["trans_fd_mm"]) for r in fd_rows]
        rot_fd_values = [float(r["rot_fd_mm"]) for r in fd_rows]

        fd_stats = summary_stats(fd_values)
        fd_ex_stats = summary_stats(fd_excluding_first)
        trans_stats = summary_stats(trans_fd_values)
        rot_stats = summary_stats(rot_fd_values)

        scan_minutes = "NA"
        if nvol_i is not None and tr_f is not None:
            scan_minutes = nvol_i * tr_f / 60.0

        motion_params_file = outdir / "motion_parameters_mcflirt_ordered.tsv"
        fd_timeseries_file = outdir / "fd_power_2012_timeseries.tsv"
        fd_vector_file = outdir / "fd_power_2012_mm.txt"
        subject_summary_file = outdir / "fd_power_qc_subject_summary.tsv"
        subject_metadata_file = outdir / "fd_power_qc_metadata.json"

        write_tsv(
            motion_params_file,
            fd_rows,
            [
                "volume_index",
                "rot_x_rad",
                "rot_y_rad",
                "rot_z_rad",
                "trans_x_mm",
                "trans_y_mm",
                "trans_z_mm",
            ],
        )

        write_tsv(fd_timeseries_file, fd_rows, fd_timeseries_fields)
        write_vector(fd_vector_file, [f"{v:.10g}" for v in fd_values])

        for t in FD_THRESHOLDS_MM:
            L = lab(t)
            mask = [1 if v > t else 0 for v in fd_values]
            keep = [1 if v <= t else 0 for v in fd_values]

            write_vector(outdir / f"fd_outlier_mask_gt_{L}mm.txt", mask)
            write_vector(outdir / f"fd_keep_mask_le_{L}mm.txt", keep)

            outlier_indices = [i for i, m in enumerate(mask) if m == 1]
            spike_file = outdir / f"fd_spike_regressors_gt_{L}mm.tsv"

            with open(spike_file, "w", newline="") as f:
                if outlier_indices:
                    cols = ["volume_index"] + [
                        f"fd_gt_{L}mm_outlier_{idx:04d}"
                        for idx in outlier_indices
                    ]
                    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
                    w.writeheader()
                    for i in range(len(fd_values)):
                        r = {"volume_index": i}
                        for idx in outlier_indices:
                            r[f"fd_gt_{L}mm_outlier_{idx:04d}"] = 1 if i == idx else 0
                        w.writerow(r)
                else:
                    cols = ["volume_index", f"no_fd_gt_{L}mm_outliers"]
                    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
                    w.writeheader()
                    for i in range(len(fd_values)):
                        w.writerow({
                            "volume_index": i,
                            f"no_fd_gt_{L}mm_outliers": 0,
                        })

        max_abs_trans = [
            max(abs(r["trans_x_mm"]) for r in fd_rows),
            max(abs(r["trans_y_mm"]) for r in fd_rows),
            max(abs(r["trans_z_mm"]) for r in fd_rows),
        ]

        max_abs_rot = [
            max(abs(r["rot_x_rad"]) for r in fd_rows),
            max(abs(r["rot_y_rad"]) for r in fd_rows),
            max(abs(r["rot_z_rad"]) for r in fd_rows),
        ]

        row.update({
            "status": "OK",
            "nvol_func": nvol_func,
            "nrow_par": len(par_rows),
            "tr_sec": tr_sec,
            "scan_minutes": fmt(scan_minutes),
            "mean_fd": fmt(fd_stats["mean"]),
            "mean_fd_excluding_first": fmt(fd_ex_stats["mean"]),
            "median_fd": fmt(fd_stats["median"]),
            "min_fd": fmt(fd_stats["min"]),
            "max_fd": fmt(fd_stats["max"]),
            "std_fd": fmt(fd_stats["std"]),
            "p75_fd": fmt(fd_stats["p75"]),
            "p90_fd": fmt(fd_stats["p90"]),
            "p95_fd": fmt(fd_stats["p95"]),
            "mean_trans_fd": fmt(trans_stats["mean"]),
            "mean_rot_fd": fmt(rot_stats["mean"]),
            "max_trans_fd": fmt(trans_stats["max"]),
            "max_rot_fd": fmt(rot_stats["max"]),
            "max_abs_translation_x": fmt(max_abs_trans[0]),
            "max_abs_translation_y": fmt(max_abs_trans[1]),
            "max_abs_translation_z": fmt(max_abs_trans[2]),
            "max_abs_rotation_x_rad": fmt(max_abs_rot[0]),
            "max_abs_rotation_y_rad": fmt(max_abs_rot[1]),
            "max_abs_rotation_z_rad": fmt(max_abs_rot[2]),
            "max_abs_rotation_x_deg": fmt(max_abs_rot[0] * 57.295779513),
            "max_abs_rotation_y_deg": fmt(max_abs_rot[1] * 57.295779513),
            "max_abs_rotation_z_deg": fmt(max_abs_rot[2] * 57.295779513),
        })

        for t in FD_THRESHOLDS_MM:
            L = lab(t)
            n_bad = sum(v > t for v in fd_values)
            pct_bad = 100.0 * n_bad / len(fd_values)
            remaining = len(fd_values) - n_bad
            remaining_min = "NA"
            if tr_f is not None:
                remaining_min = remaining * tr_f / 60.0

            row[f"fd_gt_{L}_n"] = n_bad
            row[f"fd_gt_{L}_pct"] = f"{pct_bad:.6f}"
            row[f"remaining_volumes_fd_le_{L}"] = remaining
            row[f"remaining_minutes_fd_le_{L}"] = fmt(remaining_min)

        flag_mean = safe_float(row["mean_fd"]) is not None and float(row["mean_fd"]) > MEAN_FD_REVIEW_MM
        flag_pct = safe_float(row["fd_gt_0p5_pct"]) is not None and float(row["fd_gt_0p5_pct"]) > FD05_PERCENT_REVIEW
        rem05 = safe_int(row["remaining_volumes_fd_le_0p5"])
        flag_remaining = rem05 is not None and rem05 < MIN_REMAINING_VOLUMES_REVIEW

        row["flag_mean_fd_gt_0p30"] = int(flag_mean)
        row["flag_fd05_pct_gt_20"] = int(flag_pct)
        row["flag_fd05_remaining_lt_150"] = int(flag_remaining)

        if flag_remaining:
            row["qc_status_preliminary"] = "EXCLUDE_CANDIDATE_REVIEW_REQUIRED"
        elif flag_mean or flag_pct:
            row["qc_status_preliminary"] = "HIGH_MOTION_REVIEW"
        else:
            row["qc_status_preliminary"] = "PASS_PRELIMINARY_FD"

        row["notes"] = ";".join(notes)

        write_tsv(subject_summary_file, [row], dataset_fields)

        subject_metadata = {
            "subject": subject,
            "created": datetime.now().isoformat(),
            "pipeline_name": PIPELINE_NAME,
            "pipeline_version": PIPELINE_VERSION,
            "input_func": str(func),
            "input_par": str(par),
            "outputs": {
                "motion_parameters": str(motion_params_file),
                "fd_timeseries": str(fd_timeseries_file),
                "fd_vector": str(fd_vector_file),
                "subject_summary": str(subject_summary_file),
            },
            "config": CONFIG,
            "summary": row,
        }

        with open(subject_metadata_file, "w") as f:
            json.dump(subject_metadata, f, indent=2)

        log(
            f"{subject}: {row['qc_status_preliminary']} "
            f"meanFD={row['mean_fd']} "
            f"maxFD={row['max_fd']} "
            f"FD>0.5%={row['fd_gt_0p5_pct']}"
        )

    except Exception as e:
        outdir.mkdir(parents=True, exist_ok=True)
        row["status"] = "FAILED"
        row["qc_status_preliminary"] = "FAILED"
        row["notes"] = str(e)

        with open(outdir / "fd_power_qc_FAILED.txt", "w") as f:
            f.write(f"Subject: {subject}\n")
            f.write(f"Error: {e}\n")
            f.write(f"Time: {datetime.now().isoformat()}\n")

        log(f"{subject}: FAILED {e}")

    all_rows.append(row)

write_tsv(DATASET_SUMMARY, all_rows, dataset_fields)

review_rows = [
    r for r in all_rows
    if r["qc_status_preliminary"] in [
        "HIGH_MOTION_REVIEW",
        "EXCLUDE_CANDIDATE_REVIEW_REQUIRED",
        "FAILED",
    ]
]

write_tsv(REVIEW_CANDIDATES, review_rows, dataset_fields)

status_counts = {}
for r in all_rows:
    k = r["qc_status_preliminary"]
    status_counts[k] = status_counts.get(k, 0) + 1

with open(STAGE_SUMMARY, "w") as f:
    f.write(f"Stage: {PIPELINE_NAME}\n")
    f.write(f"Version: {PIPELINE_VERSION}\n")
    f.write(f"Run ID: {RUN_ID}\n")
    f.write(f"Date: {datetime.now()}\n")
    f.write(f"Base directory: {BASE}\n")
    f.write("\nPurpose:\n")
    f.write("Compute corrected Power FD from MCFLIRT motion parameters using the correct MCFLIRT .par column order.\n")
    f.write("This stage produces preliminary FD review flags only; it does not make final motion-exclusion decisions.\n")
    f.write("\nConfiguration:\n")
    f.write(f"- Radius for rotation-to-mm conversion: {RADIUS_MM} mm\n")
    f.write(f"- FD thresholds saved: {FD_THRESHOLDS_MM} mm\n")
    f.write(f"- Preliminary review flags: mean FD > {MEAN_FD_REVIEW_MM} mm OR FD>0.5mm percent > {FD05_PERCENT_REVIEW}% OR remaining FD<=0.5mm volumes < {MIN_REMAINING_VOLUMES_REVIEW}\n")
    f.write("\nQC status counts:\n")
    for k, v in sorted(status_counts.items()):
        f.write(f"- {k}: {v}\n")
    f.write("\nFiles:\n")
    f.write(f"- Dataset summary: {DATASET_SUMMARY}\n")
    f.write(f"- Review candidates: {REVIEW_CANDIDATES}\n")
    f.write(f"- Config JSON: {CONFIG_JSON}\n")
    f.write(f"- Master log: {MASTERLOG}\n")

log(f"Dataset summary: {DATASET_SUMMARY}")
log(f"Review candidates: {REVIEW_CANDIDATES}")
log(f"Stage summary: {STAGE_SUMMARY}")
log("PIPELINE COMPLETE")
PY

chmod +x 02_fd_power_qc_v1.py
python3 -m py_compile 02_fd_power_qc_v1.py
