# HBN GAD rs-fMRI Feasibility Workflow

This repository contains early scripts for a preliminary feasibility and quality-control workflow for a planned resting-state fMRI study using the Healthy Brain Network (HBN) dataset.

The project focuses on amygdala–vmPFC functional connectivity in youth with Generalized Anxiety Disorder (GAD) compared with healthy controls.

At the current stage, we are working on the initial motion and QC scripts, including:

- `01_mcflirt_v2.sh` — MCFLIRT-based motion correction
- `02_fd_power_qc_v1.py` — Power FD calculation and motion QC summary
- `03_dvars_qc_v1.py` — DVARS-based QC summary (03_dvars_qc_v1.py computes DVARS using FSL fsl_motion_outliers and applies a subject-specific boxplot-based DVARS review flag. Final exclusion decisions are not made by DVARS alone and should be based on the combined FD + DVARS + usable-data QC stage.)

This repository is currently under development and documents the feasibility/QC stage before proceeding to final preprocessing, denoising, ROI extraction, connectivity estimation, and statistical analysis.

Raw imaging data and restricted phenotypic data are not included.
