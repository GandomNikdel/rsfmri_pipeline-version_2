# rsfmri_pipeline_version_2

Resting-state fMRI preprocessing pipeline developed for reproducible neuroimaging analysis using FSL.

## Overview

This repository contains Bash scripts for preprocessing resting-state fMRI data.

Current implemented modules:

- Motion correction (MCFLIRT)
- Framewise Displacement (FD) quality control

Additional preprocessing modules will be added as the pipeline develops.

---

## Software

- FSL
- Bash
- Python 3

---

## Pipeline

Current workflow:

1. Motion correction (MCFLIRT)
2. Framewise Displacement QC

---

## Directory structure

```
scripts/
    01_mcflirt_v2.sh
    02_fd_qc_v2.sh
```

---

## Version

Current version:

Version 2

---

## Status

This pipeline is under active development.
