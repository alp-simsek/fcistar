# Pipeline

## Overview
Backend pipeline that pulls live data, runs the FCI* estimation, and writes structured outputs consumed by the frontend.

This Python code is based on the original MATLAB replication code in `backend/refs/matlab/`.

## Folder structure
- `backend/pipeline/` — main Python scripts for data ingestion, estimation, and output writing.
- `backend/pipeline/_aux_fun/` — helper functions used by `estimate_CCS_struct.py`.
- `backend/data/raw/` — raw input data used by the pipeline; not intended for hand editing.
- `backend/data/raw/assembled/` — assembled intermediate input files used by the estimation script.
- `backend/data/output/` — estimation outputs consumed by the frontend.
- `backend/refs/matlab/` — reference MATLAB code on which this Python code is based.

## Main scripts
- `assemble_data.py` — downloads and assembles the raw input data used by the estimation pipeline.
- `estimate_CCS_struct.py` — runs the estimation and writes the frontend-facing outputs to `backend/data/output/`.
