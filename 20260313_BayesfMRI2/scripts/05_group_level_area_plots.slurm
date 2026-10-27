#!/usr/bin/env bash
# Plot the completed area-filtered NMD group analysis on Hazel.
# Arguments: fixed checkout, completed analysis directory, optional explicit plot output directory.
# Update the analysis directory's plots by default; do not create separate fix versions.
# Use afterok for the first plot job; a completed analysis needs no new fit dependency.
# Read saved models only; render through the existing private virtual display.
#SBATCH --job-name=nmd_area_plots
#SBATCH --account=statmind_cpu
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20261025_BayesfMRI_AreaFiltering/logs/plots-%j.out
#SBATCH --error=/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20261025_BayesfMRI_AreaFiltering/logs/plots-%j.err
set -euo pipefail

export NMD_CODE=$(realpath "$1")
git -C "$NMD_CODE" diff --exit-code HEAD
export NMD_CODE_COMMIT=$(git -C "$NMD_CODE" rev-parse HEAD)
export NMD_RUN_DIR=$(realpath "$2")
export NMD_PLOT_RUN_DIR=${3:-$NMD_RUN_DIR}
mkdir -p "$NMD_PLOT_RUN_DIR"
export NMD_PLOT_RUN_DIR=$(realpath "$NMD_PLOT_RUN_DIR")
source /rs1/researchers/a/afmejia/Software/ywu226/envs/nmd-bayesfmri-0.15.0-20261025/activate.sh
export TMPDIR=/share/statmind/ywu226/nmd-area-plots/job-$SLURM_JOB_ID/tmp
mkdir -p "$TMPDIR"
NMD_JOB_RECORD="$NMD_PLOT_RUN_DIR/analysis_records/plots-job-$SLURM_JOB_ID"
mkdir -p "$NMD_JOB_RECORD"
cp "$NMD_CODE/20260313_BayesfMRI2/scripts/05_group_level_area_plots.slurm" "$NMD_JOB_RECORD/submitted.slurm"
cp "$NMD_CODE/20260313_BayesfMRI2/scripts/05_group_level_area_plots.R" "$NMD_JOB_RECORD/plots.R"
printf '%s\n' "$NMD_CODE_COMMIT" > "$NMD_JOB_RECORD/code_commit.txt"
scontrol show job "$SLURM_JOB_ID" > "$NMD_JOB_RECORD/slurm_job.txt"
cat >> "$NMD_PLOT_RUN_DIR/run_record.org" <<EOF

* Plot job $SLURM_JOB_ID
Submitted checkout: $NMD_CODE
Code commit: $NMD_CODE_COMMIT
Analysis input directory: $NMD_RUN_DIR
Plot output directory: $NMD_PLOT_RUN_DIR
Started UTC: $(date -u '+%Y-%m-%d %H:%M:%S UTC')
Resources: statmind_cpu / compute / normal; 1 CPU; 8 GB; 4 hours
EOF
trap 'nmd_exit=$?; printf "Plot shell exit: %s at %s\n" "$nmd_exit" "$(date -u)" >> "$NMD_PLOT_RUN_DIR/run_record.org"' EXIT
cd "$(dirname "$TMPDIR")"
hostname
date -u
bash "$NMD_ENV_ROOT/run_headless_R.sh" "$NMD_CODE/20260313_BayesfMRI2/scripts/05_group_level_area_plots.R"
date -u
