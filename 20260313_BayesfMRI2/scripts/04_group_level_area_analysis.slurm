#!/usr/bin/env bash
# Fit all 40 participants and derive areas 20/30/40 on Hazel.
# Submit with the fixed checkout as argument 1; optional argument 2 resumes a completed fit.
# Inputs are read-only; samples, checkpoints, and records remain in persistent project storage.
#SBATCH --job-name=nmd_area_fit
#SBATCH --account=statmind_cpu
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20261025_BayesfMRI_AreaFiltering/logs/analysis-%j.out
#SBATCH --error=/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20261025_BayesfMRI_AreaFiltering/logs/analysis-%j.err
set -euo pipefail

# Use the committed checkout and frozen environment
export NMD_CODE=$(realpath "$1")
git -C "$NMD_CODE" diff --exit-code HEAD
export NMD_CODE_COMMIT=$(git -C "$NMD_CODE" rev-parse HEAD)
source /rs1/researchers/a/afmejia/Software/ywu226/envs/nmd-bayesfmri-0.15.0-20261025/activate.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
NMD_OUTPUT=/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20261025_BayesfMRI_AreaFiltering
export NMD_RUN_DIR=${2:-$NMD_OUTPUT/full/group_area_20_30_40}
export TMPDIR=/share/statmind/ywu226/nmd-area-analysis/job-$SLURM_JOB_ID/tmp
mkdir -p "$NMD_RUN_DIR/analysis_records" "$TMPDIR"
NMD_JOB_RECORD="$NMD_RUN_DIR/analysis_records/analysis-job-$SLURM_JOB_ID"
mkdir -p "$NMD_JOB_RECORD"
cp "$NMD_CODE/20260313_BayesfMRI2/scripts/04_group_level_area_analysis.slurm" "$NMD_JOB_RECORD/submitted.slurm"
cp "$NMD_CODE/20260313_BayesfMRI2/scripts/04_group_level_area_analysis.R" "$NMD_JOB_RECORD/analysis.R"
printf '%s\n' "$NMD_CODE_COMMIT" > "$NMD_JOB_RECORD/code_commit.txt"
scontrol show job "$SLURM_JOB_ID" > "$NMD_JOB_RECORD/slurm_job.txt"

# INLA defaults to up to 16 threads independently of OMP_NUM_THREADS.
# Only this job's PSOCK workers read this private profile; the main script sets its own option.
export R_PROFILE_USER="$NMD_JOB_RECORD/worker_profile.R"
cat > "$R_PROFILE_USER" <<'RS'
library(methods)
INLA::inla.setOption(num.threads = "1:1")
RS
cat >> "$NMD_RUN_DIR/run_record.org" <<EOF

* Analysis job $SLURM_JOB_ID
Submitted checkout: $NMD_CODE
Code commit: $NMD_CODE_COMMIT
Started UTC: $(date -u '+%Y-%m-%d %H:%M:%S UTC')
Node: $(hostname)
Resources: statmind_cpu / compute / normal; 4 CPUs; 96 GB; 48 hours
Output directory: $NMD_RUN_DIR
Logs: $NMD_OUTPUT/logs/analysis-$SLURM_JOB_ID.{out,err}
EOF
trap 'nmd_exit=$?; printf "Analysis shell exit: %s at %s\n" "$nmd_exit" "$(date -u)" >> "$NMD_RUN_DIR/run_record.org"' EXIT
cd "$(dirname "$TMPDIR")"
hostname
date -u
Rscript --vanilla "$NMD_CODE/20260313_BayesfMRI2/scripts/04_group_level_area_analysis.R"
date -u
