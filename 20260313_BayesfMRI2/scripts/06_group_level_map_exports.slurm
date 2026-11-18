#!/usr/bin/env bash
# Export group maps on Hazel. Pass the fixed checkout and either "first" or "remaining FIRST_JOB".
# Submit the second phase after inspecting the first contrast in a browser and FSLeyes.
#SBATCH --job-name=nmd_map_export
#SBATCH --account=statmind_cpu
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20261025_BayesfMRI_AreaFiltering/logs/exports-%j.out
#SBATCH --error=/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20261025_BayesfMRI_AreaFiltering/logs/exports-%j.err
set -euo pipefail
checkout=$(realpath "$1")
git -C "$checkout" diff --exit-code HEAD
NMD_CODE_COMMIT=$(git -C "$checkout" rev-parse HEAD)
export NMD_CODE_COMMIT
export NMD_EXPORT_PHASE=$2
case "$2" in
  first) export NMD_FIRST_EXPORT_JOB=$SLURM_JOB_ID ;;
  remaining) [[ $3 =~ ^[0-9]+$ ]]; export NMD_FIRST_EXPORT_JOB=$3 ;;
  *) exit 2 ;;
esac
source /rs1/researchers/a/afmejia/Software/ywu226/envs/nmd-bayesfmri-0.15.0-20261025/activate.sh
run=/rs1/researchers/a/afmejia/Projects/brain_stimulation_utah/20261025_BayesfMRI_AreaFiltering/full/group_area_20_30_40
record="$run/analysis_records/exports-job-$SLURM_JOB_ID"
export TMPDIR="/share/statmind/ywu226/nmd-map-exports/job-$SLURM_JOB_ID"
mkdir -p "$record" "$TMPDIR"
echo "$NMD_CODE_COMMIT" > "$record/code_commit.txt"
scontrol show job "$SLURM_JOB_ID" > "$record/slurm_job.txt"
cd "$TMPDIR"
exec 9>"$run/analysis_records/map_exports.lock"
flock -n 9
bash "$NMD_ENV_ROOT/run_headless_R.sh" "$checkout/20260313_BayesfMRI2/scripts/06_group_level_map_exports.R"
