#!/usr/bin/env bash
# Set job requirements. Number of nodes equals number of lengths to execute.
#SBATCH -N 1
#SBATCH -t 120:00:00

# Store maximum runtime, set dataset and lengths to use to use.
MAX_TIME=431700s

# Store name of this job.
job="s8_1000_10"

module load 2021
module load Python/3.9.5-GCCcore-10.3.0

# Create directories for results.
mkdir -p "${TMPDIR}/${job}/results/"

# Copy existing results if there are any.
if [[ -d $HOME/quantification_anonymous/jobs/${job}/results/ ]]; then
    cp -r "${HOME}/quantification_anonymous/jobs/${job}/results" "${TMPDIR}/${job}/"
fi

# Export scratch path for python script.
export LISA_DIR="$TMPDIR"

timeout ${MAX_TIME} python3 "${HOME}/quantification_anonymous/run_program.py" \
    dfs_bnb_parallel -j ${job} -l -le 10

# Copy results from scratch to home.
cp -r "${TMPDIR}/${job}/results" "${HOME}/quantification_anonymous/jobs/${job}/"
