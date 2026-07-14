#!/usr/bin/env bash
# Build (once) and run the ATAC-seq Snakemake workflow in Docker.
#
# Usage:
#   ./run_pipeline.sh                                 # main pipeline, 4 cores
#   ./run_pipeline.sh snakefile_ATACseq --cores 16    # main pipeline, 16 cores
#   ./run_pipeline.sh snakefile_ATAC_QC --cores 16    # QC pipeline (run AFTER the main one)
#   ./run_pipeline.sh snakefile_ATACseq -n            # dry run
#
# The current directory (code + ref/ genomes + data/ FASTQs) is mounted at
# /workflow; results/ are written back here. The pre-built conda envs live in the
# image at /opt/wf-conda and are reused via --conda-prefix.
set -euo pipefail

IMAGE="atacseq-spikein:latest"

SNAKEFILE="${1:-snakefile_ATACseq}"
if [ "$#" -gt 0 ]; then shift; fi
if [ "$#" -eq 0 ]; then set -- --cores 4; fi   # default snakemake args

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    echo ">> Building $IMAGE (first time only; pre-builds the conda envs, ~15-30 min)..."
    docker build -t "$IMAGE" .
fi

echo ">> snakemake -s ${SNAKEFILE} $*"
docker run --rm \
    -v "$(pwd)":/workflow \
    -e HOME=/tmp \
    --user "$(id -u):$(id -g)" \
    "$IMAGE" -s "${SNAKEFILE}" "$@"
