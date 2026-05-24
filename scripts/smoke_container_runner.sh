#!/usr/bin/env bash
# Stage 37: Container Runner Smoke Test
#
# Creates a temporary workspace under /tmp, copies a test profile, creates
# placeholder FASTQs matching the profile's sample sheet, and runs a
# bind-mounted dry-run.
#
# Usage:
#   bash scripts/smoke_container_runner.sh docker <image-tag>
#   bash scripts/smoke_container_runner.sh singularity <sif-path>
#
# Set DOCKER_CMD to "sudo docker" or another Docker executable if needed.
#
# Does NOT download public data.  Does NOT write outputs into the repository.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# ------------------------------------------------------------------------
# Parse mode and image
# ------------------------------------------------------------------------
MODE="${1:-}"
IMAGE="${2:-}"

if [[ -z "$MODE" || -z "$IMAGE" ]]; then
    echo "Usage: $0 docker <image-tag>"
    echo "       $0 singularity <sif-path>"
    exit 1
fi

if [[ "$MODE" != "docker" && "$MODE" != "singularity" ]]; then
    echo "ERROR: mode must be 'docker' or 'singularity', got '$MODE'"
    exit 1
fi

# Configurable Docker command (e.g., "sudo docker")
DOCKER_CMD="${DOCKER_CMD:-docker}"

# ------------------------------------------------------------------------
# Create temporary smoke workspace
# ------------------------------------------------------------------------
SMOKE_DIR="$(mktemp -d /tmp/smoke_container_runner.XXXXXX)"
trap 'rm -rf "$SMOKE_DIR"' EXIT

echo "=== Smoke workspace: $SMOKE_DIR ==="

# Copy test profile
PROFILE="$REPO_ROOT/test/profiles/chipseq_pe_noctrl"
if [[ ! -f "$PROFILE/config.yaml" ]]; then
    echo "ERROR: test profile not found at $PROFILE"
    exit 1
fi
cp "$PROFILE/config.yaml" "$SMOKE_DIR/"
cp "$PROFILE/samples.tsv" "$SMOKE_DIR/"

# Parse sample sheet to discover required FASTQ filenames
FASTQ_FILES=($(tail -n +2 "$SMOKE_DIR/samples.tsv" | awk -F'\t' '{print $2; if($3!="") print $3}' | sort -u))
if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "ERROR: no FASTQ files found in sample sheet"
    exit 1
fi

# Create placeholder FASTQs in the smoke dir
for fq in "${FASTQ_FILES[@]}"; do
    touch "$SMOKE_DIR/$fq"
    echo "  Created placeholder: $fq"
done

# Create conda cache directories
CONDA_CACHE="$SMOKE_DIR/conda_cache"
mkdir -p "$CONDA_CACHE/snakemake" "$CONDA_CACHE/home" "$CONDA_CACHE/xdg-cache"

echo "=== Smoke workspace prepared ==="

# ------------------------------------------------------------------------
# Version check
# ------------------------------------------------------------------------
echo "=== Checking Snakemake version ==="
if [[ "$MODE" == "docker" ]]; then
    $DOCKER_CMD run --rm "$IMAGE" --version
else
    singularity exec "$IMAGE" snakemake --version
fi

# ------------------------------------------------------------------------
# PyYAML import check
# ------------------------------------------------------------------------
echo "=== Checking PyYAML import ==="
if [[ "$MODE" == "docker" ]]; then
    $DOCKER_CMD run --rm --entrypoint python "$IMAGE" -c "import yaml; print('runner ok')"
else
    singularity exec "$IMAGE" python -c "import yaml; print('runner ok')"
fi

# ------------------------------------------------------------------------
# Bind-mounted dry-run
# ------------------------------------------------------------------------
echo "=== Running bind-mounted dry-run ==="
if [[ "$MODE" == "docker" ]]; then
    $DOCKER_CMD run --rm \
        -v "$SMOKE_DIR":/smoke \
        -v "$REPO_ROOT":/workspace \
        -v "$CONDA_CACHE":/conda_cache \
        -e HOME=/conda_cache/home \
        -e XDG_CACHE_HOME=/conda_cache/xdg-cache \
        -u "$(id -u):$(id -g)" \
        "$IMAGE" \
        -s /workspace/workflow/Snakefile \
        --directory /smoke \
        --configfile /smoke/config.yaml \
        --cores 1 \
        --use-conda \
        --conda-prefix /conda_cache/snakemake \
        -n
else
    singularity exec \
        --pwd /workspace \
        --bind "$SMOKE_DIR":/smoke,"$REPO_ROOT":/workspace,"$CONDA_CACHE":/conda_cache \
        --env XDG_CACHE_HOME=/conda_cache/xdg-cache \
        "$IMAGE" \
        snakemake \
        -s /workspace/workflow/Snakefile \
        --directory /smoke \
        --configfile /smoke/config.yaml \
        --cores 1 \
        --use-conda \
        --conda-prefix /conda_cache/snakemake \
        -n
fi

echo ""
echo "=== Smoke test PASSED ($MODE) ==="
echo "Workspace cleaned: $SMOKE_DIR"
