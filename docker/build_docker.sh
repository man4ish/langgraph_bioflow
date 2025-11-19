#!/bin/bash
# =============================================
# Build and push Docker images for multiple pipelines
# =============================================

set -e  # Exit immediately if a command exits with a non-zero status
set -o pipefail  # Fail if any command in a pipeline fails

# -----------------------------
# Build general genomics pipeline image
# -----------------------------
GENOMICS_DOCKERFILE="Dockerfile.genomics"
GENOMICS_IMAGE="man4ish/genomics:latest"

echo "Building Genomics Docker image..."
docker build -f ${GENOMICS_DOCKERFILE} -t ${GENOMICS_IMAGE} .
echo "Genomics Docker image built: ${GENOMICS_IMAGE}"

# -----------------------------
# Build Cell Ranger scRNA-seq image
# -----------------------------
CELLRANGER_DOCKERFILE="Dockerfile.cellranger_scRNAseq"
CELLRANGER_IMAGE="man4ish/cellranger-scrnaseq:latest"

echo "Building Cell Ranger scRNA-seq Docker image..."
docker build -f ${CELLRANGER_DOCKERFILE} -t ${CELLRANGER_IMAGE} .
echo "Cell Ranger scRNA-seq Docker image built: ${CELLRANGER_IMAGE}"

# -----------------------------
# Build RNA-seq pipeline image with specific tool versions
# -----------------------------
STAR_VERSION="2.6.0c"
KALLISTO_VERSION="v0.46.1"
SAMTOOLS_VERSION="1.9"

RNASEQ_DOCKERFILE="Dockerfile.rnaseq"
RNASEQ_IMAGE="docker.io/man4ish/rnaseq:latest"

echo "Building RNA-seq Docker image with:"
echo "  STAR: ${STAR_VERSION}"
echo "  Kallisto: ${KALLISTO_VERSION}"
echo "  Samtools: ${SAMTOOLS_VERSION}"

docker build \
    --build-arg star_version=${STAR_VERSION} \
    --build-arg kallisto_version=${KALLISTO_VERSION} \
    --build-arg samtools_version=${SAMTOOLS_VERSION} \
    -f ${RNASEQ_DOCKERFILE} \
    -t ${RNASEQ_IMAGE} .

# -----------------------------
# Push RNA-seq image if build succeeded
# -----------------------------
echo "Docker image built successfully: ${RNASEQ_IMAGE}"
echo "Pushing image to Docker Hub..."
docker push ${RNASEQ_IMAGE}
echo "Docker image pushed successfully."
