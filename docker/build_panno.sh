#!/bin/bash

# Default image name
IMAGE_NAME="aangeloo/panno:0.3.1"
PLATFORMS="linux/amd64,linux/arm64"

echo "Building multi-platform Docker image for $IMAGE_NAME..."
echo "Platforms: $PLATFORMS"

# Check if buildx is available
if ! docker buildx version > /dev/null 2>&1; then
    echo "Error: docker buildx is not installed or configured."
    exit 1
fi

# Build the image
# Note: To build for multiple platforms, you typically need to --push to a registry 
# or use a local registry. If you just want to build for the local architecture, 
# remove the --platform flag or use --load.

echo "Running build command..."
docker buildx build \
    --platform "$PLATFORMS" \
    -t "$IMAGE_NAME" \
    -f docker/Dockerfile.panno \
    --push .

if [ $? -eq 0 ]; then
    echo "Successfully built and pushed $IMAGE_NAME"
else
    echo "Build failed. Note: multi-platform builds usually require --push to a registry."
    echo "To build locally for your current architecture only, run:"
    echo "docker build -t $IMAGE_NAME -f docker/Dockerfile.panno ."
fi
