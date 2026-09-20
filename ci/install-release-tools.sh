#!/usr/bin/env bash
set -euo pipefail

APPTAINER_VERSION=1.5.3
package="apptainer_${APPTAINER_VERSION}_amd64.deb"
curl --fail --location --silent --show-error \
  "https://github.com/apptainer/apptainer/releases/download/v${APPTAINER_VERSION}/${package}" \
  --output "${RUNNER_TEMP}/${package}"
sudo apt-get update
sudo apt-get install --yes musl-tools "${RUNNER_TEMP}/${package}"
rustup target add x86_64-unknown-linux-musl
