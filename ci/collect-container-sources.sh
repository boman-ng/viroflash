#!/usr/bin/env bash
set -euo pipefail
kind=${1:?Expected oci or sif}
assets=$(realpath "${2:?Expected candidate asset directory}")
work=${3:?Expected NEW private working directory}
[[ ! -e "$work" ]] || { echo 'Refusing to overwrite source evidence' >&2; exit 1; }
umask 077
mkdir -m700 "$work"
work=$(realpath "$work")
case "$kind" in
  oci) base=$(awk '/^FROM / {print $2}' Dockerfile) ;;
  sif) base=$(awk '/^From: / {print $2}' Apptainer.def) ;;
  *) echo 'Expected oci or sif' >&2; exit 1 ;;
esac
[[ "$base" =~ @sha256:[0-9a-f]{64}$ ]]
cp "$assets/$kind-packages.tsv" "$work/expected-packages.tsv"
ca_file=$(python3 -c 'import ssl; print(ssl.get_default_verify_paths().cafile or "")')
[[ -f "$ca_file" ]] || { echo 'Host trusted CA bundle is required for HTTPS source retrieval' >&2; exit 1; }
if ! docker run --rm --network host --platform linux/amd64 --user 0 --volume "$work:/out" \
  --env http_proxy --env https_proxy --env no_proxy \
  --volume "$ca_file:/run/viroflash-ca.pem:ro" \
  --volume "$(realpath ci/download-container-sources.sh):/collect.sh:ro" \
  --entrypoint /bin/sh "$base" -c 'sh /collect.sh; status=$?; chown -R "$1:$2" /out; exit "$status"' \
  sh "$(id -u)" "$(id -g)"; then
  if [[ "$kind" == sif && -s "$work/unavailable-sources.tsv" ]]; then
    python3 ci/fetch-ubuntu-sources.py "$work"
  else
    echo 'Corresponding-source collection incomplete' >&2; exit 1
  fi
fi
python3 ci/package-container-sources.py "$work" "$assets/$kind-corresponding-source.tar.gz"
