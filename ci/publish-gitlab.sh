#!/usr/bin/env bash
set -euo pipefail

VERSION="${RELEASE_TAG#v}"
project_api="${CI_API_V4_URL}/projects/${CI_PROJECT_ID}"
package_url="${project_api}/packages/generic/viroflash/${VERSION}"

# Fail on existing releases/assets and on unexpected API errors before uploading.
require_absent() {
  local status
  status="$(curl --silent --show-error --output /dev/null --write-out '%{http_code}' \
    --header "JOB-TOKEN: ${CI_JOB_TOKEN}" "$1")"
  if [[ "$status" != 404 ]]; then
    echo "Refusing to publish: $1 returned HTTP $status (expected 404)" >&2
    exit 1
  fi
}

require_absent "${project_api}/releases/${RELEASE_TAG}"
assets=(dist/*.tar.gz dist/*.sif dist/*.sha256 dist/*.tsv dist/build.json dist/SHA256SUMS)
for asset in "${assets[@]}"; do
  require_absent "${package_url}/${asset##*/}"
done

printf '%s' "$CI_REGISTRY_PASSWORD" | docker login "$CI_REGISTRY" \
  --username "$CI_REGISTRY_USER" --password-stdin
if docker manifest inspect "${CI_REGISTRY_IMAGE}:${VERSION}" >/dev/null 2>&1; then
  echo "Versioned OCI tag ${VERSION} already exists and is immutable" >&2
  exit 1
fi
tags=("${VERSION}")
if [[ "${VERSION}" != *-* ]]; then tags+=("${VERSION%.*}" latest); fi
for tag in "${tags[@]}"; do
  docker tag viroflash:release-test "${CI_REGISTRY_IMAGE}:${tag}"
  docker push "${CI_REGISTRY_IMAGE}:${tag}"
done

links='[]'
for asset in "${assets[@]}"; do
  name="${asset##*/}"
  url="${package_url}/${name}"
  curl --fail --silent --show-error --header "JOB-TOKEN: ${CI_JOB_TOKEN}" \
    --upload-file "$asset" "$url"
  links="$(jq --arg name "$name" --arg url "$url" \
    '. + [{name: $name, url: $url, direct_asset_path: ("/" + $name), link_type: "package"}]' <<<"$links")"
done

jq -n --arg tag "$RELEASE_TAG" --arg revision "$CI_COMMIT_SHA" --argjson links "$links" \
  '{name: ("Viroflash " + $tag), tag_name: $tag,
    description: ("Linux x86_64 static binary, SIF image and SHA-256 checksums. Built from " + $revision + "."),
    assets: {links: $links}}' |
  curl --fail --silent --show-error --request POST \
    --header "JOB-TOKEN: ${CI_JOB_TOKEN}" --header 'Content-Type: application/json' \
    --data-binary @- "${project_api}/releases"
