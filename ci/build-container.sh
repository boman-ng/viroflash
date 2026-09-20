#!/usr/bin/env bash
set -euo pipefail
kind=$1 output=$2 temporary=$3 version=$4 revision=$5 source_url=$6
package="viroflash-$version-linux-x86_64"
data="$temporary/release-test"
if [[ "$kind" == oci ]]; then
  image="viroflash:candidate-$revision"
  docker build --platform linux/amd64 --file Dockerfile --tag "$image" \
    --build-arg "VERSION=$version" --build-arg "VCS_REF=$revision" --build-arg "SOURCE_URL=$source_url" "$output/$package"
  docker inspect "$image" > "$temporary/oci-inspect.json"
  jq -e --arg version "$version" --arg revision "$revision" '.[0].Config.Labels | .["org.opencontainers.image.version"] == $version and .["org.opencontainers.image.revision"] == $revision' "$temporary/oci-inspect.json"
  docker history --no-trunc "$image" > "$temporary/oci-history.txt"
  for file in LICENSE THIRD_PARTY_NOTICES.md licenses/manifest.json examples/sample.fastq; do
    docker run --rm --network none --entrypoint cat "$image" "/usr/share/viroflash/$file" | cmp - "$file"
  done
  docker run --rm --network none --user "$(id -u):$(id -g)" --volume "$data:/work" "$image" index --host-fa host.fa --target-fa target.fa --out oci-index --threads 2
  docker run --rm --network none --user "$(id -u):$(id -g)" --volume "$data:/work" "$image" run --r1 sample.fastq --index oci-index --out oci-result --threads 2
  docker save "$image" | gzip -n > "$output/viroflash-$version-oci.tar.gz"
  docker run --rm --network none --entrypoint dpkg-query "$image" -W '-f=${binary:Package}\t${Version}\t${source:Package}\t${source:Version}\n' > "$output/oci-packages.tsv"
elif [[ "$kind" == sif ]]; then
  definition=$(realpath Apptainer.def)
  image="$output/viroflash-$version-x86_64.sif"
  # Relative %files paths keep private workspace names out of embedded recipes.
  (cd "$output"; apptainer build --force --build-arg "PACKAGE=$package" --build-arg "VERSION=$version" \
    --build-arg "REVISION=$revision" --build-arg "SOURCE_URL=$source_url" "$image" "$definition")
  apptainer inspect --json "$image" > "$temporary/sif-inspect.json"
  jq -e --arg version "$version" --arg revision "$revision" '.data.attributes.labels | .["org.opencontainers.image.version"] == $version and .["org.opencontainers.image.revision"] == $revision' "$temporary/sif-inspect.json"
  for file in LICENSE THIRD_PARTY_NOTICES.md licenses/manifest.json examples/sample.fastq; do
    apptainer exec "$image" cat "/usr/share/viroflash/$file" | cmp - "$file"
  done
  apptainer run --bind "$data:/work" --pwd /work "$image" index --host-fa /work/host.fa --target-fa /work/target.fa --out /work/sif-index --threads 2
  apptainer run --bind "$data:/work" --pwd /work "$image" run --r1 /work/sample.fastq --index /work/sif-index --out /work/sif-result --threads 2
  apptainer exec "$image" dpkg-query -W '-f=${binary:Package}\t${Version}\t${source:Package}\t${source:Version}\n' > "$output/sif-packages.tsv"
else
  echo 'Expected oci or sif' >&2; exit 1
fi
for report in report.csv report.html; do cmp "$data/result/$report" "$data/$kind-result/$report"; done
jq -e '.status == "SUCCESS" and .input_fragments == 20 and .selected_fragments == 20 and .aligned_fragments == 20' "$data/$kind-result/perf.json"
