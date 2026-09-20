#!/usr/bin/env bash
# Build and validate local artifacts. This script never publishes.
set -euo pipefail

root=$(git rev-parse --show-toplevel)
cd "$root"
: "${RUNNER_TEMP:?Set a temporary directory for build evidence}"
output=${OUTPUT_DIR:-dist}
[[ ! -e "$output" ]] || { echo 'Output directory must be new' >&2; exit 1; }
version=$(cargo metadata --locked --no-deps --format-version 1 | jq -r '.packages[0].version')
revision=$(git rev-parse HEAD)
source_url=${SOURCE_URL:-https://github.com/boman-ng/viroflash}
if [[ -n "${RELEASE_TAG:-}" ]]; then
  [[ "$RELEASE_TAG" =~ ^v[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.-]+)?$ && "$RELEASE_TAG" == "v$version" ]] || {
    echo 'Release tag must match Cargo version' >&2; exit 1;
  }
  [[ "$(git rev-parse "$RELEASE_TAG^{commit}")" == "$revision" ]] || exit 1
fi
[[ -z "$(git status --porcelain)" ]] || { echo 'Candidate source must be committed and clean' >&2; exit 1; }
python3 ci/collect-licenses.py --check
cargo test --release --test smoke --locked
mkdir -p "$output" "$RUNNER_TEMP"
output=$(realpath "$output")
RUNNER_TEMP=$(realpath "$RUNNER_TEMP")
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER=${CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER:-musl-gcc}
export LIBZ_SYS_STATIC=1 ZLIB_NO_PKG_CONFIG=1
cargo_root=$(realpath "${CARGO_HOME:-$HOME/.cargo}")
rust_sysroot=$(rustc --print sysroot)
printf -v CARGO_ENCODED_RUSTFLAGS '%s\037' '-C' 'target-feature=+crt-static' '-C' 'link-arg=-Wl,--no-dynamic-linker' \
  "--remap-path-prefix=$root=/usr/src/viroflash" "--remap-path-prefix=$cargo_root=/usr/local/cargo" \
  "--remap-path-prefix=$rust_sysroot=/usr/local/rust" '-C' "link-arg=-Wl,-Map,$RUNNER_TEMP/linker.map"
export CARGO_ENCODED_RUSTFLAGS=${CARGO_ENCODED_RUSTFLAGS%$'\037'}
export CFLAGS_x86_64_unknown_linux_musl="-ffile-prefix-map=$root=/usr/src/viroflash -ffile-prefix-map=$cargo_root=/usr/local/cargo"
cargo build --release --locked --target x86_64-unknown-linux-musl
binary="${CARGO_TARGET_DIR:-target}/x86_64-unknown-linux-musl/release/viroflash"
file "$binary" | grep -q x86-64
readelf --program-headers "$binary" > "$RUNNER_TEMP/elf-program-headers.txt"
readelf --dynamic "$binary" > "$RUNNER_TEMP/elf-dynamic.txt"
if grep -q INTERP "$RUNNER_TEMP/elf-program-headers.txt" || grep -q '(NEEDED)' "$RUNNER_TEMP/elf-dynamic.txt"; then
  echo 'Executable requires a runtime loader or shared libraries' >&2; exit 1
fi

package="viroflash-$version-linux-x86_64"
mkdir "$output/$package"
install -m755 "$binary" "$output/$package/viroflash"
cp README.md LICENSE THIRD_PARTY_NOTICES.md SECURITY.md CONTRIBUTING.md CITATION.cff CHANGELOG.md "$output/$package/"
cp -R licenses examples docs "$output/$package/"
# Package permissions must not inherit a private build-evidence umask.
find "$output/$package" -type d -exec chmod 755 {} +
find "$output/$package" -type f -exec chmod 644 {} +
chmod 755 "$output/$package/viroflash"
tar --owner=0 --group=0 --numeric-owner -czf "$output/$package.tar.gz" -C "$output" "$package"
mkdir "$RUNNER_TEMP/unpacked"
tar -xzf "$output/$package.tar.gz" -C "$RUNNER_TEMP/unpacked"
installed="$RUNNER_TEMP/unpacked/$package"
data="$RUNNER_TEMP/release-test"
mkdir "$data"
cp "$installed"/examples/* "$data/"
"$installed/viroflash" index --host-fa "$data/host.fa" --target-fa "$data/target.fa" --out "$data/index" --threads 2
"$installed/viroflash" run --r1 "$data/sample.fastq" --index "$data/index" --out "$data/result" --threads 2
test "$(ls "$data/result")" = "$(printf '%s\n' perf.json report.csv report.html)"
jq -e '.status == "SUCCESS" and .input_fragments == 20 and .prescreen_passed_fragments == 20 and .selected_fragments == 20 and .aligned_fragments == 20' "$data/result/perf.json"
test "$(wc -l < "$data/result/report.csv")" -eq 2
awk -F, 'NR == 2 {exit !($3 == "target" && $5 == 20 && $6 == 100)}' "$data/result/report.csv"

for kind in oci sif; do
  bash ci/build-container.sh "$kind" "$output" "$RUNNER_TEMP" "$version" "$revision" "$source_url"
  bash ci/collect-container-sources.sh "$kind" "$output" "$RUNNER_TEMP/$kind-sources"
done
# Existing publishing scripts consume this local tag only after all checks pass.
if [[ -n "${RELEASE_TAG:-}" ]]; then
  docker tag "viroflash:candidate-$revision" viroflash:release-test
fi
jq -n --arg version "$version" --arg revision "$revision" --arg source "$source_url" \
  --arg rust "$(rustc --version)" --arg oci "$(awk '/^FROM / {print $2}' Dockerfile)" \
  --arg sif "$(awk '/^From: / {print $2}' Apptainer.def)" \
  '{version:$version, revision:$revision, source:$source, rust:$rust, oci_base:$oci, sif_base:$sif,
    local_image:("viroflash:candidate-"+$revision), publication:"not performed by builder"}' > "$output/build.json"
runtime="$rust_sysroot/lib/rustlib/x86_64-unknown-linux-musl/lib/self-contained"
(cd "$runtime"; sha256sum *.o *.a) > "$output/runtime-objects.sha256"
(cd "$output"; for asset in *.tar.gz *.sif; do sha256sum "$asset" > "$asset.sha256"; done
 sha256sum *.tar.gz *.sif *.tsv build.json runtime-objects.sha256 > SHA256SUMS
 sha256sum -c SHA256SUMS)
echo "Validated local artifacts: $output"
