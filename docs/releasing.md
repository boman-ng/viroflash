# CI and Release Strategy

## Continuous Integration

The `CI` workflow runs formatting, type checks, Clippy with warnings denied, debug and release
tests, and the Phase 0 verifier with the committed lockfile. Its Docker and Apptainer jobs build a
deterministic HOST/TARGET fixture, create a reusable index, run analysis through `--index`, and
assert that the output directory contains exactly:

```text
perf.json
report.csv
report.html
```

The Docker smoke runs as UID/GID `65532`; the mounted fixture directory must be writable by that
identity. The Apptainer smoke runs as the invoking host user.

## Tagged Releases

Releases are immutable and tag-driven. A release tag must be `vMAJOR.MINOR.PATCH`, optionally with
a SemVer prerelease suffix, and must match the package version already reviewed in `Cargo.toml`.
Changing that package version requires explicit authorization and a separate reviewed change.

The release workflow:

1. validates the tag against `Cargo.toml`;
2. builds and inspects the static Linux amd64 binary;
3. unpacks that binary and exercises the current `index` then `run --index` contract;
4. builds and smoke-tests the Apptainer image through the same contract;
5. publishes the amd64 OCI image only if its versioned tag does not already exist;
6. verifies checksums before creating the GitHub release.

Published assets are the Linux amd64 archive and checksum, the amd64 SIF and checksum, and OCI
tags for the exact version plus stable convenience tags. Do not replace versioned assets or OCI
tags. Repair a released defect with a newly authorized patch release.

## Local Docker Check

```bash
test_dir="$(mktemp -d)"
python3 .github/scripts/create-smoke-fixture.py "${test_dir}"
chmod -R a+rwX "${test_dir}"
docker build --tag viroflash:local .
docker run --rm viroflash:local --help
docker run --rm --volume "${test_dir}:/work" viroflash:local index \
  --host-fa /work/host.fa --target-fa /work/target.fa --threads 2 --out /work/index
docker run --rm --volume "${test_dir}:/work" viroflash:local run \
  --r1 /work/sample.fastq --index /work/index --threads 2 --out /work/result
find "${test_dir}/result" -mindepth 1 -maxdepth 1 -printf '%f\n' | sort
```

## Local Apptainer Check

Build the release binary first, then the image:

```bash
cargo build --release --locked
apptainer build viroflash.sif Apptainer.def
test_dir="$(mktemp -d)"
python3 .github/scripts/create-smoke-fixture.py "${test_dir}"
apptainer run viroflash.sif --help
apptainer run --bind "${test_dir}:/work" --pwd /work viroflash.sif index \
  --host-fa /work/host.fa --target-fa /work/target.fa --threads 2 --out /work/index
apptainer run --bind "${test_dir}:/work" --pwd /work viroflash.sif run \
  --r1 /work/sample.fastq --index /work/index --threads 2 --out /work/result
find "${test_dir}/result" -mindepth 1 -maxdepth 1 -printf '%f\n' | sort
```

The listed result must be exactly the three successful artifacts above. Verify downloaded release
assets with `sha256sum --check` before use.
