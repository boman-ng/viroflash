# Packaging and releases

The supported distribution is Linux x86_64: a static executable archive, OCI image
and SIF. Local builds and publication are separate steps.

## Build and verify

Install the pinned Rust toolchain, Python 3.9 or newer, jq, a musl C toolchain,
file, binutils, Docker and Apptainer. Source retrieval uses HTTPS and the host
Python installation's trusted CA bundle, mounted read-only in temporary Linux
containers. Source-download containers use host networking and inherit the host's
optional http_proxy/https_proxy/no_proxy environment; analysis runs use
the existing isolated execution settings.
Debian/Ubuntu also need build-essential,
pkg-config and zlib1g-dev for the native test build. The CI setup helper installs
musl-tools and Apptainer; local users may supply these tools themselves.

From a clean, committed checkout:

```bash
rustup target add x86_64-unknown-linux-musl
bash ci/check.sh
cargo test --release --test smoke --locked
RUNNER_TEMP=$(mktemp -d) OUTPUT_DIR="$PWD/dist/local-candidate" bash ci/build-release.sh
```

OUTPUT_DIR must not exist. The builder has no publication step. It builds the
static binary once, copies it into the archive and both containers, tests the
installed archive, compares container CSV/HTML output, verifies bundled notices
and collects exact container source packages. Failure leaves partial artifacts
for investigation; they are not a completed candidate. Rerun into new directories.

The candidate contains the binary archive, OCI archive, SIF, separate container
package inventories and corresponding-source archives, build metadata and
checksums. Sources include the distribution's .dsc and referenced original and
Debian/Ubuntu patch archives. Missing exact sources block container completion.
Base digests and toolchain identity are recorded; byte-identical builds are not
claimed. Build evidence in RUNNER_TEMP includes the linker map and container inspection.

Regenerate notices using Python 3.11+ after a reviewed dependency/toolchain change with
`python3 ci/collect-licenses.py`; `--check` validates the committed inventory.
The collector verifies cached crate archives against Cargo.lock and compares the
extracted source before copying notices. Runtime notice versions and checksums
are pinned separately. No notices are automatically relicensed.

## Exercise the archive and containers

Replace VERSION and REVISION with the values in build.json. From the asset directory:

```bash
sha256sum -c SHA256SUMS
tar -xzf viroflash-VERSION-linux-x86_64.tar.gz
cd viroflash-VERSION-linux-x86_64
./viroflash index --host-fa examples/host.fa --target-fa examples/target.fa --out demo-index --threads 2
./viroflash run --r1 examples/sample.fastq --index demo-index --out demo-results --threads 2
```

Use a writable working directory and new output names for each container:

```bash
docker load -i viroflash-VERSION-oci.tar.gz
docker run --rm --network none --user "$(id -u):$(id -g)" --volume "$PWD:/work" viroflash:candidate-REVISION index --host-fa /usr/share/viroflash/examples/host.fa --target-fa /usr/share/viroflash/examples/target.fa --out oci-index --threads 2
docker run --rm --network none --user "$(id -u):$(id -g)" --volume "$PWD:/work" viroflash:candidate-REVISION run --r1 /usr/share/viroflash/examples/sample.fastq --index oci-index --out oci-results --threads 2
apptainer run --bind "$PWD:/work" --pwd /work viroflash-VERSION-x86_64.sif index --host-fa /usr/share/viroflash/examples/host.fa --target-fa /usr/share/viroflash/examples/target.fa --out sif-index --threads 2
apptainer run --bind "$PWD:/work" --pwd /work viroflash-VERSION-x86_64.sif run --r1 /usr/share/viroflash/examples/sample.fastq --index sif-index --out sif-results --threads 2
```

All runs should report 20 input fragments and 20 supports for one target group.
Examples validate installation, not biological sensitivity or specificity.

## Authorized publication and optional GitLab mirror

Publishing, new versions, tags and ownership transfers are separate authorized
actions. Confirm the rights to distribute the code and included materials under
their stated licenses; passing technical checks does not establish those rights.
Review source, history and generated artifacts for secrets and identifying data
before publication; keep sensitive audit evidence outside distributed artifacts.
See [security and privacy reporting](../SECURITY.md).
A release tag must be vMAJOR.MINOR.PATCH (optionally with a SemVer prerelease suffix),
matching Cargo.toml.

GitHub and GitLab reuse ci/check.sh and ci/build-release.sh. A v* tag invokes release
validation before upload. Archives, SIF, OCI export, corresponding-source bundles,
package inventories, metadata and checksums are distributed together. Stable tags
also update minor/latest image tags. Uploads are not transactional: inspect partial
publication before any retry, and never overwrite existing versioned assets.

The optional GitLab mirror retains its Docker executor for ordinary CI and a
protected Ubuntu 22.04 x86_64 shell runner tagged viroflash-release for publication.
It needs Python 3.9+, Rustup, Docker, Apptainer, musl, jq and build utilities.
Enable Container and Package Registries, protect tags, disable duplicate generic
package files and exempt releases from cleanup. Built-in job credentials are used;
do not introduce personal tokens. Mirroring and pushing tags are not configured
by these scripts.
