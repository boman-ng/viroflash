# Packaging and releases

End users install the prebuilt Linux x86_64 archive or use a SIF/OCI image. Building is a maintainer task. The GitHub repository is private; users need repository access or an administrator-supplied package.

## Build a local package

Use the pinned Rust toolchain, the committed lockfile, and a musl C toolchain:

```bash
rustup target add x86_64-unknown-linux-musl
CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER=musl-gcc \
RUSTFLAGS='-C target-feature=+crt-static -C link-arg=-Wl,--no-dynamic-linker' \
cargo build --release --locked --target x86_64-unknown-linux-musl
version=$(cargo metadata --locked --no-deps --format-version 1 | jq -r '.packages[0].version')
package="viroflash-${version}-linux-x86_64"
mkdir -p "dist/${package}"
install -m755 target/x86_64-unknown-linux-musl/release/viroflash "dist/${package}/viroflash"
cp README.md "dist/${package}/README.md"
tar -czf "dist/${package}.tar.gz" -C dist "${package}"
```

The archive contains the static executable and README. Extract it, install with the system `install` command, then exercise `index` and `run` using the installed binary. Release checks verify that the binary needs no runtime loader or shared libraries. Builds and tests do not need `evaluation/`.

Build each container from the same static executable, without recompiling Rust:

```bash
docker build -f Dockerfile --build-arg VERSION="$version" -t viroflash:local "dist/$package"
apptainer build --force --build-arg VERSION="$version" "dist/viroflash-${version}-x86_64.sif" Apptainer.def
```

Docker uses the unpacked package as its build context. Apptainer reads the musl executable from `target/x86_64-unknown-linux-musl/release/`; `--force` replaces the base image's inherited version label and any existing local SIF. Container validation belongs to release preparation; ordinary CI runs Rust formatting, Clippy, and core tests in one build profile.

## Publish only when authorized

A release tag must be `vMAJOR.MINOR.PATCH` (optionally with a SemVer prerelease suffix), matching `Cargo.toml`. Package/schema version changes and publishing require explicit authorization. Never replace an existing versioned asset or image tag.

The release workflow uses one job and one Rust release build. It validates the tag, installs and tests the archive, builds SIF and OCI images from that binary, and compares their CSV/HTML against the installed executable. After those checks, it publishes the OCI tags and creates the GitHub release with archive/SIF checksums. Stable releases also update the minor-version and `latest` image tags. No intermediate artifact uploads or downloads are needed.

Local packages built from uncommitted work are unpublished artifacts. Do not present them as an existing GitHub release or point users to an older release as if it contained the current implementation.
