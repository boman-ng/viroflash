# CI and release strategy

## Continuous integration

The `CI` workflow runs for every pull request and every push to `master`. It has three required jobs:

- `Rust quality and tests` checks formatting and types, denies all Clippy warnings, runs debug and optimized tests, and builds the production binary with the committed lockfile.
- `Docker image` builds the OCI image, checks the embedded version command, and runs a positive synthetic dataset through a bind-mounted work directory as the image's unprivileged user.
- `Apptainer image` builds the Linux binary and SIF image on Ubuntu 22.04, then runs the same mounted-data smoke test as the host user.

Before merging, wait for all three jobs to pass. Configure them as required status checks if the repository plan supports protected private branches. The current private repository plan does not, so this gate is procedural rather than server-enforced. Dependabot may update pinned action commits and container base versions, but those changes should pass all three jobs before merge.

## Versioned releases

Releases are immutable and tag-driven. The package version in `Cargo.toml` is the source of truth. The workflow accepts only a `vMAJOR.MINOR.PATCH` tag (with an optional SemVer prerelease suffix) whose version exactly matches the Cargo package.

Prepare a release as follows:

1. Update the version in `Cargo.toml`, run `cargo check` to refresh `Cargo.lock`, then run `cargo check --locked` to verify it, update user-facing documentation, and merge the change through CI.
2. Create an annotated or signed tag at the reviewed commit, for example `git tag -s v0.2.0`.
3. Push only that tag after CI is green: `git push origin v0.2.0`.
4. Wait for the `Release` workflow. It publishes the GitHub release only after every binary, SIF, and OCI job succeeds.

The workflow publishes:

- `viroflash-VERSION-linux-x86_64.tar.gz`, a static Linux amd64 executable plus the README;
- `viroflash-VERSION-x86_64.sif`, an immutable Apptainer amd64 image;
- SHA-256 checksum files for both downloadable artifacts;
- `ghcr.io/boman-ng/viroflash:VERSION`, `:MAJOR.MINOR`, and, for stable versions, `:latest`, as a `linux/amd64` OCI image;

The repository deliberately does not publish to crates.io (`publish = false`). The OCI version tag and release assets are immutable. If a released build is defective, fix it in a new patch release; move only the convenience tags (`MAJOR.MINOR` and `latest`) forward. Never replace a versioned asset or OCI tag with different bytes.

The container job refuses to replace an existing versioned OCI tag, and GitHub also refuses to create the same release twice. Use the Actions UI to rerun a failed job only before it has published external state. If an interrupted release has already published its OCI version tag, inspect that partial release and create a new patch version instead of mutating it.

## Local image builds

Build and run the Docker image:

```bash
docker build --tag viroflash:local .
docker run --rm viroflash:local version
docker run --rm -v "$PWD:/work" viroflash:local run --help
```

The Docker image runs as UID/GID `65532`. A bind-mounted output directory must therefore be writable by that identity. Input data and indexes are not embedded in the image.

Build the static Linux binary and Apptainer image on Ubuntu 22.04 or an ABI-compatible system:

```bash
rustup target add x86_64-unknown-linux-musl
sudo apt-get install musl-tools
CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER=musl-gcc \
RUSTFLAGS="-C target-feature=+crt-static -C link-arg=-Wl,--no-dynamic-linker" \
  cargo build --release --locked --target x86_64-unknown-linux-musl
cp target/x86_64-unknown-linux-musl/release/viroflash target/release/viroflash
apptainer build viroflash.sif Apptainer.def
apptainer run --bind "$PWD:/work" --pwd /work viroflash.sif version
```

Apptainer runs with the invoking host UID and is the preferred image for shared HPC filesystems. Published binaries, SIF files, and OCI images are amd64-only. The standalone binary is statically linked and has no host glibc or zlib requirement.

## Verification

Verify downloaded files before use:

```bash
sha256sum --check viroflash-0.2.0-linux-x86_64.tar.gz.sha256
sha256sum --check viroflash-0.2.0-x86_64.sif.sha256
```

Prefer versioned image tags or recorded OCI digests in production and HPC workflows. `latest` is only a discovery convenience.

GitHub artifact attestations are intentionally disabled while this repository is private on a non-Enterprise plan. GitHub supports attestations for private repositories only on Enterprise Cloud. Enable them if the repository becomes public or moves to Enterprise Cloud.
