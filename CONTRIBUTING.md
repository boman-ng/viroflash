# Contributing

Maintainer: boman-ng <boman.ngs@gmail.com>. Use this repository's issues and pull
requests for questions and proposed changes. Submit a small synthetic reproduction
for bugs; do not upload patient data, internal sample names, credentials or private
logs. Report sensitive issues using [SECURITY.md](SECURITY.md).

Install the pinned Rust toolchain with Rustup and a C build toolchain (on Debian or
Ubuntu: `build-essential pkg-config zlib1g-dev`). Build with
`cargo build --release --locked`. The lockfile and toolchain file are authoritative.

Before submitting changes, run:

```bash
bash ci/check.sh
cargo test --release --test smoke --locked
git diff --check
```

Keep changes focused and use Conventional Commits. Preserve the frozen analysis
profile and documented [analysis contract](docs/analysis-contract.md). Add meaningful
Rust regression tests when behavior changes. Fixtures must be synthetic or have
documented public redistribution rights. Preserve third-party notices and record
the origin of adaptations. Submit only work you are entitled to contribute under
the applicable license; no additional CLA is imposed by this repository.

See [release preparation](docs/releasing.md) for artifact checks. Changes to versions,
publishing and ownership transfers require explicit maintainer authorization.
