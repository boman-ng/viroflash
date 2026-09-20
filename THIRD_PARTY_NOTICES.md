# Third-party notices

Viroflash's original work is Copyright (c) 2026 boman-ng, under [BSD-3-Clause](LICENSE).
This license does not replace third-party licenses. Original notice texts and a
versioned inventory are in [licenses/](licenses/).

## Code and linked components

The inventory records exact Cargo.lock versions, crate checksums, source repositories,
declared licenses and notice checksums. `release-build` includes build-time dependencies;
it does not assert that every listed component is linked into the executable.
`test-or-other-platform` components are not part of the supported Linux release.
Upstream notice files are copied without changing their legitimate attribution.

- minimap2 and minimap2-sys Rust wrappers are MIT OR Apache-2.0. The embedded minimap2
  2.30 C implementation is MIT, with Dana-Farber Cancer Institute, Broad Institute
  and additional header-specific notices. Wrapper metadata does not replace C notices.
- The local symmetric DUST implementation adapts minimap2's sdust.c. Its upstream
  MIT notice and the BSD-3-Clause notice for original Rust changes both apply.
- libz-sys bundles C zlib under the zlib license; its Rust wrapper has separate
  MIT/Apache terms. Release builds request static bundled zlib. flate2's zlib-rs
  backend is a separate component in the inventory.
- Rust's standard library includes the notices in its installed COPYRIGHT-library.html.
  The pinned musl target also includes musl, startup objects and an unwinder. Their
  original notices include the GCC Runtime Library Exception and LLVM terms.

Runtime provenance is grounded in the pinned Rust
[musl toolchain recipe](https://github.com/rust-lang/rust/blob/1.94.1/src/ci/docker/scripts/musl-toolchain.sh)
and the runtime notice lock. The build records a linker map, toolchain identity
and bundled runtime object hashes. Re-review these
materials after changing the target or toolchain.

## Containers

OCI and SIF install project notices under `/usr/share/viroflash/`. Debian/Ubuntu
packages retain their own copyright files. They are not relicensed by Viroflash.
Each image is accompanied by an installed-package inventory and an exact-version
source bundle, including distribution patch/build files in source packages.
Source descriptor versions and SHA-256 checksums must match before the image passes
distribution checks. A moving mirror link alone does not complete those checks.

## Fixtures and maintenance

The quickstart and tests use independently generated synthetic sequences; no private
panel or real sample is included. New datasets need documented source and usage rights.

The locked compile-time macro dependency paste 1.0.15 has
[RUSTSEC-2024-0436](https://rustsec.org/advisories/RUSTSEC-2024-0436.html), an
unmaintained advisory. It is tracked as a maintenance risk, not a demonstrated runtime
vulnerability. It is retained for the locked minimap2 binding; review an upstream
replacement during dependency maintenance. No advisory is globally suppressed.
