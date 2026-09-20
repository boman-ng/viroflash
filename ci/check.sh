#!/usr/bin/env bash
set -euo pipefail

cargo fmt --all -- --check
rustfmt --edition 2021 --check .github/scripts/create-smoke-fixture.rs
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --all-targets --locked
python3 ci/collect-licenses.py --check
python3 ci/check-docs.py
fixture=$(mktemp -d)
rustc --edition 2021 .github/scripts/create-smoke-fixture.rs -o "$fixture/create"
"$fixture/create" "$fixture/data"
for name in host.fa target.fa sample.fastq; do cmp "examples/$name" "$fixture/data/$name"; done
