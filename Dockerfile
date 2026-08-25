# syntax=docker/dockerfile:1

ARG RUST_VERSION=1.94.1
ARG DEBIAN_RELEASE=bookworm-slim

FROM rust:${RUST_VERSION}-bookworm AS builder
WORKDIR /src

COPY Cargo.toml Cargo.lock ./
COPY src ./src

RUN cargo build --release --locked

FROM debian:${DEBIAN_RELEASE} AS runtime

ARG VERSION=dev
ARG VCS_REF=unknown
LABEL org.opencontainers.image.title="viroflash" \
      org.opencontainers.image.description="Viral candidate detection command-line tool" \
      org.opencontainers.image.source="https://github.com/boman-ng/viroflash" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}"

RUN apt-get update \
    && apt-get install --no-install-recommends --yes libgcc-s1 zlib1g \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --gid 65532 viroflash \
    && useradd --uid 65532 --gid 65532 --no-create-home --home-dir /work --shell /usr/sbin/nologin viroflash \
    && mkdir /work \
    && chown 65532:65532 /work

COPY --from=builder /src/target/release/viroflash /usr/local/bin/viroflash

USER 65532:65532
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/viroflash"]
CMD ["--help"]
