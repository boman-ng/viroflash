# syntax=docker/dockerfile:1
FROM debian:bookworm-slim

ARG VERSION=dev
ARG VCS_REF=unknown
LABEL org.opencontainers.image.title="viroflash" \
      org.opencontainers.image.description="Reference-group fragment evidence command-line tool" \
      org.opencontainers.image.source="https://github.com/boman-ng/viroflash" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}"

RUN mkdir /work && chown 65532:65532 /work
COPY viroflash /usr/local/bin/viroflash

USER 65532:65532
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/viroflash"]
CMD ["--help"]
