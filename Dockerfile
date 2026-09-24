FROM debian:bookworm-slim@sha256:3783cc01769c7b2b1b83a5c5ad96c815348e28ed7da68e2e3687004faa906251

ARG VERSION=dev
ARG VCS_REF=unknown
ARG SOURCE_URL=https://github.com/boman-ng/viroflash
LABEL org.opencontainers.image.title="viroflash" \
      org.opencontainers.image.description="Reference-group fragment evidence command-line tool" \
      org.opencontainers.image.source="${SOURCE_URL}" \
      org.opencontainers.image.licenses="BSD-3-Clause" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.revision="${VCS_REF}"

RUN mkdir /work && chown 65532:65532 /work
COPY viroflash /usr/local/bin/viroflash
COPY LICENSE THIRD_PARTY_NOTICES.md /usr/share/viroflash/
COPY licenses /usr/share/viroflash/licenses/
COPY examples /usr/share/viroflash/examples/

USER 65532:65532
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/viroflash"]
CMD ["--help"]
