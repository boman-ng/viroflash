#!/bin/sh
# Runs only inside an ephemeral, pinned distribution base image.
set -eu
for file in /etc/apt/sources.list.d/*.sources; do
  [ -f "$file" ] || continue
  sed -i -e 's/^Types: deb$/Types: deb deb-src/' -e 's|http://|https://|g' "$file"
done
for file in /etc/apt/sources.list /etc/apt/sources.list.d/*.list; do
  [ -f "$file" ] || continue
  [ "$file" != /etc/apt/sources.list.d/viroflash-source.list ] || continue
  sed -i 's|http://|https://|g' "$file"
  sed -n 's/^deb /deb-src /p' "$file"
done > /etc/apt/sources.list.d/viroflash-source.list
dpkg-query -W '-f=${binary:Package}\t${Version}\t${source:Package}\t${source:Version}\n' > /out/base-packages.tsv
cmp /out/expected-packages.tsv /out/base-packages.tsv
apt-get -o Acquire::https::CaInfo=/run/viroflash-ca.pem -o Acquire::https::Timeout=60 -o Acquire::Retries=2 update
cut -f3,4 /out/base-packages.tsv | sort -u > /out/requested-sources.tsv
status=0
while IFS="$(printf '\t')" read -r package version; do
  mkdir -p "/out/sources/$package"
  if ! (cd "/out/sources/$package"; apt-get -o Acquire::https::CaInfo=/run/viroflash-ca.pem -o Acquire::https::Timeout=60 -o Acquire::Retries=2 source --yes --download-only --only-source "$package=$version"); then
    printf '%s\t%s\n' "$package" "$version" >> /out/unavailable-sources.tsv
    status=1
  fi
done < /out/requested-sources.tsv
exit "$status"
