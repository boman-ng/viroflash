#!/usr/bin/env python3
"""Validate Debian source descriptors against exact installed source versions."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile

work, output = map(Path, sys.argv[1:])
if output.exists():
    raise SystemExit("Refusing to overwrite corresponding-source bundle")
manifest = {"packages": [], "files": {}}
for row in (work / "requested-sources.tsv").read_text().splitlines():
    package, version = row.split("\t")
    descriptors = list((work / "sources" / package).glob("*.dsc"))
    if len(descriptors) != 1:
        raise SystemExit("Missing or ambiguous source descriptor")
    descriptor = descriptors[0]
    fields = {}
    key = None
    for line in descriptor.read_text().splitlines():
        if line[:1].isspace() and key:
            fields[key] += "\n" + line.strip()
        elif ": " in line:
            key, value = line.split(": ", 1)
            fields[key] = value
        elif line.endswith(":"):
            key = line[:-1]
            fields[key] = ""
    if fields.get("Source") != package or fields.get("Version") != version:
        raise SystemExit("Source version differs from installed package")
    checksums = fields.get("Checksums-Sha256", "").strip().splitlines()
    if not checksums:
        raise SystemExit("Descriptor has no SHA-256 source hashes")
    for entry in checksums:
        digest, size, name = entry.split()
        if Path(name).name != name:
            raise SystemExit("Unsafe source filename")
        data = (descriptor.parent / name).read_bytes()
        if len(data) != int(size) or hashlib.sha256(data).hexdigest() != digest:
            raise SystemExit("Source archive checksum mismatch")
    manifest["packages"].append({"package": package, "version": version,
                                 "descriptor": descriptor.relative_to(work).as_posix()})
for path in sorted((work / "sources").rglob("*")):
    if path.is_file():
        manifest["files"][path.relative_to(work).as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
with output.open("xb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
    with tarfile.open(fileobj=compressed, mode="w") as archive:
        items = {name: (work / name).read_bytes() for name in manifest["files"]}
        items["source-manifest.json"] = (json.dumps(manifest, indent=2) + "\n").encode()
        for name, data in sorted(items.items()):
            member = tarfile.TarInfo(name)
            member.size = len(data)
            member.mode = 0o644
            archive.addfile(member, io.BytesIO(data))
print(f"Packaged {len(manifest['packages'])} exact corresponding source versions")
