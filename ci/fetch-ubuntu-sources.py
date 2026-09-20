#!/usr/bin/env python3
"""Retrieve exact retired Ubuntu source versions from official Launchpad records."""
import hashlib
import json
from pathlib import Path
import re
import sys
import urllib.parse
import urllib.request

work = Path(sys.argv[1])
missing = work / "unavailable-sources.tsv"
if not missing.exists():
    raise SystemExit(0)

def fetch(url):
    host = urllib.parse.urlsplit(url).hostname or ""
    if not (host in {"api.launchpad.net", "launchpad.net", "launchpadlibrarian.net"} or host.endswith(".launchpadlibrarian.net")):
        raise ValueError("Unexpected source host")
    with urllib.request.urlopen(url, timeout=120) as response:
        destination = urllib.parse.urlsplit(response.url)
        if destination.scheme != "https" or not (destination.hostname in {"launchpad.net", "api.launchpad.net", "launchpadlibrarian.net"} or
                (destination.hostname or "").endswith(".launchpadlibrarian.net")):
            raise ValueError("Unexpected source redirect")
        return response.read()

for row in missing.read_text().splitlines():
    package, version = row.split("\t")
    if not re.fullmatch(r"[a-z0-9][a-z0-9+.-]+", package):
        raise ValueError("Invalid source package name")
    query = urllib.parse.urlencode({"ws.op": "getPublishedSources", "source_name": package,
                                   "version": version, "exact_match": "true"})
    publications = json.loads(fetch("https://api.launchpad.net/1.0/ubuntu/+archive/primary?" + query))["entries"]
    records = [p for p in publications if p["source_package_name"] == package and p["source_package_version"] == version]
    if not records:
        raise ValueError("Exact Ubuntu source unavailable")
    urls = json.loads(fetch(records[0]["self_link"] + "?ws.op=sourceFileUrls"))
    if not urls:
        raise ValueError("Ubuntu source record has no files")
    provenance = []
    destination = work / "sources" / package
    destination.mkdir(parents=True, exist_ok=True)
    for url in urls:
        name = urllib.parse.unquote(urllib.parse.urlsplit(url).path.rsplit("/", 1)[-1])
        if Path(name).name != name or name in ("", ".", ".."):
            raise ValueError("Unsafe source filename")
        data = fetch(url)
        path = destination / name
        if path.exists() and path.read_bytes() != data:
            raise ValueError("Conflicting source content")
        path.write_bytes(data)
        provenance.append({"source": url, "file": name, "sha256": hashlib.sha256(data).hexdigest()})
    (destination / "retrieval.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print("Retrieved exact Ubuntu source:", package, version)
# package-container-sources.py independently verifies every .dsc version and checksum.
