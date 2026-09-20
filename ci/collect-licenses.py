#!/usr/bin/env python3
"""Copy upstream notices without rewriting them; record their exact provenance."""
import argparse
import fnmatch
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile
import urllib.request


def sha(data):
    return hashlib.sha256(data).hexdigest()


def metadata(filtered=False, offline=True):
    command = ["cargo", "metadata", "--locked", "--format-version", "1"]
    if offline:
        command += ["--offline"]
    if filtered:
        command += ["--filter-platform", "x86_64-unknown-linux-musl"]
    return json.loads(subprocess.check_output(command))


def release_packages(meta):
    nodes = {n["id"]: n for n in meta["resolve"]["nodes"]}
    result, todo = set(), [meta["resolve"]["root"]]
    while todo:
        key = todo.pop()
        if key in result:
            continue
        result.add(key)
        todo.extend(d["pkg"] for d in nodes[key]["deps"]
                    if any(k["kind"] != "dev" for k in d["dep_kinds"]))
    return result


def generate(output):
    import tomllib  # Python 3.11+ is needed only when regenerating the inventory.
    meta = metadata()
    active = release_packages(metadata(True))
    lock = tomllib.loads(Path("Cargo.lock").read_text())
    checksums = {(p["name"], p["version"]): p.get("checksum") for p in lock["package"]}
    manifest = {"cargo_lock_sha256": sha(Path("Cargo.lock").read_bytes()),
                "target": "x86_64-unknown-linux-musl", "packages": [], "runtime": []}
    output.mkdir(parents=True, exist_ok=True)

    def copy(source, relative):
        data = source.read_bytes()
        data.decode("utf-8")  # A binary named LICENSE is not a reviewed notice.
        destination = output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(data)
        return {"path": relative.as_posix(), "sha256": sha(data)}

    for package in sorted(meta["packages"], key=lambda p: (p["name"], p["version"])):
        if package["id"] == meta["resolve"]["root"]:
            continue
        root = Path(package["manifest_path"]).parent
        cached = root.parent.parent.parent / "cache" / root.parent.name / (root.name + ".crate")
        if sha(cached.read_bytes()) != checksums[(package["name"], package["version"])]:
            raise ValueError("Crate archive does not match Cargo.lock")
        with tarfile.open(cached) as archive:
            for member in archive:
                if member.isfile():
                    relative = Path(member.name).relative_to(root.name)
                    if (root / relative).read_bytes() != archive.extractfile(member).read():
                        raise ValueError(f"Locally modified dependency: {package['name']}")
        prefix = Path("crates") / f"{package['name']}-{package['version']}"
        files = []
        for source in sorted(root.rglob("*")):
            if not source.is_file() or not any(fnmatch.fnmatch(source.name.lower(), p)
                    for p in ("license*", "copying*", "copyright*", "notice*")):
                continue
            if not source.resolve().is_relative_to(root.resolve()):
                raise ValueError("License symlink escapes crate")
            files.append(copy(source, prefix / source.relative_to(root)))
        # Several minimap2 headers carry notices absent from the root license.
        if package["name"] == "minimap2-sys":
            for source in sorted((root / "minimap2").glob("*")):
                if source.suffix not in (".h", ".c"):
                    continue
                text = source.read_text()
                match = re.match(r"\s*(/\*.*?\*/)", text, re.S)
                if match and re.search(r"copyright|permission|license", match[1], re.I):
                    relative = prefix / "headers" / (source.name + ".txt")
                    target = output / relative
                    target.parent.mkdir(parents=True, exist_ok=True)
                    data = (match[1] + "\n").encode()
                    target.write_bytes(data)
                    files.append({"path": relative.as_posix(), "sha256": sha(data),
                                  "source_file": source.relative_to(root).as_posix(),
                                  "source_sha256": sha(source.read_bytes())})
        if package["id"] in active and not files:
            raise ValueError(f"Missing release dependency notice: {package['name']}")
        manifest["packages"].append({
            "name": package["name"], "version": package["version"],
            "source": package["source"], "crate_sha256": checksums[(package["name"], package["version"])],
            "repository": package["repository"], "declared_license": package["license"],
            "scope": "release-build" if package["id"] in active else "test-or-other-platform",
            "modified": False, "notices": files,
        })
    sysroot = Path(subprocess.check_output(["rustc", "--print", "sysroot"], text=True).strip())
    rust_version = subprocess.check_output(["rustc", "--version"], text=True).strip()
    source = sysroot / "share/doc/rust/COPYRIGHT-library.html"
    manifest["runtime"].append({"name": "Rust standard library (including bundled third-party notices)",
        "version": rust_version, "source": "https://github.com/rust-lang/rust/tree/1.94.1/library",
        "modified": False, "notices": [copy(source, Path("runtime/Rust-COPYRIGHT-library.html"))]})
    for item in json.loads(Path("ci/runtime-licenses.lock.json").read_text()):
        target = output / item["path"]
        if not target.exists():
            with urllib.request.urlopen(item["source"], timeout=60) as response:
                target.write_bytes(response.read())
        if sha(target.read_bytes()) != item["sha256"]:
            raise ValueError("Runtime license checksum mismatch")
        manifest["runtime"].append({k: item[k] for k in ("name", "version", "source")} | {
            "modified": False, "notices": [{k: item[k] for k in ("path", "sha256")}]})
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(f"Collected notices for {len(manifest['packages'])} locked packages")


def check(output):
    data = json.loads((output / "manifest.json").read_text())
    if data["cargo_lock_sha256"] != sha(Path("Cargo.lock").read_bytes()):
        raise ValueError("Lockfile changed: re-review and regenerate notices")
    for component in data["packages"] + data["runtime"]:
        if component.get("scope") == "release-build" and not component["notices"]:
            raise ValueError("Release dependency has no notice")
        for item in component["notices"]:
            path = output / item["path"]
            if not path.resolve().is_relative_to(output.resolve()) or sha(path.read_bytes()) != item["sha256"]:
                raise ValueError("License notice missing or altered")
    actual = {(p["name"], p["version"]) for p in metadata(offline=False)["packages"] if p["source"]}
    if actual != {(p["name"], p["version"]) for p in data["packages"]}:
        raise ValueError("License inventory does not cover Cargo.lock")
    print("License inventory matches Cargo.lock and original notice checksums")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("licenses"))
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    (check if args.check else generate)(args.output)
