#!/usr/bin/env python3
"""Check repository-relative Markdown links; external URLs are reviewed separately."""
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote, urlsplit

root = Path(__file__).resolve().parent.parent
files = subprocess.check_output(["git", "ls-files", "-z", "--", "*.md"], cwd=root).decode().split("\0")
errors = []
for name in filter(None, files):
    path = root / name
    for target in re.findall(r"\[[^\]]*\]\(([^\s)]+)(?:\s+[^)]*)?\)", path.read_text()):
        link = urlsplit(target.strip("<>"))
        if link.scheme or not link.path:
            continue
        if not (path.parent / unquote(link.path)).exists():
            errors.append(f"{name}: missing link {target}")
if errors:
    raise SystemExit("\n".join(errors))
print("Repository-relative Markdown links exist")
