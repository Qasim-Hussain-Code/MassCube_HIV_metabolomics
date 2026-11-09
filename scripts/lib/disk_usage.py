"""Bytes used under each path given, printed one total per line.

Used instead of du for the checks made between stages and before each
download. On Windows, os.scandir returns file sizes from the directory
listing itself, while MSYS du opens every file; over the conda
environment's tens of thousands of files that is the difference between
seconds and half an hour with the virus scanner running. Hard links are
counted once by (device, inode) where the platform reports them; Windows
directory listings do not, so the conda package cache, whose files are
hard links into env/, is reported but expected to be emptied by the
install stage.

Usage: python disk_usage.py PATH [PATH ...]
"""

from __future__ import annotations

import os
import sys


def tree_bytes(top: str) -> int:
    total = 0
    seen: set[tuple[int, int]] = set()
    stack = [top]
    while stack:
        d = stack.pop()
        try:
            it = os.scandir(d)
        except OSError:  # locked or vanished between listing and opening
            continue
        with it:
            for e in it:
                try:
                    if e.is_dir(follow_symlinks=False):
                        stack.append(e.path)
                    elif e.is_file(follow_symlinks=False):
                        st = e.stat(follow_symlinks=False)
                        key = (st.st_dev, st.st_ino)
                        if st.st_ino and key in seen:
                            continue
                        if st.st_ino:
                            seen.add(key)
                        total += st.st_size
                except OSError:
                    continue
    return total


if __name__ == "__main__":
    for p in sys.argv[1:]:
        print(tree_bytes(p) if os.path.isdir(p) else (os.path.getsize(p) if os.path.exists(p) else 0))
