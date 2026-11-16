#!/usr/bin/env python3
# du -sb replacement: scandir gets sizes from the listing, MSYS du opens every file (1.7 s vs 40+ min over env/)

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
                        # hard links once; Windows listings report no inode, hence 01 empties the conda cache
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
