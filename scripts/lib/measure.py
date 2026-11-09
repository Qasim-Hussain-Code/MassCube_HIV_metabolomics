"""Run one command and record its elapsed time and peak resident memory.

GNU time is not available on Windows, so two figures are taken and the
larger is reported. The command's own peak is exact: on Windows the peak
working set is read from its process handle after it exits, elsewhere the
children's maximum resident set comes from getrusage. The whole process
tree is also sampled, every 50 ms for the first two seconds and every half
second after that, summing the resident set sizes of the command and all
its descendants, because MassCube parallelises through loky worker
processes that the command's own figure would not include.

Usage: python measure.py --label NAME --log logs/measure.tsv -- command [args...]

Appends one row to the log and exits with the command's exit code.
"""

from __future__ import annotations

import argparse
import datetime
import os
import subprocess
import sys
import time

import psutil


def tree_rss(proc: psutil.Process) -> int:
    try:
        members = [proc] + proc.children(recursive=True)
    except psutil.NoSuchProcess:
        return 0
    total = 0
    for p in members:
        try:
            total += p.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return total


def own_peak_bytes(child: subprocess.Popen) -> int:
    """Peak resident memory of the command itself, after it has exited."""
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes

        class Counters(ctypes.Structure):
            _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                        ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                        ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]

        c = Counters()
        c.cb = ctypes.sizeof(c)
        # The Popen object still holds the exited process's handle.
        ok = ctypes.WinDLL("psapi").GetProcessMemoryInfo(wintypes.HANDLE(int(child._handle)), ctypes.byref(c), c.cb)
        return int(c.PeakWorkingSetSize) if ok else 0
    import resource
    kb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return kb if sys.platform == "darwin" else kb * 1024  # macOS reports bytes, Linux kilobytes


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--log", required=True)
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd

    started = datetime.datetime.now().astimezone()
    t0 = time.monotonic()
    child = subprocess.Popen(cmd)
    ps = psutil.Process(child.pid)
    peak_tree = 0
    while child.poll() is None:
        peak_tree = max(peak_tree, tree_rss(ps))
        time.sleep(0.05 if time.monotonic() - t0 < 2 else 0.5)
    elapsed = time.monotonic() - t0
    peak_mb = max(peak_tree, own_peak_bytes(child)) / 1024**2

    new = not os.path.exists(a.log)
    os.makedirs(os.path.dirname(os.path.abspath(a.log)), exist_ok=True)
    with open(a.log, "a", encoding="utf-8") as fh:
        if new:
            fh.write("label\tstarted\telapsed_s\tpeak_rss_mb\texit_code\n")
        fh.write(f"{a.label}\t{started:%Y-%m-%dT%H:%M:%S%z}\t{elapsed:.1f}\t{peak_mb:.0f}\t{child.returncode}\n")
    return child.returncode


if __name__ == "__main__":
    sys.exit(main())
