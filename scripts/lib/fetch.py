# Resumable range downloads. EBI gave 3 to 16 KB/s per connection from here, so files go
# over several connections, each range in its own .partN, and nothing is kept that fails its checksum.
# stdlib only: 02 runs before env/ exists.

from __future__ import annotations

import concurrent.futures as cf
import hashlib
import http.client
import os
import time
import urllib.error
import urllib.request

class DownloadError(RuntimeError):
    pass


def _request(url: str, method: str = "GET", headers: dict | None = None):
    return urllib.request.urlopen(urllib.request.Request(url, method=method, headers=headers or {}), timeout=180)


def remote_size(url: str, tries: int = 6) -> tuple[int, str]:
    last = None
    for i in range(tries):
        try:
            with _request(url, "HEAD") as r:
                return int(r.headers["Content-Length"]), r.headers.get("Last-Modified", "")
        except (OSError, http.client.HTTPException, KeyError) as e:
            last = e
            time.sleep(min(60, 5 * 2**i))
    raise DownloadError(f"HEAD {url} failed after {tries} tries: {last}")


def _fetch_range(url: str, part: str, start: int, end: int, tries: int) -> None:
    # end is inclusive, as in the HTTP Range header
    want = end - start + 1
    for i in range(tries):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have == want:
            return
        if have > want:
            os.remove(part)
            have = 0
        try:
            with _request(url, headers={"Range": f"bytes={start + have}-{end}"}) as r:
                if r.status != 206:
                    raise DownloadError(f"server ignored the range request for {url} (status {r.status})")
                with open(part, "ab") as fh:
                    while True:
                        buf = r.read(1 << 18)
                        if not buf:
                            break
                        fh.write(buf)
        except (OSError, http.client.HTTPException) as e:
            if isinstance(e, urllib.error.HTTPError) and e.code in (403, 404, 410):
                raise DownloadError(f"{url}: HTTP {e.code}") from e
            time.sleep(min(120, 5 * 2**min(i, 5)))
    have = os.path.getsize(part) if os.path.exists(part) else 0
    if have != want:
        raise DownloadError(f"range {start}-{end} of {url} incomplete after {tries} tries ({have} of {want} bytes)")


def file_hash(path: str, algorithm: str = "sha256") -> str:
    h = hashlib.new(algorithm)
    with open(path, "rb") as fh:
        for buf in iter(lambda: fh.read(1 << 22), b""):
            h.update(buf)
    return h.hexdigest()


def download(url: str, dest: str, size: int | None = None, sha256: str | None = None,
             md5: str | None = None, connections: int = 4, tries: int = 30) -> dict:
    if os.path.exists(dest):
        got = file_hash(dest, "sha256" if sha256 else "md5") if (sha256 or md5) else None
        if got is None or got == (sha256 or md5):
            return {"size": os.path.getsize(dest), "sha256": got if sha256 else None, "seconds": 0.0, "resumed": True}
        os.remove(dest)
    if size is None:
        size, _ = remote_size(url)
    t0 = time.monotonic()
    n = 1 if size < (8 << 20) else max(1, connections)
    bounds = [(k * size // n, (k + 1) * size // n - 1) for k in range(n)]
    parts = [f"{dest}.part{k}" for k in range(n)]
    with cf.ThreadPoolExecutor(max_workers=n) as pool:
        futures = [pool.submit(_fetch_range, url, p, s, e, tries) for p, (s, e) in zip(parts, bounds)]
        for f in futures:
            f.result()
    tmp = dest + ".tmp"
    with open(tmp, "wb") as out:
        for p in parts:
            with open(p, "rb") as fh:
                for buf in iter(lambda: fh.read(1 << 22), b""):
                    out.write(buf)
    for p in parts:
        os.remove(p)
    assembled = os.path.getsize(tmp)
    if assembled != size:
        os.remove(tmp)
        raise DownloadError(f"{url}: assembled {assembled} bytes, expected {size}")
    result = {"size": size, "seconds": time.monotonic() - t0, "resumed": False}
    if sha256:
        got = file_hash(tmp, "sha256")
        if got != sha256:
            os.remove(tmp)
            raise DownloadError(f"{url}: sha256 {got} does not match the published {sha256}; file deleted")
        result["sha256"] = got
    if md5:
        got = file_hash(tmp, "md5")
        if got != md5:
            os.remove(tmp)
            raise DownloadError(f"{url}: md5 {got} does not match the published {md5}; file deleted")
        result["md5"] = got
    os.replace(tmp, dest)
    return result
