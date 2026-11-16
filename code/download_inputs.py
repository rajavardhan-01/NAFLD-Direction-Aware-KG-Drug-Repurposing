"""Step D: download every public input into data/ and record it in data/SOURCES.json.

Sources
  OptimusKG   Harvard Dataverse, doi:10.7910/DVN/IYNGEV (all files under nodes/ and edges/)
  MSigDB      2024.1.Hs KEGG legacy and Hallmark gene-symbol GMT files
  CTD         CTD_diseases_pathways.csv.gz

The combined files at the dataset root (nodes.parquet, edges.parquet and the two
largest_connected_component_* files) repeat the per-type files and are not downloaded.

Every file is checked: expected type (parquet magic bytes / GMT tab format / gzip CSV),
not empty, not an HTML error page. Dataverse MD5 checksums are verified. Nothing is
downloaded if it needs a log-in: the script stops and reports instead.

    python code/download_inputs.py            # downloads missing files, re-checks all
"""
import datetime as dt
import gzip
import hashlib
import json
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DOI = "doi:10.7910/DVN/IYNGEV"
DV = "https://dataverse.harvard.edu"
MSIGDB = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
GMTS = ["c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt", "h.all.v2024.1.Hs.symbols.gmt"]
CTD = "https://ctdbase.org/reports/CTD_diseases_pathways.csv.gz"
UA = {"User-Agent": "Paper2PublicRun-downloader/1.0"}


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def fetch(url, dest, tries=6):
    """fetch_once with retries on dropped connections (Dataverse resets long transfers)."""
    for i in range(1, tries + 1):
        try:
            return fetch_once(url, dest)
        except (ConnectionError, TimeoutError, OSError) as e:
            if i == tries:
                raise
            print(f"    retry {i}/{tries - 1} for {dest.name}: {e}", flush=True)
            time.sleep(10 * i)


def fetch_once(url, dest):
    """Download url to dest via a .partial file. Returns the final URL after redirects."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".partial")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as f:
        final = r.geturl()
        ctype = r.headers.get("Content-Type", "")
        while chunk := r.read(1 << 20):
            f.write(chunk)
    head = tmp.read_bytes()[:512].lower()
    if b"<html" in head or b"<!doctype" in head or "text/html" in ctype:
        tmp.unlink()
        sys.exit(f"STOP: {url} returned an HTML page (final URL {final}). "
                 "It may need a log-in or licence acceptance; download it by hand.")
    tmp.replace(dest)
    return final


def digest(path, algo):
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            h.update(chunk)
    return h.hexdigest()


def check_parquet(path):
    b = path.read_bytes()
    if len(b) < 12 or b[:4] != b"PAR1" or b[-4:] != b"PAR1":
        raise ValueError(f"{path.name}: not a parquet file")
    try:
        import pyarrow.parquet as pq
        md = pq.ParquetFile(path).metadata
        return {"rows": md.num_rows, "columns": md.num_columns}
    except ImportError:
        return {"rows": None, "note": "pyarrow not available for row count"}


def check_gmt(path):
    n = genes = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 3:
            raise ValueError(f"{path.name}: line {n + 1} is not GMT (name, url, genes...)")
        n += 1
        genes += len(parts) - 2
    if n == 0:
        raise ValueError(f"{path.name}: empty")
    return {"gene_sets": n, "memberships": genes}


def check_ctd(path):
    lines = data = 0
    with gzip.open(path, "rt", encoding="utf-8") as f:
        header = []
        for line in f:
            lines += 1
            if line.startswith("#"):
                header.append(line.strip())
            elif line.strip():
                data += 1
    if data == 0:
        raise ValueError(f"{path.name}: no data rows")
    version = next((h for h in header if "report created" in h.lower()), None)
    return {"lines": lines, "data_rows": data, "ctd_header": version}


def main():
    sources = []
    # ---- OptimusKG ----
    meta_url = f"{DV}/api/datasets/:persistentId/?persistentId={DOI}"
    meta = json.load(urllib.request.urlopen(urllib.request.Request(meta_url, headers=UA), timeout=120))
    v = meta["data"]["latestVersion"]
    version = f"{v['versionNumber']}.{v['versionMinorNumber']}"
    lic = v.get("license", {}).get("name") if isinstance(v.get("license"), dict) else v.get("license")
    okg = DATA / "optimuskg"
    (okg / "dataset_metadata.json").parent.mkdir(parents=True, exist_ok=True)
    (okg / "dataset_metadata.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    files = [f for f in v["files"] if f.get("directoryLabel") in ("nodes", "edges")]
    print(f"OptimusKG {DOI} version {version} ({v.get('releaseTime')}), licence {lic}: {len(files)} files")
    for f in files:
        df = f["dataFile"]
        dest = okg / f["directoryLabel"] / df["filename"]
        url = f"{DV}/api/access/datafile/{df['id']}"
        when = None
        if not dest.exists() or dest.stat().st_size != df["filesize"]:
            when = now()
            fetch(url, dest)
        md5 = df.get("md5") or (df.get("checksum") or {}).get("value")
        if md5 and digest(dest, "md5") != md5:
            raise ValueError(f"{dest.name}: MD5 mismatch with Dataverse")
        info = check_parquet(dest)
        sources.append({"source": "OptimusKG", "file": str(dest.relative_to(ROOT)).replace("\\", "/"),
                        "url": url, "doi": DOI, "version": version, "release_time": v.get("releaseTime"),
                        "licence": lic, "downloaded_utc": when, "bytes": dest.stat().st_size,
                        "md5_verified": bool(md5), "sha256": digest(dest, "sha256"), **info})
        print(f"  {dest.name:45s} {dest.stat().st_size:>12,} B  rows={info['rows']}")
    # ---- MSigDB ----
    for g in GMTS:
        dest = DATA / "msigdb" / g
        when = None
        if not dest.exists():
            when = now()
            fetch(MSIGDB + g, dest)
        info = check_gmt(dest)
        sources.append({"source": "MSigDB", "file": f"data/msigdb/{g}", "url": MSIGDB + g,
                        "version": "2024.1.Hs", "downloaded_utc": when, "bytes": dest.stat().st_size,
                        "sha256": digest(dest, "sha256"), **info})
        print(f"  {g:45s} {dest.stat().st_size:>12,} B  sets={info['gene_sets']} members={info['memberships']}")
    # ---- CTD ----
    dest = DATA / "ctd" / "CTD_diseases_pathways.csv.gz"
    when = None
    if not dest.exists():
        when = now()
        fetch(CTD, dest)
    info = check_ctd(dest)
    sources.append({"source": "CTD", "file": "data/ctd/CTD_diseases_pathways.csv.gz", "url": CTD,
                    "version": info["ctd_header"], "downloaded_utc": when, "bytes": dest.stat().st_size,
                    "sha256": digest(dest, "sha256"), **info})
    print(f"  {dest.name:45s} {dest.stat().st_size:>12,} B  rows={info['data_rows']}  {info['ctd_header']}")

    # keep first download time when re-checking
    out = DATA / "SOURCES.json"
    old = {}
    if out.exists():
        old = {s["file"]: s for s in json.loads(out.read_text(encoding="utf-8"))["files"]}
    for s in sources:
        if s["downloaded_utc"] is None:
            mtime = dt.datetime.fromtimestamp((ROOT / s["file"]).stat().st_mtime, dt.timezone.utc)
            s["downloaded_utc"] = old.get(s["file"], {}).get("downloaded_utc") or mtime.isoformat(timespec="seconds")
    out.write_text(json.dumps({"checked_utc": now(), "files": sources}, indent=2), encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)} ({len(sources)} files)")


if __name__ == "__main__":
    main()
