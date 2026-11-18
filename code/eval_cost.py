"""Runtime and peak memory of every method, each in a fresh process (Windows working-set counters).

Reported memory = peak working set of the process minus the working set right after the shared data were loaded
(the extra memory the method needs), runtime = wall time of the method including its own setup. Median of 3 runs.

    python code/eval_cost.py
"""
import ctypes
import gc
import json
import statistics
import subprocess
import sys
import time
from ctypes import wintypes
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
METHODS = ["SP", "NP", "RWR", "MPC", "DWPC", "HeteSim", "Proposed"]


class PMC(ctypes.Structure):
    _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t), ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]


def mem():
    k = ctypes.windll.kernel32
    k.GetCurrentProcess.restype = wintypes.HANDLE
    p = PMC()
    p.cb = ctypes.sizeof(p)
    ctypes.windll.psapi.GetProcessMemoryInfo(wintypes.HANDLE(k.GetCurrentProcess()), ctypes.byref(p), p.cb)
    return p.WorkingSetSize, p.PeakWorkingSetSize


def child(m):
    sys.path.insert(0, str(Path(__file__).parent))
    from evalcore import (Data, Routes, adjacency, hetesim_scores, nearest_seed_excl_self, np_scores, path_scores, pathway_links,
                          rwr_scores, sp_scores, start_links)
    D = Data()
    A = adjacency(D)
    seeds_i = D.seeds.i.to_numpy()
    linked = pathway_links(D, D.seeds)
    gc.collect()
    ws0, _ = mem()
    t0 = time.perf_counter()
    if m == "SP":
        sp_scores(D, nearest_seed_excl_self(A, seeds_i))
    elif m == "NP":
        np_scores(D, nearest_seed_excl_self(A, seeds_i), A)
    elif m == "RWR":
        rwr_scores(D, A, seeds_i, r=0.5)
    else:
        R = Routes(D, A, seeds_i, linked)
        if m == "MPC":
            path_scores(D, R, start_links(D, seeds_i, "none"))
        elif m == "DWPC":
            path_scores(D, R, start_links(D, seeds_i, "none"), alpha=0.4)
        elif m == "HeteSim":
            hetesim_scores(D, R, start_links(D, seeds_i, "none"))
        else:
            path_scores(D, R, start_links(D, seeds_i, "proposed"))
    dt = time.perf_counter() - t0
    _, peak = mem()
    print(json.dumps({"method": m, "runtime_s": dt, "extra_peak_MB": max(peak - ws0, 0) / 1e6}))


def main():
    rows = []
    for m in METHODS:
        runs = [json.loads(subprocess.run([sys.executable, "-I", __file__, "--child", m], capture_output=True, text=True, check=True).stdout.strip().splitlines()[-1])
                for _ in range(3)]
        rows.append({"method": m, "runtime_s": round(statistics.median(r["runtime_s"] for r in runs), 3),
                     "extra_peak_MB": round(statistics.median(r["extra_peak_MB"] for r in runs), 1)})
        print(rows[-1], flush=True)
    import pandas as pd
    out = ROOT / "results/Evaluation/main"
    pd.DataFrame(rows).to_csv(out / "table_B_cost.csv", index=False)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--child":
        child(sys.argv[2])
    else:
        main()
