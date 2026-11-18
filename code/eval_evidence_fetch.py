"""PubMed evidence screen for the judged pool (union of the top-20 drugs of all methods).

One fixed query per drug (drug name in title/abstract AND a liver-disease term), the same for every drug and
blind to the method and rank that selected it. Writes results/Evaluation/evidence/pubmed_hits.json (per drug: counts
and up to 8 records with title, year, publication types, MeSH flags and the abstract) for manual grading.

    python code/eval_evidence_fetch.py
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/Evaluation/evidence"
OUT.mkdir(parents=True, exist_ok=True)
UA = {"User-Agent": "Paper2PublicRun-evidence-screen/1.0"}
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DISEASE = ('(NAFLD[tiab] OR NASH[tiab] OR MASLD[tiab] OR MASH[tiab] OR "fatty liver"[tiab] OR steatohepatitis[tiab] OR '
           '"liver fibrosis"[tiab] OR "hepatic fibrosis"[tiab] OR "liver steatosis"[tiab] OR "hepatic steatosis"[tiab])')
SALT = re.compile(r"\b(DISODIUM|TETRASODIUM|SODIUM|HYDROCHLORIDE|HYDROBROMIDE|CITRATE|PHOSPHATE|CILEXETIL)\b")
MANUAL = {"F19 131I": "F19", "CHOLINE FENOFIBRATE": "fenofibrate", "PIRINIXIC ACID": "pirinixic acid OR WY-14643",
          "GW501516": "GW501516 OR GW0742 OR cardarine", "AVE8134": "AVE8134", "ASK-8007": "ASK-8007",
          "NITRIC OXIDE": "\"nitric oxide\"", "CHOLESTEROL": "cholesterol", "BARDOXOLONE METHYL": "bardoxolone", "EFATUTAZONE": "efatutazone OR CS-7017", "BALAGLITAZONE": "balaglitazone", "CIGLITAZONE": "ciglitazone", "RIVOGLITAZONE": "rivoglitazone", "NETOGLITAZONE": "netoglitazone", "LERIGLITAZONE": "leriglitazone OR MIN-102", "CLOFIBRATE": "clofibrate", "GOLIMUMAB": "golimumab", "CERTOLIZUMAB PEGOL": "certolizumab", "INT131": "INT131", "ARHALOFENATE": "arhalofenate"}


def get(url, data=None):
    for k in range(4):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, data=data, headers=UA), timeout=60).read()
        except Exception as e:  # noqa: BLE001
            time.sleep(2 + 2 * k)
    raise RuntimeError(url)


def parent(name):
    return re.sub(r"\s+", " ", SALT.sub("", name)).strip()


def term(par):
    t = MANUAL.get(par, par.lower())
    return t if " OR " in t else f'"{t}"'


def main():
    r = pd.read_csv(ROOT / "results/Evaluation/main/all_rankings.csv")
    main_m = ["SP", "NP", "RWR", "MPC", "DWPC", "HeteSim", "Proposed"]
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    names = sorted(r[r.method.isin(main_m) & (r["rank"] <= K)].drug.unique())
    tag = f"top{K}"
    if len(sys.argv) > 2:                      # extra pool (drugs tied at the cut), a CSV with a `drug` column
        names = sorted(pd.read_csv(sys.argv[2]).drug.unique())
        tag = Path(sys.argv[2]).stem
    pars = sorted({parent(n) for n in names})
    res = {}
    for p in pars:
        tm = term(p)
        q = f"(({tm})[tiab]) AND {DISEASE}" if " OR " not in tm else f"(({' OR '.join(f'({x.strip()}[tiab])' for x in tm.split(' OR '))})) AND {DISEASE}"
        j = json.loads(get(EU + "esearch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "term": q, "retmode": "json", "retmax": 300, "sort": "relevance"})))
        ids = j["esearchresult"]["idlist"]
        n = int(j["esearchresult"]["count"])
        recs = []
        if ids:
            x = get(EU + "efetch.fcgi", urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids), "retmode": "xml"}).encode())
            root = ET.fromstring(x)
            for a in root.findall(".//PubmedArticle"):
                pm = a.findtext(".//PMID")
                title = "".join(a.find(".//ArticleTitle").itertext()) if a.find(".//ArticleTitle") is not None else ""
                ab = " ".join("".join(t.itertext()) for t in a.findall(".//Abstract/AbstractText"))
                yr = a.findtext(".//PubDate/Year") or (a.findtext(".//PubDate/MedlineDate") or "")[:4]
                pt = [t.text for t in a.findall(".//PublicationType")]
                mesh = [m.text for m in a.findall(".//MeshHeading/DescriptorName")]
                recs.append({"pmid": pm, "year": yr, "title": title, "pubtypes": pt, "humans": "Humans" in mesh,
                             "animals": "Animals" in mesh, "abstract": ab})
        # drug must really appear in the title/abstract (the tiab match is already required; keep recs as returned)
        res[p] = {"query": q, "count": n, "returned": len(recs), "records": recs}
        print(p, n, len(recs), flush=True)
        time.sleep(0.4)
    (OUT / f"pubmed_hits_{tag}.json").write_text(json.dumps({"retrieved": time.strftime("%Y-%m-%d"), "parents": res}, indent=1), encoding="utf-8")
    pd.DataFrame([{"drug": p, "n_hits": v["count"]} for p, v in res.items()]).to_csv(OUT / f"pubmed_counts_{tag}.csv", index=False)


if __name__ == "__main__":
    main()
