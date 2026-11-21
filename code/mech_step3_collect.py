"""Step 3: collect independent experimental publications for every seed-gene drug-target pair of the 116 MPC drugs.

For each pair two PubMed queries are run (public E-utilities, read-only):
  (a) mechanism  : drug terms AND target terms
  (b) disease    : (a) AND a liver-disease term
"Independent" = primary research articles from PubMed, not the pharmacology databases that built the graph
(OptimusKG, DrugBank, ChEMBL, DGIdb, GtoPdb, Open Targets) and not the source signature paper; reviews, editorials,
comments, letters and retracted papers are dropped in step 4. The same queries are used for every group of drugs.

Writes results/Evaluation/mechanistic/step3_pubmed_records.json  (query, counts and records per pair).

    python code/mech_step3_collect.py
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
OUT = ROOT / "results/Evaluation/mechanistic"
UA = {"User-Agent": "Paper2PublicRun-mechanistic-screen/1.0"}
EU = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DISEASE = ('(NAFLD[tiab] OR NASH[tiab] OR MASLD[tiab] OR MASH[tiab] OR "fatty liver"[tiab] OR steatohepatitis[tiab] OR '
           '"liver fibrosis"[tiab] OR "hepatic fibrosis"[tiab] OR "hepatic steatosis"[tiab] OR "liver steatosis"[tiab])')
SALT = re.compile(r"\b(DISODIUM|TETRASODIUM|SODIUM|HYDROCHLORIDE|HYDROBROMIDE|CITRATE|PHOSPHATE|CILEXETIL|ACETATE|ISOPROPYL|METHYL)\b")
DRUG_ALIAS = {
    "BARDOXOLONE METHYL": ["bardoxolone", "CDDO-Me"], "BARDOXOLONE": ["bardoxolone", "CDDO-Me"],
    "CHOLINE FENOFIBRATE": ["fenofibrate", "fenofibric acid"], "PIRINIXIC ACID": ["pirinixic acid", "WY-14643"],
    "GW501516": ["GW501516", "GW0742", "cardarine"], "COLLAGENASE CLOSTRIDIUM HISTOLYTICUM": ["collagenase clostridium histolyticum", "clostridial collagenase"],
    "LANATOSIDE C": ["lanatoside C"], "ERITORAN TETRASODIUM": ["eritoran"], "ERITORAN": ["eritoran"], "F19 131I": ["F19"],
    "EFATUTAZONE HYDROCHLORIDE": ["efatutazone", "CS-7017"], "EFATUTAZONE": ["efatutazone", "CS-7017"], "LERIGLITAZONE": ["leriglitazone", "MIN-102"],
    "CERTOLIZUMAB PEGOL": ["certolizumab"], "LAPAQUISTAT ACETATE": ["lapaquistat", "TAK-475"], "TUROFEXORATE ISOPROPYL": ["turofexorate", "WAY-362450"],
    "ABITUZUMAB": ["abituzumab"], "OLSALAZINE SODIUM": ["olsalazine"], "BALSALAZIDE DISODIUM": ["balsalazide"],
}
GENE_ALIAS = {
    "PPARG": ["PPARG", "PPAR-gamma", "PPARgamma", "PPARγ", "peroxisome proliferator-activated receptor gamma"],
    "PPARA": ["PPARA", "PPAR-alpha", "PPARalpha", "PPARα", "peroxisome proliferator-activated receptor alpha"],
    "TNF": ["TNF", "tumor necrosis factor", "TNF-alpha", "TNFα"], "IL6": ["IL6", "IL-6", "interleukin-6"],
    "ABCB11": ["ABCB11", "BSEP", "bile salt export pump"], "MERTK": ["MERTK", "MER tyrosine kinase", "Mer kinase"],
    "TGFB1": ["TGFB1", "TGF-beta", "TGF-β", "transforming growth factor"], "ATP1A1": ["ATP1A1", "Na+/K+-ATPase", "Na,K-ATPase", "sodium-potassium ATPase"],
    "FXYD2": ["FXYD2", "Na+/K+-ATPase", "Na,K-ATPase"], "CYP2C19": ["CYP2C19"], "IL1B": ["IL1B", "IL-1beta", "IL-1β", "interleukin-1 beta"],
    "ITGB5": ["ITGB5", "integrin beta5", "alphavbeta5", "integrin"], "FAP": ["fibroblast activation protein", "FAP"],
    "NR1I2": ["NR1I2", "PXR", "pregnane X receptor"], "COL6A2": ["collagen"], "COL1A2": ["collagen"], "COL3A1": ["collagen"], "COL5A1": ["collagen"],
    "IL1A": ["IL1A", "IL-1alpha", "IL-1α", "interleukin-1 alpha"], "TLR4": ["TLR4", "toll-like receptor 4"],
    "TOP2A": ["TOP2A", "topoisomerase II", "topoisomerase 2"], "ABCC4": ["ABCC4", "MRP4"], "FDFT1": ["FDFT1", "squalene synthase"],
    "CYP2E1": ["CYP2E1"], "TXN": ["thioredoxin", "TXN"], "ADRB3": ["ADRB3", "beta3-adrenergic", "beta-3 adrenergic", "β3-adrenergic"],
    "AGTR1": ["AGTR1", "AT1 receptor", "angiotensin II type 1"], "SPP1": ["osteopontin", "SPP1"], "CA12": ["carbonic anhydrase XII", "carbonic anhydrase 12", "carbonic anhydrase"],
    "CHUK": ["CHUK", "IKKalpha", "IKK-alpha", "IKKα"], "LAMA3": ["laminin"], "LAMA2": ["laminin"], "NR1H4": ["NR1H4", "FXR", "farnesoid X receptor"],
    "CNR1": ["CNR1", "CB1 receptor", "cannabinoid receptor 1"], "MARCKS": ["MARCKS"], "LEPR": ["LEPR", "leptin receptor"],
    "CXCL12": ["CXCL12", "SDF-1", "stromal cell-derived factor"], "IRS1": ["IRS1", "insulin receptor substrate"], "FES": ["FES kinase", "FES tyrosine kinase"],
    "TXNRD1": ["TXNRD1", "thioredoxin reductase"],
}


def get(url, data=None):
    for k in range(5):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, data=data, headers=UA), timeout=60).read()
        except Exception:  # noqa: BLE001
            time.sleep(2 + 3 * k)
    raise RuntimeError(url)


def drug_terms(name):
    if name in DRUG_ALIAS:
        return DRUG_ALIAS[name]
    base = re.sub(r"\s+", " ", SALT.sub("", name)).strip()
    return [base.lower() if " " in base or base.isalpha() else base]


def tiab(terms):
    return "(" + " OR ".join(f'"{t}"[tiab]' for t in terms) + ")"


def esearch(q, retmax):
    j = json.loads(get(EU + "esearch.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "term": q, "retmode": "json", "retmax": retmax, "sort": "relevance"})))
    return int(j["esearchresult"]["count"]), j["esearchresult"]["idlist"]


def efetch(ids):
    out = {}
    for i in range(0, len(ids), 100):
        x = get(EU + "efetch.fcgi", urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids[i:i + 100]), "retmode": "xml"}).encode())
        for a in ET.fromstring(x).findall(".//PubmedArticle"):
            pm = a.findtext(".//PMID")
            t = a.find(".//ArticleTitle")
            out[pm] = {"pmid": pm, "year": a.findtext(".//PubDate/Year") or (a.findtext(".//PubDate/MedlineDate") or "")[:4],
                       "title": "".join(t.itertext()) if t is not None else "",
                       "abstract": " ".join("".join(x_.itertext()) for x_ in a.findall(".//Abstract/AbstractText")),
                       "pubtypes": [p.text for p in a.findall(".//PublicationType")]}
        time.sleep(0.4)
    return out


def main():
    pairs = pd.read_csv(OUT / "step2_drug_target_pairs.csv")
    pairs = pairs[pairs.is_seed_gene].drop_duplicates(["Drug_ID", "Gene"])
    res = {}
    for n, r in enumerate(pairs.itertuples(), 1):
        dt, gt = drug_terms(r.Drug), GENE_ALIAS.get(r.Gene, [r.Gene])
        qa = f"{tiab(dt)} AND {tiab(gt)}"
        qb = f"{qa} AND {DISEASE}"
        ca, ia = esearch(qa, 40)
        cb, ib = esearch(qb, 40) if ca else (0, [])
        recs = efetch(list(dict.fromkeys(ia + ib))) if ia else {}
        res[f"{r.Drug_ID}|{r.Gene}"] = {"drug": r.Drug, "gene": r.Gene, "drug_terms": dt, "gene_terms": gt, "query": qa, "count": ca,
                                        "count_with_disease": cb, "ids_disease": ib, "records": list(recs.values())}
        print(n, len(pairs), r.Drug, r.Gene, ca, cb, flush=True)
        time.sleep(0.4)
    (OUT / "step3_pubmed_records.json").write_text(json.dumps({"retrieved": time.strftime("%Y-%m-%d"), "pairs": res}, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
