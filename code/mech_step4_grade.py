"""Step 4: grade drug-target-disease evidence with the M0-M4 framework (rule-based, from step 3 records).

M0  no independent primary publication that mentions both the drug and the target
M1  mention only: primary publications mention both, but none reports an experiment on the drug at the target
M2  experimental drug-target engagement: binding / enzyme / receptor / reporter assay, IC50, EC50, Ki, recombinant protein
M3  experimental pharmacological action in cells, animals or humans, consistent with the annotated action
    (if the annotated direction is unknown or mixed, any reported action counts)
M4  M3 plus the same publication links the drug (through this target) to NAFLD / NASH / steatosis / liver fibrosis

Contradictory evidence (flags, not a grade): a primary record that reports the OPPOSITE action to the annotated one without
the matching action, or a negative-result statement (no effect / did not inhibit) for the drug at the target.
Reviews, editorials, comments, letters, news, meta-analyses and retracted papers are not counted.
The rules are keyword based; step 4 writes the record-level flags so every grade can be audited.

    python code/mech_step4_grade.py
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/Evaluation/mechanistic"
sys.path.insert(0, str(Path(__file__).parent))
from mech_step3_collect import DISEASE  # noqa: E402,F401

NOT_PRIMARY = {"Review", "Editorial", "Comment", "Letter", "News", "Meta-Analysis", "Systematic Review", "Retracted Publication",
               "Retraction of Publication", "Published Erratum", "Biography", "Interview", "Guideline", "Practice Guideline"}
IN_VITRO = re.compile(r"\b(binding|binds|bound|ic50|ec50|\bki\b|\bkd\b|enzym\w*|recombinant|purified|reporter|luciferase|transactivation|"
                      r"ligand|co-crystal|crystal structure|radioligand|affinity|biochemical|cell-free|in vitro assay|kinase assay)\b", re.I)
FUNCTIONAL = re.compile(r"\b(mice|mouse|murine|rats?|zebrafish|in vivo|knockout|knock-out|knockdown|sirna|shrna|xenograft|hepatocytes|hepg2|huh-?7|"
                        r"macrophages|fibroblasts|stellate|cultured|transfect\w*|cell lines?|primary cells|human subjects|healthy volunteers)\b", re.I)
DOWN = re.compile(r"\b(inhibit\w*|antagon\w*|block\w*|suppress\w*|reduc\w*|downregulat\w*|down-regulat\w*|neutraliz\w*|decreas\w*|"
                  r"attenuat\w*|degrad\w*|cleav\w*|depriv\w*|deplet\w*)\b", re.I)
UP = re.compile(r"\b(agoni\w*|activat\w*|induc\w*|upregulat\w*|up-regulat\w*|enhanc\w*|stimulat\w*|increas\w*|potentiat\w*|ligand-activated)\b", re.I)
NEG = re.compile(r"\b(no effect|did not (?:inhibit|affect|alter|activate|reduce|change|bind)|failed to|not significantly|no significant effect|"
                 r"without effect|had no|does not (?:inhibit|affect|bind|activate))\b", re.I)
DIS = re.compile(r"\b(nafld|nash|masld|mash|fatty liver|steatohepatitis|hepatic steatosis|liver steatosis|liver fibrosis|hepatic fibrosis)\b", re.I)


def mentions(text, terms):
    t = text.lower()
    return any(x.lower() in t for x in terms)


def main():
    recs = json.loads((OUT / "step3_pubmed_records.json").read_text(encoding="utf-8"))["pairs"]
    pairs = pd.read_csv(OUT / "step2_drug_target_pairs.csv")
    pairs = pairs[pairs.is_seed_gene].drop_duplicates(["Drug_ID", "Gene"]).set_index(["Drug_ID", "Gene"])
    rows, rec_rows = [], []
    for key, v in recs.items():
        did, gene = key.split("|")
        pr = pairs.loc[(did, gene)]
        ann = pr.action_direction
        n = dict(primary=0, mention_only=0, in_vitro=0, functional_match=0, opposite_only=0, negative=0, disease=0, m4_records=0)
        best = {"in_vitro": [], "functional_match": [], "m4": [], "opposite": [], "negative": []}
        for r in v["records"]:
            if set(r["pubtypes"]) & NOT_PRIMARY:
                continue
            txt = f"{r['title']}. {r['abstract']}"
            if not (mentions(txt, v["drug_terms"]) and mentions(txt, v["gene_terms"])):
                continue
            n["primary"] += 1
            sents = [x for x in re.split(r"(?<=[.;!?])\s+", txt) if mentions(x, v["drug_terms"]) and mentions(x, v["gene_terms"])]
            ctx = " ".join(sents)                        # only sentences that name BOTH the drug and the target
            iv = bool(IN_VITRO.search(ctx))
            fn = bool(FUNCTIONAL.search(txt)) and bool(sents)
            has_down, has_up = bool(DOWN.search(ctx)), bool(UP.search(ctx))
            if ann == "down":
                match, opp = has_down, has_up and not has_down
            elif ann == "up":
                match, opp = has_up, has_down and not has_up
            else:
                match, opp = has_down or has_up, False
            neg = bool(NEG.search(ctx))
            dis = bool(DIS.search(txt))
            if not (iv or fn):
                n["mention_only"] += 1
            if iv:
                n["in_vitro"] += 1; best["in_vitro"].append(r["pmid"])
            func_ok = fn and match and bool(sents)
            if func_ok:
                n["functional_match"] += 1; best["functional_match"].append(r["pmid"])
            if func_ok and dis:
                n["m4_records"] += 1; best["m4"].append(r["pmid"])
            if dis:
                n["disease"] += 1
            if opp:
                n["opposite_only"] += 1; best["opposite"].append(r["pmid"])
            if neg:
                n["negative"] += 1; best["negative"].append(r["pmid"])
            rec_rows.append({"Drug_ID": did, "Drug": v["drug"], "Gene": gene, "annotated_direction": ann, "pmid": r["pmid"], "year": r["year"],
                             "title": r["title"][:160], "in_vitro_terms": iv, "functional_terms": fn, "action_match": match, "opposite_only": opp,
                             "negative_phrase": neg, "disease_term": dis})
        if n["m4_records"]:
            m = 4
        elif n["functional_match"]:
            m = 3
        elif n["in_vitro"] or (n["primary"] and n["mention_only"] < n["primary"]):
            m = 2
        elif n["primary"]:
            m = 1
        else:
            m = 0
        contra = bool(n["opposite_only"] or n["negative"])
        rows.append({"Drug_ID": did, "Drug": v["drug"], "Gene": gene, "annotated_direction": ann, "direction_check": pr.direction_check,
                     "pubmed_hits": v["count"], **n, "M": m, "contradictory_evidence": contra,
                     "pmids_engagement": ";".join(best["in_vitro"][:5]), "pmids_action": ";".join(best["functional_match"][:5]),
                     "pmids_disease": ";".join(best["m4"][:5]), "pmids_contradiction": ";".join((best["opposite"] + best["negative"])[:5])})
    g = pd.DataFrame(rows).sort_values(["M", "Drug"], ascending=[False, True])
    g.to_csv(OUT / "step4_pair_grades.csv", index=False)
    pd.DataFrame(rec_rows).to_csv(OUT / "step4_record_flags.csv", index=False)
    print(g.M.value_counts().sort_index().to_dict(), "contradictory pairs:", int(g.contradictory_evidence.sum()))


if __name__ == "__main__":
    main()
