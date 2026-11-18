"""Export the graph tables needed by the offline evaluation (gene-gene edges, drug targets, pathways).

Read-only Neo4j access. Writes results/Evaluation/inputs/*.csv (gene_gene.csv is large and git-ignored).

    NEO4J_PASSWORD=... NEO4J_DATABASE=neo4j python code/eval_export_graph.py
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from neo4j_util import NAFLD, ROOT, driver, read_session  # noqa: E402

OUT = ROOT / "results/Evaluation/inputs"
OUT.mkdir(parents=True, exist_ok=True)
CAND = pd.read_csv(ROOT / "results/Direction_consistent_drugs/candidate_drugs_direction_summary.csv").Drug_ID.tolist()
Q = {
    "gene_gene": "MATCH (a:Gene)-[:INTERACTS_WITH]-(b:Gene) WHERE a.node_id < b.node_id RETURN DISTINCT a.node_id AS a, b.node_id AS b",
    "genes": "MATCH (g:Gene) RETURN g.node_id AS gene_id, g.symbol AS symbol",
    "drug_targets": """MATCH (d:Drug)-[r:ACTS_ON]->(g:Gene) RETURN DISTINCT d.node_id AS drug_id, d.name AS drug, g.node_id AS gene_id""",
    "pathway_genes": "MATCH (g:Gene)-[:IN_PATHWAY]->(p:Pathway) RETURN DISTINCT g.node_id AS gene_id, p.node_id AS pathway_id, p.name AS pathway, coalesce(p.source,'Reactome') AS library",
    "gene_degree": "MATCH (g:Gene) RETURN g.node_id AS gene_id, COUNT { (g)-[:INTERACTS_WITH]-() } AS degree_rel",
    "drug_liver_links": """UNWIND $ids AS id MATCH (d:Drug {node_id: id})-[r:TREATS|OFF_LABEL_FOR|CONTRAINDICATED_IN]->(x:Disease)
        RETURN DISTINCT id AS drug_id, type(r) AS rel, x.node_id AS disease_id""",
    "candidate_links_all": """UNWIND $ids AS id MATCH (d:Drug {node_id: id})-[r:ACTS_ON]->(g:Gene)
        RETURN d.node_id AS Drug_ID, g.node_id AS Gene_ID, r.layer AS layer, r.action AS okg_action,
               r.relation_assertions AS relation_assertions, r.action_labels AS ext_action_labels, r.push AS ext_push,
               r.external_action_direction AS external_action_direction""",
    "pathway_nafld": "MATCH (p:Pathway)-[:IMPLICATED_IN]->(:Disease {node_id: $n}) RETURN p.node_id AS pathway_id, p.name AS pathway",
}
drv = driver()
with read_session(drv) as s:
    for k, q in Q.items():
        df = pd.DataFrame(s.run(q, n=NAFLD, ids=CAND).data())
        df.to_csv(OUT / f"{k}.csv", index=False)
        print(k, df.shape)
