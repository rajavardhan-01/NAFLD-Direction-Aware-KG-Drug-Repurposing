"""Offline evaluation core: the proposed method, an ablation family and six baselines on one shared graph.

Everything here works on tables exported once from Neo4j (code/eval_export_graph.py) and needs no database.

Methods (all score the same universe of candidate drugs; higher score = better):
  SP        shortest path             1 / (1 + min distance(target, seed)), unsigned PPI, the target itself excluded
  NP        network proximity         -Z of the mean closest-seed distance against degree-matched random target sets
  RWR       random walk with restart  mean stationary score of the drug's targets, uniform restart on the seeds
  MPC       metapath counting         number of BB routes, no direction filter
  DWPC      degree-weighted path count, same routes as MPC, weight prod (deg_u * deg_v)^-alpha
  HeteSim   cosine of the two half-path distributions of B and BB (equal-weight mean)
  Proposed  BB route count with the direction filter (only drugs with >= 1 consistent and no inconsistent link)
The route rules (hub limit, pathway size, known liver drugs) are those of code/metapaths.py; the route counts are
computed with sparse matrix products and reproduce the Cypher result (checked in code/eval_run.py).
"""
import re
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[1]
INP = ROOT / "results/Evaluation/inputs"
HUB_MAX = 500
SIZE_MIN, SIZE_MAX = 5, 200
KNOWN_LIVER = {"EFO_0003095", "EFO_1001249", "EFO_0001422"}
MESH_NAFLD = "MESH:D065626"
STAGES = ["1_steatosis", "2_NASH", "3_fibrosis", "4_HCC"]


def bh_qvalues(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def norm(v):
    return re.sub(r"[^A-Z0-9]+", "_", str(v).upper()).strip("_")


class Data:
    """Static tables. Perturbations never change these, they build modified copies of the arrays."""

    def __init__(self):
        g = pd.read_csv(INP / "genes.csv")
        self.gene_ids = g.gene_id.to_numpy()
        self.symbol = g.symbol.to_numpy()
        self.idx = {x: i for i, x in enumerate(self.gene_ids)}
        self.N = len(self.gene_ids)
        e = pd.read_csv(INP / "gene_gene.csv")
        a = e.a.map(self.idx).to_numpy()
        b = e.b.map(self.idx).to_numpy()
        keep = a != b
        self.edges = np.unique(np.sort(np.c_[a[keep], b[keep]], axis=1), axis=0)
        self.deg_rel = pd.read_csv(INP / "gene_degree.csv").set_index("gene_id").degree_rel.reindex(self.gene_ids).to_numpy()
        seeds = pd.read_csv(ROOT / "seed_genes/seed_genes_main_Ryaboshapkina2017.csv")
        seeds = seeds[seeds.ensembl_id.isin(self.idx)].copy()
        seeds["i"] = seeds.ensembl_id.map(self.idx)
        self.seeds = seeds.reset_index(drop=True)
        # candidate drug-seed links and their direction class
        L = pd.read_csv(ROOT / "results/Direction_consistent_drugs/candidate_links_direction_check.csv")
        L["i"] = L.Gene_ID.map(self.idx)
        a_, s_ = L.Action_direction, L.seed_direction
        L["dcheck"] = np.select([((a_ == "down") & (s_ == "increasing")) | ((a_ == "up") & (s_ == "decreasing")),
                                 ((a_ == "down") & (s_ == "decreasing")) | ((a_ == "up") & (s_ == "increasing"))],
                                ["consistent", "inconsistent"], "undecided")
        assert (L.dcheck == L.direction_check.map({"direction-consistent": "consistent", "direction-inconsistent": "inconsistent",
                                                   "not decided": "undecided"})).all()
        self.links = L[["Drug", "Drug_ID", "Gene_ID", "i", "dcheck"]].drop_duplicates().reset_index(drop=True)
        self.candidates = self.links.Drug_ID.drop_duplicates().tolist()
        self.name = self.links.drop_duplicates("Drug_ID").set_index("Drug_ID").Drug.to_dict()
        # liver status
        dl = pd.read_csv(INP / "drug_liver_links.csv")
        liver = pd.read_csv(ROOT / "results/Task4/liver_disease_set.csv")
        liver_ids = set(liver.loc[~liver.excluded_non_hepatic, "disease_id"])
        self.known = set(dl[dl.rel.isin(["TREATS", "OFF_LABEL_FOR"]) & dl.disease_id.isin(KNOWN_LIVER)].drug_id)
        self.ci_liver = set(dl[(dl.rel == "CONTRAINDICATED_IN") & dl.disease_id.isin(liver_ids)].drug_id)
        # all targets of a drug (for SP / NP / RWR)
        dt = pd.read_csv(INP / "drug_targets.csv")
        dt["i"] = dt.gene_id.map(self.idx)
        dt = dt.dropna(subset=["i"]).astype({"i": int})
        self.drug_targets_all = dt
        self.drug_outdeg = dt.groupby("drug_id").size()
        self.gene_indeg = np.bincount(dt.i, minlength=self.N)
        # pathways
        pg = pd.read_csv(INP / "pathway_genes.csv")
        pg["i"] = pg.gene_id.map(self.idx)
        pg = pg.dropna(subset=["i"]).astype({"i": int})
        pg["db"] = np.where(pg.pathway_id.str.startswith("REACT:"), "Reactome",
                            np.where(pg.library == "KEGG_legacy", "KEGG", np.where(pg.library == "Hallmark", "Hallmark", "other")))
        self.memb = pg[["pathway_id", "pathway", "i", "db"]].drop_duplicates(["pathway_id", "i"]).reset_index(drop=True)
        self.pathway_meta = self.memb.drop_duplicates("pathway_id").set_index("pathway_id")[["pathway", "db"]]
        self.nafld_pathways = set(pd.read_csv(INP / "pathway_nafld.csv").pathway_id)
        self._ctd = None
        self.universe = [d for d in self.candidates if d not in self.known and d not in self.ci_liver]

    # ---------- CTD support for the pathway rule (independent of the seeds) ----------
    def ctd(self):
        if self._ctd is None:
            c = pd.read_csv(ROOT / "data/ctd/CTD_diseases_pathways.csv.gz", comment="#", header=None,
                            names=["DiseaseName", "DiseaseID", "PathwayName", "PathwayID", "Genes"])
            c = c[c.DiseaseID.eq(MESH_NAFLD)]
            kegg_nodes = [p for p in self.pathway_meta.index if self.pathway_meta.db[p] == "KEGG"]
            kegg_by_name = {norm(n[len("KEGG:KEGG_"):]): n for n in kegg_nodes}
            sup = {}
            for (pid, name), g in c.groupby(["PathwayID", "PathwayName"]):
                genes = {v for t in g.Genes.dropna() for v in re.split(r"[|;]", t)} - {""}
                if pid.startswith("REACT:"):
                    sup[pid] = len(genes)
                elif pid.startswith("KEGG:"):
                    node = kegg_by_name.get(norm(name))
                    if node:
                        sup[node] = len(genes)
            self._ctd = sup
        return self._ctd


def adjacency(D, edges=None, removed_nodes=None):
    e = D.edges if edges is None else edges
    A = sp.coo_matrix((np.ones(len(e), dtype=np.float32), (e[:, 0], e[:, 1])), shape=(D.N, D.N))
    A = (A + A.T).tocsr()
    A.data[:] = 1
    if removed_nodes is not None and len(removed_nodes):
        keep = np.ones(D.N, dtype=np.float32)
        keep[removed_nodes] = 0
        Dg = sp.diags(keep)
        A = (Dg @ A @ Dg).tocsr()
        A.eliminate_zeros()
    return A


# ------------------------------------------------------------------ pathway rule
def pathway_links(D, seeds, memb=None):
    """Step 0.7 of the pipeline (seed enrichment + CTD rule) in vectorised form. Returns the set of linked pathway ids."""
    memb = D.memb if memb is None else memb
    bgi = np.unique(memb.i.to_numpy())
    pos = {g: k for k, g in enumerate(bgi)}
    pids = memb.pathway_id.drop_duplicates().to_numpy()
    prow = {p: k for k, p in enumerate(pids)}
    M = sp.csr_matrix((np.ones(len(memb)), (memb.pathway_id.map(prow), memb.i.map(pos))), shape=(len(pids), len(bgi)))
    M.data[:] = 1
    K = np.asarray(M.sum(1)).ravel()
    db = D.pathway_meta.db.reindex(pids).to_numpy()
    s = seeds[seeds.i.isin(pos)]
    groups = {"all_histology": s[s.main_stage != "0_genetic_risk"]}
    for st in STAGES:
        sub = s[s.main_stage == st]
        groups[st] = sub
        for d in ["increasing", "decreasing"]:
            groups[f"{st}_{d}"] = sub[sub.direction == d]
    N = len(bgi)
    sup = D.ctd()
    sig = set()
    for g, sub in groups.items():
        n = len(sub)
        if n == 0:
            continue
        ind = np.zeros(len(bgi))
        ind[[pos[i] for i in sub.i]] = 1
        k = np.asarray(M @ ind).ravel().astype(int)
        p = np.ones(len(pids))
        m = k > 0
        p[m] = hypergeom.sf(k[m] - 1, N, K[m].astype(int), n)
        for dbn in ("Reactome", "KEGG", "Hallmark"):
            sel = np.where(db == dbn)[0]
            if len(sel) == 0:
                continue
            q = bh_qvalues(p[sel])
            sig.update(pids[sel[q < 0.05]])
    out = set()
    for pth in sig:
        if D.pathway_meta.db[pth] == "Hallmark" or sup.get(pth, 0) >= 2:
            out.add(pth)
    return out


def pathway_gene_stats(D, linked, memb=None):
    """For the linked pathways of size 5-200 return np[g] (eligible pathways of gene g), the size-sum used by DWPC
    and the sparse eligible membership matrix (pathway x gene)."""
    memb = D.memb if memb is None else memb
    size = memb.groupby("pathway_id").size()
    npw_all = np.bincount(memb.i, minlength=D.N)
    elig = [p for p in linked if p in size.index and SIZE_MIN <= size[p] <= SIZE_MAX]
    me = memb[memb.pathway_id.isin(elig)]
    return elig, size, npw_all, me


# ------------------------------------------------------------------ route-based scores
class Routes:
    """Vectors needed by every path-based method for one (graph, seed set, pathway set)."""

    def __init__(self, D, A, seeds_i, linked, memb=None, alpha=0.4, hub=HUB_MAX):
        self.D, self.A, self.hub = D, A, hub
        self.deg = np.asarray(A.sum(1)).ravel()
        self.S = np.zeros(D.N, dtype=np.float32)
        self.S[seeds_i] = 1
        elig, size, npw_all, me = pathway_gene_stats(D, linked, memb)
        self.elig = elig
        self.np = np.bincount(me.i, minlength=D.N).astype(np.float32)
        self.H = (self.deg <= hub).astype(np.float32)
        self.v = self.S * self.np
        w = self.A @ self.v
        self.B = np.where(self.deg <= hub, w, 0)                      # B routes of a start gene
        x = self.H * w
        nH = self.A @ self.H
        self.BB = np.where(self.deg <= hub, self.A @ x - nH * self.v, 0)
        self.me, self.size, self.npw_all = me, size, npw_all
        self._w = {}
        self.alpha = alpha

    def weighted(self, alpha):
        """DWPC vectors for the exponent alpha (relationship-specific endpoint degrees)."""
        if alpha in self._w:
            return self._w[alpha]
        D, A = self.D, self.A
        d = np.maximum(self.deg, 1)
        Acoo = A.tocoo()
        W = sp.csr_matrix(((d[Acoo.row] * d[Acoo.col]) ** (-alpha), (Acoo.row, Acoo.col)), shape=A.shape)
        me = self.me
        sz = self.size.reindex(me.pathway_id).to_numpy()
        npw = np.maximum(self.npw_all, 1)
        n_impl = max(len(self.elig), 1)
        wp = np.bincount(me.i, weights=(sz ** -alpha) * (n_impl ** -alpha), minlength=D.N)
        vw = self.S * (npw ** -alpha) * wp
        wv = W @ vw
        Bw = np.where(self.deg <= self.hub, wv, 0)
        x = self.H * wv
        # diagonal term: paths t-h-t (Gene* = start gene)
        Wh = W.multiply(self.H[None, :]).tocsr()
        diag = np.asarray(W.multiply(Wh.T).sum(1)).ravel() * vw
        BBw = np.where(self.deg <= self.hub, W @ x - diag, 0)
        self._w[alpha] = (Bw, BBw)
        return self._w[alpha]


def start_links(D, seeds_i, mode):
    """Start links (drug, gene index) of the path-based methods and the drugs allowed by the filter.
    mode: none | decidable | consistent_any | proposed"""
    L = D.links[D.links.i.isin(seeds_i)]
    L = L[L.Drug_ID.isin(D.universe)]
    g = L.groupby("Drug_ID").dcheck.agg(lambda s: (int((s == "consistent").sum()), int((s == "inconsistent").sum())))
    nc = g.map(lambda t: t[0])
    ni = g.map(lambda t: t[1])
    if mode == "none":
        use = L
    elif mode == "decidable":
        ok = (nc + ni) > 0
        use = L[L.Drug_ID.isin(ok[ok].index) & (L.dcheck != "undecided")]
    elif mode == "consistent_any":
        ok = nc > 0
        use = L[L.Drug_ID.isin(ok[ok].index) & (L.dcheck == "consistent")]
    elif mode == "proposed":
        ok = (nc > 0) & (ni == 0)
        use = L[L.Drug_ID.isin(ok[ok].index) & (L.dcheck == "consistent")]
    else:
        raise ValueError(mode)
    return use


def path_scores(D, R, use, alpha=None):
    """BB routes (or weighted routes) per drug under the pipeline rules; drugs without a kept B route score 0."""
    byd = use.groupby("Drug_ID").i.agg(list)
    out = {}
    for d, gi in byd.items():
        gi = np.array(sorted(set(gi)))
        b = R.B[gi].sum()
        if b <= 0:
            out[d] = 0.0
            continue
        if alpha is None:
            out[d] = float(R.BB[gi].sum())
        else:
            Bw, BBw = R.weighted(alpha)
            dd = D.drug_outdeg.get(d, 1)
            out[d] = float((((dd * np.maximum(D.gene_indeg[gi], 1)) ** -alpha) * BBw[gi]).sum())
    return pd.Series(out).reindex(D.universe).fillna(0.0)


def hetesim_scores(D, R, use):
    """HeteSim of the B (even) and BB (odd, edge-node midpoint) metapaths, equal-weight mean."""
    A, H, S = R.A, R.H, R.S
    drugs = sorted(use.Drug_ID.unique())
    di = {d: k for k, d in enumerate(drugs)}
    Hm = sp.diags(H)
    # drug -> start gene (hub rule on the target as in B)
    rows = [di[d] for d in use.Drug_ID]
    cols = use.i.to_numpy()
    ok = R.deg[cols] <= R.hub
    U1 = sp.csr_matrix((np.ones(ok.sum()), (np.array(rows)[ok], cols[ok])), shape=(len(drugs), D.N))
    U1.data[:] = 1
    U1 = sp.diags(1 / np.maximum(np.asarray(U1.sum(1)).ravel(), 1)) @ U1
    # disease -> pathway -> Gene*
    me = R.me
    if len(me) == 0:
        return pd.Series(0.0, index=D.universe)
    pid = {p: k for k, p in enumerate(R.elig)}
    Mp = sp.csr_matrix((np.ones(len(me)), (me.pathway_id.map(pid), me.i)), shape=(len(R.elig), D.N))
    Mp.data[:] = 1
    Mp = Mp @ sp.diags(S)                                   # only seed genes are Gene*
    Mp.eliminate_zeros()
    rown = np.asarray(Mp.sum(1)).ravel()
    keep = rown > 0
    Mn = sp.diags(np.where(keep, 1 / np.maximum(rown, 1), 0)) @ Mp
    PRg = np.asarray(Mn.sum(0)).ravel() / max(keep.sum(), 1)       # disease -> uniform pathway -> uniform member seed gene
    # ---- B: Drug -> G1 -> Gene* | Gene* <- Pathway <- NAFLD
    AS = A @ sp.diags(S)                                      # start gene -> neighbouring seed genes
    ns = np.asarray(AS.sum(1)).ravel()
    U2 = sp.diags(np.where(ns > 0, 1 / np.maximum(ns, 1), 0)) @ AS
    PL = (U1 @ U2).tocsr()
    num = PL @ PRg
    nl = np.sqrt(np.asarray(PL.multiply(PL).sum(1)).ravel())
    nr = np.sqrt((PRg ** 2).sum())
    hb = np.where((nl > 0) & (nr > 0), num / np.maximum(nl * nr, 1e-30), 0)
    # ---- BB: Drug -> G1 -> G2 -> edge(G2,Gene*) | edge <- Gene* <- Pathway <- NAFLD
    AH = A @ Hm                                               # start gene -> neighbours that may be G2 (<= hub limit)
    nh = np.asarray(AH.sum(1)).ravel()
    U12 = sp.diags(np.where(nh > 0, 1 / np.maximum(nh, 1), 0)) @ AH
    pl = (U1 @ U12).tocsr()                                   # mass on G2
    nsg = ns                                                  # seed neighbours of every G2 candidate
    inv_ns = np.where(nsg > 0, 1 / np.maximum(nsg, 1), 0)
    nl2 = np.sqrt(np.asarray(pl.multiply(pl) @ inv_ns).ravel())
    nHg = np.asarray(A @ H).ravel()                           # G2 neighbours of a Gene*
    xg = np.where(nHg > 0, PRg / np.maximum(nHg, 1), 0)
    nr2 = np.sqrt((np.where(nHg > 0, PRg ** 2 / np.maximum(nHg, 1), 0)).sum())
    AHS = (Hm @ A @ sp.diags(S)).tocsr()                      # G2 -> Gene*
    t = AHS @ xg                                              # per G2: sum over its seed neighbours
    dot = np.asarray(pl @ (inv_ns * t)).ravel()
    hbb = np.where((nl2 > 0) & (nr2 > 0), dot / np.maximum(nl2 * nr2, 1e-30), 0)
    s = pd.Series((hb + hbb) / 2, index=drugs)
    return s.reindex(D.universe).fillna(0.0)


# ------------------------------------------------------------------ distance / propagation baselines
def nearest_seed_excl_self(A, seeds_i, limit=8):
    """Distance from every node to the nearest seed other than itself (inf if none within `limit`/unreachable)."""
    seeds_i = np.asarray(seeds_i)
    D2 = dijkstra(A, directed=False, indices=seeds_i, unweighted=True, limit=limit)    # (n_seeds, N)
    D2[np.arange(len(seeds_i)), seeds_i] = np.inf          # a seed is not its own nearest seed
    return D2.min(0)


def sp_scores(D, nd):
    t = D.drug_targets_all[D.drug_targets_all.drug_id.isin(D.universe)]
    dist = nd[t.i.to_numpy()]
    s = pd.Series(dist, index=t.drug_id.to_numpy()).groupby(level=0).min()
    sc = 1 / (1 + s)
    sc[~np.isfinite(s)] = 0
    return sc.reindex(D.universe).fillna(0.0)


def np_scores(D, nd, A, draws=1000, seed=0):
    """-Z of the mean closest-seed distance; random target sets are degree-matched (same size, same degree bin);
    unreachable targets are dropped from the drug and from the random pool."""
    rng = np.random.default_rng(seed)
    deg = np.asarray(A.sum(1)).ravel()
    reach = np.isfinite(nd)
    pool = np.where(reach & (deg > 0))[0]
    bins = np.minimum(np.floor(np.log2(deg[pool] + 1) * 2).astype(int), 40)
    by_bin = {b: pool[bins == b] for b in np.unique(bins)}
    t = D.drug_targets_all[D.drug_targets_all.drug_id.isin(D.universe)]
    out, zero_var = {}, 0
    for d, g in t.groupby("drug_id"):
        gi = g.i.to_numpy()
        gi = gi[reach[gi] & (deg[gi] > 0)]
        if len(gi) == 0:
            continue
        dc = nd[gi].mean()
        bb = np.minimum(np.floor(np.log2(deg[gi] + 1) * 2).astype(int), 40)
        samp = np.empty((draws, len(gi)))
        for k, b in enumerate(bb):
            cand = by_bin.get(b)
            if cand is None or len(cand) == 0:
                cand = pool
            samp[:, k] = nd[rng.choice(cand, size=draws)]
        m = samp.mean(1)
        sd = m.std()
        if sd == 0:
            zero_var += 1
            continue
        out[d] = -(dc - m.mean()) / sd
    s = pd.Series(out).reindex(D.universe)
    s.attrs["zero_variance"] = zero_var
    s.attrs["no_reachable_target"] = int(s.isna().sum())
    return s.fillna(-1e9)


def rwr_scores(D, A, seeds_i, r=0.5, tol=1e-8, maxit=500):
    deg = np.asarray(A.sum(0)).ravel()
    W = A @ sp.diags(np.where(deg > 0, 1 / np.maximum(deg, 1), 0))        # column-normalised
    p0 = np.zeros(D.N)
    p0[seeds_i] = 1 / len(seeds_i)
    p = p0.copy()
    for _ in range(maxit):
        pn = (1 - r) * (W @ p) + r * p0
        if np.abs(pn - p).sum() < tol:
            p = pn
            break
        p = pn
    t = D.drug_targets_all[D.drug_targets_all.drug_id.isin(D.universe)]
    s = pd.Series(p[t.i.to_numpy()], index=t.drug_id.to_numpy()).groupby(level=0).mean()
    return s.reindex(D.universe).fillna(0.0)
