"""Shared Neo4j helpers for this run.

Credentials come only from the environment: NEO4J_PASSWORD (required),
NEO4J_URI (default bolt://127.0.0.1:7687), NEO4J_USERNAME (default neo4j).
Every write goes through `write_session()`, which refuses any database other
than TARGET_DB ('paper2public').
"""
import datetime as dt
import json
import os
from pathlib import Path

from neo4j import GraphDatabase

ROOT = Path(__file__).resolve().parents[1]
TARGET_DB = os.environ.get("NEO4J_DATABASE", "paper2public")   # set NEO4J_DATABASE if the graph lives in another database (e.g. the default 'neo4j')
URI = os.environ.get("NEO4J_URI", "bolt://127.0.0.1:7687")
USER = os.environ.get("NEO4J_USERNAME", "neo4j")
NAFLD = "EFO_0003095"


def driver():
    pw = os.environ.get("NEO4J_PASSWORD")
    if not pw:
        raise SystemExit("Set NEO4J_PASSWORD in the environment.")
    return GraphDatabase.driver(URI, auth=(USER, pw))


def check_target(drv, database=TARGET_DB):
    if database != TARGET_DB:
        raise SystemExit(f"Refusing to write to '{database}': only '{TARGET_DB}' is allowed.")
    with drv.session(database="system") as s:
        row = s.run("SHOW DATABASES YIELD name, currentStatus WHERE name = $n RETURN currentStatus",
                    n=database).single()
    if row is None or row["currentStatus"] != "online":
        raise SystemExit(f"Database '{database}' is not online.")


def write_session(drv, database=TARGET_DB):
    check_target(drv, database)
    return drv.session(database=database)


def read_session(drv, database=TARGET_DB):
    return drv.session(database=database)


def counts(sess):
    """Node counts per label and relationship counts per type."""
    labels = {r["l"]: r["n"] for r in sess.run(
        "CALL db.labels() YIELD label AS l CALL (l) { MATCH (n) WHERE l IN labels(n) RETURN count(n) AS n } RETURN l, n")}
    rels = {r["t"]: r["n"] for r in sess.run(
        "CALL db.relationshipTypes() YIELD relationshipType AS t "
        "CALL (t) { MATCH ()-[r]->() WHERE type(r) = t RETURN count(r) AS n } RETURN t, n")}
    total_n = sess.run("MATCH (n) RETURN count(n) AS n").single()["n"]
    total_r = sess.run("MATCH ()-[r]->() RETURN count(r) AS n").single()["n"]
    return {"nodes_total": total_n, "relationships_total": total_r,
            "nodes": dict(sorted(labels.items())), "relationships": dict(sorted(rels.items()))}


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def write_log(path, log):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(log, indent=2, default=str), encoding="utf-8")


def run_batches(sess, query, rows, batch=10000, **params):
    """Run `query` (using $rows) over rows in batches, each batch in its own transaction."""
    total = 0
    for i in range(0, len(rows), batch):
        res = sess.execute_write(lambda tx, b: tx.run(query, rows=b, **params).single(), rows[i:i + batch])
        total += res[0] if res is not None else 0
    return total
