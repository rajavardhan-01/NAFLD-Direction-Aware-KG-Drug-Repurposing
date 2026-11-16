// Undo step 0.4 (or run: python code/write_new_links_to_graph.py --undo)
MATCH ()-[r:ACTS_ON {layer: 'added curated drug-target links'}]->() DELETE r;
MATCH (d:Drug {source: 'added curated drug-target links'}) WHERE NOT (d)--() DELETE d;
