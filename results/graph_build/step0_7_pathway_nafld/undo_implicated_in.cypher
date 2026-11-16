// Undo step 0.7
MATCH (:Pathway)-[r:IMPLICATED_IN {source: 'seed enrichment + CTD'}]->(:Disease) DELETE r;
