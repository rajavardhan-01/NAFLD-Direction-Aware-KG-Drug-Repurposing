// Undo step 0.6
MATCH (p:Pathway) WHERE p.source IN ['KEGG_legacy','Hallmark'] DETACH DELETE p;
