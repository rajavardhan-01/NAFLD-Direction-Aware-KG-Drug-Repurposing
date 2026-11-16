// Undo step 0.2: remove only the properties ChEMBL filled (OptimusKG values were never overwritten)
MATCH (d:Drug) WHERE d.smiles_source STARTS WITH 'ChEMBL' REMOVE d.smiles, d.smiles_source;
MATCH (d:Drug) WHERE d.molecule_type_source STARTS WITH 'ChEMBL' REMOVE d.molecule_type, d.molecule_type_source;
MATCH (d:Drug) WHERE d.max_phase_source STARTS WITH 'ChEMBL' REMOVE d.max_phase, d.max_phase_source;
MATCH (d:Drug) WHERE d.recovered_name_source STARTS WITH 'ChEMBL' REMOVE d.recovered_name, d.recovered_name_source;
