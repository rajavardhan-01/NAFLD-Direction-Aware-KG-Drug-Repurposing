// Undo step 0.0: the database was empty before the import, so undo = delete everything.
// Run only against paper2public.
MATCH (n) CALL (n) { DETACH DELETE n } IN TRANSACTIONS OF 10000 ROWS;
