// Removes ONLY the properties added in step 0.5
MATCH ()-[r:ACTS_ON]->() WHERE r.external_action_direction IS NOT NULL
REMOVE r.external_action_direction, r.external_action_label, r.external_action_sources, r.external_action_release, r.external_action_added
RETURN count(r) AS removed;
