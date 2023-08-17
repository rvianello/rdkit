CREATE INDEX fpidx ON pgbfp USING gin (f);
CREATE INDEX maccsfpidx ON pgbfp USING gin (maccsf);

SET rdkit.tanimoto_threshold = 0.5;


SET enable_indexscan=off;
SET enable_bitmapscan=off;
SET enable_seqscan=on;

SELECT
    id, tanimoto_sml(rdkit_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol), f) AS sml
FROM
	pgbfp
WHERE rdkit_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol) % f
ORDER BY sml DESC, id limit 10;

SELECT
    id, tanimoto_sml(maccs_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol), maccsf) AS sml
FROM
	pgbfp
WHERE maccs_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol) % maccsf
ORDER BY sml DESC, id limit 10;


SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=off;

SELECT
    id, tanimoto_sml(rdkit_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol), f) AS sml
FROM
	pgbfp
WHERE rdkit_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol) % f
ORDER BY sml DESC, id limit 10;

SELECT
    id, tanimoto_sml(maccs_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol), maccsf) AS sml
FROM
	pgbfp
WHERE maccs_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol) % maccsf
ORDER BY sml DESC, id limit 10;

SET enable_indexscan=on;
SET enable_bitmapscan=on;
SET enable_seqscan=on;

SELECT
    id, tanimoto_sml(rdkit_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol), f) AS sml
FROM
	pgbfp
WHERE rdkit_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol) % f
ORDER BY sml DESC, id limit 10;

SELECT
    id, tanimoto_sml(maccs_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol), maccsf) AS sml
FROM
	pgbfp
WHERE maccs_fp('CC(=NNC(=O)C1=CC(=C(C=C1)OC)OC)C2=CC=CC=N2'::mol) % maccsf
ORDER BY sml DESC, id limit 10;

DROP INDEX fpidx;
DROP INDEX maccsfpidx;
