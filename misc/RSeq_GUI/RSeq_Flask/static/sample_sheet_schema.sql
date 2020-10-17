CREATE TABLE IF NOT EXISTS samples (
    sample_name TEXT PRIMARY KEY,
    rseq_params TEXT,
    workflow_id INT NOT NULL,
    dag TEXT,
    run_snakemake TEXT,
    progress TEXT,
    log TEXT,
    status TEXT,
    args TEXT,
    targets TEXT,
    rule_info TEXT,
    resources TEXT
);

