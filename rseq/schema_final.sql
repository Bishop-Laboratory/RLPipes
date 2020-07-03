

CREATE TABLE IF NOT EXISTS runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    run_name TEXT UNIQUE NOT NULL,
    owner TEXT NOT NULL,
    pipeline_id INT NOT NULL,
    current_stage TEXT NOT NULL,
    messages TEXT,
    errors TEXT,
    results TEXT,
    last_update DATETIME
);

CREATE TABLE IF NOT EXISTS pipelines (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    pipeline_name TEXT UNIQUE NOT NULL,
    created DATETIME NOT NULL,
    description TEXT,
    pipeline_author TEXT,
    pipeline TEXT NOT NULL
);
