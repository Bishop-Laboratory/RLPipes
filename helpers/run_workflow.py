import snakemake as snk
import json
import sys
import io
import os
import time
import pathlib
from contextlib import redirect_stdout


def make_snakes(config_file):

    # config_file = "tests/config.json"
    config = json.load(open(config_file))
    config['helpers_dir'] = sys.path[0]
    snake_path = os.path.join(sys.path[0], "rseq_workflow.smk")
    outdir = config['outdir'][0]
    snake_args = config['snake_args']
    threads = config['threads'][0]

    # Make DAG
    pathlib.Path(outdir + "/dags/").mkdir(parents=True, exist_ok=True)
    out = io.StringIO()
    with redirect_stdout(out):
        snk.snakemake(snake_path, printdag=True, config=config)
        out = out.getvalue()
        out_file = outdir + '/dags/dag.gv'

        if os.path.exists(out_file):
            os.remove(out_file)

        with open(out_file, 'a') as stdout_log:
            stdout_log.writelines(out)

        out_svg = outdir + '/dags/dag.png'
        os.system('cat ' + out_file + ' | dot -Tpng -o ' + out_svg)
        os.remove(out_file)

    # Run pipeline
    kwargs = {
        "printdag": True,
        "dryrun": True
    }
    snk.snakemake(snake_path, config=config, cores=threads,
                  **kwargs)


if __name__ == "__main__":
    configs = json.load(open(sys.argv[1]))
    make_snakes(config_file=sys.argv[1])
