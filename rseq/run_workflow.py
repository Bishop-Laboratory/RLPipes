import snakemake as snk
import json
import sys
import io
import time
import os
import pathlib
from contextlib import redirect_stdout
import warnings


def make_snakes(config_file, snake_args, verify=True):
    # config_file = "rseq_out/config.json"
    config = json.load(open(config_file))
    snake_path = os.path.join(config['helpers_dir'][0], "rseq_workflow.smk")
    outdir = config['outdir'][0]
    threads = config['threads'][0]
    
    # Set use_conda as True if not supplied
    if 'use_conda' not in snake_args.keys():
        warnings.warn("'use_conda' is set to False by user. It is expected that user has all dependencies installed.")
        time.sleep(4)
    elif not snake_args['use_conda']:
        warnings.warn("'use_conda' is set to False by user. It is expected that user has all dependencies installed.")
        time.sleep(4)

    # Set conda_frontend='mamba' if not supplied
    if 'conda_frontend' not in snake_args.keys():
        snake_args['conda_frontend'] = 'mamba'
    elif kwargs['conda_frontend'] == 'conda':
        warnings.warn("'conda_frontend' is set to 'conda' by user. This will lead to slower operations, consider using "
              "'mamba' instead.")
        time.sleep(4)

    if verify:
        # Make DAG and perform dry-run
        good_exit = snk.snakemake(snake_path, config=config, cores=threads, dryrun=True, **snake_args)
    
        pathlib.Path(outdir + "/dags/").mkdir(parents=True, exist_ok=True)
        out = io.StringIO()
        with redirect_stdout(out):
            snk.snakemake(snake_path, printdag=True, config=config, **snake_args)
            out = out.getvalue()
            out_file = outdir + '/dags/dag.gv'
    
            if os.path.exists(out_file):
                os.remove(out_file)
    
            with open(out_file, 'a') as stdout_log:
                stdout_log.writelines(out)
    
            out_svg = outdir + '/dags/dag.png'
            os.system('cat ' + out_file + ' | dot -Tpng -o ' + out_svg)
            os.remove(out_file)
    else:
        # Run snakemake
        good_exit = snk.snakemake(snake_path, config=config, cores=threads, **snake_args)

    # Check exit status
    if good_exit:
        sys.exit(0)
    else:
        sys.exit(1)
