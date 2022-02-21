import snakemake as snk
import json
import pandas as pd
import sys
import io
import time
import click
import os
import pathlib
from contextlib import redirect_stdout
import warnings


def make_snakes(
    run_dir, snake_args, src_dir, bwamem2=False, macs3=False, threads=1, groupby=None,
    noexp=False, debug=False, noreport=False, tsv=False, useaws=False, verify=True
):
    if tsv:
        config = pd.read_csv(os.path.join(run_dir, "config.tsv"), sep = "\t").fillna("").to_dict("list")
    else:
        config = json.load(open(os.path.join(run_dir, "config.json")))
        
    config["debug"] = debug
    config["run_dir"] = run_dir
    config["src"] = src_dir
    config["bwamem2"] = bwamem2
    config["macs3"] = macs3
    config["groupby"] = groupby
    config['noexp'] = noexp
    config['noreport'] = noreport
    config['useaws'] = useaws
    snake_path = os.path.join(config["src"], "rlpipes.smk")
    
    # Check groupby
    if groupby is not None:
        def testgroup(group):
            assert group in config.keys()
        try:
            [testgroup(group) for group in groupby.split(" ")]
        except AssertionError:
            raise click.BadParameter("groupby '" + groupby + "' invalid. All groups must be found in config and separated with a space. E.g., '-G group1 group2'")

    if debug:
        snake_args["notemp"] = True

    # Set use_conda as True if not supplied
    if "use_conda" not in snake_args.keys():
        snake_args["use_conda"] = True
        time.sleep(4)
    elif not snake_args["use_conda"]:
        warnings.warn(
            "'use_conda' is set to False by user. It is expected that user has all dependencies installed."
        )
        time.sleep(4)

    # Set conda_frontend='mamba' if not supplied
    if "conda_frontend" not in snake_args.keys():
        snake_args["conda_frontend"] = "mamba"
    elif snake_args["conda_frontend"] == "conda":
        warnings.warn(
            "'conda_frontend' is set to 'conda' by user. This will lead to slower operations, consider using "
            "'mamba' instead."
        )
        time.sleep(4)

    if verify:
        # Make DAG and perform dry-run
        good_exit = snk.snakemake(
            snake_path,
            workdir=run_dir,
            config=config,
            cores=threads,
            dryrun=True,
            **snake_args
        )

        pathlib.Path(run_dir).mkdir(parents=True, exist_ok=True)
        out = io.StringIO()
        with redirect_stdout(out):
            snk.snakemake(
                snake_path,
                workdir=run_dir,
                printdag=True,
                config=config,
                cores=threads,
                **snake_args
            )
            out = out.getvalue()
            out_file = run_dir + "/dag.gv"

            if os.path.exists(out_file):
                os.remove(out_file)

            with open(out_file, "a") as stdout_log:
                stdout_log.writelines(out)

            out_png = run_dir + "/dag.png"
            os.system("cat " + out_file + " | dot -Tpng -o " + out_png)
            os.remove(out_file)
        if good_exit:
            return out_png
        else:
            sys.exit(1)
    else:
        # Run snakemake
        good_exit = snk.snakemake(
            snake_path, workdir=run_dir, config=config, cores=threads, **snake_args
        )

        # Check exit status
        if good_exit:
            sys.exit(0)
        else:
            sys.exit(1)
