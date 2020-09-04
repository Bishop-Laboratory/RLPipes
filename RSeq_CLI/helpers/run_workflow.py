import snakemake as snk
import json
import sys
import io
import os
import pathlib
from contextlib import redirect_stdout


def make_snakes(config_file, run_path, snake_path, dryrun=False, dag=False, force=False):

    if dryrun:
        print("Dry run")
    elif dag:
        print("Dag run")

    configs = json.load(open(config_file))

    for sample_name, config in configs.items():

        if not dryrun and not dag:
            # Unlock any previous runs TODO: this should be removed probably... We need unique dirs for multiple runs!
            snk.snakemake(snake_path, unlock=True, dryrun=dryrun, config=config, forceall=force)

        cores = config['cores'][0]

        pathlib.Path(run_path + '/' + config['sample_name'][0]).mkdir(parents=True, exist_ok=True)

        if dag:
            out = io.StringIO()
            with redirect_stdout(out):
                snk.snakemake(snake_path, force_incomplete=True, dryrun=dryrun, printdag=dag,
                              config=config, cores=cores, forceall=force)
                out = out.getvalue()

                out_file = run_path + '/' + config['sample_name'][0] + '/dag.gv'

                if os.path.exists(out_file):
                    os.remove(out_file)

                with open(out_file, 'a') as stdout_log:
                    stdout_log.writelines(out)

                out_svg = run_path + '/' + config['sample_name'][0] + '/dag.gv.svg'
                os.system('cat ' + out_file + ' | dot -Tsvg -o ' + out_svg)
        else:
            snk.snakemake(snake_path, force_incomplete=True, dryrun=dryrun, printdag=dag,
                          config=config, cores=cores, forceall=force)

    print("OUT OF SNAKE")


if __name__ == "__main__":
    print(sys.argv)
    make_snakes(config_file=sys.argv[1], run_path=sys.argv[2], snake_path=sys.argv[3],
                dryrun=(sys.argv[4] == 'True'), dag=(sys.argv[5] == 'True'),
                force=(sys.argv[6] == 'True'))
