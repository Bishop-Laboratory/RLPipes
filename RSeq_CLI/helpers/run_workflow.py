import snakemake as snk
import json
import sys
import io
import os
import pathlib
from contextlib import redirect_stdout
import tempfile


def make_snakes(config_file, run_path, snake_path, dryrun=False, dag=False, force=False, notemp=False):

    if dryrun:
        print("Dry run")
    elif dag:
        print("Dag run")

    configs = json.load(open(config_file))

    for sample_name, config in configs.items():

        if not dryrun and not dag:
            # Unlock any previous runs TODO: this should be removed probably... We need unique dirs for multiple runs!
            snk.snakemake(snake_path, unlock=True, dryrun=dryrun, config=config, forceall=force, notemp=notemp)

        cores = config['cores'][0]

        pathlib.Path(run_path + "/dags/").mkdir(parents=True, exist_ok=True)

        if dag:
            out = io.StringIO()
            with redirect_stdout(out):
                snk.snakemake(snake_path, dryrun=dryrun, printdag=dag,
                              config=config, cores=cores, forceall=force, notemp=notemp)
                out = out.getvalue()

                out_file = run_path + '/dags/' + sample_name + '.dag.gv'

                if os.path.exists(out_file):
                    os.remove(out_file)

                with open(out_file, 'a') as stdout_log:
                    stdout_log.writelines(out)

                out_svg = run_path + '/dags/' + sample_name + '.dag.png'
                os.system('cat ' + out_file + ' | dot -Tpng -o ' + out_svg)
                os.remove(out_file)
        else:
            print(force)
            snk.snakemake(snake_path, dryrun=dryrun, printdag=dag, printreason=True,
                          config=config, cores=cores, forceall=force, notemp=notemp)

    print("OUT OF SNAKE")


if __name__ == "__main__":
    print(sys.argv)
    make_snakes(config_file=sys.argv[1], run_path=sys.argv[2], snake_path=sys.argv[3],
                dryrun=(sys.argv[4] == 'True'), dag=(sys.argv[5] == 'True'),
                force=(sys.argv[6] == 'True'), notemp=(sys.argv[7] == '--keepTmp'))
