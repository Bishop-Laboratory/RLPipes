import snakemake as snk
import json
import sys
import io
import os
import time
import pathlib
from contextlib import redirect_stdout


def make_snakes(config_file, run_path, snake_path, dryrun=False, dag=False, force=False, notemp=False, reason=False):

    configs = json.load(open(config_file))
    for sample_name, config in configs.items():
        if sample_name in ['dryrun', 'keepTmp', 'dag', 'force', 'reason', 'cores']:
            continue

        print("\nCurrent sample: '" + sample_name + "'\n")
        if len(config.keys()) < 20:
            time.sleep(.5)

        if not dryrun and not dag:
            # Unlock any previous runs? TODO: this should be removed probably... We need unique dirs for multiple runs!
            cores = int(config['cores'][0])
            snk.snakemake(snake_path, unlock=True, dryrun=dryrun, config=config, forceall=force,
                          printreason=reason, notemp=notemp, cores=cores)
            snk.snakemake(snake_path, unlock=False, dryrun=dryrun, config=config, forceall=force,
                          printreason=reason, notemp=notemp, cores=cores)
        else:
            pathlib.Path(run_path + "/dags/").mkdir(parents=True, exist_ok=True)

            if dag:
                out = io.StringIO()
                with redirect_stdout(out):
                    snk.snakemake(snake_path, dryrun=dryrun, printdag=dag, printreason=reason,
                                  config=config, forceall=force, notemp=notemp)
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
                snk.snakemake(snake_path, dryrun=dryrun, printdag=dag, printreason=reason,
                              config=config, forceall=force, notemp=notemp)


if __name__ == "__main__":
    configs = json.load(open(sys.argv[1]))
    notemp = (not configs['keepTmp'][0])
    dryrun = configs['dryrun'][0]
    dag = configs['dag'][0]
    force = configs['force'][0]
    reason = configs['reason'][0]

    make_snakes(config_file=sys.argv[1], run_path=os.path.dirname(sys.argv[1]), snake_path=sys.argv[2],
                dryrun=dryrun, dag=dag, force=force, notemp=notemp, reason=reason)
