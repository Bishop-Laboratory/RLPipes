import os
import json
import random
import string
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, current_app as app
)
from werkzeug.utils import secure_filename
from rseq.models import Run, Pipeline
from rseq.database import db_session
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
import snakemake as snk
import sqlite3
from contextlib import redirect_stdout
import io
import threading
LOCK = threading.Lock()


bp = Blueprint('runs', __name__, url_prefix='/runs')
ALLOWED_EXTENSIONS = {'txt', 'csv', 'tsv', 'tab', 'xls', 'xlsx'}


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


def random_string(string_length=8):
    # From https://pynative.com/python-generate-random-string/
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(string_length))


def get_run(run_id):
    """
    Get the info for an RSeq run.

    :param run_id: id of the run
    :return: run information from database
    """
    run = db_session.query(Run).filter_by(id=run_id).first()
    return run


@bp.route('/create', methods=('GET', 'POST'))
def create_run():
    if request.method == 'POST':
        run_name = request.form['run_name']
        mode = request.form['mode']
        file = request.files['sample_sheet']
        error = None
        if not run_name:
            error = 'run_name is required.'
        elif not file:
            error = "file is required"
        elif not allowed_file(file.filename):
            error = '{} is not an allowed filename'.format(file.filename)

        if file and allowed_file(file.filename) and error is None:
            run_str = random_string(50)
            run_path = os.path.join(app.config['UPLOAD_FOLDER'], run_str)
            os.makedirs(run_path, exist_ok=True)
            sample_sheet_path = os.path.join(run_path, "samples.csv")
            file.save(sample_sheet_path)
            run_now = Run(run_name=run_name,
                          run_path=run_path,
                          mode=mode,
                          status='<strong class="text-muted">uninitialized</strong>')
            db_session.add(run_now)
            db_session.commit()
            message = "Successfully created new run " + run_name
            flash(message, category="success")

            return redirect("/")

        flash(error, category="danger")

    return render_template('create_run.html')


@bp.route('/<int:run_id>/edit', methods=('GET', 'POST'))
def edit_run(run_id):

    run = get_run(run_id)
    if request.method == 'POST':
        run_name = request.form['run_name']
        file = request.files['sample_sheet']
        error = None

        if not run_name:
            error = 'run_name is required.'
        elif not file:
            error = "file is required"
        elif not allowed_file(file.filename):
            error = '{} is not an allowed filename'.format(file.filename)

        if file and allowed_file(file.filename) and error is None:
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))

            # Update run object
            run.run_name = run_name
            db_session.commit()
            message = "Successfully edited run " + run_name
            flash(message, category="success")

            return redirect("/")

        flash(error, category="danger")

    return render_template('edit_run.html', run=run)


@bp.route('/<int:run_id>/delete', methods=('GET',))
def delete_run(run_id):
    run = get_run(run_id)
    msg = "Deleted run " + run.run_name
    db_session.delete(run)
    db_session.commit()
    flash(msg, category="warning")
    return redirect("/")


@bp.route('/<int:run_id>/sample_sheet', methods=['GET', 'POST'])
def sample_sheet_edit(run_id):
    run = get_run(run_id)
    sample_sheet_path = run.run_path + "/samples.csv"
    table_data = pd.read_csv(sample_sheet_path)
    table_data = table_data.replace(np.nan, '', regex=True)
    table_data = table_data.replace('[]', '')
    table_data_t = table_data.transpose()
    return render_template('sample_sheet.html', data=table_data_t.to_dict(),
                           colnames=table_data.columns.values.tolist(), run_id=run_id)


@bp.route('/<int:run_id>/sample_sheet/save', methods=['GET', 'POST'])
def sample_sheet_save(run_id):
    # Get the table from save
    # TODO: Issue with save function -- conflict with empty lists
    json_data = request.get_json()
    new_sample_sheet = pd.DataFrame(json_data[1:], columns=json_data[0])
    new_sample_sheet = new_sample_sheet.replace('[]', '')
    # Overwrite original
    run = get_run(run_id)
    table_path = run.run_path + "/samples.csv"
    new_sample_sheet.to_csv(table_path, index=False)

    return json.dumps(True)


@bp.route('/<int:run_id>/initialize')
def initialize_run(run_id):
    # Takes the sample sheet and uses R to build the final run info
    RSeq = importr('RSeq')
    run = get_run(run_id)
    sample_sheet_path = run.run_path + "/samples.csv"
    sample_sheet_path_initialized = run.run_path + "/samples.init.json"
    RSeq.initialize_run(samples=sample_sheet_path, mode=run.mode,
                        output_json=sample_sheet_path_initialized)
    # # Fix issue with formatting -- only for multi-sample mode
    # samples = json.load(open(run.sample_sheet_path_initialized))
    # samples = json.loads(samples[0])
    # json.dump(samples, open(run.sample_sheet_path_initialized, "w"))
    # finish init and return to dashboard
    run.status = '<strong class="text-info">Initialized. Ready for pre-flight ☑</strong>'
    db_session.commit()
    message = "Successfully edited run " + run.run_name
    flash(message, category="success")
    return redirect('/runs/' + str(run_id) + '/monitor')
    # Get the config file

    # Initialize sample db
    # run.sample_database = os.path.splitext(run.sample_sheet_path)[0] + ".init.db"
    # conn = sqlite3.connect(run.sample_database)
    # # sample_db = '/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/instance/uploads/lzjuoygzegsnaqcxfkhw__xXxXx__sample_sheet_example.init.db'
    # # conn = sqlite3.connect(sample_db)
    # c = conn.cursor()
    # c.execute('''
    #     CREATE TABLE IF NOT EXISTS samples (
    #     sample_name TEXT PRIMARY KEY,
    #     rseq_params TEXT,
    #     dag TEXT,
    #     run_snakemake TEXT,
    #     progress TEXT,
    #     log TEXT,
    #     status TEXT,
    #     args TEXT,
    #     targets TEXT,
    #     rule_info TEXT,
    #     resources TEXT
    # );
    # '''.replace("\n", ""))

    # # Populate
    # samples = json.load(open(run.sample_sheet_path_initialized))
    # # samples = json.load(open('/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/instance/uploads/gdroxmmczhctzjlurnut__xXxXx__sample_sheet_example_-_Copy.init.json'))
    # samples = json.loads(samples[0])


    # sample_tuples = []
    # for sample_name, rseq_params in samples.items():
    #     sample_tuple = (sample_name, json.dumps(rseq_params), None, None, "", "", '{"running": false}', None, "", "", "")
    #     sample_tuples.append(sample_tuple)
    #
    # c.executemany('INSERT INTO samples VALUES (?,?,?,?,?,?,?,?,?,?,?)', sample_tuples)
    # conn.commit()
    # conn.close()




@bp.route('/<int:run_id>/terminate')
def terminate_run(run_id):
    # TODO: Terminate run function
    return json.dumps(True)





# def register(run_snakemake, args, sample):
#
#     sample["run_snakemake"] = run_snakemake
#     sample["args"] = json.dumps(dict(
#         targets=args.target,
#         cluster=args.cluster,
#         workdir=args.directory,
#         touch=args.touch,
#         forcetargets=args.force,
#         forceall=args.forceall,
#         forcerun=args.forcerun,
#         prioritytargets=args.prioritize,
#         stats=args.stats,
#         keepgoing=args.keep_going,
#         jobname=args.jobname,
#         immediate_submit=args.immediate_submit,
#         ignore_ambiguity=args.allow_ambiguity,
#         lock=not args.nolock,
#         force_incomplete=args.rerun_incomplete,
#         ignore_incomplete=args.ignore_incomplete,
#         jobscript=args.jobscript,
#         notemp=args.notemp,
#         latency_wait=args.latency_wait,
#     ))
#
#     target_rules = []
#
#     def log_handler(msg):
#         if msg["level"] == "rule_info":
#             target_rules.append(msg["name"])
#
#     run_snakemake(list_target_rules=True, log_handler=log_handler)
#     for target in args.target:
#         target_rules.remove(target)
#     sample["targets"] = args.target + target_rules
#
#     resources = []
#
#     def log_handler(msg):
#         if msg["level"] == "info":
#             resources.append(msg["msg"])
#
#     run_snakemake(list_resources=True, log_handler=log_handler)
#     sample["resources"] = resources
#     sample["snakefilepath"] = os.path.abspath(args.snakefile)
#
#
# def run_snakemake(**kwargs):
#     args = dict(sample["args"])
#     args.update(kwargs)
#     sample["run_snakemake"](**args)


@bp.route('/<int:run_id>/monitor')
def monitor(run_id):
    # run = get_run(run_id)
    # conn = sqlite3.connect(run.sample_database)
    # sample_db = '/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/instance/uploads/lzjuoygzegsnaqcxfkhw__xXxXx__sample_sheet_example.init.db'
    # conn = sqlite3.connect(sample_db)
    # conn.row_factory = sqlite3.Row
    # conn.row_factory = lambda cursor, row: row[0]
    # c = conn.cursor()
    # targets = c.execute("SELECT targets FROM samples").fetchall()

    # conn.close()
    return render_template('run_info_page.html',
                           run_id=run_id)


# @bp.route("/<int:run_id>/sample/<string:sample_name>/dag")
# def dag(run_id, sample_name):
#     run = get_run(run_id)
#     samples = json.load(open(run.sample_sheet_path_running))
#     sample = samples[sample_name]
#     if sample["dag"] is None:
#         def record(msg):
#             if msg["level"] == "d3dag":
#                 sample["dag"] = msg
#             elif msg["level"] in ("error", "info"):
#                 sample["log"].append(msg)
#             samples[sample_name] = sample
#
#         run_snakemake(printd3dag=True, log_handler=record)
#
#     return json.dumps(sample["dag"])
#
#
# @bp.route("/<int:run_id>/sample/<string:sample_name>/log/<int:log_id>")
# def log(run_id, sample_name, log_id):
#     run = get_run(run_id)
#     samples = json.load(open(run.sample_sheet_path_running))
#     sample = samples[sample_name]
#     logs = sample["log"][log_id:]
#     return json.dumps(logs)
#
#
# @bp.route("/<int:run_id>/sample/<string:sample_name>/progress")
# def progress(run_id, sample_name):
#     run = get_run(run_id)
#     samples = json.load(open(run.sample_sheet_path_running))
#     sample = samples[sample_name]
#     return json.dumps(sample["progress"])


def check_logs(run_id):

    run = get_run(run_id)
    log_file = run.log_file


    with sqlite3.connect(run.sample_database) as conn:
        # sample_db = '/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/instance/uploads/lzjuoygzegsnaqcxfkhw__xXxXx__sample_sheet_example.init.db'
        # conn = sqlite3.connect(sample_db)
        c = conn.cursor()
        if level == "progress":
            print("Progress")
            c.execute("UPDATE samples SET progress = ? WHERE sample_name = ?", (msg, sample))
        elif level in ("info", "error", "job_info", "job_finished"):
            print("LOginn")
            prev_log = c.execute("SELECT log FROM samples WHERE sample_name = ?", (sample,)).fetchone()
            prev_log = prev_log[0]
            print(prev_log)
            if prev_log == "":
                current_log = {0: json.dumps(msg)}
            else:
                top_ind = prev_log.keys().max
                current_ind = top_ind + 1
                current_log = prev_log
                current_log[current_ind] = json.dumps(msg)

            print(current_log)
            # sample["log"].append(msg)
            c.execute("UPDATE samples SET log = ? WHERE sample_name = ?", (json.dumps(current_log), sample))
        conn.commit()
        conn.close()


def _execute(run_id, dryrun=False, force=False, dag=False):
    run = get_run(run_id)
    run_sheet = run.run_path + '/samples.init.json'
    config_dict = json.load(open(run_sheet))
    for sample, config in config_dict.items():
        sample_path = os.path.join(run.run_path, sample)
        os.makedirs(sample_path, exist_ok=True)
        config_file = os.path.join(sample_path, 'config.json')
        dag_file = os.path.join(sample_path, 'dag.svg')
        json.dump(config, open(config_file, 'w'))
        snake_file = 'modules/rseq.smk'
        cmd = "snakemake -s " + snake_file + " --configfile " + config_file
        if dryrun:
            cmd += " -np"
            if force:
                cmd += " -R output"
        elif dag:
            cmd += " --dag | dot -Tsvg > " + dag_file
        os.system(cmd)

    return json.dumps(True)

# For each sample in run, set up the log handler and execute snakemake
    # for sample, params in samples.items():
    #     print(sample)
    #     print(sample)
    #     # Clear previous progress and logs if they exist
    #     print("resetting logs")
    #     conn = sqlite3.connect(run.sample_database)
    #     # sample_db = '/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/instance/uploads/lzjuoygzegsnaqcxfkhw__xXxXx__sample_sheet_example.init.db'
    #     # conn = sqlite3.connect(sample_db)
    #     # conn.row_factory = sqlite3.Row
    #     # conn.row_factory = lambda cursor, row: row[0]
    #     c = conn.cursor()
    #     c.execute("UPDATE samples SET progress = ? WHERE sample_name = ?", ('', sample))
    #     c.execute("UPDATE samples SET log = ? WHERE sample_name = ?", ('', sample))
    #     conn.commit()
    #     conn.close()




        # snake_file = 'modules/Snakefile_multiInput'
        # cmd = "snakemake -s " + snake_file
        # if dryrun:
        #     cmd += " -np"
        #
        # if force:
        #     cmd += " -R output"
        # print(cmd)
        # os.system(cmd)
        # snk.snakemake(snakefile='modules/Snakefile_multiInput', dryrun=dryrun)

        # return json.dumps(True)



    # with LOCK:
    #     sample["status"]["running"] = True
    # run_snakemake(dryrun=dryrun)
    # with LOCK:
    #     sample["status"]["running"] = False


@bp.route('/<int:run_id>/execute_run')
def execute_run(run_id):
    # snk.snakemake(snakefile='modules/Snakefile_multiInput', dryrun=dryrun)

    _execute(run_id)
    return json.dumps(True)


@bp.route('/<int:run_id>/execute_dryrun')
def execute_dryrun(run_id):
    # snk.snakemake(snakefile='modules/Snakefile_multiInput', dryrun=True)

    _execute(run_id, dryrun=True, force=True)
    return json.dumps(True)

@bp.route('/<int:run_id>/execute_dag')
def execute_dag(run_id):
    # snk.snakemake(snakefile='modules/Snakefile_multiInput', dryrun=True)
    _execute(run_id, dryrun=False, force=False, dag=True)
    return json.dumps(True)


@bp.route("/<int:run_id>/sample/<string:sample_name>/status")
def status(run_id):
    with LOCK:
        return json.dumps(sample["status"])


@bp.route("/<int:run_id>/sample/<string:sample_name>/targets")
def targets(run_id):
    return json.dumps(sample["targets"])


@bp.route("/<int:run_id>/sample/<string:sample_name>/get_args")
def get_args(run_id):
    return json.dumps(sample["args"])


@bp.route("/<int:run_id>/sample/<string:sample_name>/set_args", methods=["POST"])
def set_args():
    sample["args"].update(
        {name: value for name, value in request.form.items() if not name.endswith("[]")}
    )
    targets = request.form.getlist("targets[]")
    if targets != sample["args"]["targets"]:
        sample["dag"] = None
    sample["args"]["targets"] = targets
    return ""


def make_snake_config(run_id):
    run = get_run(run_id)
    # sample_json = '/home/UTHSCSA/millerh1/Bishop.lab/Projects/RSeq/instance/uploads/lzjuoygzegsnaqcxfkhw__xXxXx__sample_sheet_example.init.json'
    # samples = json.load(open(sample_json))
    samples = json.load(open(run.sample_sheet_path_initialized))
    samp_df = pd.DataFrame.from_dict(samples)
    samp_df = samp_df.transpose()
    samp_df = samp_df.applymap(lambda x: x if not isinstance(x, list) else x[0] if len(x) else '')
    samp_dict = pd.DataFrame.to_dict(samp_df)
    config_path = os.path.splitext(run.sample_sheet_path)[0] + ".config.json"
    json.dump(samp_dict, config_path)
    return config_path




