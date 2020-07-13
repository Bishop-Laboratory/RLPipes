import os
import json
import random
import string
from flask import (
    Blueprint, flash, redirect, render_template, session, url_for, request, current_app
)
from werkzeug.utils import secure_filename
from rseq.models import Run
from rseq.database import db_session
import pandas as pd
import numpy as np
from rpy2.robjects.packages import importr
import subprocess

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
            run_path = os.path.join(current_app.config['UPLOAD_FOLDER'], run_str)
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
            file.save(os.path.join(current_app.config['UPLOAD_FOLDER'], filename))

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
    rseq_r = importr('RSeq')
    run = get_run(run_id)
    sample_sheet_path = run.run_path + "/samples.csv"
    sample_sheet_path_initialized = run.run_path + "/samples.init.json"
    rseq_r.initialize_run(samples=sample_sheet_path, mode=run.mode,
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


@bp.route('/<int:run_id>/monitor')
def monitor(run_id):
    return render_template('run_info_page.html',
                           run_id=run_id)


def _execute(run_id, dryrun=False, force=False, dag=False, cores=-1):
    run = get_run(run_id)
    run_sheet = run.run_path + '/samples.init.json'
    config_dict = json.load(open(run_sheet))
    if cores == -1:
        import multiprocessing
        cores = multiprocessing.cpu_count()

    cmds = []
    config_files = []
    for sample, config in config_dict.items():
        sample_path = os.path.join(run.run_path, sample)
        os.makedirs(sample_path, exist_ok=True)
        config_file = os.path.join(sample_path, 'config.json')
        print(config_file)
        config_files.append(config_file)
        dag_file = os.path.join(sample_path, 'dag.svg')

        json.dump(config, open(config_file, 'w'))
        snake_file = 'rseq/rseq.smk'
        cmd = "snakemake -s " + snake_file + " --configfile " + config_file
        if dryrun:
            cmd += " -np"
            if force:
                cmd += " -R output"
        elif dag:
            cmd += " --dag | dot -Tsvg > " + dag_file

        cmds.append(cmd)

    cmd_final = " && ".join(cmds)

    if dryrun or dag:
        os.system(cmd_final)
    else:
        args = ['python', 'rseq/make_snakes.py']
        args.extend(config_files)
        proc = subprocess.Popen(args=args)
        session['run_data'] = {run_id: proc.pid}

    # TODO: Get the module names, log files, and status from the dry-run with parser function
    # with subprocess.Popen(cmd_final, shell=True, stdout=subprocess.PIPE) as proc:
    #     stdout_log = open(run.run_path + '/output_log.txt', 'wb')
    #     error_log = open(run.run_path + '/error_log.txt', 'wb')
    #     stdout_log.write(proc.stdout.read())
    #     # error_log.write(proc.stderr.read())

    return json.dumps(True)


@bp.route('/<int:run_id>/terminate_run')
def terminate_run(run_id):
    """Kills an active snakemake run

    When a run is created, a PID is assigned to the global session.
    This function takes the PID, gets the session ID, and uses that to send a kill signal
    to the whole process tree underneath the PID.

    REQUIRES that subprocess.Popen not be used in shell mode (works best with a python script)
    e.g., subprocess.Popen(['python', 'make_snakes.py', 'config_file.json'])

    :param run_id: Id of run to kill
    :return: redirects to monitor page
    """
    # TODO: Add fast/unsafe kill method
    pid = session['run_data'][str(run_id)]
    cmd = 'kill $(ps -s $$ -o pid={pid})'.format(pid=pid)
    os.system(cmd)
    return redirect(url_for('runs.monitor', run_id=run_id))


@bp.route('/<int:run_id>/execute_run')
def execute_run(run_id):
    print("EXECURTIN")
    _execute(run_id)
    return redirect(url_for('runs.monitor', run_id=run_id))


@bp.route('/<int:run_id>/execute_dryrun')
def execute_dryrun(run_id):
    _execute(run_id, dryrun=True, force=True)
    return redirect(url_for('runs.monitor', run_id=run_id))


@bp.route('/<int:run_id>/execute_dag')
def execute_dag(run_id):
    _execute(run_id, dryrun=False, force=False, dag=True)
    return redirect(url_for('runs.monitor', run_id=run_id))
