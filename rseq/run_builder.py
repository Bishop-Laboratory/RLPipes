import os
import json
import random
import string
import requests
import time
import ast
import re
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
        manual = False
        error = None
        if 'manual_toggle' in request.form.keys() and request.form['manual_toggle'] == 'on':
            manual = True
            to_parse = request.form['hiddeninput'].replace('null', '""')
            res = ast.literal_eval(to_parse)
            print(res)
            sample_sheet = pd.DataFrame(res, columns=['sample_name', 'experiment', 'control', 'genome']).dropna()
            exp_lst = [exp for exp in sample_sheet['experiment'] if exp != ""]
            if not exp_lst:
                error = "'experiment' column must not be empty"
            # Drop columns with empty column names
            sample_sheet = sample_sheet[[col for col in sample_sheet.columns if col != ""]]
            print(sample_sheet)
        else:
            file = request.files['sample_sheet']
            if not run_name:
                error = 'run_name is required.'
            elif not allowed_file(file.filename):
                error = '{} is not an allowed filename'.format(file.filename)

        if error is None:
            run_str = random_string(50)
            run_path = os.path.join(current_app.config['UPLOAD_FOLDER'], run_str)
            os.makedirs(run_path, exist_ok=True)
            sample_sheet_path = os.path.join(run_path, "samples.csv")
            if not manual:
                file.save(sample_sheet_path)
            else:
                sample_sheet.to_csv(sample_sheet_path, index=False)

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
    return redirect('/')


@bp.route('/<int:run_id>/preflight')
def pre_flight(run_id):
    # Modifier allows other requests to use pre_flight without side-effects like db changes and redirection
    modifier = request.args.getlist('modify')

    run = get_run(run_id)
    _execute(run_id, dryrun=True)
    time.sleep(1)
    error_file = run.run_path + '/logs-dryrun/' + 'error.log'
    if os.path.exists(error_file):
        error = json.load(open(error_file))
        msg = "Error encountered. Check logs for more info.\n" + error['0']['msg']
        flash(msg, category="danger")
    else:
        _execute(run_id, dryrun=False, dag=True)
        if not modifier:
            run.status = '<strong class="text-info">Ready for takeoff 🛫</strong>'
            db_session.commit()
            msg = "preflight checks for run: '" + run.run_name + "' succeeded! Ready for take off!."
            flash(msg, category="success")

    if not modifier:
        return redirect("/")
    else:
        return json.dumps(True)


def modify_digraph(digraph, completed_jobs):
    """
    Modifies a graphviz digraph object to reflected completed snakemake jobs
    :param digraph: digraph to modify
    :param completed_jobs: numeric list with the completed jobs inside
    :return: Updated digraph
    """
    for job in completed_jobs:
        digraph = re.sub(r'(\b%s\[label = "[a-zA-Z0-9: _\\n]+", color = "[0-9\\. ]+", style=")rounded("\];)' % str(job),
                          r"\1rounded,dashed\2", digraph)
    return digraph


@bp.route('/<int:run_id>/track_flight')
def track_flight(run_id):
    """Checks snakemake logs and updates run status"""

    # TODO: Also collect the unfinished jobs (CASE: someone deletes part of the run's file system)
    # TODO: Perhaps also rerun dryrun periodically

    run = get_run(run_id)
    samples = request.args.getlist('samples[]')
    msg_now = ""
    fail = False
    complete = False
    status_dict = {}
    for sample in samples:

        dag_file = run.run_path + "/" + sample + "/dag.gv"

        with open(dag_file) as f:
            gv = f.read().replace('\n', '')

        job_info_log = run.run_path + "/logs/" + sample + "/job_info.log"
        if not os.path.exists(job_info_log):
            # CASE: only dry-run so far
            print("Only dry")
            # Get the first entry in the run_info (the first step for when the run starts)
            job_info_log = run.run_path + "/logs-dryrun/" + sample + "/job_info.log"
            job_info = json.load(open(job_info_log))
            key_int = min([int(x) for x in job_info.keys()])
            msg_now = "Ready for launch!"
        else:
            # CASE: Some progress already made
            job_info = json.load(open(job_info_log))
            # Get the last entry in the run_info (the ind for the current step)
            key_int = max([int(x) for x in job_info.keys()])

        # Get the run info
        run_info_log = run.run_path + "/logs/" + sample + "/run_info.log"
        if os.path.exists(run_info_log):
            run_info_log = json.load(open(run_info_log))

        else:
            run_info_log = {}

        print(run_info_log)

        # Find the finished jobs and modify the digraph
        job_finished_log = run.run_path + "/logs/" + sample + "/job_finished.log"
        if os.path.exists(job_finished_log):
            job_finished_log = json.load(open(job_finished_log))
            finished_jobs = [int(val['jobid']) for key, val in job_finished_log.items()]
            gv = modify_digraph(gv, completed_jobs=finished_jobs)

        print(run.status)
        if run.status == '<strong class="text-info">Ready for takeoff 🛫</strong>':
            msg_now = "ready for takeoff"
            complete = False
        elif run.status == '<strong class="text-info">Safely aborting</strong>':
            print("Aborting")
            print(run_info_log.keys())
            reset = False
            complete = False
            last_entry = ""
            print(len(run_info_log.keys()))
            if len(run_info_log.keys()) == 0:
                # CASE: mNo run info log keys
                reset = True
            else:
                last_entry = run_info_log[str(max([int(x) for x in run_info_log.keys()]))]['msg']

            if last_entry[:13] == 'Complete log:':
                print("Aborted")

                # Sends a get request to update the DAG
                requests.get(request.url_root + '/runs/' + str(run_id) + '/preflight?modify=true')

                run.status = '<strong class="text-info">Ready for takeoff 🛫</strong>'
                db_session.commit()

                msg = "Run '" + run.run_name + "' successfully aborted."
                flash(msg, 'info')
                msg_now = "aborted"
            elif reset and msg_now != "aborting":
                # CASE: True need for reset
                print("Needs reset")
            else:
                print("Still aborting")
                msg_now = "aborting"
        elif run.status not in ['<strong class="text-danger">Failed! ❌</strong>', '<strong class="text-info">Complete! ✔️</strong>']:
            # CASE: either the run is still in flight, its finished, or it failed
            completed_log = run.run_path + "/logs/" + sample + "/completed.log"  # This will exist if run has finished
            failed_log = run.run_path + "/logs/" + sample + "/failed.log"  # This will exist if run failed with non-zero exit code
            if os.path.exists(failed_log):
                # CASE: it failed
                msg = "Run '" + run.run_name + "' failed. Please see error logs for " \
                      + "additional details. Re-run after addressing errors. Please also contact package maintainer."
                flash(msg, 'danger')
                fail = True
                complete = False
                msg_now = "failed"
            elif os.path.exists(completed_log):
                # CASE: it finished
                msg = "Run '" + run.run_name + "' complete!"
                flash(msg, 'success')
                complete = True
                msg_now = "completed"
            else:
                # CASE: still flying
                complete = False
                msg_now = "In flight"

        # Get current step
        run_now = job_info[str(key_int)]['name']

        # Gather into a response
        status_dict[sample] = {
            "current_step": run_now,
            "current_step_ind": key_int,
            "DAG": gv
        }

    tracker = {'run_status': run.status,
               'status_dict': status_dict,
               'msg': msg_now}

    if fail:
        run.status = '<strong class="text-danger">Failed! ❌</strong>'
        db_session.commit()
    elif complete:
        run.status = '<strong class="text-info">Complete! ✔️</strong>'
        db_session.commit()

    return json.dumps(tracker)


@bp.route('/<int:run_id>/monitor')
def monitor(run_id):
    run = get_run(run_id)
    run_sheet = run.run_path + '/samples.init.json'
    config_dict = json.load(open(run_sheet))
    # Check on status of each sample
    for sample, config in config_dict.items():
        dry_run = run.run_path + "/" + sample + "/dryrun.log"

        # Add some additional info
        with open(dry_run) as f_dryrun:
            dry_run_lines = f_dryrun.readlines()

        if 'Nothing to be done.\n' in dry_run_lines:
            config_dict[sample]['current_step'] = 'done'
            config_dict[sample]['status'] = 'Complete'
            config_dict[sample]['card_color'] = 'success'
            libstr = ""
            if config['paired_end']:
                libstr += "Paired-end"
            else:
                libstr += "Single-end"
            libstr += " (" + str(config['read_length']) + ")"
            if config['strand_specific']:
                libstr += ", Strand-specific"
            else:
                libstr += ", Unstranded"
            config_dict[sample]['library'] = "Strand-specific"
        else:
            # Get number of completed jobs
            r_dryrun = re.compile('^\t\d+\n$')
            num_complete = list(filter(r_dryrun.match, dry_run_lines))
            num_complete = re.sub(pattern='^\t(\d+)\n$', string=num_complete[0], repl='\\1')
            if num_complete == 0:
                config_dict[sample]['status'] = 'Ready'
                config_dict[sample]['card_color'] = 'primary'
                config_dict[sample]['progress'] = "0%"
                config_dict[sample]['current_step'] = 'None'
            else:
                dag_file = run.run_path + "/" + sample + "/dag.gv"
                with open(dag_file) as f_dag:
                    dag_lines = f_dag.readlines()
                r_dag = re.compile(r'.+\d+\[label.+')
                num_total = list(filter(r_dag.match, dag_lines))
                num_total = max([int(re.sub(pattern=r'\t([\d]+)\[label.+\n', repl="\\1", string=x)) for x in num_total])
                config_dict[sample]['current_step'] = dry_run_lines[4][3:]
                config_dict[sample]['status'] = 'In progress'
                config_dict[sample]['card_color'] = 'info'
                print(num_complete)
                print(num_total)
                prog = 100 * int(num_complete)/num_total
                config_dict[sample]['progress'] = str(prog.__round__()) + "%"
    from flask import jsonify
    return render_template('run_info_page.html', run_id=run_id, sample_dict=config_dict,
                           run_name=run.run_name, run_status=run.status)


@bp.route('/<int:run_id>/takeoff')
def execute_run(run_id):
    # check for re-run request
    print(request.args)
    force = False
    if 'rerun' in request.args:
        force = True
    _execute(run_id, force=force)
    run = get_run(run_id)
    run.status = '<strong class="text-info">In flight ✈️</strong>'
    db_session.commit()
    return json.dumps({'run_status': run.status})


@bp.route('/<int:run_id>/abort')
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
    run = get_run(run_id)
    run.status = '<strong class="text-info">Safely aborting</strong>'
    db_session.commit()
    return json.dumps({'run_status': run.status})


def _execute(run_id, dryrun=False, dag=False, force=False):
    run = get_run(run_id)
    run_sheet = run.run_path + '/samples.init.json'
    config_dict = json.load(open(run_sheet))
    config_files = []
    for sample, config in config_dict.items():
        sample_path = os.path.join(run.run_path, sample)
        os.makedirs(sample_path, exist_ok=True)
        config_file = os.path.join(sample_path, 'config.json')
        json.dump(config, open(config_file, 'w'))
        config_files.append(config_file)

    args = ['python', 'rseq/make_snakes.py']
    args.extend(config_files)
    args.extend([run.run_path, str(dryrun), str(dag), str(force)])
    print(args)
    proc = subprocess.Popen(args=args)
    session['run_data'] = {run_id: proc.pid}

    return json.dumps(True)
