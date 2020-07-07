import os
import random
import string
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, current_app as app
)
from werkzeug.utils import secure_filename
from rseq.models import Run, Pipeline
from rseq.database import db_session
import pandas as pd

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
        file = request.files['sample_sheet']
        print(file.filename)
        error = None

        if not run_name:
            error = 'run_name is required.'
        elif not file:
            error = "file is required"
        elif not allowed_file(file.filename):
            error = '{} is not an allowed filename'.format(file.filename)

        if file and allowed_file(file.filename) and error is None:
            filename = random_string(20) + "__xXxXx__" + secure_filename(file.filename)
            sample_sheet_path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            print(sample_sheet_path)
            file.save(sample_sheet_path)
            run_now = Run(run_name=run_name,
                          sample_sheet=file.filename,
                          sample_sheet_path=sample_sheet_path)
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
    print(run)
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


@bp.route('/<int:run_id>/sample_sheet', methods=['GET', 'POST'])
def sample_sheet(run_id):
    run = get_run(run_id)
    print(run)

    table_path = run.sample_sheet_path
    table_data = pd.read_csv(table_path)
    return render_template('sample_sheet.html', data=table_data.to_json())


