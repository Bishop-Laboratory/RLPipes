import os
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for, current_app as app
)
from werkzeug.utils import secure_filename
from rseq.db import get_db

bp = Blueprint('runs', __name__, url_prefix='/runs')
ALLOWED_EXTENSIONS = {'txt', 'csv', 'tsv', 'tab', 'xls', 'xlsx'}


def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@bp.route('/run_builder', methods=('GET', 'POST'))
def run_builder():
    if request.method == 'POST':
        pipeline_name = request.form['pipeline_name']
        run_name = request.form['run_name']
        file = request.files['file']

        # if 'file' not in request.files:
        #     flash('No file part')
        #     return redirect(request.url)
        # file = request.files['file']
        # if file.filename == '':
        #     flash('No selected file')
        #     return redirect(request.url)

        db = get_db()
        error = None

        pipeline_id = db.execute(
            'SELECT id FROM pipelines WHERE pipeline_name = ?', (pipeline_name,)
        ).fetchone()

        if not pipeline_name:
            error = 'runs is required.'
        elif not run_name:
            error = 'run_name is required.'
        elif pipeline_id is None:
            error = 'Pipeline {} does not exist.'.format(pipeline_name)
        elif not file:
            error = "file is required"
        elif not allowed_file(file.filename):
            error = '{} is not an allowed filename'.format(file.filename)

        if file and allowed_file(file.filename) and error is None:
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            db.execute(
                'INSERT INTO runs (run_name, pipeline_id) VALUES (?, ?)',
                (run_name, pipeline_id[0])
            )
            db.commit()
            message = "Successfully created new run" + run_name
            flash(message, category="success")

            return redirect(url_for('runs.run_builder',
                                    run_name=run_name,
                                    pipeline_name=pipeline_name,
                                    filename=filename))

        flash(error, category="danger")

    return render_template('runs/run_builder.html')


