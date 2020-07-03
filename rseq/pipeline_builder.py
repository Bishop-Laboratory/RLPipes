import functools
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for
)
from rseq.db import get_db

bp = Blueprint('pipelines', __name__, url_prefix='/pipelines')


@bp.route('/pipeline_builder', methods=('GET', 'POST'))
def pipeline_builder():
    if request.method == 'POST':
        pipeline_name = request.form['pipeline_name']
        pipeline_author = request.form['pipeline_author']
        db = get_db()
        error = None

        pipeline_id = db.execute(
            'SELECT id FROM pipelines WHERE pipeline_name = ?', (pipeline_name,)
        ).fetchone()

        if pipeline_id is not None:
            error = 'Pipeline {} already exists.'.format(pipeline_name)
        elif not pipeline_author:
            error = 'pipeline_author is required.'

        if error is None:
            db.execute(
                'INSERT INTO pipelines (pipeline_name, pipeline_author) VALUES (?, ?)',
                (pipeline_name, pipeline_author)
            )
            db.commit()
            message = "Pipeline " + pipeline_name + " successfully created!"
            flash(message, category="success")
            return redirect(url_for('pipelines.pipeline_builder'))

        flash(error, category="danger")

    return render_template('pipelines/pipeline_builder.html')


