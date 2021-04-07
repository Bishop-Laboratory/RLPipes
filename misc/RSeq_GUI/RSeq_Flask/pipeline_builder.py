import functools
from flask import (
    Blueprint, flash, g, redirect, render_template, request, session, url_for
)
from rseq.models import Run, Pipeline
from rseq.database import db_session
bp = Blueprint('pipelines', __name__, url_prefix='/pipelines')


@bp.route('/pipeline_builder', methods=('GET', 'POST'))
def pipeline_builder():
    if request.method == 'POST':
        pipeline_name = request.form['pipeline_name']
        print(request.form)
        error = None

        pipeline_id = db_session.query(Pipeline).filter_by(pipeline_name=pipeline_name).first()

        if pipeline_id is not None:
            error = 'Pipeline {} already exists.'.format(pipeline_name)
        # elif not pipeline_author:
        #     error = 'pipeline_author is required.'

        if error is None:
            pipe_now = Pipeline(pipeline_name=pipeline_name)
            db_session.add(pipe_now)
            db_session.commit()
            message = "Pipeline " + pipeline_name + " successfully created!"
            flash(message, category="success")
            return redirect(url_for('pipelines.pipeline_builder'))

        flash(error, category="danger")

    return render_template('pipeline_builder.html')


