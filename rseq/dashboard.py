from flask import (
    Blueprint, flash, g, redirect, render_template, request, url_for
)

from rseq.db import get_db

bp = Blueprint('dashboard', __name__)


@bp.route('/')
def index():
    db = get_db()
    pipelines = db.execute(
        'SELECT pipeline_name'
        ' FROM pipelines'
    ).fetchall()

    runs = db.execute(
        'SELECT run_name'
        ' FROM runs'
    ).fetchall()

    # reports = db.execute(
    #     'SELECT report_name'
    #     ' FROM reports'
    # ).fetchall()

    return render_template('dashboard/index.html', pipelines=pipelines, runs=runs)
