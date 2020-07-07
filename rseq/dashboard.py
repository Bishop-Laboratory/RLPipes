import os
from flask import (
    Blueprint, flash, g, redirect, render_template, jsonify, request, session, url_for, current_app as app
)
import pandas as pd
from rseq.database import db_session

bp = Blueprint('dashboard', __name__)


@bp.route('/')
def dashboard():
    return render_template('dashboard.html')


@bp.route('/_get_table')
def get_table():
    runs = pd.read_sql_table('Runs', con=app.config['SQLALCHEMY_DATABASE_URI'])
    runs['actions'] = runs['id'].apply(
        lambda x: '''\
        <a class="btn btn-primary" role="button" href="/runs/{0}/edit">Edit</a>\
        <a class="btn btn-danger" role="button" href="/runs/{0}/delete">Delete</a>\
        '''.format(x))
    runs['sample_sheet'] = runs[['id', 'sample_sheet']].apply(
        lambda x: '<a class="btn btn-primary" role="button" href="/runs/{0}/sample_sheet">{1}</a>'.format(
            x[0], x[1]), axis=1)
    print(runs)
    return jsonify(run_table=runs.to_html(classes='table table-striped" id = "run_table',
                                          index=False, border=0, escape=False))






