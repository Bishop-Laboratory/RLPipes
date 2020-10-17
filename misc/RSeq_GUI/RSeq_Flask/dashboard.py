import os
from flask import (
    Blueprint, flash, g, redirect, render_template, jsonify, request, session, url_for, current_app as app
)
import pandas as pd
from rseq.database import db_session

bp = Blueprint('dashboard', __name__)


@bp.route('/')
def dashboard():
    runs = pd.read_sql_table('Runs', con=app.config['SQLALCHEMY_DATABASE_URI'])
    run_exists = True
    if len(runs['run_name']) == 0:
        run_exists = False
    return render_template('dashboard.html', run_exists=run_exists)


@bp.route('/_get_table')
def get_table():
    runs = pd.read_sql_table('Runs', con=app.config['SQLALCHEMY_DATABASE_URI'])
    runs['edit'] = runs['id'].apply(
        lambda x: '<a class="btn btn-primary" role="button" href="/runs/{0}/edit">Edit</a>'.format(x))
    runs['delete'] = runs[['id', 'run_name']].apply(
        lambda x: '''
            <!-- Button to Open the Modal -->
            <button type="button" class="btn btn-danger" data-toggle="modal" data-target="#delete_run_{0}">
                Delete
            </button>
            <!-- The Modal -->
            <div class="modal fade" id="delete_run_{0}">
            <div class="modal-dialog modal-dialog-centered">
              <div class="modal-content">
            
                <!-- Modal Header -->
                <div class="modal-header">
                  <h4 class="modal-title">Confirm Delete</h4>
                  <button type="button" class="close" data-dismiss="modal">&times;</button>
                </div>
            
                <!-- Modal body -->
                <div class="modal-body">
                  Confirm deletion of {1}. <strong>This action cannot be undone</strong>. 
                  If run is currently processing, it will be
                  terminated immediately. Any changes made to the sample sheet will not be saved.
                  <br>
                  Are you sure?
                </div>
            
                <!-- Modal footer -->
                <div class="modal-footer">
                  <a class="btn btn-danger" role="button" href="/runs/{0}/delete">Delete</a>
                  <button type="button" class="btn btn-secondary" data-dismiss="modal">Close</button>
                </div>
            
              </div>
            </div>
            </div>
            '''.replace('\n', '').format(x[0], x[1]), axis=1)
    runs['sample_sheet'] = runs['run_path'] + '/samples.csv'
    runs['sample_sheet'] = runs[['id', 'sample_sheet']].apply(
        lambda x: '<a class="btn btn-primary" role="button" href="/runs/{0}/sample_sheet">{1}</a>'.format(
            x[0], 'samples'), axis=1)
    action = []
    for ind, row in runs.iterrows():
        if row['status'] == '<strong class="text-muted">uninitialized</strong>':
            action.append(
                '<a class="btn btn-primary" role="button" href="/runs/{0}/initialize">Initialize</a>'.format(row['id']))
        elif row['status'] == '<strong class="text-info">Initialized. Ready for pre-flight ☑</strong>':
            action.append(
                '<a class="btn btn-info" role="button" href="/runs/{0}/preflight">Preflight</a>'.format(row['id']))
        elif row['status'] in ('<strong class="text-info">In flight ✈️</strong>',
                               '<strong class="text-info">Safely aborting</strong>',
                               '<strong class="text-info">Ready for takeoff 🛫</strong>',
                               '<strong class="text-info">Complete! ✔️</strong>'):
            action.append("""
                <a class="btn btn-info" role="button" href="/runs/{0}/monitor">Dashboard</a>
            """.replace("\n", "").format(row['id']))
        else:
            action.append(
                '<a class="btn btn-primary" role="button" href="/runs/{0}/something">something</a>'.format(row['id']))
    runs['action'] = action
    runs.head()
    cols_selected = ['time_created', 'run_name', 'sample_sheet', 'outdir', 'status', 'edit', 'delete', 'action']
    runs = runs[cols_selected]
    return jsonify(run_table=runs.to_html(classes='table table-striped" id = "run_table',
                                          index=False, border=0, escape=False))

