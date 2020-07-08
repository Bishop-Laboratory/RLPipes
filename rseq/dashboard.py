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
    runs['edit'] = runs['id'].apply(
        lambda x: '''\
        <a class="btn btn-primary" role="button" href="/runs/{0}/edit">Edit</a>\
        '''.format(x))
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
    runs['sample_sheet'] = runs[['id', 'sample_sheet']].apply(
        lambda x: '<a class="btn btn-primary" role="button" href="/runs/{0}/sample_sheet">{1}</a>'.format(
            x[0], x[1]), axis=1)
    cols_selected = ['time_created', 'run_name', 'sample_sheet', 'outdir', 'status', 'edit', 'delete']
    runs = runs[cols_selected]
    print(runs)
    return jsonify(run_table=runs.to_html(classes='table table-striped" id = "run_table',
                                          index=False, border=0, escape=False))






