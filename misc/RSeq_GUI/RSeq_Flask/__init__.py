import os
from flask import Flask


def create_app(test_config=None):
    # create and configure the app
    app = Flask(__name__, instance_relative_config=True)
    app.config['UPLOAD_FOLDER'] = os.path.join(app.instance_path, 'uploads')
    app.config['SQLALCHEMY_DATABASE_URI'] = "sqlite:///" + os.path.join(app.instance_path, "rseq.db")

    # Import blueprints
    from . import run_builder
    app.register_blueprint(run_builder.bp)

    from . import pipeline_builder
    app.register_blueprint(pipeline_builder.bp)

    from . import dashboard
    app.register_blueprint(dashboard.bp)

    if test_config is None:
        # load the instance config, if it exists, when not testing
        app.config.from_pyfile('config.py', silent=True)
    else:
        # load the test config if passed in
        app.config.from_mapping(test_config)

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    # ensure the uploads folder exists
    try:
        os.makedirs(app.config['UPLOAD_FOLDER'])
    except OSError:
        pass

    # Set up the database
    from rseq.database import init_db, init_engine
    init_engine(app.config['SQLALCHEMY_DATABASE_URI'])
    init_db()
    app.secret_key = 'super secret key'
    app.config['SESSION_TYPE'] = 'filesystem'
    return app
