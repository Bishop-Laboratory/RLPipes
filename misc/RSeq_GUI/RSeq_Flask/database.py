from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, create_session
from sqlalchemy.ext.declarative import declarative_base

engine = None
db_session = scoped_session(lambda: create_session(autocommit=False, autoflush=False, bind=engine))
Base = declarative_base()


def init_engine(uri, **kwargs):
    global engine
    engine = create_engine(uri, **kwargs)
    return engine


def init_db():
    from rseq.models import Run
    from rseq.models import Pipeline
    Base.metadata.create_all(bind=engine)
