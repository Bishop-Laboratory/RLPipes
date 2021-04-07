import datetime
from sqlalchemy import Column, Integer, String, DateTime
from sqlalchemy.sql import func
from rseq.database import Base


class Run(Base):
    __tablename__ = 'Runs'
    id = Column(Integer, primary_key=True)
    run_name = Column(String(1000), unique=False, nullable=False)
    author = Column(String(1000), nullable=True)
    email_address = Column(String(1000))
    outdir = Column(String(1000))
    mode = Column(String(1000))
    run_path = Column(String(1000))
    # config_file = Column(String(1000))
    # log_file = Column(String(1000))
    time_created = Column(DateTime(timezone=True), server_default=func.now())
    status = Column(String(1000))
    steps = Column(String(10000))
    time_updated = Column(DateTime(timezone=True), onupdate=func.now())

    def __repr__(self):
        if self.email_address is not None:
            title_str = '<Run: %s; notify: %s>' % (
                self.run_name, self.email_address)
        else:
            title_str = '<Run: %s>' % (
                self.run_name)
        return title_str


class Pipeline(Base):
    __tablename__ = 'Pipelines'
    id = Column(Integer, primary_key=True)
    pipeline_name = Column(String(1000), unique=True, nullable=False)
    author = Column(String(1000))
    email_address = Column(String(1000))
    snake_file = Column(String(1000))
    steps = Column(String(10000))
    time_created = Column(DateTime(timezone=True), server_default=func.now())
    time_updated = Column(DateTime(timezone=True), onupdate=func.now())

    def __repr__(self):
        return '<Pipeline %r>' % self.pipeline_name


