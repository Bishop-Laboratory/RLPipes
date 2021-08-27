import click
import ast
import requests
import pandas as pd
import numpy as np
from pysradb.sraweb import SRAweb
import requests
from .run_workflow import make_snakes
import os
import re
import pkg_resources  # part of setuptools
import json

# Init SRAdb
db = SRAweb(os.environ.get("NCBI_API_KEY", None))

# Allows passing strings to CLI and eval as python objects
# From https://stackoverflow.com/questions/47631914/how-to-pass-several-list-of-arguments-to-click-option
# Had to add the "str(value)" part or default would throw ValueErrors.
class PythonLiteralOption(click.Option):
  def type_cast_value(self, ctx, value):
    try:
      return ast.literal_eval(str(value))
    except ValueError:
      raise click.BadParameter(value)
    
    
# Constants
AVAILABLE_MODES=["DRIP", "DRIPc", "qDRIP", "sDRIP", "ssDRIP", "R-ChIP", 
                 "RR-ChIP", "RDIP", "S1-DRIP", "DRIVE", "RNH-CnR", "MapR"]
SRA_URL = "https://www.ncbi.nlm.nih.gov/sra/"
redict = {
  "fastq": '^.+\.f[ast]q$|^.+f[ast]q\w.+f[ast]q$',
  "bam": '^.+\.bam$',
  "public": '^GSM[0-9]+$|^SRX[0-9]+$'
}

this_dir, this_filename = os.path.split(__file__)
DATA_PATH = os.path.join(os.path.dirname(this_dir), "src", "data", "available_genomes.tsv.xz")
__version__ = pkg_resources.require("rseq")[0].version

# Help text
snakeHelp = """
  A dctionary of arguments to pass to the snakemake python API. Default: "{'use_conda': True}".
  Please visit the snakemake API reference for the full list of available options.
"""
modeHelp = """
  Specify the type of sequencing (e.g., "DRIP"). The available options are currently:
  DRIP, DRIPc, qDRIP, sDRIP, ssDRIP, R-ChIP, RR-ChIP, RDIP, S1-DRIP, DRIVE, RNH-CnR, and MapR
"""

# Get the shared options
# From https://stackoverflow.com/questions/40182157/shared-options-and-flags-between-commands
verify_run_options = [
    click.option("--smargs", "-s", cls=PythonLiteralOption, help=snakeHelp, default="{'use_conda': True}")
]


# Function for addint these options to the click command
def add_options(options):
  def _add_options(func):
      for option in reversed(options):
          func = option(func)
      return func
  return _add_options
  

def validate_genome(ctx, param, value):
  """Validate genome input"""
  if value is not None:
    available_genomes = pd.read_table(DATA_PATH)
    try:
      assert value in available_genomes.UCSC_orgID.to_list()
    except AssertionError:
      raise click.BadParameter("'" + value + "' is not a valid UCSC genome ID (e.g., 'hg38' is valid).")
  return(value)


def validate_mode(ctx, param, value):
  """Validate mode input"""
  if value is not None:
    try:
      assert value in AVAILABLE_MODES
    except AssertionError:
      raise click.BadParameter("'" + value + "' is not a valid mode. (RSeqCLI build --help for more info)")
  return(value)


def validate_samples(ctx, param, value):
  """Validate mode input"""
  params=ctx.params
  samps = pd.read_csv(value)
  
  # First, check for matching pattern
  exp = samps.experiment[0]
  samptype = [key for key, val in redict.items() if re.match(val, exp)][0]
  
  # Wrangle controls if provided
  if "control" in samps.columns:
    # This line adds all controls an independent "experiment" samples
    # This eases the process of querying sample info and running the pipeline later.
    samps = pd.concat([
      samps,
      samps.assign(experiment=samps.control).assign(control=pd.NA).dropna(subset=['experiment'])
    ])
    samps=samps.assign(control=samps.control.apply(lambda x: pd.NA if pd.isna(x) else x ))
  
  
  if samptype == "public":
    # Query the SRAdb and wrangle with original data
    SRA_COLS=['experiment_title', 'experiment_accession', "organism_taxid ", "run_accession"]
    newSamps = pd.concat(samps.experiment.apply(lambda x: db.sra_metadata(x)).values.tolist())[SRA_COLS]
    
    # Get the latest genomes. From https://stackoverflow.com/questions/15705630/get-the-rows-which-have-the-max-value-in-groups-using-groupby
    available_genome = pd.read_table(DATA_PATH)
    latest_genomes = available_genome[available_genome.groupby(axis=0, by=['taxId'])['year'].transform(max) == available_genome['year']]
    latest_genomes = latest_genomes.rename(columns={'taxId': 'organism_taxid '})
    newSamps['organism_taxid '] = newSamps['organism_taxid '].astype(np.int64)
    newSamps=newSamps.set_index('organism_taxid ')
    latest_genomes=latest_genomes.set_index('organism_taxid ')
    newSamps=newSamps.join(latest_genomes, how='left')
    newSamps=newSamps[['experiment_title', 'experiment_accession', 'run_accession', 'UCSC_orgID']].rename(columns={
      "experiment_title": "name", "experiment_accession": "experiment", "run_accession": "run", "UCSC_orgID": "genome"
    }).set_index('experiment')
    
    if "genome" in samps.columns:
      newSamps=newSamps.drop('genome', axis=1)
    
    # Finally, join the dataframes 
    samps=samps.set_index('experiment').join(newSamps).reset_index(level=0)
    
  return(samps)


@click.group()
@click.version_option(__version__)
@click.pass_context
def cli(ctx, **kwargs):
  """
  RSeq: An R-loop mapping pipeline with built-in QC.
  """
  pass


@cli.command("build")
@click.argument(
  "samples", type=click.File('rb'), callback=validate_samples
)
@click.option(
  "--mode", "-m", 
  callback=validate_mode,
  help=modeHelp
)
@click.option(
  "--genome", "-g", callback=validate_genome,
  help="UCSC genome for samples (e.g., 'hg38'). Not required if providing public data accessions."
)
@click.option("--outdir", "-t", help="Output folder. Default: 'rseq_out/'", default="rseq_out/")
@click.option("--threads", "-o", help="Number of threads to use. Default: 1", default=1)
@click.option("--name", "-n", help="Sample names for use in output report. By default, inferred from inputs.")
@click.option("--bwamem2", help="Align with BWA-MEM2 instead of BWA. BWA MEM2 Needs > 70GB RAM avaialble to build index, but shows > 3x speed increase. Default: False.", default=False)
@click.pass_context
def build(ctx, samples, mode, genome, outdir, threads, name, bwamem2):
    """
    Configure an RSeq workflow.
    
    SAMPLES: A CSV file with at least one column "experiment" that provides the
    path to either local fastq files, bam files, or public sample accessions (SRX or GSM). 
    Input controls should be in the "control" column. \n
    If providing paired-end fastq files, enter: "exp_1.fastq~exp_2.fastq".\n
    Columns may also include "genome" and "mode" columns. These will 
    override the -g, -m, and -n  options.\n
    "genome" (-g/--genome) is not required if providing public data accessions.\n
     \n
    Example #1: "RSeqCLI build -m DRIP samples.csv"; where "samples.csv" is:\n
    
    experiment\n
    SRX113812\n
    SRX113813\n 
     \n
    
    Example #2: "RSeqCLI build samples.csv"; where "samples.csv" is:\n
    
    experiment,control,genome,mode\n
    qDRIP_siGL3_1.fq~qDRIP_siGL3_2.fq,,hg38,qDRIP\n
    DRIPc_3T3.fq,Input_3T3.fq,mm10,DRIPc
    """
    # Add in potentially, missing columns
    if not "mode" in samples.columns:
      samples['mode'] = mode
    if not 'genome' in samples.columns:
      samples['genome'] = genome
    if not 'name' in samples.columns:
      samples['name'] = name
    
    # Create output folder if it doesn't already exist
    if not os.path.exists(outdir):
      os.mkdir(outdir)
    
    # Compile for snakemake
    outjson=os.path.join(outdir, "config.json")
    samples.fillna('').reset_index().drop('index', axis=1).to_json(path_or_buf=outjson)
    
    print("Success! The config file is located here: " + outjson)
    

@cli.command("check")
@click.argument("config")
@add_options(verify_run_options)
def check(config, **kwargs):
    """
    Validate an RSeq workflow.
    
    CONFIG: required. Filepath to config.json generated by RSeqCLI build.
    """
    click.echo(f"{config}")
    smargs=kwargs['smargs']
    click.echo("smargs, type: {}  value: {}".format(
        type(smargs), smargs))
    make_snakes(config, smargs, verify=True)


@cli.command("run")
@click.argument("config")
@add_options(verify_run_options)
def run(config, **kwargs):
    """
    Execute an RSeq workflow.
    
    CONFIG: required. Filepath to config.json generated by RSeqCLI build.
    """
    click.echo(f"{config}")
    smargs=kwargs['smargs']
    click.echo("smargs, type: {}  value: {}".format(
        type(smargs), smargs))

