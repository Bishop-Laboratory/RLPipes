import click
import ast
import requests
import pandas as pd
import numpy as np
from pysradb.sraweb import SRAweb
import requests
import pysam
import re
import pkg_resources  # part of setuptools
import json
import warnings
import sys
import os
from rseq.run_workflow import make_snakes

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
SRA_COLS=[
  'experiment_title', 'experiment_accession', 
  "organism_taxid ", "run_accession",
  "library_layout", "run_total_bases",
  "run_total_spots"
]
redict = {
  "fastq": '^.+\\.f[ast]q$|^.+f[ast]q~.+f[ast]q$',
  "bam": '^.+\\.bam$',
  "public": '^GSM[0-9]+$|^SRX[0-9]+$'
}
# __file__=os.path.abspath("../../rseq/cli.py")
this_dir = os.path.dirname(__file__)
DATA_PATH = os.path.abspath(os.path.join(this_dir, "src", "data", "available_genomes.tsv.xz"))
GENSIZE_PATH = os.path.abspath(os.path.join(this_dir, "src", "data", "eff_gen_size.tsv.xz"))
SRC_DIR = os.path.abspath(os.path.join(this_dir, "src"))
N_BAM_READS_CHECK=1000

# Set verion
__version__ = pkg_resources.require("rseq")[0].version

# Help text
snakeHelp = """
  Dict of arguments passed to the snakemake python API. Default: "{'use_conda': True}".
  Read the snakemake API reference for the full list of options.
"""
modeHelp = """
  The type of sequencing (e.g., "DRIP"). The available options are currently:
  DRIP, DRIPc, qDRIP, sDRIP, ssDRIP, R-ChIP, RR-ChIP, RDIP, S1-DRIP, DRIVE, RNH-CnR, and MapR
"""

# Get the shared options
# From https://stackoverflow.com/questions/40182157/shared-options-and-flags-between-commands
verify_run_options = [
    click.option("--smargs", "-s", cls=PythonLiteralOption, help=snakeHelp, default="{'use_conda': True}"),
    click.option("--threads", "-t", help="Number of threads to use. Default: 1", default=1),
    click.option(
      "--bwamem2", is_flag=True,
      help="Align with BWA-MEM2 instead of BWA. BWA MEM2 Needs > 70GB RAM avaialble to build index, but shows > 3x speed increase. Default: False.",
      default=False
    ),
    click.option(
      "--macs2", help="Call peaks using macs2 instead of macs2", is_flag=True, default=False
    ),
    click.option("--debug", is_flag=True, help="Run pipeline on subsampled number of reads (for testing).", default=False)
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


def validate_run_dir(ctx, param, value):
  try:
      os.makedirs(value, exist_ok=True)  
  except FileNotFoundError:
    raise click.BadParameter("'" + value + "' could not be created using `os.makedirs(" + value + ", exist_ok=True)` please re-check this path.")
  except FileExistsError:
    raise click.BadParameter("RUN_DIR must be a directory. User supplied '" + value + "'")
  return(os.path.abspath(value))

def validate_run_dir_prepped(ctx, param, value):
  try:
      assert os.path.exists(value) and os.path.exists(os.path.join(value, "config.json")) 
  except AssertionError:
    raise click.BadParameter("Configuration file '" + os.path.join(value, "config.json") + "' is not found. Have you run 'RSeqCLI build' yet?")
  return(os.path.abspath(value))


def bam_info(bamfile, n_bam_reads_check=1000):
    """Tests whether bam file is paired end and checks read length. Requires pysam."""
    save = pysam.set_verbosity(0)
    samfile=pysam.AlignmentFile(bamfile, "rb")
    pysam.set_verbosity(save)
    numPair = sum([x.is_paired for x in samfile.head(n=n_bam_reads_check)])
    read_len=sum([x.infer_read_length() for x in samfile.head(n=n_bam_reads_check)])//n_bam_reads_check
    return({"paired_end": numPair > n_bam_reads_check/2,
            "read_len": read_len})
    

def validate_samples(ctx, param, value):
  """Validate and wrangle sampels input"""
  params=ctx.params
  samps = pd.read_csv(value)
  
  # First, check for matching pattern
  exp = samps.experiment[0]
  samptype = [key for key, val in redict.items() if re.match(val, exp)][0]
  samps['file_type'] = samptype
  
  # Wrangle controls if provided
  if "control" in samps.columns:
      controls=True
      samps = pd.concat([
        samps,
        samps.assign(experiment=samps.control).assign(control=pd.NA).dropna(subset=['experiment'])
      ])
      samps=samps.assign(control=samps.control.apply(lambda x: pd.NA if pd.isna(x) else x ))
  else:
    controls=False
    samps['control']=""
  
  if samptype == "public":
    
    # Init SRAdb
    db = SRAweb(os.environ.get("NCBI_API_KEY", None))
    
    # Query the SRAdb and wrangle with original data
    newSamps = pd.concat(samps.experiment.apply(lambda x: db.sra_metadata(x)).values.tolist())[SRA_COLS]
    
    # Get the read length
    newSamps=newSamps.astype({'run_total_bases': "int64", 'run_total_spots': 'int64'})
    newSamps['read_length'] = newSamps.run_total_bases // newSamps.run_total_spots
    
    # Get the latest genomes. 
    # From https://stackoverflow.com/questions/15705630/get-the-rows-which-have-the-max-value-in-groups-using-groupby
    available_genome = pd.read_table(DATA_PATH)
    latest_genomes = available_genome[available_genome.groupby(axis=0, by=['taxId'])['year'].transform(max) == available_genome['year']]
    latest_genomes = latest_genomes.rename(columns={'taxId': 'organism_taxid '})
    newSamps['organism_taxid '] = newSamps['organism_taxid '].astype(np.int64)
    
    # Necessary to avoid taxid conflict between sacCer2/3
    newSamps.loc[newSamps['organism_taxid '] == 4932, 'organism_taxid '] = 559292
    newSamps=newSamps.set_index('organism_taxid ')
    latest_genomes=latest_genomes.set_index('organism_taxid ')
    newSamps=newSamps.join(latest_genomes, how='left')
    newSamps=newSamps[
      ['experiment_title', 'experiment_accession', 'run_accession', 
       'library_layout', 'UCSC_orgID', 'read_length']
    ].rename(columns={
      "experiment_title": "name", "library_layout": "paired_end",
      "experiment_accession": "experiment", 
      "run_accession": "run", "UCSC_orgID": "genome"
    }).set_index('experiment')
    
    # Set paired end
    newSamps['paired_end'] = newSamps['paired_end'] == "PAIRED"
    
    if "genome" in samps.columns:
      newSamps=newSamps.drop('genome', axis=1)
    
    # Finally, join the dataframes 
    samps=samps.set_index('experiment').join(newSamps).reset_index(level=0)
  elif samptype == "bam":
    # Check which are paired-end
    samps['paired_end'] = [bam_info(bam, N_BAM_READS_CHECK)['paired_end'] for bam in samps['experiment']]
    samps['read_length'] = [bam_info(bam, N_BAM_READS_CHECK)['read_len'] for bam in samps['experiment']]
    samps['run'] = [os.path.abspath(bam) for bam in samps['experiment']]
    samps['name'] = samps['experiment'] = [os.path.splitext(os.path.basename(exp))[0] for exp in samps['experiment']]
    if controls:
      samps['control'] = [
        os.path.splitext(os.path.basename(exp))[0] if not pd.isna(exp) else exp for exp in samps['control']
      ]
    else:
      samps['control']=""
  
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
  "run_dir", type=click.Path(), callback=validate_run_dir
)
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
@click.option("--name", "-n", help="Sample names for use in output report. By default, inferred from inputs.")
@click.pass_context
def build(ctx, samples, mode, genome, run_dir, name):
    """
    Configure an RSeq workflow.
    
    RUN_DIR: Directory for RSeq Execution. Will be created if it does not exist.
    
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
      try: 
        assert genome is not None
      except AssertionError:
        raise click.BadParameter(message="Genome cannot be missing when running local files.")
      samples['genome'] = genome
    if not 'name' in samples.columns:
      samples['name'] = name
      
    # Get the effective genome sizes
    eff_gen_size=pd.read_table(GENSIZE_PATH)
    sizes=eff_gen_size['read_length'].unique()
    samples['read_length'] = [min(sizes, key=lambda x:abs(x-sizeCheck)) for sizeCheck in samples['read_length']]
    samples = pd.merge(samples, eff_gen_size.rename(columns={"UCSC_orgID": "genome"}), how='inner')
  
      
    # Compile to json for snakemake
    outjson = os.path.join(run_dir, "config.json")
    outdict = samples.fillna('').reset_index().drop('index', axis=1).to_dict("list")
    with open(outjson, 'w') as f:
      json.dump(outdict, f, ensure_ascii=False)
    
    print("\nSuccess! RSeq has been initialized at the specified directory: " + run_dir)
    print("\nRun 'RSeqCLI check " + run_dir + "' to verify the configuration.\n")
    

@cli.command("check")
@click.argument(
  "run_dir", type=click.Path(), callback=validate_run_dir_prepped
)
@add_options(verify_run_options)
def check(run_dir, threads, debug, bwamem2, macs2, **kwargs):
    """
    Validate an RSeq workflow.
    
    RUN_DIR: Directory configured with `RSeqCLI build` and ready for checking and execution.
    """
    smargs = kwargs['smargs']
    dagfile = make_snakes(
      snake_args=smargs, run_dir=run_dir, 
      src_dir=SRC_DIR, threads=threads, 
      bwamem2=bwamem2, macs2=macs2,
      debug=debug, verify=True
    )
    print("\nSuccess! The DAG has been generated successfully. You can view it here: " + dagfile)
    print("\nRun 'RSeqCLI run " + run_dir + "' to execute the workflow.\n")
    

@cli.command("run")
@click.argument(
  "run_dir", type=click.Path(), callback=validate_run_dir_prepped
)
@add_options(verify_run_options)
def run(run_dir, threads, debug, bwamem2, macs2, **kwargs):
    """
    Execute an RSeq workflow.
    
    RUN_DIR: Directory configured with `RSeqCLI build` and ready for checking and execution.
    """
    smargs = kwargs['smargs']
    exitcode = make_snakes(
      snake_args=smargs, run_dir=run_dir, 
      src_dir=SRC_DIR, threads=threads, 
      bwamem2=bwamem2, macs2=macs2,
      debug=debug, verify=False
    )
    print(exitcode)
    print("Success! RSeqCLI will now close.")

