# pytest script for checking RSeqCLI build against user inputs
from click.testing import CliRunner
from rseq.cli import build, check, run
import os
import shutil

RSEQ_OUT_BAM1='rseq_out_bams/'
BAMSAMPS='test_data/bam_test_samples_1.csv'
RSEQ_OUT_BAM2='rseq_out_bams_macs2/'

def test_run_1():
  if os.path.exists(RSEQ_OUT_BAM1):
    shutil.rmtree(RSEQ_OUT_BAM1)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_BAM1, BAMSAMPS, "-g", "hg38"])
  checkres = runner.invoke(check, [RSEQ_OUT_BAM1])
  runres = runner.invoke(run, [RSEQ_OUT_BAM1, "--debug"])
  assert runres.exit_code == 0


def test_run_2():
  if os.path.exists(RSEQ_OUT_BAM2):
    shutil.rmtree(RSEQ_OUT_BAM2)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_BAM2, BAMSAMPS, "-g", "hg38"])
  checkres = runner.invoke(check, [RSEQ_OUT_BAM2])
  runres = runner.invoke(run, [RSEQ_OUT_BAM2, "--debug", "--macs2"])
  assert runres.exit_code == 0
  
