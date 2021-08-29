# pytest script for checking RSeqCLI build against user inputs
from click.testing import CliRunner
from rseq.cli import build, check, run
import os
import shutil

RSEQ_OUT_BAM1='tests/rseq_out_bams/'
BAMSAMPS='tests/test_data/bam_test_samples_1.csv'
RSEQ_OUT_BAM2='tests/rseq_out_bams_macs2/'

def test_run_1():
  if os.path.exists(RSEQ_OUT_BAM1):
    shutil.rmtree(RSEQ_OUT_BAM1)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_BAM1, BAMSAMPS, "-g", "hg38"])
  checkres = runner.invoke(check, [RSEQ_OUT_BAM1])
  # runres = runner.invoke(run, [RSEQ_OUT_BAM1, "--debug"])
  assert checkres.exit_code == 0
