# pytest script for checking RSeqCLI build against user inputs
from click.testing import CliRunner
from rseq.cli import build, check, run
import os
import shutil

RSEQ_OUT_BAM1='tests/rseq_out_bams/'
BAMSAMPS='tests/test_data/bam_test_samples_1.csv'
RSEQ_OUT_FQ2='tests/rseq_out_fqs2/'
FQSAMPS2='tests/test_data/fq_test_samples_2.csv'

def test_run_1():
  if os.path.exists(RSEQ_OUT_BAM1):
    shutil.rmtree(RSEQ_OUT_BAM1)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_BAM1, BAMSAMPS])
  checkres = runner.invoke(check, [RSEQ_OUT_BAM1])
  runres = runner.invoke(run, [RSEQ_OUT_BAM1, "--debug"])
  assert runres.exit_code == 0


# def test_run_2():
#   if os.path.exists(RSEQ_OUT_FQ2):
#     shutil.rmtree(RSEQ_OUT_FQ2)
#   runner = CliRunner()
#   buildres = runner.invoke(build, [RSEQ_OUT_FQ2, FQSAMPS2])
#   checkres = runner.invoke(check, [RSEQ_OUT_FQ2])
#   runres = runner.invoke(run, [RSEQ_OUT_FQ2, "--debug", "--bwamem2"])
#   assert runres.exit_code == 0

