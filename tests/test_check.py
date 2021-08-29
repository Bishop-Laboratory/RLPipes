# pytest script for checking RSeqCLI build against user inputs
from click.testing import CliRunner
from rseq.cli import build, check
import os
import shutil

RSEQ_OUT_PUBLIC='rseq_out_public/'
PUBSAMPS='test_data/samples.csv'
RSEQ_OUT_BAM='rseq_out_bams/'
BAMSAMPS='test_data/bam_test_samples_1.csv'

def test_check():
  if os.path.exists(RSEQ_OUT_PUBLIC):
    shutil.rmtree(RSEQ_OUT_PUBLIC)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_PUBLIC, PUBSAMPS])
  checkres = runner.invoke(check, [RSEQ_OUT_PUBLIC])
  assert checkres.exit_code == 0


def test_check_bam():
  if os.path.exists(RSEQ_OUT_BAM):
    shutil.rmtree(RSEQ_OUT_BAM)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_BAM, BAMSAMPS, "-g", "hg38"])
  checkres = runner.invoke(check, [RSEQ_OUT_BAM])
  assert checkres.exit_code == 0

