# pytest script for checking RSeqCLI build against user inputs
from click.testing import CliRunner
from rlpipes.cli import build
import os
import shutil

RSEQ_OUT_PUBLIC='tests/rseq_out_public/'
PUBSAMPS='tests/test_data/samples.csv'
RSEQ_OUT_BAM='tests/rseq_out_bams/'
BAMSAMPS='tests/test_data/bam_test_samples_1.csv'
RSEQ_OUT_FQ='tests/rseq_out_fqs/'
FQSAMPS='tests/test_data/fq_test_samples_3.csv'
FQSAMPSGZ='tests/test_data/fq_test_samples_1.csv'
FQSAMPSGZ2='tests/test_data/fq_test_samples_4.csv'


def test_build():
  if os.path.exists(RSEQ_OUT_PUBLIC):
    shutil.rmtree(RSEQ_OUT_PUBLIC)
  runner = CliRunner()
  result = runner.invoke(build, [RSEQ_OUT_PUBLIC, PUBSAMPS])
  assert result.exit_code == 0

def test_build_fail_mode():
  if os.path.exists(RSEQ_OUT_PUBLIC):
    shutil.rmtree(RSEQ_OUT_PUBLIC)
  runner = CliRunner()
  result = runner.invoke(build, [RSEQ_OUT_PUBLIC, PUBSAMPS, "-m", "bisDRIP"])
  assert result.exit_code == 2

def test_build_fail_mkdir():
  runner = CliRunner()
  result = runner.invoke(build, ['/asd/as/d/sad/asd/rseqOutdir/', PUBSAMPS, "-m", "bisDRIP"])
  assert result.exit_code == 2
  
def test_build_bam():
  if os.path.exists(RSEQ_OUT_BAM):
    shutil.rmtree(RSEQ_OUT_BAM)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_BAM, BAMSAMPS, "-g", "hg38"])
  assert buildres.exit_code == 0


def test_build_fq():
  if os.path.exists(RSEQ_OUT_FQ):
    shutil.rmtree(RSEQ_OUT_FQ)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_FQ, FQSAMPS, "-g", "hg38"])
  assert buildres.exit_code == 0


def test_build_fq_gz():
  if os.path.exists(RSEQ_OUT_FQ):
    shutil.rmtree(RSEQ_OUT_FQ)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_FQ, FQSAMPSGZ, "-g", "mm10"])
  assert buildres.exit_code == 0


def test_build_fq_gz2():
  if os.path.exists(RSEQ_OUT_FQ):
    shutil.rmtree(RSEQ_OUT_FQ)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_FQ, FQSAMPSGZ2, "-g", "hg38"])
  assert buildres.exit_code == 0
