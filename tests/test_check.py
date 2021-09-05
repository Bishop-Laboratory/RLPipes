# pytest script for checking RSeqCLI build against user inputs
from click.testing import CliRunner
from rlpipes.cli import build, check
import os
import shutil

RSEQ_OUT_PUBLIC='tests/rseq_out_public/'
RSEQ_OUT_PUBLICRNA='tests/rseq_out_public_rna/'
PUBSAMPS='tests/test_data/samples.csv'
PUBRNASAMPS="tests/test_data/samples_rna.csv"
RSEQ_OUT_BAM='tests/rseq_out_bams/'
BAMSAMPS='tests/test_data/bam_test_samples_1.csv'
RSEQ_OUT_FQ='tests/rseq_out_fqs/'
FQSAMPS='tests/test_data/fq_test_samples_1.csv'

def test_check():
  if os.path.exists(RSEQ_OUT_PUBLIC):
    shutil.rmtree(RSEQ_OUT_PUBLIC)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_PUBLIC, PUBSAMPS])
  checkres = runner.invoke(check, [RSEQ_OUT_PUBLIC])
  assert checkres.exit_code == 0
  

def test_check_tsv():
  if os.path.exists(RSEQ_OUT_PUBLIC):
    shutil.rmtree(RSEQ_OUT_PUBLIC)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_PUBLIC, PUBSAMPS])
  checkres = runner.invoke(check, [RSEQ_OUT_PUBLIC, "--tsv"])
  assert checkres.exit_code == 0


def test_check_failgroupby():
  if os.path.exists(RSEQ_OUT_PUBLIC):
    shutil.rmtree(RSEQ_OUT_PUBLIC)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_PUBLIC, PUBSAMPS])
  checkres = runner.invoke(check, [RSEQ_OUT_PUBLIC, "-G", "mode asd"])
  assert checkres.exit_code == 2
  

def test_check_passgroupby():
  if os.path.exists(RSEQ_OUT_PUBLIC):
    shutil.rmtree(RSEQ_OUT_PUBLIC)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_PUBLIC, PUBSAMPS])
  checkres = runner.invoke(check, [RSEQ_OUT_PUBLIC, "-G", "mode"])
  assert checkres.exit_code == 0

def test_check_rna():
  if os.path.exists(RSEQ_OUT_PUBLICRNA):
    shutil.rmtree(RSEQ_OUT_PUBLICRNA)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_PUBLICRNA, PUBRNASAMPS])
  checkres = runner.invoke(check, [RSEQ_OUT_PUBLICRNA, "-G", "mode"])
  assert checkres.exit_code == 0

def test_check_bam():
  if os.path.exists(RSEQ_OUT_BAM):
    shutil.rmtree(RSEQ_OUT_BAM)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_BAM, BAMSAMPS, "-g", "hg38"])
  checkres = runner.invoke(check, [RSEQ_OUT_BAM])
  assert checkres.exit_code == 0
  

def test_check_fq():
  if os.path.exists(RSEQ_OUT_FQ):
    shutil.rmtree(RSEQ_OUT_FQ)
  runner = CliRunner()
  buildres = runner.invoke(build, [RSEQ_OUT_FQ, FQSAMPS, "-g", "hg38"])
  checkres = runner.invoke(check, [RSEQ_OUT_FQ])
  assert checkres.exit_code == 0

